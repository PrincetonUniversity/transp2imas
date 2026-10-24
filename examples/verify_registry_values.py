#!/usr/bin/python
"""Read back registered mappings from HDF5 and compare them with values computed from the CDF.

One real conversion of ITER 38530P94 (5 time slices) per Data Dictionary (4.1.1 and 3.39.0) through
the actual converter; each ``test_<registry id>`` recomputes the expected IMAS values with numpy
straight from the CDF variables (unit factor, sum, interpolation as stated in the registry) and
compares them with what is stored.  ``mapping_registry.yaml`` lists this script in ``tests`` for
exactly the ids in ``COVERED`` (checked both ways), so a registry entry cannot claim this test
without a test body, and a test cannot exist unregistered.

Sign conventions (ip_sign, q_sign, toroidal_rotation_sign) are tested as a constant sign and an
exact magnitude; the choice of sign itself is covered by the COCOS verifiers.

    python examples/verify_registry_values.py
"""

import re
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.io import netcdf_file

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
ITER_CASE = ROOT / 'examples' / 'data' / 'iter' / '38530P94.CDF'
T16_CASE = ROOT / 'examples' / 'data' / 'd3d' / '147634T16' / '147634T16.CDF'
# case label -> (CDF, strict).  Every registered mapping is verified on the strict ITER case (a missing
# variable or IDS there is a failure).  Other cases run the same tests where their data allow and skip,
# with the reason, what they do not contain; DIII-D 147634T16 is the bundled MMM07 reference.
CASES = {'ITER': (ITER_CASE, True), 'D3D_T16': (T16_CASE, False)}
IDS_NAMES = ('wall', 'tf', 'ic_antennas', 'ec_launchers', 'nbi', 'radiation', 'core_sources', 'core_profiles',
             'core_transport', 'edge_sources', 'equilibrium', 'summary', 'distribution_sources', 'waves')
DD_VERSIONS = ('4.1.1', '3.39.0')
RTOL = 1e-6

try:
    import imas
except ImportError:      # pragma: no cover - environment dependent
    imas = None

COVERED = set()


def covers(*ids):
    def mark(function):
        COVERED.update(ids)
        return function
    return mark


def convert_and_open(cdf, dd_version, outdir):
    result = subprocess.run(
        [sys.executable, str(ROOT / 'convert.py'), '-i', str(cdf), '-b', 'HDF5', '-u', str(outdir),
         '--dd-version', dd_version, '-new', '--ids', ','.join(IDS_NAMES)],
        cwd=ROOT, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(result.stdout[-1500:] + result.stderr[-1500:])
    match = re.search(r'ids stored in (\S+)', result.stdout)
    if not match:
        raise RuntimeError('convert.py did not report the stored entry:\n' + result.stdout[-1500:])
    path = Path(match.group(1))
    database, _major, shot, run = path.parts[-4:]
    entry = imas.DBEntry(imas.ids_defs.HDF5_BACKEND, database, int(shot), int(run), str(path.parents[3]),
                         dd_version=dd_version)
    entry.open()
    return entry


def arr(node):
    return np.asarray(node, dtype=float)


def stack(slices, getter):
    """Per-time arrays stacked to (nt, nx); ``None`` when any slice is empty."""
    rows = [arr(getter(s)) for s in slices]
    if not rows or any(r.size == 0 for r in rows):
        return None
    return np.array(rows)


def constant_sign(test, values):
    """Magnitude check helper: returns +1 or -1 if ``values`` is one sign everywhere (zeros ignored)."""
    nonzero = values[values != 0]
    test.assertTrue(nonzero.size, 'all-zero array cannot show a sign')
    sign = np.sign(nonzero[0])
    test.assertTrue(np.all(np.sign(nonzero) == sign), 'sign changes within the profile')
    return sign


class RegistryValueTests:
    """Mixin; the concrete classes below bind one Data Dictionary version each."""

    dd_version = None
    case_path = None
    strict = True
    entry = None
    cdf = None
    tmp = None

    @classmethod
    def setUpClass(cls):
        if imas is None:
            raise unittest.SkipTest('IMAS-Python is not available')
        if not cls.case_path.is_file():
            raise unittest.SkipTest(f'{cls.case_path.name} is not installed')
        cls.tmp = tempfile.TemporaryDirectory()
        cls.entry = convert_and_open(cls.case_path, cls.dd_version, cls.tmp.name)
        with netcdf_file(str(cls.case_path), 'r', mmap=False) as nc:
            cls.cdf = {name: np.array(var.data, dtype=float) for name, var in nc.variables.items()
                       if var.data.dtype.kind in 'fiu'}
        cls.ids = {}

    @classmethod
    def tearDownClass(cls):
        if cls.entry is not None:
            cls.entry.close()
        if cls.tmp is not None:
            cls.tmp.cleanup()

    def get(self, name):
        if name not in self.ids:
            try:
                self.ids[name] = self.entry.get(name)
            except Exception as exc:                  # the case did not produce this IDS
                if self.strict:
                    raise
                raise unittest.SkipTest(f'{name} was not written for {self.case_path.name} ({type(exc).__name__})')
        return self.ids[name]

    def var(self, name):
        if name not in self.cdf and not self.strict:
            raise unittest.SkipTest(f'{name} is not in {self.case_path.name}')
        self.assertIn(name, self.cdf, f'{name} missing from {self.case_path.name}')
        return self.cdf[name]

    def names(self, pattern):
        regex = re.compile(pattern)
        found = sorted(n for n in self.cdf if regex.fullmatch(n))
        if not found and not self.strict:
            raise unittest.SkipTest(f'no variable matching {pattern} in {self.case_path.name}')
        self.assertTrue(found, f'no CDF variable matches {pattern}')
        return found

    def close(self, actual, expected, label):
        actual, expected = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
        self.assertEqual(actual.shape, expected.shape, f'{label}: shape')
        np.testing.assert_allclose(actual, expected, rtol=RTOL, atol=1e-12 * max(1.0, np.abs(expected).max()),
                                   err_msg=label)

    # ---- helpers over IDSs with several sources / processes / ions --------------------------

    def source_arrays(self, getter):
        out = []
        for source in self.get('core_sources').source:
            values = stack(source.profiles_1d, getter)
            if values is not None:
                out.append(values)
        return out

    def assert_some_source(self, getter, expected, label):
        for values in self.source_arrays(getter):
            if values.shape == expected.shape and np.allclose(values, expected, rtol=RTOL, atol=1e-12 * max(1.0, np.abs(expected).max())):
                return
        self.fail(f'{label}: no core_sources source carries the expected profile')

    def ion_arrays(self, getter):
        out = []
        for ion_index in range(max(len(s.ion) for s in self.get('core_profiles').profiles_1d)):
            rows = []
            for s in self.get('core_profiles').profiles_1d:
                rows.append(arr(getter(s.ion[ion_index])) if ion_index < len(s.ion) else np.array([]))
            if rows and all(r.size for r in rows):
                out.append(np.array(rows))
        return out

    def assert_some_ion(self, getter, expected, label):
        for values in self.ion_arrays(getter):
            if values.shape == expected.shape and np.allclose(values, expected, rtol=RTOL, atol=1e-12 * max(1.0, np.abs(expected).max())):
                return
        self.fail(f'{label}: no ion carries the expected profile')

    def expected_betas(self):
        """The betas by the Data Dictionary definitions, recomputed with numpy straight from the CDF.

        beta_tor = 2 mu0 <PTOWB>_V / b0^2, beta_pol = 4 int PTOWB dV / (r0 mu0 Ip^2), beta_tor_norm =
        100 beta_tor a b0 / Ip[MA], with r0 the middle of the limiter, b0 = BZXR/r0 and a the boundary
        minor radius rebuilt from the RMCB/RMSB moments; magnitudes only (the callers compare up to the
        global sign).  Independent of ``plasma_betas``.
        """
        mu0 = 4.0e-7 * np.pi
        rlim = self.var('RLIM')[0] * 0.01
        rlim = rlim[np.isfinite(rlim) & (rlim > 0)]
        r0 = 0.5 * (rlim.max() + rlim.min())
        time, time3 = self.var('TIME'), self.var('TIME3')
        on_time3 = lambda series: np.interp(time3, time, series)
        b0 = on_time3(np.abs(self.var('BZXR')) * 0.01) / r0
        ip = on_time3(np.abs(self.var('PCURC')))
        dvol = self.var('DVOL') * 1e-6
        pdv = np.sum(self.var('PTOWB') * dvol, axis=1)
        beta_tor = 2 * mu0 * pdv / np.sum(dvol, axis=1) / b0 ** 2
        beta_pol = 4 * pdv / (r0 * mu0 * ip ** 2)
        theta = np.linspace(0.0, 2 * np.pi, 720, endpoint=False)
        self.var('RMCB0')
        radius, k = np.zeros((len(time3), theta.size)), 0
        while f'RMCB{k}' in self.cdf:
            radius += on_time3(self.cdf[f'RMCB{k}'] * 0.01)[:, None] * (np.cos(k * theta)[None, :] if k else 1.0)
            if k and f'RMSB{k}' in self.cdf:
                radius += on_time3(self.cdf[f'RMSB{k}'] * 0.01)[:, None] * np.sin(k * theta)[None, :]
            k += 1
        minor = 0.5 * (radius.max(axis=1) - radius.min(axis=1))
        return {'beta_tor': beta_tor, 'beta_pol': beta_pol, 'beta_tor_norm': 100 * beta_tor * minor * b0 / (ip / 1e6)}

    # ---- hardware ---------------------------------------------------------------------------

    @covers('wall_limiter_outline')
    def test_wall_limiter_outline(self):
        unit = self.get('wall').description_2d[0].limiter.unit[0]
        r, z = self.var('RLIM')[0] * 0.01, self.var('YLIM')[0] * 0.01
        if self.dd_version.startswith('3.'):
            self.assertEqual(int(unit.closed), 1)
            r, z = r[:-1], z[:-1]
        elif not (np.isclose(r[0], r[-1]) and np.isclose(z[0], z[-1])):
            r, z = np.append(r, r[0]), np.append(z, z[0])
        self.close(unit.outline.r, r, 'outline.r')
        self.close(unit.outline.z, z, 'outline.z')

    @covers('tf_vacuum_field')
    def test_tf_vacuum_field(self):
        tf = self.get('tf')
        node = tf.b_field_phi_vacuum_r if hasattr(tf, 'b_field_phi_vacuum_r') else tf.b_field_tor_vacuum_r
        actual = arr(node.data)
        self.close(np.abs(actual), np.abs(self.var('BZXR')) * 0.01, 'b_field_vacuum_r')
        constant_sign(self, actual)

    @covers('ic_antenna_frequency', 'ic_antenna_power', 'ic_total_power')
    def test_ic_antennas(self):
        ic = self.get('ic_antennas')
        powers = self.names(r'PICHA\d')
        self.assertEqual(len(ic.antenna), len(powers))
        for i, name in enumerate(powers):
            number = name[len('PICHA'):]
            self.close(ic.antenna[i].frequency.data, self.var(f'FREQA{number}'), f'FREQA{number}')
            self.close(ic.antenna[i].power_launched.data, self.var(name), name)
        self.close(ic.power_launched.data, self.var('PICHTOT'), 'PICHTOT')

    @covers('ec_beam_power', 'ec_steering_angles')
    def test_ec_beams(self):
        ec = self.get('ec_launchers')
        powers = self.names(r'PECIN\d')
        self.assertEqual(len(ec.beam), len(powers))
        for i, name in enumerate(powers):
            number = name[len('PECIN'):]
            beam = ec.beam[i]
            self.close(beam.power_launched.data, self.var(name), name)
            phi = np.radians(self.var(f'PHAIECH{number}'))
            theta = np.radians(self.var(f'THETECH{number}'))
            # Prater: TRANSP lab-frame aiming angles -> angles in the launcher frame
            beta = np.arcsin(np.sin(theta) * np.sin(phi - np.pi))
            cos_beta = -np.cos(beta)
            alpha = np.where(np.abs(cos_beta) > 1e-10, np.arcsin(np.cos(np.pi - theta) / np.where(
                np.abs(cos_beta) > 1e-10, cos_beta, 1.0)), 0.0)
            self.close(beam.steering_angle_tor, beta, f'steering_angle_tor {number}')
            self.close(beam.steering_angle_pol, alpha, f'steering_angle_pol {number}')

    @covers('nbi_unit_power', 'nbi_unit_energy')
    def test_nbi_units(self):
        nbi = self.get('nbi')
        powers = self.names(r'PINJ\d\d')
        self.assertEqual(len(nbi.unit), len(powers))
        for i, name in enumerate(powers):
            number = name[len('PINJ'):]
            self.close(nbi.unit[i].power_launched.data, self.var(name), name)
            self.close(nbi.unit[i].energy.data, self.var(f'EINJ{number}_E1'), f'EINJ{number}_E1')

    # ---- radiation --------------------------------------------------------------------------

    @covers('radiation_electron_emissivity')
    def test_radiation_electron_emissivity(self):
        processes = [stack(p.profiles_1d, lambda s: s.electrons.emissivity) for p in self.get('radiation').process]
        for name in ('PRAD_BR', 'PRAD_CY'):
            expected = self.var(name) * 1e6
            self.assertTrue(any(v is not None and v.shape == expected.shape and np.allclose(v, expected, rtol=RTOL)
                                for v in processes), f'no radiation process carries {name}')

    @covers('radiation_line_emissivity')
    def test_radiation_line_emissivity(self):
        processes = [stack(p.profiles_1d, lambda s: s.emissivity_ion_total) for p in self.get('radiation').process]
        expected = self.var('PRAD_LI') * 1e6
        self.assertTrue(any(v is not None and v.shape == expected.shape and np.allclose(v, expected, rtol=RTOL)
                            for v in processes), 'no radiation process carries PRAD_LI as emissivity_ion_total')

    def assert_some_line_ion(self, expected, label, states=False):
        found = []
        for process in self.get('radiation').process:
            slices = process.profiles_1d
            width = max((len(s.ion) for s in slices), default=0)
            for i in range(width):
                ions = [s for s in slices if i < len(s.ion)]
                if len(ions) != len(slices):
                    continue
                if not states:
                    found.append(stack(ions, lambda s: s.ion[i].emissivity))
                    continue
                for j in range(max(len(s.ion[i].state) for s in ions)):
                    if all(j < len(s.ion[i].state) for s in ions):
                        found.append(stack(ions, lambda s: s.ion[i].state[j].emissivity))
        self.assertTrue(any(v is not None and v.shape == expected.shape and np.allclose(v, expected, rtol=RTOL)
                            for v in found), f'{label}: no line_radiation ion carries the expected profile')

    @covers('radiation_line_ion_emissivity')
    def test_radiation_line_ion_emissivity(self):
        for name in self.names(r'PRLS_(?!TOK$|SINGL$)[A-Z0-9]+'):
            self.assert_some_line_ion(self.var(name) * 1e6, name)

    @covers('radiation_line_state_emissivity')
    def test_radiation_line_state_emissivity(self):
        for name in self.names(r'PRLX_(?!SINGL)[A-Z0-9]+_[0-9]+'):
            self.assert_some_line_ion(self.var(name) * 1e6, name, states=True)

    # ---- core_sources -----------------------------------------------------------------------

    @covers('core_sources_nbi_electron_energy')
    def test_nbi_electron_energy(self):
        self.assert_some_source(lambda s: s.electrons.energy, self.var('PBE') * 1e6, 'PBE')

    @covers('core_sources_nbi_ion_energy')
    def test_nbi_ion_energy(self):
        expected = (self.var('PBI') + self.var('PBTH')) * 1e6
        self.assert_some_source(lambda s: s.total_ion_energy, expected, 'PBI+PBTH')

    @covers('core_sources_nbi_electron_particles')
    def test_nbi_electron_particles(self):
        self.assert_some_source(lambda s: s.electrons.particles, self.var('SBE') * 1e6, 'SBE')

    @covers('core_sources_ohmic_energy')
    def test_ohmic_energy(self):
        self.assert_some_source(lambda s: s.electrons.energy, self.var('POH') * 1e6, 'POH')

    @covers('core_sources_bootstrap_current')
    def test_bootstrap_current(self):
        cs = self.get('core_sources')
        b0 = np.abs(arr(cs.vacuum_toroidal_field.b0))
        expected = self.var('CURBS') * 1e4 * self.var('CONPLJB')
        for values in self.source_arrays(lambda s: s.j_parallel):
            scaled = values * b0[:, None]
            if scaled.shape == expected.shape and (np.allclose(scaled, expected, rtol=1e-6, atol=1e-6 * np.abs(expected).max())
                                                   or np.allclose(scaled, -expected, rtol=1e-6, atol=1e-6 * np.abs(expected).max())):
                return
        self.fail('no core_sources source carries CURBS*CONPLJB/b0 as j_parallel')

    @covers('core_sources_radiation_sink')
    def test_radiation_sink(self):
        for name in ('PRAD_BR', 'PRAD_CY', 'PRAD_LI'):
            self.assert_some_source(lambda s: s.electrons.energy, -self.var(name) * 1e6, f'-{name}')

    def assert_some_source_ion(self, getter, expected, label):
        """Some core_sources source carries ``expected`` as ion[j].<quantity> for one species index j."""
        for source in self.get('core_sources').source:
            slices = source.profiles_1d
            for j in range(max((len(s.ion) for s in slices), default=0)):
                if all(j < len(s.ion) for s in slices):
                    values = stack(slices, lambda s: getter(s.ion[j]))
                    if values is not None and values.shape == expected.shape and np.allclose(values, expected, rtol=RTOL):
                        return
        self.fail(f'{label}: no core_sources ion carries the expected profile')

    def species_sum(self, pattern, factor=1e6):
        names = self.names(pattern)
        return sum(self.cdf[n] for n in names) * factor

    def skip_before_dd4(self, what):
        if self.dd_version.startswith('3.'):
            raise unittest.SkipTest(f'{what} does not exist in DD {self.dd_version}')

    @covers('core_sources_ec_electron_energy')
    def test_ec_electron_energy(self):
        self.assert_some_source(lambda s: s.electrons.energy, self.var('PEECH') * 1e6, 'PEECH')

    @covers('core_sources_ic_electron_energy')
    def test_ic_electron_energy(self):
        self.assert_some_source(lambda s: s.electrons.energy, self.var('PEICH') * 1e6, 'PEICH')

    @covers('core_sources_ic_ion_energy')
    def test_ic_ion_energy(self):
        self.assert_some_source(lambda s: s.total_ion_energy, self.var('PIICH') * 1e6, 'PIICH')

    @covers('core_sources_fusion_electron_energy')
    def test_fusion_electron_energy(self):
        self.assert_some_source(lambda s: s.electrons.energy, self.species_sum(r'PFE_[A-Z0-9]+'), 'sum PFE_*')

    @covers('core_sources_fusion_ion_energy')
    def test_fusion_ion_energy(self):
        expected = self.species_sum(r'PFI_[A-Z0-9]+') + self.species_sum(r'PFTH_[A-Z0-9]+')
        self.assert_some_source(lambda s: s.total_ion_energy, expected, 'sum PFI_* + PFTH_*')

    @covers('core_sources_fusion_ion_particles')
    def test_fusion_ion_particles(self):
        self.skip_before_dd4('total_ion_particles')
        self.assert_some_source(lambda s: s.total_ion_particles, self.species_sum(r'SFTH_[A-Z0-9]+'), 'sum SFTH_*')

    @covers('core_sources_fusion_ion_species')
    def test_fusion_ion_species(self):
        for name in self.names(r'SFTH_[A-Z0-9]+'):
            self.assert_some_source_ion(lambda i: i.particles, self.var(name) * 1e6, name)
            self.assert_some_source_ion(lambda i: i.energy, self.var('PFTH_' + name[5:]) * 1e6, 'PFTH_' + name[5:])

    @covers('core_sources_nbi_ion_particles')
    def test_nbi_ion_particles(self):
        self.skip_before_dd4('total_ion_particles')
        self.assert_some_source(lambda s: s.total_ion_particles, self.var('SBTH') * 1e6, 'SBTH')

    @covers('core_sources_nbi_ion_species')
    def test_nbi_ion_species(self):
        for name in self.names(r'SBTH_[A-Z0-9]'):
            self.assert_some_source_ion(lambda i: i.particles, self.var(name) * 1e6, name)
            self.assert_some_source_ion(lambda i: i.energy, self.var('PBTH_' + name[5:]) * 1e6, 'PBTH_' + name[5:])

    @covers('core_sources_equipartition')
    def test_equipartition(self):
        qie = self.var('QIE') * 1e6
        self.assert_some_source(lambda s: s.electrons.energy, -qie, '-QIE')
        self.assert_some_source(lambda s: s.total_ion_energy, qie, '+QIE')

    @covers('core_sources_cold_neutral_particles')
    def test_cold_neutral_particles(self):
        expected = (self.var('SCEW') + self.var('SCEV')) * 1e6
        self.assert_some_source(lambda s: s.electrons.particles, expected, 'SCEW+SCEV')

    @covers('core_sources_cold_neutral_energy')
    def test_cold_neutral_energy(self):
        self.assert_some_source(lambda s: s.electrons.energy, -self.var('PION') * 1e6, '-PION')

    def grid_poloidal_flux(self, slices_of_ids, label):
        """Every 1D profile set of an IDS must carry rho_pol_norm and psi built from PLFLX and PSI0_TR."""
        x, xb, plflx = self.var('X'), self.var('XB'), self.var('PLFLX') * 2 * np.pi
        axis = self.var('PSI0_TR') * 2 * np.pi
        rel = np.array([np.interp(x[i], np.concatenate(([0.0], xb[i])), np.concatenate(([0.0], plflx[i])))
                        for i in range(len(x))])
        rho_pol = np.sqrt(np.clip(rel / plflx[:, -1:], 0.0, None))
        checked = 0
        for slices in slices_of_ids:
            values = stack(slices, lambda s: s.grid.rho_pol_norm)
            if values is None:
                continue
            self.close(values, rho_pol, f'{label} rho_pol_norm')
            self.close_global_sign(stack(slices, lambda s: s.grid.psi), rel + axis[:, None], f'{label} psi')
            checked += 1
        self.assertTrue(checked, f'{label}: no profile set carries the poloidal-flux grid')

    @covers('grid_poloidal_flux_core_sources')
    def test_grid_poloidal_flux_core_sources(self):
        self.grid_poloidal_flux([s.profiles_1d for s in self.get('core_sources').source], 'core_sources')

    @covers('grid_poloidal_flux_radiation')
    def test_grid_poloidal_flux_radiation(self):
        self.grid_poloidal_flux([p.profiles_1d for p in self.get('radiation').process], 'radiation')

    @covers('grid_poloidal_flux_waves')
    def test_grid_poloidal_flux_waves(self):
        self.grid_poloidal_flux([w.profiles_1d for w in self.get('waves').coherent_wave], 'waves')

    @covers('grid_poloidal_flux_distribution_sources')
    def test_grid_poloidal_flux_distribution_sources(self):
        self.grid_poloidal_flux([s.profiles_1d for s in self.get('distribution_sources').source], 'distribution_sources')

    # ---- core_profiles ----------------------------------------------------------------------

    def profiles(self, getter):
        return stack(self.get('core_profiles').profiles_1d, getter)

    @covers('cp_electron_temperature')
    def test_te(self):
        self.close(self.profiles(lambda s: s.electrons.temperature), self.var('TE'), 'TE')

    @covers('cp_electron_density')
    def test_ne(self):
        self.close(self.profiles(lambda s: s.electrons.density), self.var('NE') * 1e6, 'NE')

    @covers('cp_ion_average_temperature')
    def test_t_i_average(self):
        self.close(self.profiles(lambda s: s.t_i_average), self.var('TI'), 'TI')

    @covers('cp_ion_temperature')
    def test_ion_temperature(self):
        checked, temperatures = 0, set()
        for index in range(max(len(s.ion) for s in self.get('core_profiles').profiles_1d)):
            slices = [s for s in self.get('core_profiles').profiles_1d if index < len(s.ion)]
            if not all(len(s.ion[index].density_thermal) for s in slices):
                continue                    # fast-ion populations: see cp_fast_ion_temperature
            actual = stack(slices, lambda s: s.ion[index].temperature)
            # thermal species carry TI, impurities carry TX (the same in some CDFs, different in others)
            matches = [name for name in ('TI', 'TX') if name in self.cdf and actual.shape == self.cdf[name].shape
                       and np.allclose(actual, self.cdf[name], rtol=RTOL, atol=1e-9)]
            self.assertTrue(matches, f'ion {index} temperature equals neither TI nor TX')
            temperatures.update(matches)
            checked += 1
        self.assertTrue(checked)
        self.assertIn('TI', temperatures, 'no thermal species carries TI')

    @covers('cp_fast_ion_temperature')
    def test_fast_ion_temperature(self):
        for name in ('EBEAM_D', 'EFUSN_T', 'EFUSN_3', 'EFUSN_4'):
            self.assert_some_ion(lambda ion: ion.temperature, self.var(name), name)

    @covers('cp_ion_density_thermal')
    def test_ion_density_thermal(self):
        for name in ('ND', 'NT', 'NHE3', 'NHE4'):
            if name in self.cdf:
                self.assert_some_ion(lambda ion: ion.density_thermal, self.var(name) * 1e6, name)

    @covers('cp_ion_radial_velocity')
    def test_ion_radial_velocity(self):
        for name in ('VNDNC_D', 'VNDNC_T', 'VNDNC_HE3', 'VNDNC_HE4'):
            if name not in self.cdf:
                continue
            if np.any(self.var(name)):
                self.assert_some_ion(lambda ion: ion.velocity.radial, self.var(name) * 0.01, name)

    @covers('cp_ion_rotation_frequency')
    def test_ion_rotation_frequency(self):
        omega = self.var('OMEGA')
        sign = constant_sign(self, self.ion_arrays(lambda ion: ion.rotation_frequency_tor)[0])
        flip = sign * np.sign(omega[omega != 0][0])      # +1 when the stored sign equals OMEGA's
        self.assert_some_ion(lambda ion: ion.rotation_frequency_tor, flip * omega, 'OMEGA')

    @covers('cp_ion_to_electron_density')
    def test_ni_over_ne(self):
        self.close(self.profiles(lambda s: s.n_i_total_over_n_e), self.var('NI') / self.var('NE'), 'NI/NE')

    @covers('cp_pressure_thermal')
    def test_pressure_thermal(self):
        self.close(self.profiles(lambda s: s.pressure_thermal), self.var('PPLAS'), 'PPLAS')

    @covers('cp_pressure_electron')
    def test_pressure_electron(self):
        self.close(self.profiles(lambda s: s.electrons.pressure_thermal), self.var('UE') * 1e6 * 2 / 3, 'UE')

    @covers('cp_pressure_ion_total')
    def test_pressure_ion_total(self):
        self.close(self.profiles(lambda s: s.pressure_ion_total), self.var('UI') * 1e6 * 2 / 3, 'UI')

    @covers('cp_pressure_anisotropy')
    def test_pressure_anisotropy(self):
        thermal = self.var('PPLAS')
        self.close(self.profiles(lambda s: s.pressure_perpendicular), thermal + self.var('UFASTPP') * 1e6, 'perp')
        self.close(self.profiles(lambda s: s.pressure_parallel), thermal + 2 * self.var('UFASTPA') * 1e6, 'par')

    @covers('cp_q_profile')
    def test_q_profile(self):
        x, xb, q, q0 = self.var('X'), self.var('XB'), self.var('Q'), self.var('Q0')
        expected = np.array([np.interp(x[i], np.concatenate(([0.0], xb[i])), np.concatenate(([q0[i]], q[i])))
                             for i in range(len(x))])
        actual = self.profiles(lambda s: s.q)
        constant_sign(self, actual)
        self.close(np.abs(actual), np.abs(expected), '|q|')

    @covers('cp_zeff')
    def test_zeff(self):
        self.close(self.profiles(lambda s: s.zeff), self.var('ZEFFP'), 'ZEFFP')

    @covers('cp_magnetic_shear')
    def test_magnetic_shear(self):
        self.close(self.profiles(lambda s: s.magnetic_shear), self.var('SHAT'), 'SHAT')

    @covers('cp_parallel_conductivity')
    def test_parallel_conductivity(self):
        eta = self.var('ETA_USE') * 0.01
        expected = np.where(np.abs(eta) > 0, 1.0 / np.where(eta == 0, 1.0, eta), 0.0)
        self.close(self.profiles(lambda s: s.conductivity_parallel), expected, 'conductivity_parallel')

    @covers('cp_phi_potential')
    def test_phi_potential(self):
        x, xb, v = self.var('X'), self.var('XB'), self.var('VRPOT')
        expected = np.array([np.interp(x[i], np.concatenate(([0.0], xb[i])) if xb[i][0] > 0 else xb[i],
                                       np.concatenate(([v[i][0]], v[i])) if xb[i][0] > 0 else v[i])
                             for i in range(len(x))])
        self.close(self.profiles(lambda s: s.phi_potential), expected, 'phi_potential')


    # ---- core_transport ---------------------------------------------------------------------

    def model(self, name):
        for model in self.get('core_transport').model:
            if str(model.identifier.name) == name:
                return model
        self.fail(f'no core_transport model named {name}')

    @staticmethod
    def momentum(profile):
        return profile.momentum_phi if hasattr(profile, 'momentum_phi') else profile.momentum_tor

    def model_stack(self, name, getter):
        return stack(self.model(name).profiles_1d, getter)

    def surface_flux(self, div_name):
        """Unit-analytic flux of a per-cm^3 divergence: sum(div*DVOL[cm^3]) / SURF[cm^2] -> per m^2."""
        integrated = np.cumsum(self.var(div_name) * self.var('DVOL'), axis=1) * 1e4
        surf = self.var('SURF')
        return np.divide(integrated, surf, out=np.zeros_like(integrated), where=np.abs(surf) > 0)

    def assert_some_model_ion(self, model_name, getter, expected, label):
        model = self.model(model_name)
        for index in range(max(len(s.ion) for s in model.profiles_1d)):
            values = stack([s for s in model.profiles_1d if index < len(s.ion)], lambda s: getter(s.ion[index]))
            if values is not None and values.shape == expected.shape and np.allclose(
                    values, expected, rtol=RTOL, atol=1e-9 * max(1.0, np.abs(expected).max())):
                return
        self.fail(f'{label}: no ion of model {model_name} carries the expected profile')

    @covers('ct_combined_electron_particle_dv')
    def test_combined_electron_particle_dv(self):
        self.close(self.model_stack('combined', lambda s: s.electrons.particles.d), self.var('DIFFE') * 1e-4, 'DIFFE')
        self.close(self.model_stack('combined', lambda s: s.electrons.particles.v), self.var('VELE') * 0.01, 'VELE')

    @covers('ct_combined_electron_heat_dv')
    def test_combined_electron_heat_dv(self):
        self.close(self.model_stack('combined', lambda s: s.electrons.energy.d), self.var('CONDE') * 1e-4, 'CONDE')
        self.close(self.model_stack('combined', lambda s: s.electrons.energy.v), self.var('VEL_TE') * 0.01, 'VEL_TE')

    @covers('ct_combined_ion_heat_dv')
    def test_combined_ion_heat_dv(self):
        self.close(self.model_stack('combined', lambda s: s.total_ion_energy.d), self.var('CONDI') * 1e-4, 'CONDI')
        self.close(self.model_stack('combined', lambda s: s.total_ion_energy.v), self.var('VEL_TI') * 0.01, 'VEL_TI')

    @covers('ct_combined_momentum_dv')
    def test_combined_momentum_dv(self):
        self.close(self.model_stack('combined', lambda s: self.momentum(s).d), self.var('CHPHI') * 1e-4, 'CHPHI')
        self.close(self.model_stack('combined', lambda s: self.momentum(s).v), self.var('VMO_NET') * 0.01, 'VMO_NET')

    @covers('ct_combined_conductivity')
    def test_combined_conductivity(self):
        x, xb, eta = self.var('X'), self.var('XB'), self.var('ETA_USE') * 0.01
        eta_b = np.array([np.interp(xb[i], x[i], eta[i]) for i in range(len(x))])
        expected = np.divide(1.0, eta_b, out=np.zeros_like(eta_b), where=np.abs(eta_b) > 0)
        self.close(self.model_stack('combined', lambda s: s.conductivity_parallel), expected, 'conductivity_parallel')

    @covers('ct_combined_species_dv')
    def test_combined_species_dv(self):
        for sp, dname, vname in (('D', 'DIFFD', 'VELD'), ('T', 'DIFFT', 'VELT'), ('HE3', 'DIFFHE3', 'VELHE3'),
                                 ('HE4', 'DIFFHE4', 'VELHE4')):
            if dname in self.cdf:
                self.assert_some_model_ion('combined', lambda i: i.particles.d, self.var(dname) * 1e-4, dname)
            if vname in self.cdf:
                self.assert_some_model_ion('combined', lambda i: i.particles.v, self.var(vname) * 0.01, vname)

    @covers('ct_observed_fluxes')
    def test_observed_fluxes(self):
        self.close(self.model_stack('combined', lambda s: s.electrons.particles.flux), self.surface_flux('EPTR_OBS'), 'EPTR_OBS')
        self.close(self.model_stack('combined', lambda s: s.electrons.energy.flux), self.surface_flux('EETR_OBS'), 'EETR_OBS')
        self.close(self.model_stack('combined', lambda s: s.total_ion_energy.flux), self.surface_flux('IETR_OBS'), 'IETR_OBS')

    @covers('momentum_flux_observed')
    def test_observed_momentum_flux(self):
        actual = self.model_stack('combined', lambda s: self.momentum(s).flux)
        expected = self.surface_flux('AMTR_OBS')
        # the flux changes sign with radius; only one global sign (ip_sign) may differ from the CDF
        sign = 1.0 if np.allclose(actual, expected, rtol=RTOL, atol=1e-9) else -1.0
        self.close(actual, sign * expected, 'AMTR_OBS flux up to one global sign')

    @covers('ct_observed_species_flux')
    def test_observed_species_flux(self):
        for name in ('PTRH_OBS', 'PTRD_OBS', 'PTRT_OBS', 'PTRHE3_OBS', 'PTRHE4_OBS'):
            if name in self.cdf:
                self.assert_some_model_ion('combined', lambda i: i.particles.flux, self.surface_flux(name), name)

    @covers('ct_solver_fluxes')
    def test_solver_fluxes(self):
        channels = ((lambda s: s.electrons.particles.flux, 'EPTR_MOD'),
                    (lambda s: s.electrons.energy.flux, 'EETR_MOD'),
                    (lambda s: s.total_ion_energy.flux, 'IETR_MOD'))
        active = [(getter, name) for getter, name in channels if np.any(self.var(name))]
        solver = [m for m in self.get('core_transport').model if str(m.identifier.name) == 'transport_solver']
        if not active and not any(np.any(self.cdf[n]) for n in self.cdf if re.fullmatch(r'(AMTR|PTR[A-Z0-9]+)_MOD', n)):
            self.assertEqual(solver, [], 'a solver model whose channels are all zero must not be written')
            return
        for getter, name in active:
            self.close(self.model_stack('transport_solver', getter), self.surface_flux(name), name)
        for getter, name in channels:
            if (getter, name) not in active:        # all zero: the solver did not solve this channel
                self.assertIsNone(self.model_stack('transport_solver', getter),
                                  f'{name} is identically zero in the CDF and must not be written as a zero flux')

    @covers('ct_solver_species_flux')
    def test_solver_species_flux(self):
        for name in ('PTRH_MOD', 'PTRD_MOD', 'PTRT_MOD', 'PTRHE3_MOD', 'PTRHE4_MOD'):
            if name in self.cdf and np.any(self.cdf[name]):
                self.assert_some_model_ion('transport_solver', lambda i: i.particles.flux, self.surface_flux(name), name)

    def neoclassical_check(self, model_name, prefix_particles, prefix_heat):
        names = {}
        for prefix in (prefix_particles, prefix_heat):
            if prefix:
                names[prefix] = [n for n in self.cdf if re.fullmatch(prefix + r'[A-Z0-9]+', n)
                                 and n[len(prefix):] not in ('X', 'TOK')]
        checked = 0
        for name in names.get(prefix_particles, ()):
            token = name[len(prefix_particles):]
            expected = self.surface_flux(name)
            if token == 'E':
                self.close(self.model_stack(model_name, lambda s: s.electrons.particles.flux), expected, name)
            elif token != 'I':           # the aggregate ion particle channel has no destination
                self.assert_some_model_ion(model_name, lambda i: i.particles.flux, expected, name)
            checked += token != 'I'
        for name in names.get(prefix_heat, ()):
            token = name[len(prefix_heat):]
            expected = self.surface_flux(name)
            if token == 'E':
                self.close(self.model_stack(model_name, lambda s: s.electrons.energy.flux), expected, name)
            elif token == 'I':
                self.close(self.model_stack(model_name, lambda s: s.total_ion_energy.flux), expected, name)
            else:
                self.assert_some_model_ion(model_name, lambda i: i.energy.flux, expected, name)
            checked += 1
        self.assertTrue(checked, f'nothing checked for {model_name}')

    @covers('ct_neoclassical_fluxes')
    def test_neoclassical_fluxes(self):
        self.neoclassical_check('neoclassical', 'GFLNC_', r'QFLNC_')

    @covers('ct_neoclassical_classical_heat')
    def test_neoclassical_classical_heat(self):
        self.neoclassical_check('neoclassical_classical_heat', None, r'QFLNCC_')

    @covers('ct_glf23_diffusivities')
    def test_glf23_diffusivities(self):
        names = ('XKEGLF23', 'XKIGLF23', 'ETPHIGLF')
        if not any(n in self.cdf and np.any(self.cdf[n]) for n in names):
            # all-zero GLF23 channels are placeholders of a run that did not use GLF23: no model is invented
            self.assertNotIn('glf23', [str(m.identifier.name) for m in self.get('core_transport').model])
            return
        self.close(self.model_stack('glf23', lambda s: s.electrons.energy.d), self.var('XKEGLF23') * 1e-4, 'XKEGLF23')
        self.close(self.model_stack('glf23', lambda s: s.total_ion_energy.d), self.var('XKIGLF23') * 1e-4, 'XKIGLF23')
        self.close(self.model_stack('glf23', lambda s: self.momentum(s).d), self.var('ETPHIGLF') * 1e-4, 'ETPHIGLF')


    @staticmethod
    def dotted(profile, path):
        """``profile.electrons.energy.d``-style attribute path (DD3 spells momentum_phi as momentum_tor)."""
        node = profile
        for part in path.split('.'):
            if part == 'momentum_phi' and not hasattr(node, 'momentum_phi'):
                part = 'momentum_tor'
            node = getattr(node, part)
        return node

    @covers('ct_mmm07_diffusivities')
    def test_mmm07_models(self):
        from theory_transport_mapping import MODEL_SPECS
        scale = {'diffusivity': 1e-4, 'velocity': 0.01}
        def active(spec):
            return any(n in self.cdf and np.any(self.cdf[n]) for n in (*spec['coefficients'], *spec['auxiliary']))
        written = [str(m.identifier.name) for m in self.get('core_transport').model]
        present = [spec for spec in MODEL_SPECS if active(spec)]
        if not present:
            if self.case_path == T16_CASE:
                self.fail('DIII-D 147634T16 is the MMM07 reference case but has no MMM07 output')
            self.skipTest(f'no MMM07 output in {self.case_path.name}')
        for spec in MODEL_SPECS:
            if spec not in present:           # all-zero channels are placeholders: no model is invented
                self.assertNotIn(spec['name'], written)
        for spec in present:
            self.assertIn(spec['name'], written)
            for source, (destination, kind) in spec['coefficients'].items():
                if source not in self.cdf:
                    continue
                actual = self.model_stack(spec['name'], lambda s: self.dotted(s, destination))
                self.assertIsNotNone(actual, f"{spec['name']}: {source} -> {destination} was not written")
                # XB-grid boundary coefficients are copied (unit conversion only), never volume-integrated
                self.close(actual, self.cdf[source] * scale[kind], f"{spec['name']} {source}")

    @covers('ct_theory_auxiliary_xml')
    def test_theory_auxiliary_xml(self):
        from theory_transport_mapping import AUXILIARY_SPECS, GLF23_AUXILIARY_SPECS, GLF23_MODEL_SPEC, MODEL_SPECS
        multiplier = {name: spec[3] for name, spec in {**AUXILIARY_SPECS, **GLF23_AUXILIARY_SPECS}.items()}
        checked = 0
        for spec in (*MODEL_SPECS, GLF23_MODEL_SPEC):
            for name in spec['auxiliary']:
                if name not in self.cdf or not np.any(self.cdf[name]):
                    continue
                model = self.model(spec['name'])
                root = ET.fromstring(str(model.code.parameters))
                # MMM models list them as <quantity>, the GLF23 model as <supplementary_quantity>
                quantity = {q.get('source'): q for tag in ('quantity', 'supplementary_quantity') for q in root.iter(tag)}
                self.assertIn(name, quantity, f"{spec['name']}: auxiliary {name} is missing from the model XML")
                node = quantity[name]
                shape = tuple(int(n) for n in node.get('shape').split())
                values = np.array(node.text.split(), dtype=float).reshape(shape)
                self.close(values, self.cdf[name] * multiplier[name], f"{spec['name']} {name} (XML)")
                checked += 1
        if not checked:
            self.skipTest(f'no auxiliary theory-model channel with data in {self.case_path.name}')

    # ---- edge_sources (GGD on the X grid; one source per neutral origin and process) ---------

    def edge_sources(self, process):
        return [src for src in self.get('edge_sources').source if str(src.identifier.name) == process]

    @staticmethod
    def ggd_stack(source, getter):
        """(nt, nx) values of a one-entry GGD quantity; ``None`` when any time slice lacks it."""
        rows = []
        for ggd in source.ggd:
            node = getter(ggd)
            if len(node) == 0 or len(node[0].values) == 0:
                return None
            rows.append(arr(node[0].values))
        return np.array(rows) if rows else None

    def assert_some_edge_source(self, process, getter, expected, label):
        for source in self.edge_sources(process):
            values = self.ggd_stack(source, getter)
            if values is not None and values.shape == expected.shape and np.allclose(
                    values, expected, rtol=RTOL, atol=1e-12 * max(1.0, np.abs(expected).max())):
                return
        self.fail(f'{label}: no {process} source carries the expected profile')

    def edge_variables(self, prefix):
        found = sorted(n for n in self.cdf if re.fullmatch(prefix + r'[A-Z0-9]', n))
        self.assertTrue(found, f'no CDF variable {prefix}?')
        return found

    def each_edge_variable(self, prefixes, process, targets, factor_by_target, recombination=None):
        names = [n for prefix in prefixes for n in self.edge_variables(prefix)]
        if recombination in self.cdf:            # re-ionisation / CX of recombination neutrals (one resolved species)
            names.append(recombination)
        for name in names:
            expected = self.var(name) * 1e6
            for getter, sign in targets:
                self.assert_some_edge_source(process, getter, sign * expected, f'{name} -> {sign:+d}')

    @covers('es_ionization_particles')
    def test_edge_ionization_particles(self):
        self.each_edge_variable(('SEGF_', 'SERC_'), 'atomic_ionization',
                                [(lambda g: g.electrons.particles, 1), (lambda g: g.ion[0].particles, 1),
                                 (lambda g: g.neutral[0].particles, -1)], None, recombination='SERECO')

    @covers('es_ionization_energy')
    def test_edge_ionization_energy(self):
        self.each_edge_variable(('PIGF_', 'PIRC_'), 'atomic_ionization',
                                [(lambda g: g.ion[0].energy, 1), (lambda g: g.neutral[0].energy, -1)], None, recombination='PIRECO')

    @covers('es_ionization_momentum')
    def test_edge_ionization_momentum(self):
        if self.dd_version.startswith('3.'):
            # DD 3.39 ties momentum/toroidal to radial, so the writer records a skip instead of
            # writing a radial component that TRANSP does not provide
            for source in self.edge_sources('atomic_ionization'):
                for ggd in source.ggd:
                    self.assertEqual(len(ggd.ion[0].momentum), 0)
            return
        names = [n for prefix in ('TQIGF_', 'TQIRC_') for n in self.edge_variables(prefix)] + [n for n in ('TQIRECO',) if n in self.cdf]
        for name in names:
            expected = self.var(name) * 1e6
            for getter, sign in ((lambda g: g.ion[0].momentum, 1), (lambda g: g.neutral[0].momentum, -1)):
                for source in self.edge_sources('atomic_ionization'):
                    rows = [arr(getter(g)[0].phi) for g in source.ggd if len(getter(g)) and len(getter(g)[0].phi)]
                    if rows and len(rows) == len(source.ggd) and np.allclose(np.array(rows), sign * expected, rtol=RTOL, atol=1e-9):
                        break
                else:
                    self.fail(f'{name}: no atomic_ionization source carries it with sign {sign:+d}')

    @covers('es_cx_particles')
    def test_edge_cx_particles(self):
        self.each_edge_variable(('SFCXGF_', 'SFCXRC_'), 'charge_exchange',
                                [(lambda g: g.ion[0].particles, 1), (lambda g: g.neutral[0].particles, -1)], None, recombination='SFCXRECO')

    @covers('es_cx_energy')
    def test_edge_cx_energy(self):
        self.each_edge_variable(('PCXGF_', 'PCXRC_'), 'charge_exchange',
                                [(lambda g: g.total_ion_energy, -1), (lambda g: g.neutral[0].energy, 1)], None, recombination='PCXRECO')

    @covers('es_cx_momentum')
    def test_edge_cx_momentum(self):
        if self.dd_version.startswith('3.'):
            for source in self.edge_sources('charge_exchange'):
                for ggd in source.ggd:
                    self.assertEqual(len(ggd.momentum), 0)
            return
        names = [n for prefix in ('TQCXGF_', 'TQCXRC_') for n in self.edge_variables(prefix)] + [n for n in ('TQCXRECO',) if n in self.cdf]
        for name in names:
            expected = self.var(name) * 1e6
            for getter, sign in ((lambda g: g.momentum, -1), (lambda g: g.neutral[0].momentum, 1)):
                for source in self.edge_sources('charge_exchange'):
                    rows = [arr(getter(g)[0].phi) for g in source.ggd if len(getter(g)) and len(getter(g)[0].phi)]
                    if rows and len(rows) == len(source.ggd) and np.allclose(np.array(rows), sign * expected, rtol=RTOL, atol=1e-9):
                        break
                else:
                    self.fail(f'{name}: no charge_exchange source carries it with sign {sign:+d}')

    @covers('es_recombination_particles')
    def test_edge_recombination_particles(self):
        self.names_if_present(r'S0RECO')
        singles = sorted(n for n in self.cdf if re.fullmatch(r'S0V0[A-Z0-9]+', n))
        name = singles[0] if len(singles) == 1 else 'S0RECO'
        expected = self.var(name) * 1e6
        self.assert_some_edge_source('recombination', lambda g: g.neutral[0].particles, expected, name)
        if len(singles) == 1:
            for getter in (lambda g: g.electrons.particles, lambda g: g.ion[0].particles):
                self.assert_some_edge_source('recombination', getter, -expected, f'-{name}')

    @covers('es_recombination_energy')
    def test_edge_recombination_energy(self):
        self.names_if_present(r'P0RECO')
        self.assert_some_edge_source('recombination', lambda g: g.neutral[0].energy, self.var('P0RECO') * 1e6, 'P0RECO')

    @covers('es_recombination_momentum')
    def test_edge_recombination_momentum(self):
        self.names_if_present(r'TQ0RECO')
        if self.dd_version.startswith('3.'):
            for source in self.edge_sources('recombination'):
                for ggd in source.ggd:
                    self.assertEqual(len(ggd.neutral[0].momentum), 0)
            return
        expected = self.var('TQ0RECO') * 1e6
        for source in self.edge_sources('recombination'):
            rows = [arr(g.neutral[0].momentum[0].phi) for g in source.ggd
                    if len(g.neutral) and len(g.neutral[0].momentum) and len(g.neutral[0].momentum[0].phi)]
            if rows and len(rows) == len(source.ggd) and np.allclose(np.array(rows), expected, rtol=RTOL, atol=1e-9):
                return
        self.fail('TQ0RECO: no recombination source carries it as neutral momentum')


    # ---- equilibrium ------------------------------------------------------------------------

    def eq_scalar(self, getter):
        return np.array([float(getter(ts)) for ts in self.get('equilibrium').time_slice])

    def eq_profile(self, getter):
        return stack(self.get('equilibrium').time_slice, lambda ts: getter(ts.profiles_1d))

    def close_global_sign(self, actual, expected, label):
        """Equal to ``expected`` up to one global sign (ip_sign, q_sign, psi_sign ... choose it)."""
        actual, expected = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
        self.assertEqual(actual.shape, expected.shape, f'{label}: shape')
        atol = 1e-9 * max(1.0, np.abs(expected).max())
        self.assertTrue(np.allclose(actual, expected, rtol=RTOL, atol=atol) or
                        np.allclose(actual, -expected, rtol=RTOL, atol=atol),
                        f'{label}: not equal to the CDF values up to one global sign')

    def x_to_xb(self, name, scale=1.0, axis_value=False):
        """TRANSP zone-centre (X) profile on the boundary (XB) grid, as documented in the registry."""
        x, xb, values = self.var('X'), self.var('XB'), self.var(name) * scale
        rows = []
        for i in range(len(x)):
            xx, yy = x[i], values[i]
            if axis_value and xx[0] > 0:
                xx, yy = np.concatenate(([0.0], xx)), np.concatenate(([yy[0]], yy))
            rows.append(np.interp(xb[i], xx, yy))
        return np.array(rows)

    @covers('eq_plasma_current')
    def test_eq_plasma_current(self):
        self.close_global_sign(self.eq_scalar(lambda ts: ts.global_quantities.ip), self.var('PCURC'), 'ip')

    @covers('eq_magnetic_axis')
    def test_eq_magnetic_axis(self):
        self.close(self.eq_scalar(lambda ts: ts.global_quantities.magnetic_axis.r), self.var('RAXIS') * 0.01, 'axis r')
        self.close(self.eq_scalar(lambda ts: ts.global_quantities.magnetic_axis.z), self.var('YAXIS') * 0.01, 'axis z')

    @covers('eq_psi_magnetic_axis')
    def test_eq_psi_magnetic_axis(self):
        name = 'psi_axis' if self.dd_version.startswith('3.') else 'psi_magnetic_axis'
        self.close_global_sign(self.eq_scalar(lambda ts: getattr(ts.global_quantities, name)),
                               self.var('PSI0_TR') * 2 * np.pi, name)

    @covers('eq_q_scalars')
    def test_eq_q_scalars(self):
        self.close_global_sign(self.eq_scalar(lambda ts: ts.global_quantities.q_95), self.var('Q95'), 'q_95')
        self.close_global_sign(self.eq_scalar(lambda ts: ts.global_quantities.q_axis), self.var('Q0'), 'q_axis')

    @covers('eq_li3')
    def test_eq_li3(self):
        self.close(self.eq_scalar(lambda ts: ts.global_quantities.li_3), self.var('LI_3'), 'li_3')

    @covers('eq_betas')
    def test_eq_betas(self):
        expected = self.expected_betas()
        norm = 'beta_normal' if self.dd_version.startswith('3.') else 'beta_tor_norm'
        self.close(self.eq_scalar(lambda ts: ts.global_quantities.beta_pol), expected['beta_pol'], 'beta_pol')
        self.close(self.eq_scalar(lambda ts: ts.global_quantities.beta_tor), expected['beta_tor'], 'beta_tor')
        self.close_global_sign(self.eq_scalar(lambda ts: getattr(ts.global_quantities, norm)),
                               expected['beta_tor_norm'], norm)

    @covers('eq_plasma_volume')
    def test_eq_plasma_volume(self):
        self.close(self.eq_scalar(lambda ts: ts.global_quantities.volume), self.var('PVOL') * 1e-6, 'volume')

    @covers('eq_boundary_shape')
    def test_eq_boundary_shape(self):
        for name, getter in (('ELONG', lambda ts: ts.boundary.elongation), ('TRIANG', lambda ts: ts.boundary.triangularity),
                             ('TRIANGU', lambda ts: ts.boundary.triangularity_upper),
                             ('TRIANGL', lambda ts: ts.boundary.triangularity_lower)):
            self.close(self.eq_scalar(getter), self.var(name)[:, -1], f'boundary {name} (last XB surface)')

    @covers('eq_profile_shaping')
    def test_eq_profile_shaping(self):
        for name, getter in (('ELONG', lambda p: p.elongation), ('TRIANGU', lambda p: p.triangularity_upper),
                             ('TRIANGL', lambda p: p.triangularity_lower),
                             ('SQUARE_UO', lambda p: p.squareness_upper_outer),
                             ('SQUARE_LO', lambda p: p.squareness_lower_outer)):
            self.close(self.eq_profile(getter), self.var(name), name)

    @covers('eq_energy_mhd')
    def test_eq_energy_mhd(self):
        self.close(self.eq_scalar(lambda ts: ts.global_quantities.energy_mhd),
                   1.5 * np.sum(self.var('PMHD') * self.var('DVOL') * 1e-6, axis=1), 'energy_mhd')

    @covers('eq_profile_geometric_axis')
    def test_eq_profile_geometric_axis(self):
        self.close(self.eq_profile(lambda p: p.geometric_axis.r), self.var('RMC00') * 0.01, 'geometric_axis.r')
        self.close(self.eq_profile(lambda p: p.geometric_axis.z), self.var('YMC00') * 0.01, 'geometric_axis.z')

    @covers('eq_boundary_geometric_axis')
    def test_eq_boundary_geometric_axis(self):
        slices = self.get('equilibrium').time_slice
        for label, getter in (('boundary', lambda ts: ts.boundary), ('global_quantities', lambda ts: ts.global_quantities)):
            if not hasattr(getter(slices[0]), 'geometric_axis'):
                continue
            self.close([float(getter(ts).geometric_axis.r) for ts in slices], self.var('RMCB0') * 0.01, f'{label} r')
            self.close([float(getter(ts).geometric_axis.z) for ts in slices], self.var('YMCB0') * 0.01, f'{label} z')

    @covers('eq_profile_2d_b_field')
    def test_eq_profile_2d_b_field(self):
        field = self.var('B_FIELD')
        phi_name = 'b_field_tor' if self.dd_version.startswith('3.') else 'b_field_phi'
        for it, ts in enumerate(self.get('equilibrium').time_slice):
            prof = ts.profiles_2d[0]
            nr, nz = arr(prof.b_field_r).shape
            cube = field[it].reshape((3, nz, nr))
            for index, name in enumerate(('b_field_r', phi_name, 'b_field_z')):
                self.close_global_sign(arr(getattr(prof, name)), cube[index].T, f'{name}[{it}]')

    @covers('eq_profile_psi')
    def test_eq_profile_psi(self):
        self.close_global_sign(self.eq_profile(lambda p: p.psi),
                               self.var('PLFLX2PI') + self.var('PSI0_TR')[:, None] * 2 * np.pi, 'psi')

    @covers('eq_profile_phi')
    def test_eq_profile_phi(self):
        self.close_global_sign(self.eq_profile(lambda p: p.phi), self.var('TRFLX'), 'phi')

    @covers('eq_profile_q')
    def test_eq_profile_q(self):
        self.close_global_sign(self.eq_profile(lambda p: p.q), self.var('Q'), 'q')

    @covers('eq_profile_f')
    def test_eq_profile_f(self):
        self.close_global_sign(self.eq_profile(lambda p: p.f), self.var('GFUN') * (self.var('BZXR') * 0.01)[:, None], 'f')

    @covers('eq_profile_pressure')
    def test_eq_profile_pressure(self):
        self.close(self.eq_profile(lambda p: p.pressure), self.x_to_xb('PTOWB'), 'pressure')

    @covers('eq_profile_geometry')
    def test_eq_profile_geometry(self):
        self.close(self.eq_profile(lambda p: p.volume), np.cumsum(self.var('DVOL') * 1e-6, axis=1), 'volume')
        self.close(self.eq_profile(lambda p: p.area), np.cumsum(self.var('DAREA') * 1e-4, axis=1), 'area')
        self.close(self.eq_profile(lambda p: p.surface), self.var('SURF') * 1e-4, 'surface')

    @covers('eq_profile_field_averages')
    def test_eq_profile_field_averages(self):
        self.close(self.eq_profile(lambda p: p.b_field_average), self.x_to_xb('GB1'), 'GB1')
        self.close(self.eq_profile(lambda p: p.b_field_max), self.x_to_xb('BMAX'), 'BMAX')
        self.close(self.eq_profile(lambda p: p.b_field_min), self.x_to_xb('BMIN'), 'BMIN')
        self.close(self.eq_profile(lambda p: p.gm5), self.x_to_xb('GB2'), 'gm5 (GB2)')
        self.close(self.eq_profile(lambda p: p.gm1), self.x_to_xb('GR2I', 1e4), 'gm1 (GR2I)')
        self.close(self.eq_profile(lambda p: p.gm9), self.x_to_xb('GRI', 100), 'gm9 (GRI)')
        gm7 = self.var('GXI') * 10 * self.var('DRAV') * 0.1 / self.var('DRAVFAC')
        expected = np.array([np.interp(self.var('XB')[i], self.var('X')[i], gm7[i]) for i in range(len(gm7))])
        self.close(self.eq_profile(lambda p: p.gm7), expected, 'gm7 (GXI*DRAV/DRAVFAC)')

    @covers('eq_profile_j_tor')
    def test_eq_profile_j_tor(self):
        name = 'j_tor' if self.dd_version.startswith('3.') else 'j_phi'
        self.close_global_sign(self.eq_profile(lambda p: getattr(p, name)),
                               self.x_to_xb('CUR', 1e4, axis_value=True), name)

    @covers('eq_profile_shear_mass')
    def test_eq_profile_shear_mass(self):
        self.close(self.eq_profile(lambda p: p.magnetic_shear), self.x_to_xb('SHAT', axis_value=True), 'magnetic_shear')
        self.close(self.eq_profile(lambda p: p.mass_density), self.x_to_xb('PMASS', 1e3, axis_value=True), 'mass_density')


    # ---- summary (global_quantities/<name>/value on the TIME grid) -----------------------------

    def sm(self, name, *path):
        node = self.get('summary').global_quantities
        for part in (name,) + path:
            node = getattr(node, part)
        return arr(node.value)

    def sm_has(self, name):
        return hasattr(self.get('summary').global_quantities, name)

    def volume_integral(self, name, scale=1.0):
        """Sum of a per-cm^3 profile times the zone volume in cm^3."""
        return np.sum(self.var(name) * self.var('DVOL'), axis=1) * scale

    @covers('sm_plasma_current')
    def test_sm_plasma_current(self):
        self.close_global_sign(self.sm('ip'), self.var('PCURC'), 'ip')

    @covers('sm_direct_scalars')
    def test_sm_direct_scalars(self):
        li = 'li' if self.dd_version.startswith('3.') else 'li_3'
        self.close(self.sm(li), self.var('LI_3'), li)
        self.close(self.sm('volume'), self.var('PVOL') * 1e-6, 'volume')
        self.close_global_sign(self.sm('q_95'), self.var('Q95'), 'q_95')
        self.close_global_sign(self.sm('v_loop'), self.var('VSURC'), 'v_loop')
        self.close(self.sm('tau_energy'), self.var('TAUA1'), 'tau_energy')
        self.close(self.sm('h_98'), self.var('H98Y2'), 'h_98')
        self.close(self.sm('tau_energy_98'), self.var('TAUE98Y2'), 'tau_energy_98')
        self.close(self.sm('h_mode'), np.rint(self.var('LHMODE')), 'h_mode')

    @covers('sm_vacuum_field')
    def test_sm_vacuum_field(self):
        r0 = float(np.asarray(self.get('summary').global_quantities.r0.value))
        b0 = self.sm('b0')
        # one constant r0, the middle of the limiter, and b0 = BZXR/r0 up to the toroidal-field sign
        rlim = self.var('RLIM')[0] * 0.01
        rlim = rlim[np.isfinite(rlim) & (rlim > 0)]
        self.assertAlmostEqual(r0, 0.5 * (rlim.max() + rlim.min()), places=9)
        self.close_global_sign(b0, self.var('BZXR') * 0.01 / r0, 'b0')

    @covers('sm_betas')
    def test_sm_betas(self):
        expected = {name: np.interp(self.var('TIME'), self.var('TIME3'), values)
                    for name, values in self.expected_betas().items()}
        self.close(self.sm('beta_pol'), expected['beta_pol'], 'beta_pol')
        self.close(self.sm('beta_tor'), expected['beta_tor'], 'beta_tor')
        self.close_global_sign(self.sm('beta_tor_norm'), expected['beta_tor_norm'], 'beta_tor_norm')   # same name in DD 3.39

    @covers('sm_energies')
    def test_sm_energies(self):
        for name, cdf in (('energy_total', 'UTOTL'), ('energy_thermal', 'UTHRM'),
                          ('energy_electrons_thermal', 'UE'), ('energy_ion_total_thermal', 'UI')):
            self.close(self.sm(name), self.volume_integral(cdf), name)
        dvol_m3 = self.var('DVOL') * 1e-6
        self.close(self.sm('energy_mhd'), 1.5 * np.sum(self.var('PTOWB') * dvol_m3, axis=1), 'energy_mhd')
        self.close(self.sm('energy_diamagnetic'),
                   1.5 * np.sum((self.var('PPLAS') + self.var('UFASTPP') * 1e6) * dvol_m3, axis=1), 'energy_diamagnetic')
        self.close(self.sm('energy_fast_perpendicular'), 1.5 * self.volume_integral('UFASTPP'), 'energy_fast_perpendicular')
        self.close(self.sm('energy_fast_parallel'), 3.0 * self.volume_integral('UFASTPA'), 'energy_fast_parallel')

    @covers('sm_poloidal_field_energy')
    def test_sm_poloidal_field_energy(self):
        x, xb, bpol = self.var('X'), self.var('XB'), self.var('BPOL')
        rows = []
        for i in range(len(x)):
            xx, yy = (np.concatenate(([0.0], xb[i])), np.concatenate(([0.0], bpol[i]))) if xb[i][0] > 0 else (xb[i], bpol[i])
            rows.append(np.interp(x[i], xx, yy))
        bpol_x = np.array(rows)
        mu0 = 4e-7 * np.pi
        expected = np.sum(bpol_x ** 2 / (2 * mu0) * self.var('DVOL') * 1e-6, axis=1)
        self.close(self.sm('energy_b_field_pol'), expected, 'energy_b_field_pol')

    @covers('sm_currents')
    def test_sm_currents(self):
        darea = self.var('DAREA')
        self.close_global_sign(self.sm('current_ohm'), np.sum(self.var('CUROH') * darea, axis=1), 'current_ohm')
        self.close_global_sign(self.sm('current_non_inductive'),
                               np.sum((self.var('CUR') - self.var('CUROH')) * darea, axis=1), 'current_non_inductive')
        self.close_global_sign(self.sm('current_bootstrap'), np.sum(self.var('CURBS') * darea, axis=1), 'current_bootstrap')

    @covers('sm_powers')
    def test_sm_powers(self):
        self.close(self.sm('power_ohm'), self.var('POHT'), 'power_ohm')
        self.close(self.sm('power_radiated'), self.volume_integral('PRAD'), 'power_radiated')
        self.close(self.sm('power_radiated_inside_lcfs'), self.volume_integral('PRAD'), 'power_radiated_inside_lcfs')
        self.close(self.sm('power_bremsstrahlung'), self.volume_integral('PRAD_BR'), 'power_bremsstrahlung')
        self.close(self.sm('power_line'), self.volume_integral('PRAD_LI'), 'power_line')
        self.close(self.sm('power_synchrotron'), self.volume_integral('PRAD_CY'), 'power_synchrotron')
        hcd = self.get('summary').heating_current_drive
        self.close(arr(hcd.power_launched_nbi.value), self.var('PINJ'), 'power_launched_nbi')

    @covers('sm_neutron_rates')
    def test_sm_neutron_rates(self):
        rates = self.get('summary').fusion.neutron_rates
        self.close(arr(rates.total.value), self.var('NEUTT'), 'neutron total')
        self.close(arr(rates.thermal.value), self.var('NEUTX'), 'neutron thermal')
        self.close(arr(rates.dd.thermal.value), self.var('NEUTX_DD'), 'dd thermal')
        self.close(arr(rates.dd.beam_thermal.value), self.var('BTNTS_DD'), 'dd beam-thermal')
        self.close(arr(rates.dd.beam_beam.value), self.var('BBNTS_DD'), 'dd beam-beam')
        self.close(arr(rates.dd.total.value), self.var('NEUTX_DD') + self.var('BTNTS_DD') + self.var('BBNTS_DD'), 'dd total')


    # ---- distribution_sources (NBI birth profiles per unit with non-zero deposition) ----------

    def beam_units(self):
        """Two-digit unit numbers whose BDEPnn_TOT is not identically zero."""
        units = []
        for name in sorted(self.cdf):
            match = re.fullmatch(r'BDEP(\d\d)_TOT', name)
            if match and np.any(self.cdf[name]):
                units.append(match.group(1))
        if not units and not self.strict:
            raise unittest.SkipTest(f'no NBI unit with deposition in {self.case_path.name}')
        self.assertTrue(units, 'no NBI unit with deposition in the case')
        return units

    def source_for(self, particles):
        """The distribution source whose particle profile equals ``particles`` (nt, nx)."""
        for source in self.get('distribution_sources').source:
            values = stack(source.profiles_1d, lambda p: p.particles)
            if values is not None and values.shape == particles.shape and np.allclose(
                    values, particles, rtol=RTOL, atol=1e-12 * max(1.0, np.abs(particles).max())):
                return source
        self.fail('no distribution source carries the expected birth profile')

    @covers('ds_birth_particles')
    def test_ds_birth_particles(self):
        for unit in self.beam_units():
            self.source_for(self.var(f'BDEP{unit}_TOT') * 1e6)

    @covers('ds_birth_energy')
    def test_ds_birth_energy(self):
        e = 1.602176634e-19
        for unit in self.beam_units():
            source = self.source_for(self.var(f'BDEP{unit}_TOT') * 1e6)
            expected = sum(self.var(f'BDEP{unit}_E{k}') * 1e6 * self.var(f'EINJ{unit}_E{k}')[:, None] * e
                           for k in (1, 2, 3))
            self.close(stack(source.profiles_1d, lambda p: p.energy), expected, f'unit {unit} energy')

    @covers('ds_birth_momentum')
    def test_ds_birth_momentum(self):
        for unit in self.beam_units():
            source = self.source_for(self.var(f'BDEP{unit}_TOT') * 1e6)
            momentum = (lambda p: p.momentum_tor) if self.dd_version.startswith('3.') else (lambda p: p.momentum_phi)
            self.close_global_sign(stack(source.profiles_1d, momentum), self.var(f'TQTOT{unit}') * 1e6, f'unit {unit} torque')

    @covers('ds_global_quantities')
    def test_ds_global_quantities(self):
        e = 1.602176634e-19
        dvol_m3 = self.var('DVOL') * 1e-6
        for unit in self.beam_units():
            source = self.source_for(self.var(f'BDEP{unit}_TOT') * 1e6)
            glob = lambda name: np.array([float(getattr(g, name)) for g in source.global_quantities])
            self.close(glob('particles'), np.sum(self.var(f'BDEP{unit}_TOT') * 1e6 * dvol_m3, axis=1), 'particles')
            power = sum(self.var(f'BDEP{unit}_E{k}') * 1e6 * self.var(f'EINJ{unit}_E{k}')[:, None] * e for k in (1, 2, 3))
            self.close(glob('power'), np.sum(power * dvol_m3, axis=1), 'power')
            torque = 'torque_tor' if self.dd_version.startswith('3.') else 'torque_phi'
            self.close_global_sign(glob(torque), np.sum(self.var(f'TQTOT{unit}') * 1e6 * dvol_m3, axis=1), torque)

    @covers('ds_shinethrough')
    def test_ds_shinethrough(self):
        e = 1.602176634e-19
        for unit in self.beam_units():
            source = self.source_for(self.var(f'BDEP{unit}_TOT') * 1e6)
            power = sum(self.var(f'BPSH{unit}_E{k}') for k in (1, 2, 3) if f'BPSH{unit}_E{k}' in self.cdf)
            particles = sum(np.divide(self.var(f'BPSH{unit}_E{k}'), self.var(f'EINJ{unit}_E{k}') * e,
                                      out=np.zeros(len(power)), where=self.var(f'EINJ{unit}_E{k}') > 0)
                            for k in (1, 2, 3) if f'BPSH{unit}_E{k}' in self.cdf)
            self.close([float(g.shinethrough.power) for g in source.global_quantities], power, 'shinethrough power')
            self.close([float(g.shinethrough.particles) for g in source.global_quantities], particles, 'shinethrough particles')


    @covers('eq_boundary_psi')
    def test_eq_boundary_psi(self):
        edge = self.eq_profile(lambda p: p.psi)[:, -1]
        boundary = self.eq_scalar(lambda ts: ts.boundary.psi)
        expected = self.var('PLFLX2PI')[:, -1] + self.var('PSI0_TR') * 2 * np.pi
        self.close_global_sign(boundary, expected, 'boundary/psi')
        # one reference: the boundary flux is the edge of profiles_1d/psi, in every node that names it
        self.close(boundary, edge, 'boundary/psi against profiles_1d/psi[-1]')
        self.close(self.eq_scalar(lambda ts: ts.global_quantities.psi_boundary), edge, 'global psi_boundary')

    @covers('eq_boundary_separatrix_psi')
    def test_eq_boundary_separatrix_psi(self):
        if not self.dd_version.startswith('3.'):
            self.assertFalse(hasattr(self.get('equilibrium').time_slice[0], 'boundary_separatrix'))
            return
        expected = (self.var('PBPMHD') + self.var('PSI0_TR')) * 2 * np.pi
        self.close_global_sign(self.eq_scalar(lambda ts: ts.boundary_separatrix.psi), expected, 'boundary_separatrix/psi')

    @covers('eq_profile_2d_psi')
    def test_eq_profile_2d_psi(self):
        eq = self.get('equilibrium')
        axis_name = 'psi_axis' if self.dd_version.startswith('3.') else 'psi_magnetic_axis'
        sign = None
        for it, ts in enumerate(eq.time_slice):
            nr, nz = len(self.var('RGRID')[it]), len(self.var('ZGRID')[it])
            actual = arr(ts.profiles_2d[0].psi)
            expected = (2 * np.pi * (self.var('PSIRZ')[it] + self.var('PSI0_TR')[it])).reshape(nz, nr).T
            self.assertEqual(actual.shape, expected.shape)
            if sign is None:
                sign = 1.0 if np.allclose(actual, expected, rtol=RTOL, atol=1e-9) else -1.0
            np.testing.assert_allclose(actual, sign * expected, rtol=RTOL, atol=1e-9, err_msg=f'psi_2d slice {it}')
            span = float(actual.max() - actual.min())
            # the same reference as the axis value: the field is at the axis value at the magnetic axis
            axis = float(getattr(ts.global_quantities, axis_name))
            self.assertLess(np.min(np.abs(actual - axis)), 1e-3 * span, f'2D psi never reaches the axis flux at slice {it}')
            # one reference for the outline, the 2-D flux and boundary/psi: the outline is a contour of the
            # 2-D flux at boundary/psi.  ITER holds this to 1e-6 of the span.  DIII-D 147634T16's PSIRZ is flat
            # outside the plasma at a value 0.1-0.2% below the transport-edge flux (a property of that CDF),
            # so the tolerance there is 1% of the span: still far below the 40% offset of mixed references.
            outline = np.column_stack((arr(ts.boundary.outline.r), arr(ts.boundary.outline.z)))
            self.assertGreater(len(outline), 3, f'no boundary outline at slice {it}')
            on_outline = RegularGridInterpolator((arr(ts.profiles_2d[0].grid.dim1), arr(ts.profiles_2d[0].grid.dim2)),
                                                 actual, bounds_error=False)(outline)
            tolerance = (1e-6 if self.strict else 1e-2) * span
            self.assertFalse(np.isnan(on_outline).any(), f'outline leaves the 2-D grid at slice {it}')
            self.assertLessEqual(float(np.ptp(on_outline)), tolerance, f'outline is not a flux contour at slice {it}')
            self.assertLessEqual(abs(float(np.mean(on_outline)) - float(ts.boundary.psi)), tolerance,
                                 f'flux on the outline differs from boundary/psi at slice {it}')


    # ---- waves (coherent_wave per EC gyrotron and per IC antenna, ITER has 5 EC and 2 IC) -------

    def wave(self, antenna_name):
        for wave in self.get('waves').coherent_wave:
            if str(wave.identifier.antenna_name) == antenna_name:
                return wave
        self.fail(f'no coherent_wave for {antenna_name}')

    def wave_numbers(self, prefix):
        return [name[len(prefix):] for name in self.names(prefix + r'\d')]

    @staticmethod
    def wave_profile(wave, getter):
        return stack(wave.profiles_1d, getter)

    @staticmethod
    def wave_global(wave, getter):
        return np.array([float(getter(g)) for g in wave.global_quantities])

    @covers('wv_ec_power')
    def test_wv_ec_power(self):
        dvol_m3 = self.var('DVOL') * 1e-6
        for n in self.wave_numbers('PEECH'):
            wave = self.wave(f'beam{n}')
            density = self.var(f'PEECH{n}') * 1e6
            inside = np.cumsum(density * dvol_m3, axis=1)
            self.close(self.wave_profile(wave, lambda p: p.power_density), density, f'PEECH{n} power_density')
            self.close(self.wave_profile(wave, lambda p: p.power_inside), inside, f'PEECH{n} power_inside')
            self.close(self.wave_profile(wave, lambda p: p.electrons.power_density_thermal), density, f'PEECH{n} electrons')
            self.close(self.wave_global(wave, lambda g: g.power), inside[:, -1], f'PEECH{n} global power')
            self.close(self.wave_global(wave, lambda g: g.electrons.power_thermal), inside[:, -1], f'PEECH{n} electron power')

    @covers('wv_ec_current')
    def test_wv_ec_current(self):
        inside_name = 'current_tor_inside' if self.dd_version.startswith('3.') else 'current_phi_inside'
        global_name = 'current_tor' if self.dd_version.startswith('3.') else 'current_phi'
        for n in self.wave_numbers('ECCUR'):
            wave = self.wave(f'beam{n}')
            inside = np.cumsum(self.var(f'ECCUR{n}') * self.var('DAREA'), axis=1)
            self.close_global_sign(self.wave_profile(wave, lambda p: getattr(p, inside_name)), inside, f'ECCUR{n} inside')
            self.close_global_sign(self.wave_global(wave, lambda g: getattr(g, global_name)), inside[:, -1], f'ECCUR{n} global')

    @covers('wv_ic_power')
    def test_wv_ic_power(self):
        for n in self.wave_numbers('RFPWRAB'):
            wave = self.wave(f'antenna{n}')
            self.close(self.wave_global(wave, lambda g: g.power), self.var(f'RFPWRAB{n}'), f'RFPWRAB{n}')
            if f'PIC_F{n}' in self.cdf:
                density = self.var(f'PIC_F{n}') * 1e6
                self.close(self.wave_profile(wave, lambda p: p.power_density), density, f'PIC_F{n}')
                self.close(self.wave_profile(wave, lambda p: p.power_inside),
                           np.cumsum(density * self.var('DVOL') * 1e-6, axis=1), f'PIC_F{n} power_inside')

    @covers('wv_ic_frequency')
    def test_wv_ic_frequency(self):
        for n in self.wave_numbers('RFPWRAB'):
            self.close(self.wave_global(self.wave(f'antenna{n}'), lambda g: g.frequency), self.var(f'FREQA{n}'), f'FREQA{n}')


    # ---- core_profiles: species-resolved, global and remaining profile quantities ------------

    def ion_states(self, ion_index):
        rows = [s for s in self.get('core_profiles').profiles_1d if ion_index < len(s.ion)]
        return rows

    def assert_some_ion_global_sign(self, getter, expected, label):
        for values in self.ion_arrays(getter):
            if values.shape == expected.shape:
                atol = 1e-9 * max(1.0, np.abs(expected).max())
                if np.allclose(values, expected, rtol=RTOL, atol=atol) or np.allclose(values, -expected, rtol=RTOL, atol=atol):
                    return
        self.fail(f'{label}: no ion carries the expected profile (either sign)')

    def cp_series(self, name):
        """A core_profiles global quantity as a time series (the IDS stores one value per TIME3 slice)."""
        return arr(getattr(self.get('core_profiles').global_quantities, name))

    @covers('cp_ion_density_total')
    def test_cp_ion_density_total(self):
        checked = 0
        for sym, short in (('D', 'D'), ('T', 'T'), ('HE3', '3'), ('HE4', '4'), ('H', 'H')):
            if f'N{sym}' not in self.cdf:
                continue
            parts = [f'N{sym}'] + [n for n in (f'BDENS_{short}', f'FDENS_{short}', f'NMINI_{short}') if n in self.cdf]
            expected = sum(self.var(n) for n in parts) * 1e6
            self.assert_some_ion(lambda ion: ion.density, expected, f'total {sym} = ' + '+'.join(parts))
            checked += 1
        self.assertTrue(checked)

    @covers('cp_fast_ion_density')
    def test_cp_fast_ion_density(self):
        for name in ('BDENS_D', 'FDENS_T', 'FDENS_3', 'FDENS_4'):
            if name in self.cdf:
                self.assert_some_ion(lambda ion: ion.density_fast, self.var(name) * 1e6, name)

    @covers('cp_fast_ion_rotation')
    def test_cp_fast_ion_rotation(self):
        for name in ('OMEGB_D', 'OMEGF_T', 'OMEGF_3', 'OMEGF_4'):
            if name in self.cdf:
                self.assert_some_ion_global_sign(lambda ion: ion.rotation_frequency_tor, self.var(name), name)

    @covers('cp_impurity_density')
    def test_cp_impurity_density(self):
        for name in self.names(r'NIMPS_[A-Z0-9]+'):
            # the aggregate impurity (TOK) is read from NIMP; NIMPS_TOK is an all-zero placeholder in some CDFs
            expected = self.var('NIMP' if name == 'NIMPS_TOK' else name) * 1e6
            self.assert_some_ion(lambda ion: ion.density, expected, name)
            self.assert_some_ion(lambda ion: ion.density_thermal, expected, name)

    @covers('cp_impurity_temperature_rotation')
    def test_cp_impurity_temperature_rotation(self):
        self.assert_some_ion(lambda ion: ion.temperature, self.var('TX'), 'TX')
        self.assert_some_ion_global_sign(lambda ion: ion.rotation_frequency_tor, self.var('OMGX'), 'OMGX')

    @covers('cp_impurity_charge_states')
    def test_cp_impurity_charge_states(self):
        profiles = self.get('core_profiles').profiles_1d
        n_ions = max(len(s.ion) for s in profiles)
        for name in self.names(r'NIMP_[A-Z0-9]+_\d+'):
            expected = self.var(name) * 1e6
            found = False
            for ion_index in range(n_ions):
                n_states = max((len(s.ion[ion_index].state) for s in profiles if ion_index < len(s.ion)), default=0)
                for state_index in range(n_states):
                    values = stack(profiles, lambda s: s.ion[ion_index].state[state_index].density)
                    if values is not None and values.shape == expected.shape and np.allclose(
                            values, expected, rtol=RTOL, atol=1e-9 * max(1.0, np.abs(expected).max())):
                        found = True
            self.assertTrue(found, f'{name}: no ion charge state carries the expected density')

    @covers('cp_neutral_profiles')
    def test_cp_neutral_profiles(self):
        profiles = self.get('core_profiles').profiles_1d

        def some_neutral(getter, expected, label):
            for j in range(max(len(s.neutral) for s in profiles)):
                values = stack([s for s in profiles if j < len(s.neutral)], lambda s: getter(s.neutral[j]))
                if values is not None and values.shape == expected.shape and np.allclose(
                        values, expected, rtol=RTOL, atol=1e-9 * max(1.0, np.abs(expected).max())):
                    return
            self.fail(f'{label}: no neutral population carries the expected profile')

        for name in self.names(r'DN0[VW][A-Z0-9]+'):
            some_neutral(lambda n: n.density, self.var(name) * 1e6, name)
        for name in self.names(r'T0[VW][A-Z0-9]+'):
            some_neutral(lambda n: n.temperature, self.var(name), name)

    def names_if_present(self, pattern):
        """Like ``names`` but skips, even on the strict case, when the case simply has no such variable."""
        found = sorted(n for n in self.cdf if re.fullmatch(pattern, n))
        if not found:
            raise unittest.SkipTest(f'{self.case_path.name} has no variable matching {pattern}')
        return found

    def some_neutral(self, getter, expected, label):
        profiles = self.get('core_profiles').profiles_1d
        for j in range(max(len(s.neutral) for s in profiles)):
            values = stack([s for s in profiles if j < len(s.neutral)], lambda s: getter(s.neutral[j]))
            if values is not None and values.shape == expected.shape and np.allclose(
                    values, expected, rtol=RTOL, atol=1e-9 * max(1.0, np.abs(expected).max())):
                return
        self.fail(f'{label}: no neutral population carries the expected profile')

    @covers('cp_neutral_recombination')
    def test_cp_neutral_recombination(self):
        for name in self.names_if_present(r'N0V0_[A-Z0-9]+'):
            self.some_neutral(lambda n: n.density, self.var(name) * 1e6, name)
            temperature = 'T0V0_' + name[5:]
            if temperature in self.cdf:
                self.some_neutral(lambda n: n.temperature, self.var(temperature), temperature)

    @covers('cp_neutral_rotation_xml')
    def test_cp_neutral_rotation_xml(self):
        names = self.names_if_present(r'OM0[VW](?!0_)[A-Z0-9]+')
        root = ET.fromstring(str(self.get('core_profiles').code.parameters))
        records = {p.get('source_variable'): p for p in root.iter('profile')}
        for name in names:
            self.assertIn(name, records, f'{name} is missing from the neutral_toroidal_rotation XML')
            node = records[name]
            shape = (int(node.get('rows')), int(node.get('cols')))
            values = np.array(node.text.split(), dtype=float).reshape(shape)
            self.close(values, float(node.get('source_sign_factor')) * self.var(name), f'{name} (XML)')
            self.assertEqual(abs(float(node.get('source_sign_factor'))), 1.0)

    def single_impurity_ions(self):
        """Impurity ions that carry a per-time element (the NIMP/TOK path), as (time, a, z_n, z_ion) arrays."""
        profiles = self.get('core_profiles').profiles_1d
        times = arr(self.get('core_profiles').time)
        out = []
        for j in range(max(len(s.ion) for s in profiles)):
            if all(j < len(s.ion) and len(s.ion[j].element) for s in profiles):
                out.append((times, [float(s.ion[j].element[0].a) for s in profiles],
                            [float(s.ion[j].element[0].z_n) for s in profiles], [float(s.ion[j].z_ion) for s in profiles]))
        return out

    @covers('cp_single_impurity_mass')
    def test_cp_single_impurity_mass(self):
        self.names_if_present(r'NIMPS_TOK')            # AIMP only identifies the lumped single impurity
        expected = np.interp(arr(self.get('core_profiles').time), self.var('TIME'), self.var('AIMP'))
        self.assertTrue(any(np.allclose(a, expected, rtol=RTOL) for _t, a, _z, _zi in self.single_impurity_ions()),
                        'no ion carries AIMP as element.a')

    @covers('cp_single_impurity_charge')
    def test_cp_single_impurity_charge(self):
        self.names_if_present(r'NIMPS_TOK')
        expected = np.interp(arr(self.get('core_profiles').time), self.var('TIME'), self.var('XZIMP'))
        self.assertTrue(any(np.allclose(z, expected, rtol=RTOL) and np.allclose(zi, expected, rtol=RTOL)
                            for _t, _a, z, zi in self.single_impurity_ions()),
                        'no ion carries XZIMP as element.z_n and z_ion')

    @covers('cp_momentum_phi')
    def test_cp_momentum_phi(self):
        name = 'momentum_tor' if self.dd_version.startswith('3.') else 'momentum_phi'
        self.close_global_sign(self.profiles(lambda s: getattr(s, name)), self.var('PPHI') * 1e6, name)

    @covers('cp_global_scalars')
    def test_cp_global_scalars(self):
        self.close_global_sign(self.cp_series('ip'), self.var('PCURC'), 'ip')
        self.close_global_sign(self.cp_series('v_loop'), self.var('VSURC'), 'v_loop')
        self.close(self.cp_series('li_3'), self.var('LI_3'), 'li_3')

    @covers('cp_betas')
    def test_cp_betas(self):
        expected = self.expected_betas()
        self.close(self.cp_series('beta_tor'), expected['beta_tor'], 'beta_tor')
        self.close(self.cp_series('beta_pol'), expected['beta_pol'], 'beta_pol')
        self.close_global_sign(self.cp_series('beta_tor_norm'), expected['beta_tor_norm'], 'beta_tor_norm')

    @covers('cp_global_currents')
    def test_cp_global_currents(self):
        darea = self.var('DAREA')
        self.close_global_sign(self.cp_series('current_non_inductive'),
                               np.sum((self.var('CUR') - self.var('CUROH')) * darea, axis=1), 'current_non_inductive')
        self.close_global_sign(self.cp_series('current_bootstrap'), np.sum(self.var('CURBS') * darea, axis=1), 'current_bootstrap')

    @covers('cp_energy_diamagnetic')
    def test_cp_energy_diamagnetic(self):
        self.close(self.cp_series('energy_diamagnetic'),
                   1.5 * np.sum(self.var('PTOWB') * self.var('DVOL') * 1e-6, axis=1), 'energy_diamagnetic')

    @covers('cp_volume_averages')
    def test_cp_volume_averages(self):
        gvol = self.var('DVOL') * 1e-6
        total = np.sum(gvol, axis=1)
        average = lambda values: np.sum(values * gvol, axis=1) / total
        gq = self.get('core_profiles').global_quantities
        self.close(self.cp_series('t_e_volume_average'), average(self.var('TE')), 't_e_volume_average')
        self.close(self.cp_series('n_e_volume_average'), average(self.var('NE') * 1e6), 'n_e_volume_average')
        self.close(arr(gq.ion[0].t_i_volume_average), average(self.var('TI')), 't_i_volume_average')
        self.close(arr(gq.ion[0].n_i_volume_average), average(self.var('NI') * 1e6), 'n_i_volume_average')

    @covers('cp_peaking')
    def test_cp_peaking(self):
        gvol = self.var('DVOL') * 1e-6
        total = np.sum(gvol, axis=1)
        self.close(self.cp_series('t_e_peaking'), self.var('TE0') / (np.sum(self.var('TE') * gvol, axis=1) / total), 't_e_peaking')
        self.close(self.cp_series('t_i_average_peaking'), self.var('TI0') / (np.sum(self.var('TI') * gvol, axis=1) / total),
                   't_i_average_peaking')

    @covers('cp_j_bootstrap_and_non_inductive')
    def test_cp_j_bootstrap_and_non_inductive(self):
        b0 = np.abs(arr(self.get('core_profiles').vacuum_toroidal_field.b0))[:, None]
        self.close_global_sign(self.profiles(lambda s: s.j_bootstrap), self.var('CURBS') * 1e4 * self.var('CONPLJB') / b0, 'j_bootstrap')
        self.close_global_sign(self.profiles(lambda s: s.j_non_inductive),
                               (self.var('PLJB') - self.var('PLJBOH')) * 1e4 / b0, 'j_non_inductive')

    @covers('cp_profile_j_tor')
    def test_cp_profile_j_tor(self):
        name = 'j_tor' if self.dd_version.startswith('3.') else 'j_phi'
        self.close_global_sign(self.profiles(lambda s: getattr(s, name)), self.var('CUR') * 1e4, name)


class RegistryCoverageTests(unittest.TestCase):
    def test_registry_and_tests_list_the_same_ids(self):
        import mapping_registry as reg
        claimed = {e['id'] for e in reg.load_registry()
                   if any(Path(t).name == Path(__file__).name for t in e.get('tests') or ())}
        # COVERED is filled when RegistryValueTests is defined (import time)
        self.assertEqual(claimed, COVERED, 'registry `tests` and the covers() decorators disagree')


for _label, (_path, _strict) in CASES.items():
    for _version in DD_VERSIONS:
        _name = f'RegistryValues_{_label}_DD' + _version.replace('.', '_')
        globals()[_name] = type(_name, (RegistryValueTests, unittest.TestCase),
                                {'dd_version': _version, 'case_path': _path, 'strict': _strict})

if __name__ == '__main__':
    unittest.main(verbosity=1)
