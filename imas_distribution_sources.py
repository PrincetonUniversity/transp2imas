#!/usr/bin/python
"""Populate NBI distribution_sources from TRANSP/NUBEAM birth and deposition moments."""

import re

import numpy as np
from scipy import constants

from nc2imas import nc2imas
from utilities import append_successful_entries


class imas_distribution_sources:
    """Export NUBEAM fast-ion birth/deposition sources to distribution_sources.

    TRANSP's BDEPnn_TOT and BDEPnn_E1/E2/E3 profiles are particle birth
    source rates for individual NBI units on the normal TRANSP flux grid.
    They are therefore a substantially better match to distribution_sources
    than to the full distributions IDS.  This writer deliberately does not
    synthesize a phase-space distribution when FBM/EBA2* arrays are absent.

    One ``source`` is written per NBI unit.  Each source carries the summed
    accelerated-energy process (nbi_energy/index=0), radial particle source,
    kinetic birth-energy source, volume-integrated totals, and shinethrough
    losses when BPSHnn_E1/E2/E3 are present.
    """

    NBI_IDENTIFIER_INDEX = 1

    def __init__(self, ncimas: nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger
        self.time = np.asarray(self.im.get_ncvar('TIME3'), dtype=np.float64)
        self.dvol = np.asarray(self.im.get_ncvar('DVOL', self.im.vconv), dtype=np.float64)
        self.nt = len(self.time)
        self.create_ids()

    @staticmethod
    def _beam_number(name):
        match = re.match(r'^BDEP(\d{2})_TOT$', str(name))
        return int(match.group(1)) if match else None

    def _beam_numbers(self):
        names = self.im.get_var_names_re_('^BDEP\\d{2}_TOT$')
        nums = sorted({self._beam_number(name) for name in names})
        return [value for value in nums if value is not None and value > 0]

    def _get(self, name, mult=None):
        values = self.im.get_ncvar(name)
        if not self.im.is_valid(values):
            return np.array([])
        values = np.asarray(values, dtype=np.float64)
        if mult is not None:
            values = values * float(mult)
        return values

    def _nml_float(self, name, ib):
        """Return a per-NBI-unit namelist value as float.

        TRANSP beam/unit numbers are 1-based, while nc2imas.get_namelist_var()
        uses a zero-based ``ind`` argument and transparently handles either
        NAME = a,b,... or NAME(1)=..., NAME(2)=... namelist spellings.
        Missing/malformed values are optional metadata and must never abort the
        distribution_sources writer.
        """
        target = f'source[{int(ib)}]/species'
        try:
            value = self.im.get_namelist_var(name, float, ind=int(ib) - 1)
        except (KeyError, ValueError, TypeError, IndexError) as exc:
            self.im.mapping.skipped(f'namelist value unavailable: {type(exc).__name__}: {exc}',
                                    source=name, target=target)
            return None
        if value is None:
            self.im.mapping.skipped('namelist value not set', source=name, target=target)
            return None
        try:
            # get_namelist_var normally returns a scalar for ind>=0, but be
            # defensive in case an older parser returns a length-one array.
            arr = np.asarray(value)
            if arr.size != 1:
                self.im.mapping.skipped(f'namelist value has {arr.size} elements, expected 1',
                                        source=name, target=target)
                return None
            value = float(arr.reshape(-1)[0])
            if not np.isfinite(value):
                self.im.mapping.skipped('namelist value is not finite', source=name, target=target)
                return None
            return value
        except (ValueError, TypeError) as exc:
            self.im.mapping.skipped(f'namelist value is not numeric: {exc}', source=name, target=target)
            return None

    @staticmethod
    def _species_name(a, z_n):
        if a is None or z_n is None:
            return 'NBI ion'
        zi = int(round(z_n))
        ai = int(round(a))
        if zi == 1:
            return {1: 'H', 2: 'D', 3: 'T'}.get(ai, f'H{ai}')
        if zi == 2:
            return {3: 'He3', 4: 'He4'}.get(ai, f'He{ai}')
        return f'Z{zi}A{ai}'

    def _fill_species(self, source, ib):
        source.species.type.name = 'ion'
        source.species.type.index = 2
        source.species.type.description = 'Ion species in a single/average state'

        a = self._nml_float('ABEAMA', ib)
        z_n = self._nml_float('XZBEAMA', ib)
        ion = source.species.ion
        if a is not None or z_n is not None:
            ion.element.resize(1)
            element = ion.element[0]
            if a is not None:
                element.a = float(a)
            if z_n is not None:
                self.im.set_nuclear_charge(element, z_n)
            element.atoms_n = 1
        if z_n is not None:
            # The kinetic source represented by BDEP is the ion born after
            # ionization of the injected neutral.  For hydrogenic NBI the
            # born fast ion is singly charged.  For other nuclei, do not guess
            # a charge state from nuclear charge alone.
            if int(round(z_n)) == 1:
                ion.z_ion = 1.0
        self.im.set_species_name(ion, self._species_name(a, z_n))

    @staticmethod
    def _set_process(process, ib):
        process.type.name = 'NBI'
        process.type.index = imas_distribution_sources.NBI_IDENTIFIER_INDEX
        process.type.description = 'Neutral Beam Injection'
        process.nbi_energy.name = 'total'
        process.nbi_energy.index = 0
        process.nbi_energy.description = 'Sum over full-, half-, and third-energy NBI components'
        process.nbi_unit = int(ib)
        process.nbi_beamlets_group = 0

    def _energy_source(self, ib):
        """Return sum_k BDEP_Ek * Einj_Ek as W/m3.

        BDEP is in particles cm^-3 s^-1 and EINJ in eV.  Preserve the CDF
        profile as supplied (including small Monte-Carlo/smoothing negatives)
        instead of clipping, so volume integrals remain reproducible.
        """
        stn = f'{ib:02d}'
        terms = []
        for k in (1, 2, 3):
            bdep = self._get(f'BDEP{stn}_E{k}', self.im.vol_inv)
            energy = self._get(f'EINJ{stn}_E{k}')
            if bdep.ndim == 2 and energy.ndim == 1 and bdep.shape[0] == len(energy):
                terms.append(bdep * energy[:, None] * constants.e)
        if not terms:
            return np.array([])
        result = np.zeros_like(terms[0], dtype=np.float64)
        for term in terms:
            result += term
        return result

    def _shinethrough(self, ib):
        stn = f'{ib:02d}'
        power = np.zeros(self.nt, dtype=np.float64)
        particles = np.zeros(self.nt, dtype=np.float64)
        have_power = False
        have_particles = False
        for k in (1, 2, 3):
            p = self._get(f'BPSH{stn}_E{k}')
            e = self._get(f'EINJ{stn}_E{k}')
            if p.ndim == 1 and len(p) == self.nt:
                power += p
                have_power = True
                if e.ndim == 1 and len(e) == self.nt:
                    denom = e * constants.e
                    term = np.divide(
                        p, denom,
                        out=np.zeros_like(p, dtype=np.float64),
                        where=np.isfinite(denom) & (denom > 0.0),
                    )
                    particles += term
                    have_particles = True
        return (power if have_power else np.array([]),
                particles if have_particles else np.array([]))

    def _fill_source(self, source, ib):
        stn = f'{ib:02d}'
        particle_source = self._get(f'BDEP{stn}_TOT', self.im.vol_inv)
        if particle_source.ndim != 2 or particle_source.shape[0] != self.nt:
            return False
        if not np.any(np.isfinite(particle_source) & (particle_source != 0.0)):
            # The CDF commonly allocates BDEP arrays for configured but unused
            # beams.  Keep those hardware units in nbi, but do not create an
            # all-zero kinetic source here.
            return False

        source.process.resize(1)
        self._set_process(source.process[0], ib)
        self._fill_species(source, ib)
        # BDEP documentation in this CDF does not state unambiguously whether
        # the spatial coordinate is the actual birth point or gyro-centre.
        # Leave gyro_type unset instead of making that distinction up.

        source.profiles_1d.resize(self.nt)
        source.global_quantities.resize(self.nt)
        self.im.fill_profiles1d_grid(source)

        energy_source = self._energy_source(ib)
        have_energy = (
            energy_source.ndim == 2 and energy_source.shape == particle_source.shape
        )
        # TQTOTnn is the total toroidal beam torque density for the NBI unit
        # in N*m/cm^3. Convert to N/m^2 and orient with the IMAS phi sign.
        torque_source = self._get(f'TQTOT{stn}', self.im.vol_inv)
        if torque_source.ndim == 2 and torque_source.shape == particle_source.shape:
            torque_source = self.im.toroidal_rotation_sign * torque_source
        else:
            torque_source = np.array([])
        particle_global = np.sum(particle_source * self.dvol, axis=1)
        power_global = np.sum(energy_source * self.dvol, axis=1) if have_energy else np.array([])
        # dvol is in m^3 here; torque_source is N/m^2 = N*m/m^3.
        torque_global = np.sum(torque_source * self.dvol, axis=1) if torque_source.size else np.array([])
        shine_power, shine_particles = self._shinethrough(ib)

        nx = particle_source.shape[1]
        for it in range(self.nt):
            profile = source.profiles_1d[it]
            profile.time = float(self.time[it])
            self.im.fill_data(profile.particles, particle_source, ind=it, chk=nx)
            if have_energy:
                self.im.fill_data(profile.energy, energy_source, ind=it, chk=nx)
            if torque_source.size:
                mom = profile.momentum_phi if self.im.dd_at_least(3, 42, 0) else profile.momentum_tor
                self.im.fill_data(mom, torque_source, ind=it, chk=nx)

            glob = source.global_quantities[it]
            glob.time = float(self.time[it])
            glob.particles = float(particle_global[it])
            if power_global.size:
                glob.power = float(power_global[it])
            if torque_global.size:
                if self.im.dd_at_least(3, 42, 0):
                    glob.torque_phi = float(torque_global[it])
                else:
                    glob.torque_tor = float(torque_global[it])
            if shine_power.size:
                glob.shinethrough.power = float(shine_power[it])
            if shine_particles.size:
                glob.shinethrough.particles = float(shine_particles[it])
        return True

    def create_ids(self):
        self.logger.info('\ndistribution_sources')
        beam_numbers = self._beam_numbers()
        if not beam_numbers:
            self.logger.info('distribution_sources: no BDEPnn_TOT NBI birth profiles; skipping')
            return
        if self.dvol.ndim != 2 or self.dvol.shape[0] != self.nt:
            self.logger.warning('distribution_sources: DVOL/TIME3 mismatch; skipping')
            return

        try:
            ids = self.im.new_ids('distribution_sources')
        except ValueError as exc:          # imas.exception.IDSNameError: not in this DD
            self.im.mapping.unavailable_in_dd('distribution_sources', str(exc))
            self.logger.warning(f'distribution_sources is unavailable in target DD: {exc}')
            return

        ids = self.im.fill_properties(ids)
        self.im.fill_code_properties(ids)
        if hasattr(ids.ids_properties, 'homogeneous_time'):
            ids.ids_properties.homogeneous_time = 1
        self.im.fill_vac_tor_field(ids)
        if hasattr(ids, 'time'):
            self.im.fill_data(ids.time, self.time)

        used = append_successful_entries(
            lambda: self.im.new_ids('distribution_sources'), ids, 'source',
            [lambda source, ib=ib: self._fill_source(source, ib) for ib in beam_numbers])

        if used == 0:
            self.logger.info('distribution_sources: no supported NBI sources after validation; skipping')
            return
        self.im.putdb_imas(ids)
        self.logger.info(f'distribution_sources: wrote {used} NBI unit source(s)')


if __name__ == '__main__':
    print('Please use: python convert.py')
