#!/usr/bin/python
"""All-example checks for neoclassical particle/heat transport and neutrals.

Reads real CDF variable dimensions/units, independently reconstructs
neoclassical boundary fluxes, verifies neutral densities/temperatures/rotation,
and tests missing/invalid data. No IMAS installation required.
"""

import argparse
import ast
import json
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest import mock
import xml.etree.ElementTree as ET

import numpy as np
from scipy.io import netcdf_file

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from neoclassical_neutral_mapping import (
    attach_rotation_xml, discover_neoclassical_channels,
    discover_neutral_channels, divergence_to_surface_flux,
    validate_profile,
)


def require(check, reason):
    if not check:
        raise AssertionError(reason)


def verify_case(nc, case):
    volume = np.asarray(nc.variables['DVOL'][:], dtype=np.float64)
    surface = np.asarray(nc.variables['SURF'][:], dtype=np.float64)
    x = np.asarray(nc.variables['X'][:], dtype=np.float64)
    time = np.asarray(nc.variables['TIME3'][:], dtype=np.float64)
    require(volume.shape == surface.shape == x.shape,
            f'{case}: incompatible TRANSP X/XB finite-volume grids')
    require(volume.shape[0] == time.size, f'{case}: wrong profile time axis')

    neo = discover_neoclassical_channels(nc.variables)
    nu = discover_neutral_channels(nc.variables)
    result = {'case': str(case), 'neoclassical': {}, 'neutral_populations': [],
              'native_time3': len(time), 'radial_points': int(x.shape[1])}
    for kind, channels in neo.items():
        checked = []
        for species, name in channels.items():
            div = validate_profile(nc, name, kind, required_shape=x.shape)
            flux = divergence_to_surface_flux(div, volume, surface, kind)
            # In SI, flux [1/m2/s or W/m2] * area [cm2] / 1e4 equals
            # enclosed CDF source [N/s or W]. Check each volume element.
            enclosed = flux * surface / 1.0e4
            reconstructed = np.diff(
                np.concatenate([np.zeros((div.shape[0],1)),enclosed],axis=1),axis=1)
            expected = div*volume
            mask = np.abs(surface)>0
            # Zero area means the source has no defined flux in an IDS;
            # compare only boundaries having a meaningful area.
            good = mask & np.isfinite(reconstructed)
            if np.any(good):
                den = np.maximum(1.0, np.abs(expected[good]))
                err = float(np.max(np.abs(reconstructed[good]-expected[good])/den))
                require(err < 2e-6, f'{case}: {name} NC flux integration residual {err:.3g}')
            require(np.isfinite(flux).all(), f'{case}: {name} non-finite neoclassical flux')
            checked.append(name)
        result['neoclassical'][kind] = checked
    require(any(result['neoclassical'].values()), f'{case}: no NC source groups')

    for spec in nu:
        entry = {'origin': spec['origin'], 'species': spec['species'],
                 'sources': {}, 'rotation_preserved': False}
        density = validate_profile(nc, spec['density'], 'density', required_shape=x.shape)
        require(np.all(density>=0), f'{case}: negative neutral density in {spec["density"]}')
        entry['sources']['density'] = spec['density']
        if spec['temperature']:
            temp = validate_profile(nc, spec['temperature'], 'temperature', x.shape)
            require(np.all(temp>=0), f'{case}: negative neutral temperature in {spec["temperature"]}')
            entry['sources']['temperature'] = spec['temperature']
        if spec['rotation']:
            validate_profile(nc, spec['rotation'], 'rotation', x.shape)
            entry['sources']['rotation'] = spec['rotation']
            entry['rotation_preserved'] = True
        result['neutral_populations'].append(entry)
    require(nu, f'{case}: no neutral density channels')
    return result


class TestMappingContracts(unittest.TestCase):
    def test_aggregate_channels_are_distinct(self):
        chan = discover_neoclassical_channels([
            'GFLNC_I','GFLNC_E','GFLNC_D','GFLNC_X','QFLNC_I',
            'QFLNC_D','QFLNCC_E'])
        self.assertIn('I',chan['particles'])
        self.assertIn('X',chan['particles'])
        self.assertEqual(chan['heat_nc_classical']['E'],'QFLNCC_E')
        source = (ROOT/'imas_core_transport.py').read_text()
        self.assertIn("- {'E', 'I', 'X', 'TOK'}", source)
        self.assertIn('p.total_ion_energy.flux', source)
        self.assertIn("'neoclassical_classical_heat'", source)

    def test_neutral_populations_are_not_merged(self):
        names = ['DN0VD','T0VD','OM0VD','DN0WD','T0WD','OM0WD',
                 'N0V0_D','T0V0_D','OM0V0_D','DN0VHE4','T0VHE4']
        profiles=discover_neutral_channels(names)
        self.assertEqual([s['origin'] for s in profiles],
                         ['wall','volume','volume','recombination'])
        self.assertEqual(profiles[2]['species'],'HE4')
        self.assertEqual(profiles[-1]['density'],'N0V0_D')

    def test_flux_integrates_divergence_with_si_conversion(self):
        div=np.array([[2.,-1.,3.],[0.,4.,2.]])
        dv=np.array([[10.,20.,30.],[10.,20.,30.]])
        surf=np.array([[100.,120.,150.],[100.,120.,150.]])
        got=divergence_to_surface_flux(div,dv,surf,'heat')
        self.assertTrue(np.allclose(got[:,1],
                        np.cumsum(div*dv,axis=1)[:,1]*1e4/surf[:,1]))
        with self.assertRaisesRegex(ValueError,'incompatible'):
            divergence_to_surface_flux(div,dv[:,:2],surf,'particles')

    def test_neutral_omega_is_reversible_in_code_xml(self):
        neutral=types.SimpleNamespace(
            label='D_VOLUME',rotation_source='OM0VD',rotation_sign=-1.0,
            omega=np.array([[2.,-4.],[1.,8.]],dtype=np.float64))
        xml=attach_rotation_xml('<transp2imas_parameters/>',[neutral],np.array([0.,1.]))
        node=ET.fromstring(xml).find('neutral_toroidal_rotation/profile')
        self.assertEqual(node.attrib['source_variable'],'OM0VD')
        self.assertEqual(node.attrib['source_sign_factor'],'-1.0')
        self.assertEqual(np.fromstring(node.text,sep=' ').reshape((2,2)).tolist(),
                         neutral.omega.tolist())
        self.assertEqual(attach_rotation_xml('<root/>',[],[0.]),'<root/>')

    def test_neutral_omega_handles_imas_string_0d(self):
        """Real IMAS code.parameters is IDSString0D, not a Python str."""
        class IDSString0D:
            def __init__(self, value):
                self.value = value
            def __str__(self):
                raise AssertionError('must read IDSString0D.value, not str(wrapper)')

        neutral = types.SimpleNamespace(
            label='D_VOLUME', rotation_source='OM0VD', rotation_sign=1.0,
            omega=np.array([[1., 2.], [3., 4.]]))
        original = '<transp2imas_parameters><input_files complete="true"/></transp2imas_parameters>'
        xml = attach_rotation_xml(IDSString0D(original), [neutral], [0.0, 1.0])
        self.assertIsInstance(xml, str)
        root = ET.fromstring(xml)
        self.assertIsNotNone(root.find('input_files'))
        self.assertEqual(root.find('neutral_toroidal_rotation/profile').get('source_variable'),
                         'OM0VD')
        self.assertEqual(attach_rotation_xml(IDSString0D(original), [], [0.0]), original)

        # The input XML may legitimately be empty before provenance is loaded.
        empty = ET.fromstring(attach_rotation_xml(
            IDSString0D(''), [neutral], [0.0, 1.0]))
        self.assertEqual(empty.tag, 'transp2imas_parameters')
        self.assertIsNotNone(empty.find('neutral_toroidal_rotation'))
        self.assertEqual(attach_rotation_xml(IDSString0D(''), [], []), '')
        self.assertEqual(attach_rotation_xml(IDSString0D(np.str_(original)), [], []),
                         original)
        self.assertEqual(attach_rotation_xml(IDSString0D(original.encode('utf-8')),
                                             [], []), original)
        with self.assertRaisesRegex(TypeError, 'XML string'):
            attach_rotation_xml(IDSString0D(123), [neutral], [0.0, 1.0])

    def test_coverage_audit_distinguishes_dynamic_fields_from_aggregates(self):
        from audit_mapping_coverage import _dynamic_source_references
        names = ['GFLNC_D', 'GFLNC_E', 'GFLNC_I', 'GFLNC_X',
                 'QFLNC_I', 'QFLNCC_X', 'DN0VD', 'T0VD', 'OM0VD']
        mapped = _dynamic_source_references(names)
        self.assertTrue({'GFLNC_D', 'GFLNC_E', 'QFLNC_I',
                         'DN0VD', 'T0VD', 'OM0VD'}.issubset(mapped))
        self.assertNotIn('GFLNC_I', mapped)
        self.assertNotIn('GFLNC_X', mapped)
        self.assertNotIn('QFLNCC_X', mapped)

    def test_source_contract_rejects_wrong_dimension_or_units(self):
        class FakeVar:
            def __init__(self,dimensions,units):
                self.dimensions=dimensions
                self.units=units
            def __getitem__(self,item):
                return np.ones((2,3))
        dataset=types.SimpleNamespace(variables={'BAD':FakeVar(('TIME','X'),'EV')})
        with self.assertRaisesRegex(ValueError,'expected'):
            validate_profile(dataset,'BAD','density')
        dataset.variables['BAD']=FakeVar(('TIME3','X'),'EV')
        with self.assertRaisesRegex(ValueError,'unexpected'):
            validate_profile(dataset,'BAD','density')
        self.assertEqual(validate_profile(dataset,'BAD','temperature').shape,(2,3))

    def test_writer_neutral_temperature_has_no_density_scaling(self):
        src=(ROOT/'imas_core_profiles.py').read_text()
        self.assertIn("spec['temperature'], 'temperature', shape)", src)
        self.assertIn("spec['density'], 'density', shape)*1.0e6", src)
        self.assertNotIn("get_ncvar(f'{tstr}{symb}', 1.0e6)", src)
        self.assertIn('attach_rotation_xml(', src)


class FakeNode:
    """Small schema-shaped object to exercise actual writer assignments."""

    def __init__(self):
        self._items = []

    def __getattr__(self, name):
        if name.startswith('_'):
            raise AttributeError(name)
        node = FakeNode()
        setattr(self, name, node)
        return node

    def resize(self, count, **kwargs):
        length = int(count)
        self._items = (self._items[:length] +
                       [FakeNode() for _ in range(max(0, length-len(self._items)))])

    def __len__(self):
        return len(self._items)

    def __getitem__(self, index):
        return self._items[index]


class FakeConverter:
    def __init__(self, nc):
        self.nc = nc
        self.vars_ = list(nc.variables)
        self.cm2m = .01
        self.aconv = 1e-4
        self.vconv = 1e-6
        self.vol_inv = 1e6
        self.area_inv = 1e4
        self.ip_sign = 1.0
        self.toroidal_rotation_sign = -1.0
        self.logger = types.SimpleNamespace(info=lambda *a, **k: None,
                                            warning=lambda *a, **k: None)
        self.written = None

    def get_ncvar(self, name, scale=1., **kw):
        if name not in self.nc.variables:
            return np.array([])
        value = np.asarray(self.nc.variables[name].data, dtype=np.float64)
        if value.ndim == 0:
            value = value.reshape(1)
        return value * scale

    def is_valid(self, *args, **kwargs):
        return all(isinstance(x, np.ndarray) and x.size>0 for x in args)

    def check_var_name(self, name):
        return name in self.nc.variables

    def new_ids(self, name):
        return FakeNode()

    def fill_properties(self, ids):
        return ids

    def fill_code_properties(self, ids):
        ids.code.parameters = '<transp2imas_parameters/>'

    def fill_external_model_inputs(self, *args):
        return []

    def fill_vac_tor_field(self, ids):
        return None

    def fill_data(self, node, values, ind=-1, chk=-1, **kwargs):
        if not self.is_valid(values):
            return False
        data=np.asarray(values, dtype=np.float64)
        if ind >= 0:
            data=data[ind]
        if chk > 0 and len(data) != chk:
            raise AssertionError(f'writer assigned wrong radial count {len(data)} vs {chk}')
        node.values=np.asarray(data).copy()
        return True

    def cumulative_zone_quantity(self, values):
        return np.cumsum(values,axis=1)

    def interp_radial_profiles(self, values, xfrom, xto, axis_value=None):
        return np.array([np.interp(x, xf, y) for xf, y, x in zip(xfrom,values,xto)])

    def dd_at_least(self,major,minor=0,patch=0):
        return (4,1,1)>=(major,minor,patch)

    def set_species_name(self, node, name):
        node.name=name

    def set_nuclear_charge(self,node,charge):
        node.z_n=int(round(charge))

    def putdb_imas(self, ids):
        self.written=ids


class TestWriterIntegration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # The regular tests use no IMAS runtime. A stub permits importing the
        # real writer methods and inspecting the IDS-shaped objects they fill.
        with mock.patch.dict(sys.modules, {
            'imas': types.ModuleType('imas'),
            'netCDF4': types.ModuleType('netCDF4'),
        }):
            sys.modules['netCDF4'].Dataset = object
            from imas_core_transport import imas_core_transport
            from imas_core_profiles import imas_core_profiles
        cls.Transport = imas_core_transport
        cls.Profiles = imas_core_profiles

    def test_actual_neoclassical_writer_and_species(self):
        path=ROOT/'examples/data/d3d/147634T13/147634T13.CDF'
        if not path.is_file():
            self.skipTest('DIII-D T13 local fixture not present')
        with netcdf_file(path,'r',mmap=False) as nc:
            im=FakeConverter(nc)
            writer=self.Transport.__new__(self.Transport)
            writer.im=im
            writer.logger=im.logger
            writer.time=im.get_ncvar('TIME3')
            writer.x=im.get_ncvar('X')
            writer.xb=im.get_ncvar('XB')
            writer.dvol=im.get_ncvar('DVOL',im.vconv)
            writer.surf=im.get_ncvar('SURF',im.aconv)
            writer.nt=len(writer.time)
            writer.nx=writer.xb.shape[1]
            writer.create_ids()
            ct=im.written
            labels=[ct.model[i].identifier.name for i in range(len(ct.model))]
            self.assertIn('neoclassical',labels)
            self.assertIn('neoclassical_classical_heat',labels)
            ncmodel=ct.model[labels.index('neoclassical')]
            profile=ncmodel.profiles_1d[0]
            # Aggregate ion energy channel is not invented as an ion species.
            self.assertEqual([profile.ion[j].name for j in range(len(profile.ion))], ['D'])
            particle=profile.ion[0].particles.flux.values
            expected=divergence_to_surface_flux(
                validate_profile(nc,'GFLNC_D','particles'),
                im.get_ncvar('DVOL'), im.get_ncvar('SURF'),'particles')
            np.testing.assert_allclose(particle, expected[0],rtol=3e-14,atol=1e-14)
            self.assertTrue(np.isfinite(profile.electrons.energy.flux.values).all())
            self.assertTrue(np.isfinite(profile.total_ion_energy.flux.values).all())
            ncclass=ct.model[labels.index('neoclassical_classical_heat')]
            self.assertEqual(len(ncclass.profiles_1d[0].ion),1)
            self.assertNotIn('GFLNC_I',ncmodel.comment.split('Source channels: ')[1])

    def test_actual_transport_writer_on_every_available_cdf(self):
        """Exercise the actual transport writer against all bundled examples."""
        examples = sorted((ROOT / 'examples/data').rglob('*.CDF'))
        if not examples:
            self.skipTest('no local TRANSP CDF fixtures')
        for path in examples:
            with self.subTest(case=path.name), netcdf_file(path, 'r', mmap=False) as nc:
                im = FakeConverter(nc)
                writer = self.Transport.__new__(self.Transport)
                writer.im = im
                writer.logger = im.logger
                writer.time = im.get_ncvar('TIME3')
                writer.x = im.get_ncvar('X')
                writer.xb = im.get_ncvar('XB')
                writer.dvol = im.get_ncvar('DVOL', im.vconv)
                writer.surf = im.get_ncvar('SURF', im.aconv)
                writer.nt = len(writer.time)
                writer.nx = writer.xb.shape[1]
                writer.create_ids()
                ids = im.written
                models = {ids.model[i].identifier.name: ids.model[i]
                          for i in range(len(ids.model))}
                self.assertIn('neoclassical', models)
                self.assertIn('neoclassical_classical_heat', models)
                for model in models.values():
                    self.assertEqual(len(model.profiles_1d), writer.nt)
                profile = models['neoclassical'].profiles_1d[0]
                species = [profile.ion[j].name for j in range(len(profile.ion))]
                self.assertNotIn('I', species)
                self.assertNotIn('X', species)
                self.assertNotIn('TOK', species)
                self.assertEqual(len(species), len(set(species)))
                neutral_writer = self.Profiles.__new__(self.Profiles)
                neutral_writer.im = im
                neutral_writer.logger = im.logger
                neutrals = neutral_writer.get_neutral_data()
                self.assertTrue(neutrals)
                self.assertTrue(all(n.dens.shape == writer.x.shape for n in neutrals))
                xml = attach_rotation_xml('<root/>', neutrals, writer.time)
                self.assertEqual(
                    len(ET.fromstring(xml).findall('neutral_toroidal_rotation/profile')),
                    sum(int(n.omega.size > 0) for n in neutrals))

    def test_actual_neutral_units_and_rotation(self):
        path=ROOT/'examples/data/d3d/147634T13/147634T13.CDF'
        if not path.is_file():
            self.skipTest('DIII-D T13 local fixture not present')
        with netcdf_file(path,'r',mmap=False) as nc:
            im=FakeConverter(nc)
            writer=self.Profiles.__new__(self.Profiles)
            writer.im=im
            writer.logger=im.logger
            neutrals=writer.get_neutral_data()
            self.assertEqual([n.origin for n in neutrals],['wall','volume','recombination'])
            vol=next(n for n in neutrals if n.origin=='volume')
            np.testing.assert_allclose(vol.dens,
                validate_profile(nc,'DN0VD','density')*1e6)
            np.testing.assert_allclose(vol.temp,
                validate_profile(nc,'T0VD','temperature'))
            np.testing.assert_allclose(vol.omega,
                -validate_profile(nc,'OM0VD','rotation'))
            xml=attach_rotation_xml('<root/>',neutrals,im.get_ncvar('TIME3'))
            group=ET.fromstring(xml).find('neutral_toroidal_rotation')
            self.assertEqual(len(group.findall('profile')),3)
            self.assertTrue(writer._neutral_midplane_radius().shape==vol.dens.shape)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-root',type=Path,default=ROOT/'examples/data')
    p.add_argument('--report-dir',type=Path,default=ROOT/'build/neoclassical-neutrals')
    p.add_argument('--quick',action='store_true',help='synthetic test contracts only')
    ns=p.parse_args(argv)
    loader = unittest.TestLoader()
    suite=unittest.TestSuite((loader.loadTestsFromTestCase(TestMappingContracts),
                              loader.loadTestsFromTestCase(TestWriterIntegration)))
    outcome=unittest.TextTestRunner(verbosity=1).run(suite)
    if not outcome.wasSuccessful():
        return 1
    if ns.quick:
        return 0
    files=sorted(ns.data_root.rglob('*.CDF'))
    if not files:
        p.error('no bundled TRANSP CDF examples')
    results=[]
    errors=[]
    for path in files:
        try:
            with netcdf_file(path,'r',mmap=False) as f:
                item=verify_case(f,path.relative_to(ns.data_root))
            results.append(item)
            count=sum(len(v) for v in item['neoclassical'].values())
            neut=len(item['neutral_populations'])
            print(f'PASS {item["case"]}: neoclassical channels={count}, neutral populations={neut}, '
                  f'rotation channels={sum(int(n["rotation_preserved"]) for n in item["neutral_populations"])}')
        except Exception as exc:
            errors.append(f'{path}: {type(exc).__name__}: {exc}')
            print('FAIL',errors[-1])
    ns.report_dir.mkdir(parents=True,exist_ok=True)
    report={'checked':len(files),'passed':len(results),'failed':len(errors),
            'cases':results,'errors':errors}
    (ns.report_dir/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    (ns.report_dir/'summary.md').write_text(
        f'# Neoclassical/neutral TRANSP mapping source verification\n\n'
        f'Cases: {len(results)}/{len(files)} passed.\n\n' + '\n'.join(
            f"- `{r['case']}`: {sum(map(len,r['neoclassical'].values()))} NC channels, "
            f"{len(r['neutral_populations'])} neutral populations" for r in results)+'\n')
    return 1 if errors else 0


if __name__=='__main__':
    raise SystemExit(main())
