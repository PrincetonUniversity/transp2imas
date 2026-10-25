#!/usr/bin/python
"""Verify the services extracted from ``nc2imas``: namelist, time/radial grid, IDS identity, provenance.

Pure tests (no IMAS, no CDF) of ``namelist_reader``, ``name_patterns``, ``grid_services``, ``ids_identity``
and ``conversion_provenance``.  The same behaviour end to end is covered by ``compare_conversions.py``
(identical IDS output before and after the extraction) and by the verifiers that drive ``nc2imas``.

    python examples/verify_conversion_services.py
"""

import logging
from pathlib import Path
import sys
import tempfile
import types
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import grid_services as gs
import ids_identity as ii
import name_patterns as npat
import fortran_literals as fl
import namelist_reader as nl
import plasma_betas as pb
from conversion_plan import ConversionJournal
from conversion_provenance import ProvenanceRecorder
from mapping_status import MappingError, MappingLedger, RecordingDataset, set_strict

LOG = logging.getLogger('verify-conversion-services')


class NamelistTests(unittest.TestCase):
    def test_logicals_parse_case_insensitively_and_unknown_is_none(self):
        for text, expected in (('T', True), ('.true.', True), ('no', False), ('.F.', False), (1, True), (0, False)):
            self.assertIs(nl.parse_bool_value(text), expected, text)
        self.assertIsNone(nl.parse_bool_value('maybe'))
        self.assertIsNone(nl.parse_bool_value(None))

    def test_get_bool_finds_a_key_whatever_its_case(self):
        self.assertTrue(nl.get_bool({'nljccw': '.T.'}, 'NLJCCW'))
        self.assertEqual(nl.get_bool({}, 'NLJCCW', default='dflt'), 'dflt')
        self.assertEqual(nl.get_bool({'A': 'maybe'}, 'A', default=7), 7)

    def test_arithmetic_in_values(self):
        self.assertEqual(nl.parse_number('2*3', float), 6.0)
        self.assertEqual(nl.parse_number('1.0E-3', float), 1.0e-3)       # the sign after E is not an operator
        self.assertEqual(nl.parse_number('1.5E+3', float), 1500.0)
        self.assertEqual(nl.parse_number('2*3+1', float), 7.0)
        self.assertEqual(nl.parse_number(5, int), 5)

    def test_load_casts_numbers_and_skips_comments(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'TR.DAT'
            path.write_text("! comment\n&INNML\nNLJCCW = T\nAIMP = 12.0d0 ! carbon\nNLBCCW=1\nFREQ(1) = 1.5E9\n/\n",
                            encoding='ISO-8859-1')
            data = nl.load_namelist(str(path), LOG)
        self.assertEqual(data['AIMP'], 12.0)
        self.assertEqual(data['NLBCCW'], 1)
        self.assertEqual(data['NLJCCW'], 'T')
        self.assertEqual(data['FREQ(1)'], 1.5e9)
        self.assertNotIn('&INNML', data)

    def test_fortran_d_exponents_are_read_not_turned_into_zero(self):
        for text, expected in (('1.0D-3', 1.0e-3), ('2.5d0', 2.5), ('1D3', 1000.0), ('9.205d-2', 0.09205),
                               ('1.142d3', 1142.0), ('-4.D1', -40.0)):
            self.assertEqual(nl.parse_number(text, float), expected, text)
        self.assertEqual(nl.parse_number('1.0D-3*2', float), 2.0e-3)           # inside an expression too
        self.assertEqual(nl.fortran_number('5D-3'), 5.0e-3)
        self.assertEqual(nl.fortran_number('7'), 7)
        self.assertIsNone(nl.fortran_number('abc'))
        self.assertIsNone(nl.fortran_number("'1.0D-3'"))                       # a quoted string stays a string

    def test_names_are_not_mistaken_for_exponents(self):
        self.assertEqual(nl.normalize_exponents('ADD3 + 1D2'), 'ADD3 + 1E2')

    def test_load_reads_d_exponents_and_keeps_a_bang_inside_quotes(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'TR.DAT'
            path.write_text("TBONA(1) = 9.205d-2\nXLBTNA(3) = 1.142d3\nFILE='beam!config.nml' ! the comment\nN=7 ! seven\n",
                            encoding='ISO-8859-1')
            data = nl.load_namelist(str(path), LOG)
        self.assertEqual(data['TBONA(1)'], 0.09205)
        self.assertEqual(data['XLBTNA(3)'], 1142.0)
        self.assertEqual(data['FILE'], "'beam!config.nml'")
        self.assertEqual(data['N'], 7)

    def test_an_unparseable_number_raises_and_is_unknown_to_the_reader_not_zero(self):
        with self.assertRaises(ValueError):
            nl.parse_number('abc', float)
        with self.assertRaises(ValueError):
            nl.parse_number('10**3', float)                                    # unsupported form: never a plausible 0
        data = {'E': 'abc', 'A': '1.0D-3'}
        with self.assertLogs('verify-conversion-services', level='WARNING') as logged:
            self.assertIsNone(nl.get_var_(data, LOG, 'E', float))
        self.assertIn('treated as unknown', logged.output[0])
        self.assertEqual(nl.get_var_(data, LOG, 'A', float), 1.0e-3)
        self.assertIsNone(nl.get_var_({'L': ['1.0', 'x']}, LOG, 'L', float))   # one bad element makes the list unknown

    def test_the_required_literal_cases_from_the_review(self):
        cases = (('1.0D-3', 0.001), ('-2.5d+2', -250.0), ('2.5D0', 2.5), ('+.5e1', 5.0), ('-4.', -4.0), ('12', 12))
        for text, expected in cases:
            self.assertEqual(fl.parse_number(text), expected, text)
            self.assertEqual(nl.parse_number(text, float), float(expected), text)

    def test_malformed_numbers_are_unknown_never_zero_false_or_a_previous_value(self):
        for text in ('1.0D', 'D3', '1.0.0', '1e', '--1', '1 2', 'nan', 'inf', '1.0q3', '', '1.0D-3x', '1/0'):
            self.assertIsNone(fl.parse_number(text), text)
            self.assertIsNone(fl.parse_scalar(text), text)
            with self.assertRaises(ValueError):
                nl.parse_number(text, float)
        self.assertIs(fl.parse_scalar('.false.'), False)             # a real False is not the same as unknown
        self.assertIsNone(fl.parse_scalar('maybe'))

    def test_a_later_unparseable_assignment_does_not_keep_the_earlier_value(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'TR.DAT'
            path.write_text('A = 1.0\nA = 1.0.0\nB = 2.0\n', encoding='ISO-8859-1')
            records = nl.read_namelist(str(path), LOG)
        self.assertEqual(records.raw['A'], '1.0.0')
        self.assertEqual(records.lines['A'], 2)                       # the assignment that supplied the value
        self.assertIsNone(nl.get_var_(records.values, LOG, 'A', float))
        self.assertEqual(nl.get_var_(records.values, LOG, 'B', float), 2.0)

    def test_records_keep_raw_text_source_lines_and_gaps(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'TR.DAT'
            path.write_text("&TRIN LPREDICT_NE = 1\nTE = 1.5D-3 ! keV\nFILE='beam!x.nml'\nA=1, B=2\ninclude 'o.nml'\n/\n",
                            encoding='ISO-8859-1')
            records = nl.read_namelist(str(path), LOG)
            dictionary_view = nl.load_namelist(str(path), LOG)
        self.assertEqual(records.values['LPREDICT_NE'], 1)             # an assignment on the group-header line
        self.assertEqual((records.raw['TE'], records.lines['TE'], records.values['TE']), ('1.5D-3', 2, 0.0015))
        self.assertEqual(records.values['FILE'], "'beam!x.nml'")
        reasons = [(gap['line'], gap['reason'].split(';')[0]) for gap in records.gaps]
        self.assertIn((4, 'several assignments on one line'), reasons)
        self.assertIn((5, 'statement without an assignment'), reasons)
        self.assertEqual(dictionary_view, records.values)             # the dictionary view is unchanged

    def test_the_converter_and_curation_read_the_same_scalars(self):
        from curation.inputs import literal
        for text in ('1.0D-3', '-2.5d+2', '7', '.TRUE.', 'f', "'a!b'", "'it''s'", '1.0D', '1.0.0', 'abc', 'nan', '[1,2]'):
            shared = fl.parse_scalar(text)
            curated = literal(text)
            if isinstance(shared, bool) or shared is None or isinstance(shared, str):
                self.assertEqual(curated, shared, text)
            else:
                self.assertEqual(curated, float(shared), text)

    def test_missing_file_gives_an_empty_namelist_and_no_path_gives_none(self):
        self.assertEqual(nl.load_namelist('/nonexistent/TR.DAT', LOG), {})
        self.assertIsNone(nl.load_namelist(None, LOG))

    def test_scalar_array_and_indexed_values(self):
        data = {'A': 2.5, 'B(1)': 1.0, 'B(2)': 2.0, 'C': '1.0 2.0 3.0'}
        self.assertEqual(nl.get_vars(data, LOG, 'A', float), 2.5)
        np.testing.assert_array_equal(nl.get_vars(data, LOG, 'B', float), [1.0, 2.0])
        self.assertEqual(nl.get_vars(data, LOG, 'B', float, ind=1), 2.0)
        np.testing.assert_array_equal(nl.get_var_(data, LOG, 'C', float), [1.0, 2.0, 3.0])
        self.assertIsNone(nl.get_var_(data, LOG, 'MISSING'))
        self.assertEqual(nl.count_vars(data, LOG, 'B?'), 2)
        self.assertEqual(nl.count_vars(None, LOG, 'B?'), 0)

    def test_scale_and_resize(self):
        self.assertEqual(nl.get_var_({'A': 2.0}, LOG, 'A', float, scale=3.0), 6.0)
        np.testing.assert_array_equal(nl.get_vars({'A': 2.0}, LOG, 'A', float, size=3, resize=True), [2.0, 2.0, 2.0])
        with self.assertRaises(Exception):
            nl.get_vars({'B(1)': 1.0, 'B(2)': 2.0}, LOG, 'B', float, size=3)

    def test_keys_by_pattern_and_text_between_markers(self):
        keys = nl.keys_re({'XECECH1': 1, 'XECECH12': 2, 'XECECH': 3, 'OTHER': 4}, 'XECECH?', nd1=1, nd2=2)
        self.assertEqual(sorted(keys), ['XECECH1', 'XECECH12'])
        # the text starts one character after the first marker and ends at the next separator after it
        self.assertEqual(nl.extract_str_between(LOG, 'AA(Bxyz,rest)', '(', ',', 'dflt'), 'xyz')
        self.assertIsNone(nl.extract_str_between(LOG, 'no markers', '(', ',', 'dflt'))


class PatternTests(unittest.TestCase):
    def test_a_name_needs_a_placeholder(self):
        with self.assertRaises(Exception):
            npat.make_re_pattern('NOQUEST')

    def test_digit_ranges_and_two_placeholders(self):
        self.assertEqual(npat.make_re_pattern('VV_PF?', nd=1, nch1=1, nch2=2), r'VV_PF\d{1,1}[A-z]{1,2}$')
        self.assertEqual(npat.make_re_pattern_quest('BSTH_?4_P?'), r'BSTH_\w+4_P\w+$')
        self.assertEqual(npat.make_namelist_re_pattern(r'NAME\(?,?\)', nd1=1, nd2=3),
                         r'NAME\(\d{1,3}[A-z]{0,0},\d{1,3}[A-z]{0,0}\)$')


class GridTests(unittest.TestCase):
    @staticmethod
    def valid(*values, ind=0, name=''):
        return all(v is not None and np.size(v) > ind for v in values)

    def test_time_bases_need_two_samples_of_each(self):
        data = {'TIME': np.array([0.0, 1.0, 2.0]), 'TIME3': np.array([0.0, 1.0, 2.0]), 'X': np.zeros((3, 4))}
        get = lambda name: data[name]
        timesc, ntimes, timep, ntimep = gs.read_time_bases(get, LOG)
        self.assertEqual((ntimes, ntimep), (3, 3))
        data['TIME3'] = np.array([0.0])
        with self.assertRaises(Exception):
            gs.read_time_bases(get, LOG)

    def test_interp_time_compares_coordinates_not_lengths(self):
        time = np.array([0.0, 1.0, 2.0])
        values = np.array([0.0, 10.0, 20.0])
        np.testing.assert_allclose(gs.interp_time(self.valid, np.array([0.0, 0.5, 2.0]), time, values), [0.0, 5.0, 20.0])
        self.assertIs(gs.interp_time(self.valid, time, time, values), values)       # same coordinates: untouched
        self.assertEqual(gs.interp_time(self.valid, np.array([]), time, values), [])

    def test_history_on_time_uses_the_declared_dimension(self):
        var = lambda dim: types.SimpleNamespace(dimensions=(dim,))
        nc = types.SimpleNamespace(variables={'A': var('TIME'), 'B': var('TIME3'), 'P': types.SimpleNamespace(dimensions=('TIME3', 'X'))})
        timesc, timep = np.array([0.0, 1.0, 2.0]), np.array([0.0, 0.5, 2.0])
        get = lambda name, scale=1.0: scale * np.array([0.0, 10.0, 20.0])
        np.testing.assert_allclose(gs.history_on_time(nc, get, timesc, timep, 'A'), [0.0, 5.0, 20.0])
        np.testing.assert_allclose(gs.history_on_time(nc, get, timesc, timep, 'B', target='TIME'), [0.0, 40.0 / 3.0, 20.0])
        self.assertEqual(gs.history_on_time(nc, get, timesc, timep, 'MISSING').size, 0)
        with self.assertRaises(ValueError):
            gs.history_on_time(nc, get, timesc, timep, 'P')                 # a profile is not a history
        with self.assertRaises(ValueError):
            gs.history_on_time(nc, get, timesc, timep, 'A', target='OTHER')

    def test_zone_quantities_accumulate_to_enclosed_values(self):
        np.testing.assert_array_equal(gs.cumulative_zone_quantity([[1.0, 2.0, 3.0]]), [[1.0, 3.0, 6.0]])
        self.assertEqual(gs.cumulative_zone_quantity([]).size, 0)

    def test_radial_gradient_of_a_parabola_and_of_unusable_input(self):
        x = np.tile(np.linspace(0.0, 1.0, 21), (2, 1))
        grad = gs.radial_gradient(x ** 2, x)
        np.testing.assert_allclose(grad, 2.0 * x, atol=1e-9)
        self.assertEqual(gs.radial_gradient(x, x[:1]).size, 0)              # shape mismatch
        bad = x.copy()
        bad[0] = np.nan
        self.assertTrue(np.all(np.isnan(gs.radial_gradient(x ** 2, bad)[0])))


class Node:
    pass


class IdentityTests(unittest.TestCase):
    @staticmethod
    def at(version):
        return lambda major, minor=0, patch=0: version >= (major, minor, patch)

    def test_species_text_is_label_before_3_42_and_name_after(self):
        old, new = Node(), Node()
        ii.set_species_name(self.at((3, 39, 0)), old, 'D')
        ii.set_species_name(self.at((4, 1, 1)), new, 'D')
        self.assertEqual((old.label, new.name), ('D', 'D'))
        self.assertFalse(hasattr(old, 'name') or hasattr(new, 'label'))

    def test_z_n_is_integer_from_dd4_and_float_before(self):
        old, new = Node(), Node()
        ii.set_nuclear_charge(self.at((3, 39, 0)), old, 5.0)
        ii.set_nuclear_charge(self.at((4, 1, 1)), new, 5.0)
        self.assertIsInstance(old.z_n, float)
        self.assertIsInstance(new.z_n, int)

    def test_beam_and_unit_identity_swap_roles_between_dd3_and_dd4(self):
        for setter, label in ((ii.set_ec_beam_identity, 'beam'), (ii.set_nbi_unit_identity, 'unit')):
            dd3, dd4 = Node(), Node()
            setter(self.at((3, 39, 0)), dd3, 'ID', 'long text')
            setter(self.at((4, 1, 1)), dd4, 'ID', 'long text')
            self.assertEqual((dd3.identifier, dd3.name), ('ID', 'long text'), label)
            self.assertEqual((dd4.name, dd4.description), ('ID', 'long text'), label)

    def test_ec_mode_by_dd_version(self):
        dd3, dd4 = Node(), Node()
        dd4.polarization = Node()
        ii.set_ec_mode(self.at((3, 39, 0)), dd3, 0.8, [0.0, 1.0])
        self.assertEqual(dd3.mode, 1)
        ii.set_ec_mode(self.at((3, 39, 0)), dd3, 0.2, [0.0, 1.0])
        self.assertEqual(dd3.mode, -1)
        ii.set_ec_mode(self.at((4, 1, 1)), dd4, 1.7, [0.0, 1.0])            # clipped to 1
        np.testing.assert_array_equal(dd4.polarization.o_mode_fraction, [1.0, 1.0])

    def test_helper_binds_the_predicate(self):
        node = Node()
        ii.IdentityHelper(self.at((4, 1, 1))).set_species_name(node, 'T')
        self.assertEqual(node.name, 'T')


class ProvenanceTests(unittest.TestCase):
    def source(self):
        return types.SimpleNamespace(nml_name=None, ncfile='x.CDF', device='ITER', shot=1, run=2, occur=0, dd_version='4.1.1',
                                     signs=None, logger=LOG, mapping=None,
                                     get_ncattr=lambda key: None, get_ncattr_text=lambda key, default='': default)

    def test_the_inventory_is_cached_until_reset(self):
        recorder = ProvenanceRecorder(self.source(), '0.0.0')
        first = recorder.input_inventory()
        self.assertIs(recorder.input_inventory(), first)
        recorder.reset()
        self.assertIsNot(recorder.input_inventory(), first)

    def test_code_properties_name_the_converter_and_skip_ids_without_a_code_node(self):
        recorder = ProvenanceRecorder(self.source(), '9.9.9')
        recorder._code_parameters_xml = '<transp2imas_parameters/>'
        ids = types.SimpleNamespace(code=Node())
        ids.code.description = ids.code.version = ids.code.repository = ids.code.parameters = ''   # the optional nodes
        recorder.fill_code_properties(ids)
        self.assertEqual((ids.code.name, ids.code.version, ids.code.parameters),
                         ('transp2imas', '9.9.9', '<transp2imas_parameters/>'))
        bare = Node()
        self.assertIs(recorder.fill_code_properties(bare), bare)
        self.assertEqual(recorder.fill_external_model_inputs(bare, 'toric'), [])


class BetaTests(unittest.TestCase):
    """The IMAS beta definitions on cases with a closed form."""

    nt, nx = 3, 5

    def case(self, pressure=2.0e5, dvol=2.0, b0=5.0, r0=6.0, ip=1.0e7, a=2.0):
        return (np.full((self.nt, self.nx), pressure), np.full((self.nt, self.nx), dvol),
                np.full(self.nt, b0), r0, np.full(self.nt, ip), np.full(self.nt, a))

    def test_beta_tor_is_twice_mu0_p_over_b0_squared_for_uniform_pressure(self):
        result = pb.imas_betas(*self.case())
        np.testing.assert_allclose(result['beta_tor'], 2 * pb.MU0 * 2.0e5 / 5.0 ** 2, rtol=1e-14)

    def test_beta_pol_uses_the_integral_of_p_and_r0(self):
        result = pb.imas_betas(*self.case())
        integral = 2.0e5 * 2.0 * self.nx                                  # int p dV in J
        np.testing.assert_allclose(result['beta_pol'], 4 * integral / (6.0 * pb.MU0 * (1.0e7) ** 2), rtol=1e-14)

    def test_normalised_beta_follows_the_dd_formula(self):
        result = pb.imas_betas(*self.case())
        np.testing.assert_allclose(result['beta_tor_norm'], 100 * result['beta_tor'] * 2.0 * 5.0 / 10.0, rtol=1e-14)

    def test_beta_pol_and_beta_tor_ignore_the_signs_but_the_normalised_beta_keeps_them(self):
        base = pb.imas_betas(*self.case())
        flipped = pb.imas_betas(*self.case(b0=-5.0, ip=-1.0e7))
        np.testing.assert_allclose(flipped['beta_tor'], base['beta_tor'])
        np.testing.assert_allclose(flipped['beta_pol'], base['beta_pol'])
        np.testing.assert_allclose(flipped['beta_tor_norm'], base['beta_tor_norm'])       # B0/Ip keeps its sign
        counter = pb.imas_betas(*self.case(b0=-5.0))
        np.testing.assert_allclose(counter['beta_tor_norm'], -base['beta_tor_norm'])

    def test_a_fixed_radius_scales_the_betas_as_the_definitions_say(self):
        # same BZXR = r0*b0: beta_tor goes as r0^2 (b0 = BZXR/r0), beta_pol as 1/r0
        a = pb.imas_betas(*self.case(b0=36.0 / 6.0, r0=6.0))
        b = pb.imas_betas(*self.case(b0=36.0 / 7.0, r0=7.0))
        np.testing.assert_allclose(b['beta_tor'] / a['beta_tor'], (7.0 / 6.0) ** 2, rtol=1e-13)
        np.testing.assert_allclose(b['beta_pol'] / a['beta_pol'], 6.0 / 7.0, rtol=1e-13)

    def test_inputs_that_do_not_fit_give_nothing(self):
        pressure, dvol, b0, r0, ip, a = self.case()
        self.assertEqual(pb.imas_betas(pressure, dvol[:, :-1], b0, r0, ip, a), {})
        self.assertEqual(pb.imas_betas(pressure, dvol, b0[:-1], r0, ip, a), {})
        self.assertEqual(pb.imas_betas(pressure[0], dvol[0], b0, r0, ip, a), {})

    def test_minor_radius_of_a_shifted_ellipse_is_its_half_width(self):
        moments = np.array([[1.7, 1.8], [0.5, 0.6]])                       # R = R0 + a cos(theta), two times
        np.testing.assert_allclose(pb.boundary_minor_radius(moments), [0.5, 0.6], rtol=1e-4)
        np.testing.assert_allclose(pb.boundary_minor_radius(moments + np.array([[3.0], [0.0]])), [0.5, 0.6], rtol=1e-4)

    def test_minor_radius_with_a_triangular_sin_harmonic(self):
        # R = R0 + a cos(theta + delta sin(theta)) has half-width a exactly; its first Fourier moments approximate it
        theta, a, delta = np.linspace(0, 2 * np.pi, 4000, endpoint=False), 0.6, 0.3
        radius = 1.7 + a * np.cos(theta + delta * np.sin(theta))
        cos = np.array([[np.mean(radius * np.cos(k * theta)) * (1 if k == 0 else 2)] for k in range(8)])
        sin = np.array([[np.mean(radius * np.sin(k * theta)) * 2] for k in range(8)])
        self.assertAlmostEqual(float(pb.boundary_minor_radius(cos, sin)[0]), a, places=3)


try:
    import netCDF4
except ImportError:      # pragma: no cover - environment dependent
    netCDF4 = None


@unittest.skipIf(netCDF4 is None, 'netCDF4 is required (it is what produces masked arrays)')
class MaskedSourceTests(unittest.TestCase):
    """A masked sample is a named broken mapping, never a fill value passed on as a physical number."""

    def connection(self, directory, values, fill=-999.0):
        from netcdf_class import netcdf_connection
        path = Path(directory) / 'masked.nc'
        with netCDF4.Dataset(path, 'w') as nc:
            nc.createDimension('TIME', len(values))
            var = nc.createVariable('TE', 'f8', ('TIME',), fill_value=fill)
            var[:] = np.ma.masked_values(values, fill)
            nc.createVariable('TI', 'f8', ('TIME',))[:] = np.arange(len(values), dtype=float)
        connection = netcdf_connection.__new__(netcdf_connection)
        connection.mapping = MappingLedger()
        connection.mapping.writer = 'test'
        connection.nc = RecordingDataset(netCDF4.Dataset(path), connection.mapping)
        return connection

    def test_masked_samples_are_a_broken_mapping_under_the_strict_policy(self):
        with tempfile.TemporaryDirectory() as tmp:
            connection = self.connection(tmp, [1.0, -999.0, 3.0])
            with self.assertRaises(MappingError) as caught:
                connection.get_ncvar('TE')
            self.assertIn('masked', str(caught.exception))
            self.assertIn('TE', str(caught.exception))

    def test_under_the_lenient_policy_the_variable_is_skipped_with_the_reason(self):
        with tempfile.TemporaryDirectory() as tmp:
            connection = self.connection(tmp, [1.0, -999.0, 3.0])
            set_strict(False)
            try:
                values = connection.get_ncvar('TE')
            finally:
                set_strict(True)
            self.assertEqual(values.size, 0)                                  # not [1, -999, 3]
            skipped = connection.mapping.entries('test')
            self.assertTrue(any(row['source'] == 'TE' and 'masked' in row['reason'] for row in skipped), skipped)

    def test_an_entirely_missing_array_is_a_broken_mapping_too(self):
        with tempfile.TemporaryDirectory() as tmp:
            connection = self.connection(tmp, [-999.0, -999.0, -999.0])
            with self.assertRaises(MappingError):
                connection.get_ncvar('TE')

    def test_legitimate_zeros_are_not_missing_for_floats_and_integers(self):
        with tempfile.TemporaryDirectory() as tmp:
            connection = self.connection(tmp, [0.0, 0.0, 0.0])
            np.testing.assert_array_equal(connection.get_ncvar('TE'), [0.0, 0.0, 0.0])
            self.assertEqual(connection.get_ncvar('TI', ind=0).tolist(), [0.0])      # integer-valued zero

    def test_the_returned_array_is_a_copy_the_caller_may_change(self):
        with tempfile.TemporaryDirectory() as tmp:
            connection = self.connection(tmp, [1.0, 2.0, 3.0])
            first = connection.get_ncvar('TE')
            first[:] = -1.0
            np.testing.assert_array_equal(connection.get_ncvar('TE'), [1.0, 2.0, 3.0])

    def test_a_variable_without_masked_samples_is_read_unchanged(self):
        with tempfile.TemporaryDirectory() as tmp:
            connection = self.connection(tmp, [1.0, 2.0, 3.0])
            np.testing.assert_array_equal(connection.get_ncvar('TE'), [1.0, 2.0, 3.0])
            np.testing.assert_array_equal(connection.get_ncvar('TI', scale=2.0), [0.0, 2.0, 4.0])
            np.testing.assert_array_equal(connection.get_ncvar('TI', ind=1), [1.0])


try:
    import nc2imas as nc2imas_module
except ImportError:      # pragma: no cover - environment dependent
    nc2imas_module = None


@unittest.skipIf(nc2imas_module is None, 'the converter module needs IMAS-Python')
class ProfileShapeTests(unittest.TestCase):
    """``nc2imas.check_profile_shape``: empty or exactly (nt, nx); everything else is a named broken mapping."""

    def check(self, array, nt=3, nx=3, strict=False):
        me = types.SimpleNamespace(mapping=MappingLedger())
        me.mapping.writer = 'test'
        set_strict(strict)
        try:
            result = nc2imas_module.nc2imas.check_profile_shape(me, 'VAR', array, nt, nx)
        finally:
            set_strict(True)
        return result, me.mapping.entries('test')

    def test_the_declared_shape_an_empty_array_and_histories_pass_unchanged(self):
        full = np.zeros((3, 3))
        self.assertIs(self.check(full)[0], full)
        result, entries = self.check(np.array([]))
        self.assertEqual((result.size, entries), (0, []))
        # a scalar or 1-D time history is not a profile; its length is that of TIME, which need not be TIME3's
        for history in (np.zeros(3), np.zeros(7), np.float64(2.0)):
            result, entries = self.check(history)
            self.assertIs(result, history)
            self.assertEqual(entries, [])

    def test_each_wrong_profile_shape_is_rejected_with_its_own_reason(self):
        for array, text in ((np.zeros((2, 3)), 'time length mismatch 2 != 3'),
                            (np.zeros((3, 4)), 'radial size mismatch 4 != 3'),
                            (np.zeros((3, 3, 2)), '3-dimensional')):
            result, entries = self.check(array)
            self.assertEqual(result.size, 0, text)
            self.assertTrue(any(text in row['reason'] and row['source'] == 'VAR' for row in entries), (text, entries))

    def check_declared(self, array, dimensions, nt=3, nx=3, grid=None):
        me = types.SimpleNamespace(mapping=MappingLedger())
        me.mapping.writer = 'test'
        set_strict(False)
        try:
            result = nc2imas_module.nc2imas.check_profile_shape(me, 'VAR', array, nt, nx, dimensions=dimensions, grid=grid)
        finally:
            set_strict(True)
        return result, me.mapping.entries('test')

    def test_equal_lengths_with_different_coordinates_are_rejected_naming_the_dimensions(self):
        full = np.zeros((3, 3))
        for dimensions, grid, text in ((('TIME', 'X'), None, 'declared time dimension TIME where TIME3 is expected'),
                                       (('TIME3', 'RMAJM'), None, 'declared radial dimension RMAJM is neither X nor XB'),
                                       (('TIME3', 'XB'), 'X', 'profile is on XB but the writer expects X'),
                                       (('TIME3', 'X'), 'XB', 'profile is on X but the writer expects XB')):
            result, entries = self.check_declared(full, dimensions, grid=grid)
            self.assertEqual(result.size, 0, text)
            reason = ' | '.join(row['reason'] for row in entries)
            self.assertIn(text, reason)
            self.assertIn(str(tuple(dimensions)), reason)                  # source dimensions named
            self.assertTrue(all(row['source'] == 'VAR' and 'profile (3 times x 3 radial points)' in row['target']
                                for row in entries), entries)                # variable and destination named

    def test_declared_coordinates_that_fit_pass_unchanged(self):
        full = np.zeros((3, 3))
        for dimensions, grid in ((('TIME3', 'X'), 'X'), (('TIME3', 'XB'), 'XB'), (('TIME3', 'X'), None),
                                 (('TIME3', 'XB'), None), (None, 'X')):
            result, entries = self.check_declared(full, dimensions, grid=grid)
            self.assertIs(result, full)
            self.assertEqual(entries, [])

    def test_a_valid_single_slice_passes(self):
        single = np.zeros((1, 4))
        result, entries = self.check_declared(single, ('TIME3', 'X'), nt=1, nx=4, grid='X')
        self.assertIs(result, single)
        self.assertEqual(entries, [])
        wrong, _ = self.check_declared(np.zeros((2, 4)), ('TIME3', 'X'), nt=1, nx=4)
        self.assertEqual(wrong.size, 0)

    def test_the_shape_is_judged_before_the_coordinates(self):
        _result, entries = self.check_declared(np.zeros((3, 4)), ('TIME', 'X'))
        self.assertIn('radial size mismatch 4 != 3', entries[0]['reason'])     # the length problem is the first reported

    def test_the_explicit_ways_across_are_the_documented_transformations(self):
        # TIME -> TIME3 for a 1-D history and X <-> XB for a profile are done by named methods, never by the check
        self.assertTrue(callable(nc2imas_module.nc2imas.get_ncvar_on_time))
        self.assertTrue(callable(nc2imas_module.nc2imas.interp_radial_profiles))
        x, xb = np.array([[0.25, 0.75]]), np.array([[0.5, 1.0]])
        moved = nc2imas_module.nc2imas.interp_radial_profiles(np.array([[1.0, 3.0]]), xb, x, axis_value=0.0)
        self.assertEqual(moved.shape, (1, 2))
        np.testing.assert_allclose(moved, [[0.5, 2.0]])      # (0 -> 0, 0.5 -> 1, 1.0 -> 3) sampled at 0.25 and 0.75

    def test_the_strict_policy_raises_naming_the_variable(self):
        for array in (np.zeros((2, 3)), np.zeros((3, 3, 2))):
            with self.assertRaises(MappingError) as caught:
                self.check(array, strict=True)
            self.assertIn('VAR', str(caught.exception))


@unittest.skipIf(netCDF4 is None, 'netCDF4 is required (the backend whose automatic handling is compared)')
class MissingDataContractTests(unittest.TestCase):
    """The converter (netCDF4's automatic mask and scale) and curation (explicit decoding) give one answer."""

    def write(self, directory):
        path = Path(directory) / 'contract.nc'
        with netCDF4.Dataset(path, 'w', format='NETCDF3_CLASSIC') as nc:
            nc.createDimension('T', 4)

            def var(name, dtype, data, **attrs):
                fill = attrs.pop('_FillValue', None)
                variable = nc.createVariable(name, dtype, ('T',), fill_value=fill)
                variable.set_auto_maskandscale(False)
                variable[:] = data
                for key, value in attrs.items():
                    variable.setncattr(key, value)

            var('fill', 'f8', [1, -999, 3, 0], _FillValue=-999.0)
            var('missing', 'f8', [1, -999, 3, 0], missing_value=-999.0)
            var('valid_range', 'f8', [1, 50, 3, 0], valid_range=np.array([0.0, 10.0]))
            var('valid_min', 'f8', [1, -5, 3, 0], valid_min=0.0)
            var('valid_max', 'f8', [1, 50, 3, 0], valid_max=10.0)
            var('packed', 'i2', [100, 200, 300, 0], scale_factor=0.5, add_offset=10.0)
            var('packed_fill', 'i2', [100, -32768, 300, 0], _FillValue=np.int16(-32768), scale_factor=0.5, add_offset=10.0)
            var('packed_range', 'i2', [100, 5000, 300, 0], valid_range=np.array([0, 1000], dtype='i2'),
                scale_factor=0.5, add_offset=10.0)
            var('plain_int', 'i4', [1, 2, 3, 0])
            var('all_zero', 'f8', [0, 0, 0, 0])
            var('all_missing', 'f8', [-999, -999, -999, -999], _FillValue=-999.0)
            var('default_fill', 'f8', [1.0, netCDF4.default_fillvals['f8'], 3.0, 0.0])
            var('default_fill_int', 'i4', [1, netCDF4.default_fillvals['i4'], 3, 0])
        return path

    @staticmethod
    def backend_view(path, name):
        """What the backend's own automatic masking and scaling returns, with missing samples as NaN."""
        with netCDF4.Dataset(path) as nc:
            data = nc.variables[name][...]
            return np.ma.filled(np.ma.asarray(data).astype(float), np.nan), np.ma.getmaskarray(data)

    def test_curation_decodes_every_case_as_the_backend_does_on_both_readers(self):
        from curation.source import NativeSource
        from scipy.io import netcdf_file
        with tempfile.TemporaryDirectory() as tmp:
            path = self.write(tmp)
            with netCDF4.Dataset(path) as nc4, netcdf_file(str(path), 'r', mmap=False) as nc3:
                nc4.set_auto_maskandscale(False)
                readers = {'netCDF4 (explicit)': NativeSource(nc4, packed=True),
                           'scipy classic': NativeSource(nc3, packed=True)}
                for name in nc4.variables:
                    expected, missing = self.backend_view(path, name)
                    for label, source in readers.items():
                        decoded = source.array(name)
                        np.testing.assert_array_equal(np.isnan(decoded), missing, f'{name} missingness, {label}')
                        np.testing.assert_allclose(decoded, expected, rtol=1e-12, err_msg=f'{name} values, {label}')

    def test_packing_is_applied_once_and_sentinels_are_judged_in_the_packed_representation(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self.write(tmp)
            values, mask = self.backend_view(path, 'packed')
            np.testing.assert_array_equal(values, [60.0, 110.0, 160.0, 10.0])    # raw*0.5 + 10, once
            values, mask = self.backend_view(path, 'packed_fill')
            self.assertEqual(mask.tolist(), [False, True, False, False])         # -32768 is the raw sentinel
            values, mask = self.backend_view(path, 'packed_range')
            self.assertEqual(mask.tolist(), [False, True, False, False])         # valid_range [0, 1000] bounds raw values
            np.testing.assert_array_equal(values[[0, 2, 3]], [60.0, 160.0, 10.0])
            values, mask = self.backend_view(path, 'all_zero')
            self.assertFalse(mask.any())                                         # legitimate zeros stay valid
            np.testing.assert_array_equal(values, 0.0)

    def test_the_converter_reads_packed_data_without_applying_the_packing_again(self):
        from netcdf_class import netcdf_connection
        with tempfile.TemporaryDirectory() as tmp:
            path = self.write(tmp)
            connection = netcdf_connection.__new__(netcdf_connection)
            connection.mapping = MappingLedger()
            connection.mapping.writer = 'test'
            connection.nc = RecordingDataset(netCDF4.Dataset(path), connection.mapping)
            np.testing.assert_array_equal(connection.get_ncvar('packed'), [60.0, 110.0, 160.0, 10.0])
            np.testing.assert_array_equal(connection.get_ncvar('plain_int'), [1, 2, 3, 0])
            np.testing.assert_array_equal(connection.get_ncvar('all_zero'), 0.0)
            set_strict(False)
            try:
                for name in ('packed_fill', 'packed_range', 'default_fill', 'default_fill_int', 'all_missing', 'fill'):
                    self.assertEqual(connection.get_ncvar(name).size, 0, name)     # never a sentinel passed on as a number
            finally:
                set_strict(True)


@unittest.skipIf(nc2imas_module is None, 'the converter module needs IMAS-Python')
class ConverterExitTests(unittest.TestCase):
    """``with nc2imas(...)`` closes the entry; only a run that ended without an exception is "completed"."""

    def exit_with(self, exc_type):
        from unittest import mock
        converter = nc2imas_module.nc2imas.__new__(nc2imas_module.nc2imas)
        converter.logger = LOG
        with mock.patch.object(nc2imas_module.nc2imas, 'close_imas_database') as close:
            result = converter.__exit__(exc_type, exc_type() if exc_type else None, None)
        self.assertFalse(result)                                    # exceptions are never swallowed
        return close.call_args.kwargs

    def test_a_clean_exit_is_completed_and_an_exception_is_not(self):
        self.assertEqual(self.exit_with(None), {'completed': True})
        self.assertEqual(self.exit_with(RuntimeError), {'completed': False})


class ReadCountTests(unittest.TestCase):
    def test_reads_are_counted_per_writer_and_the_set_of_names_is_unchanged(self):
        ledger = MappingLedger()
        ledger.writer = 'a'
        for _ in range(30):
            ledger.consume('DAREA')
        ledger.consume('X')
        ledger.consume('X')
        ledger.writer = 'b'
        ledger.consume('TE')
        self.assertEqual(ledger.read_summary('a'), {'total': 32, 'variables': 2, 'most_repeated': {'DAREA': 30}})
        self.assertEqual(ledger.read_summary('b'), {'total': 1, 'variables': 1, 'most_repeated': {}})
        self.assertEqual(ledger.read_summary('absent'), {})
        self.assertEqual(ledger.consumed('a'), ['DAREA', 'X'])

    def test_only_variables_read_often_are_listed_most_first_and_capped(self):
        ledger = MappingLedger()
        ledger.writer = 'w'
        for k in range(8):
            for _ in range(11 + k):
                ledger.consume(f'V{k}')
        listed = ledger.read_summary('w', listed=3)['most_repeated']
        self.assertEqual(list(listed), ['V7', 'V6', 'V5'])

    def test_the_journal_carries_reads_only_when_there_are_some(self):
        journal = ConversionJournal.__new__(ConversionJournal)
        journal.data = {'writers': [], 'written_ids': []}
        journal.record_writer('a', [], 0.1, reads={'total': 3, 'variables': 1, 'most_repeated': {}})
        journal.record_writer('b', [], 0.1, reads={})
        self.assertEqual(journal.data['writers'][0]['reads']['total'], 3)
        self.assertNotIn('reads', journal.data['writers'][1])


if __name__ == '__main__':
    unittest.main(verbosity=1)
