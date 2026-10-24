#!/usr/bin/python
"""Source-only checks of XML input preservation and per-model IDS metadata."""

from pathlib import Path
import importlib.util
import logging
import sys
import tempfile
import types
import xml.etree.ElementTree as ET
from unittest import mock
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from input_provenance import (InputInventory, assignments, fill_model_libraries,
                              parse_namelist, references_from_text, xml_safe_text)


class Array:
    def __init__(self):
        self.items = []

    def resize(self, size):
        self.items = [Library() for _ in range(size)]

    def __getitem__(self, index):
        return self.items[index]


class Library:
    def __init__(self):
        self.name = ''
        self.parameters = ''
        self.description = ''


class Code:
    def __init__(self):
        self.library = Array()


def assert_full_input_xml(inv, source_path):
    assert inv.records, source_path
    main = inv.records[0]
    assert main.status == 'embedded'
    assert main.sha256 and len(main.sha256) == 64
    # XML must preserve full source (not just a lossy namelist parse).
    out = ET.fromstring(inv.scoped_xml('TORIC'))
    files = out.find('input_files')
    assert files is not None
    source = files.find("file[@role='TRANSP_INPUT']")
    assert source is not None
    assert source.findtext('content') == main.text
    assert source.find('parsed_namelist') is not None
    assert files.attrib['complete'] == ('false' if inv.missing or inv.truncated else 'true')


def test_bundled():
    """All CDF examples get an input inventory, even without TR.DAT."""
    base = HERE / 'data'
    cdfs = sorted(base.rglob('*.CDF'))
    assert cdfs, 'No source examples found'
    expected_adjacent = {
        '147634B13': 'ech_inp.nml',
        '147634G14': 'genray_ech.in',
        '147634T13': 'toray.in',
        '147634T16': 'mmm.nml',
        '141711Y05': 'y05_pt.nml',
    }
    for cdf in cdfs:
        runid = cdf.stem
        trdat = cdf.with_name(f'{runid}TR.DAT')
        inv = InputInventory(trdat, runid)
        assert inv.records and inv.records[0].path == trdat.resolve(), runid
        xml = ET.fromstring(inv.scoped_xml(next(iter(inv.models()), 'TRANSP')))
        assert xml.find('input_files') is not None
        if trdat.exists():
            assert_full_input_xml(inv, trdat)
        else:
            assert inv.records[0].status == 'missing', runid
            assert xml.find('input_files').attrib['complete'] == 'false'
        # Every embedded file must retain its exact decoded source content,
        # a SHA256 digest and a searchable XML assignment index.
        for record in inv.records:
            if record.status != 'embedded':
                continue
            assert record.sha256 and len(record.sha256) == 64
            raw = record.path.read_bytes()
            try:
                decoded = raw.decode('utf-8')
            except UnicodeDecodeError:
                decoded = raw.decode('ISO-8859-1')
            assert record.text == xml_safe_text(decoded), record.path
            assert record.depth <= 12, record
            assert record.path.is_file(), record.path
        if runid in expected_adjacent:
            matched = [r for r in inv.records
                       if r.path.name == expected_adjacent[runid]]
            assert len(matched) == 1 and matched[0].status == 'embedded', runid
            assert matched[0].role in ('TEMPLATE', 'ADJACENT_UNREFERENCED'), runid
            if matched[0].role == 'ADJACENT_UNREFERENCED':
                assert 'association' in ET.fromstring(inv.scoped_xml(
                    matched[0].model)).find("input_files/file[@role='ADJACENT_UNREFERENCED']").attrib
        if runid == '38530P94':
            pt = next(r for r in inv.records if r.key == 'PT_TEMPLATE')
            assert pt.text == (trdat.parent / 'p94_pt.nml').read_text()
        code = Code()
        names = fill_model_libraries(code, inv, 'waves')
        assert names == inv.models('waves')
        for obj in code.library.items:
            scoped = ET.fromstring(obj.parameters)
            assert scoped.attrib['model'] == obj.name
        rel = cdf.relative_to(base)
        print(f'PASS {rel}: source inputs={len(inv.records)}, '
              f'models={inv.models()}, missing={[r.path.name for r in inv.missing]}')
    print(f'INPUT INVENTORY COVERAGE: {len(cdfs)}/{len(cdfs)} CDF examples')


def test_missing_relative_top_level_trdat_not_doubled():
    """Regression for -i examples/data/.../T16.CDF without TR.DAT.

    The first inventory entry must be relative to the CWD, while subsequent
    explicit references are relative to the referring namelist directory.
    """
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        directory = root / 'examples/data/d3d/147634T16'
        directory.mkdir(parents=True)
        (directory / 'mmm.nml').write_text('&mmm\n model=1\n/\n')
        relative = Path(__import__('os').path.relpath(
            directory / '147634T16TR.DAT', Path.cwd()))
        inv = InputInventory(relative, '147634T16')
        assert inv.records[0].status == 'missing'
        assert inv.records[0].path == (directory / '147634T16TR.DAT').resolve()
        assert '147634T16/examples/' not in str(inv.records[0].path)
        assert len(inv.records) == 2
        assert inv.records[1].path == (directory / 'mmm.nml').resolve()
        root_xml = ET.fromstring(inv.scoped_xml('MMM'))
        missing = root_xml.find(".//file[@status='missing']")
        assert missing is not None
        assert missing.get('resolved_path') == str(inv.records[0].path)
        # A present TR.DAT may point to a nested relative file; only those
        # child references are resolved from the TR.DAT parent directory.
        (directory / '147634T16TR.DAT').write_text(
            "MMM_NML = 'mmm.nml'\n")
        present = InputInventory(relative, '147634T16')
        assert present.records[0].status == 'embedded'
        assert any(item.path == (directory / 'mmm.nml').resolve()
                   and item.status == 'embedded' for item in present.records)


def test_nested_and_missing():
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        tr = d / 'sampleTR.DAT'
        toric = d / 'toric.nml'
        child = d / 'nested.nml'
        toric.write_text("&toric\n control='literal!notcomment'\n /\nINCLUDE 'nested.nml'\n")
        child.write_text("&physics\n a=5,6\n/\n")
        tr.write_text("! CQL3D_TEMPLATE='unused.nml'\n"
                      "TORIC_TEMPLATE = 'toric.nml'  ! used\n"
                      "CQL3D_NML = 'missing.nml'\n"
                      "TORIC_INFILE='toric.nml'\n"
                      "IGNORE = 'text!not-comment'\n")
        inv = InputInventory(tr, 'sample')
        assert len(inv.records) == 4, [(x.key, x.reference) for x in inv.records]
        assert set(inv.models('waves')) == {'TORIC', 'CQL3D'}
        nested = next(r for r in inv.records if r.reference == 'nested.nml')
        assert nested.model == 'TORIC'
        missing = inv.missing
        assert len(missing) == 1 and missing[0].reference == 'missing.nml'
        assert missing[0].status == 'missing'
        toric_xml = ET.fromstring(inv.scoped_xml('TORIC'))
        assert toric_xml.find(".//file[@key='INCLUDE']") is not None
        assert toric_xml.find(".//file[@key='CQL3D_NML']") is None
        cql_xml = ET.fromstring(inv.scoped_xml('CQL3D'))
        file = cql_xml.find(".//file[@key='CQL3D_NML']")
        assert file is not None and file.attrib['status'] == 'missing'
        assert cql_xml.find('input_files').attrib['complete'] == 'false'
        assert 'CQL3D_NML' in [key for _, key, _ in assignments(tr.read_text())]
        assert ('CQL3D_TEMPLATE', 'unused.nml') not in list(references_from_text(tr.read_text()))
        code = Code()
        names = fill_model_libraries(code, inv, 'waves')
        assert names == ['CQL3D', 'TORIC']
        assert [item.name for item in code.library.items] == names
        print('PASS nested include / missing CQL3D / explicit input-only model XML')


def test_actual_nc2imas_xml_methods_without_imas():
    """Exercise actual converter XML paths with only IMAS/NetCDF import shims."""
    source = HERE / 'data' / 'iter' / '38530P94TR.DAT'
    if not source.is_file():
        return
    netcdf_mod = types.ModuleType('netcdf_class')
    imas_mod = types.ModuleType('imas_class')
    imas_package = types.ModuleType('imas')
    netcdf_mod.netcdf_connection = type('FakeNetcdf', (), {})
    imas_mod.imas_connection = type('FakeImasConnection', (), {})
    with mock.patch.dict(sys.modules, {
        'netcdf_class': netcdf_mod, 'imas_class': imas_mod,
        'imas': imas_package,
    }):
        spec = importlib.util.spec_from_file_location(
            '_source_test_nc2imas', HERE.parent / 'nc2imas.py')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

    converter = module.nc2imas.__new__(module.nc2imas)
    converter.nml_name = str(source)
    converter.ncfile = str(source.parent / '38530P94.CDF')
    converter.device = 'ITER'
    converter.shot = 38530
    converter.run = 1694
    converter.occur = 0
    converter.dd_version = '4.1.1'
    from cocos_signs import resolve_signs
    converter.signs = resolve_signs('cw', 'cw', None, dd_major=4)      # Ip and Bphi clockwise, COCOS 17
    assert (converter.signs.psi_sign, converter.signs.phi_sign, converter.signs.q_sign) == (1.0, -1.0, 1.0)
    converter._code_parameters_xml_cache = None
    converter._input_inventory_cache = None
    converter.logger = logging.getLogger('converter-input-xml')
    converter.get_ncattr = lambda name: None
    converter.get_ncattr_text = lambda name, default='': (
        '38530P94' if name == 'Runid' else default)

    xml = converter.build_code_parameters_xml()
    root = ET.fromstring(xml)
    assert root.find("input_files/file[@key='TRANSP_NAMELIST']") is not None
    pt = root.find("input_files/file[@key='PT_TEMPLATE']")
    assert pt is not None and pt.attrib['status'] == 'embedded'
    assert pt.find('parsed_namelist') is not None
    assert converter.build_code_parameters_xml() == xml
    ids = types.SimpleNamespace(code=Code())
    assert converter.fill_external_model_inputs(ids, 'waves') == ['TORBEAM', 'TORIC']
    assert [c.name for c in ids.code.library.items] == ['TORBEAM', 'TORIC']
    print('PASS actual nc2imas XML/provenance/IDS code-library methods')


def test_workflow_writer_without_imas():
    """Exercise the actual workflow writer using an IMAS-shaped fake IDS."""
    tr = HERE / 'data' / 'iter' / '38530P94TR.DAT'
    if not tr.is_file():
        return
    inventory = InputInventory(tr, '38530P94')
    placeholder = types.ModuleType('nc2imas')
    placeholder.nc2imas = object
    with mock.patch.dict(sys.modules, {'nc2imas': placeholder}):
        spec = importlib.util.spec_from_file_location(
            '_test_imas_workflow', HERE.parent / 'imas_workflow.py')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

    class Fake:
        def __init__(self):
            self.logger = logging.getLogger('workflow-provenance-regression')
            self.timesc = np.array([0.0, 1.0])
            self.output = None

        def new_ids(self, name):
            assert name == 'workflow'
            return types.SimpleNamespace(
                ids_properties=types.SimpleNamespace(homogeneous_time=0),
                time_loop=types.SimpleNamespace(component=Array(), time_end=0.0),
            )

        def fill_properties(self, ids):
            return ids

        def fill_code_properties(self, ids):
            return ids

        def input_inventory(self):
            return inventory

        def build_code_parameters_xml(self):
            return inventory.scoped_xml('TORIC')

        def get_ncattr_text(self, name, default=''):
            return default

        def putdb_imas(self, ids):
            self.output = ids

    fake = Fake()
    writer = module.imas_workflow(fake)
    assert fake.output is None, 'constructing the writer must not write'
    writer.run()
    comps = fake.output.time_loop.component.items
    assert [c.name for c in comps] == ['TRANSP'] + inventory.models()
    for comp in comps[1:]:
        root = ET.fromstring(comp.parameters)
        assert root.attrib['model'] == comp.name
        assert root.find("input_files/file[@role='TRANSP_INPUT']") is not None
    assert fake.output.time_loop.time_end == 1.0
    print('PASS workflow writer: one TRANSP + separate model components')


def test_parser_preserves_source():
    raw = "&input\n X=3 ! comment\n MESSAGE='bang!in-string'\n Y=1,\n   2,3\n/\n"
    root = ET.Element('file')
    parse_namelist(root, raw)
    params = {p.attrib['name']: p.attrib['value'] for p in root.findall('.//parameter')}
    assert params['MESSAGE'] == "'bang!in-string'"
    assert params['Y'] == '1, 2,3'
    assert params['X'] == '3'
    print('PASS Fortran comment stripping, groups and continuation')


if __name__ == '__main__':
    test_parser_preserves_source()
    test_nested_and_missing()
    test_missing_relative_top_level_trdat_not_doubled()
    test_bundled()
    test_actual_nc2imas_xml_methods_without_imas()
    test_workflow_writer_without_imas()
    print('INPUT PROVENANCE MAPPING: PASS')
