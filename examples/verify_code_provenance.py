#!/usr/bin/python
"""Verify the provenance helpers (no IMAS, NetCDF or database needed).

* source-file records: unique by resolved path, capped, named ``role:file:<basename>#sha256=<digest>``
  without any absolute path, with the file's modification time; a missing file still gets a record,
  without a digest or a timestamp,
* attaching them to an IDS in either DD layout (``reference`` from DD 3.42, ``sources`` before),
* the ``ids_properties/comment`` text and the ``code.parameters`` XML structure.

    python examples/verify_code_provenance.py
"""

import hashlib
import os
from pathlib import Path
import re
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import code_provenance as cp
from cocos_signs import resolve_signs


class FileRecordTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, name, data=b'content'):
        path = self.dir / name
        path.write_bytes(data)
        return path

    def test_name_digest_and_timestamp(self):
        path = self.write('run.CDF', b'abc')
        os.utime(path, (1_700_000_000, 1_700_000_000))
        [(name, stamp)] = cp.source_file_records([('TRANSP CDF', str(path))])
        self.assertEqual(name, f'TRANSP CDF:file:run.CDF#sha256={hashlib.sha256(b"abc").hexdigest()}')
        self.assertEqual(stamp, '2023-11-14T22:13:20Z')
        self.assertNotIn(str(self.dir), name)                       # no machine-specific path

    def test_unique_by_resolved_path_and_empty_names_skipped(self):
        path = self.write('a.dat')
        link = self.dir / 'link.dat'
        link.symlink_to(path)
        records = cp.source_file_records([('A', str(path)), ('B', str(link)), ('C', ''), ('D', None),
                                          ('E', str(self.dir / '.' / 'a.dat'))])
        self.assertEqual(len(records), 1)
        self.assertTrue(records[0][0].startswith('A:file:a.dat#sha256='))

    def test_cap(self):
        paths = [self.write(f'f{i}.dat', bytes([i])) for i in range(15)]
        self.assertEqual(len(cp.source_file_records([('R', str(p)) for p in paths])), 10)
        self.assertEqual(len(cp.source_file_records([('R', str(p)) for p in paths], cap=3)), 3)

    def test_missing_file_keeps_a_record_without_digest_or_timestamp(self):
        [(name, stamp)] = cp.source_file_records([('TRANSP TR.DAT', str(self.dir / 'absent.DAT'))])
        self.assertEqual((name, stamp), ('TRANSP TR.DAT:file:absent.DAT', ''))
        self.assertEqual(cp.file_sha256(self.dir / 'absent.DAT'), '')
        self.assertEqual(cp.file_timestamp(self.dir / 'absent.DAT'), '')

    def test_timestamp_format(self):
        path = self.write('t.dat')
        self.assertRegex(cp.file_timestamp(path), r'^\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ$')


class Node:
    """Attribute bag standing in for an IDS structure; ``resize`` creates the array entries."""

    def __init__(self, **fields):
        self.__dict__.update(fields)

    def resize(self, n):
        self.items = [Node(**self.factory()) for _ in range(n)]

    def __getitem__(self, i):
        return self.items[i]

    def __iter__(self):
        return iter(self.items)

    def __len__(self):
        return len(self.items)


def reference_node():
    node = Node(path='x', factory=lambda: {})
    node.reference = Node(factory=lambda: {'name': '', 'timestamp': ''})
    return node


def sources_node():
    return Node(path='x', sources=[], factory=lambda: {})


def fake_ids(node):
    provenance = Node(factory=lambda: {})
    provenance.node = Node(factory=lambda: {})
    provenance.node.items = [node]
    provenance.node.resize = lambda n: None
    return Node(ids_properties=Node(provenance=provenance))


class AttachTests(unittest.TestCase):
    RECORDS = [('TRANSP CDF:file:a.CDF#sha256=aa', '2024-01-02T03:04:05Z'), ('TRANSP TR.DAT:file:a.DAT', '')]

    def test_reference_layout_from_dd_3_42(self):
        node = reference_node()
        cp.attach_source_provenance(fake_ids(node), self.RECORDS)
        self.assertEqual(node.path, '')
        self.assertEqual([(r.name, r.timestamp) for r in node.reference],
                         [('TRANSP CDF:file:a.CDF#sha256=aa', '2024-01-02T03:04:05Z'),
                          ('TRANSP TR.DAT:file:a.DAT', '')])        # an empty timestamp is not written

    def test_sources_layout_before_dd_3_42(self):
        node = sources_node()
        cp.attach_source_provenance(fake_ids(node), self.RECORDS)
        self.assertEqual(node.sources, [name for name, _ in self.RECORDS])

    def test_an_ids_without_provenance_raises_for_the_caller_to_record(self):
        with self.assertRaises(AttributeError):
            cp.attach_source_provenance(Node(), self.RECORDS)


class TextTests(unittest.TestCase):
    def test_comment(self):
        signs = resolve_signs('ccw', 'cw', None, dd_major=4)
        self.assertEqual(
            cp.properties_comment('38530P94', 'ITER', 2023, '01-Jan-2024', '02-Jan-2024', 'abc', 'def', '4.1.1', signs),
            'TRANSP results for run 38530P94 ITER year: 2023 CDF date: 01-Jan-2024 build: 02-Jan-2024 '
            'tglf: abc glf23: def DD=4.1.1; COCOS=17; NLJCCW=True; NLBCCW=False; Ip_sign=+1; Bphi_sign=-1; '
            'q_sign=-1; psi_sign=-1')


class Inventory:
    def __init__(self, records=(), missing=(), truncated=False):
        self.records, self.missing, self.truncated = list(records), list(missing), truncated


class XmlTests(unittest.TestCase):
    def build(self, **overrides):
        args = dict(dd_version='4.1.1', signs=resolve_signs('cw', 'cw', None, dd_major=4), device='ITER',
                    shot=38530, imas_run=1694, occurrence=0, source_cdf='/data/run.CDF',
                    converter_version='1.2.3', metadata={'Runid': b'38530P94 ', 'CDF_date': None, 'TRANSP_DOI': 'doi'},
                    inventory=Inventory())
        args.update(overrides)
        return ET.fromstring(cp.build_code_parameters_xml(**args))

    def test_structure_and_values(self):
        root = self.build()
        self.assertEqual((root.tag, root.get('format_version'), root.get('target_dd'), root.get('target_cocos')),
                         ('transp2imas_parameters', '1', '4.1.1', '17'))
        run = root.find('run')
        self.assertEqual((run.get('device'), run.get('shot'), run.get('imas_run'), run.get('occurrence'),
                          run.get('source_cdf')), ('ITER', '38530', '1694', '0', '/data/run.CDF'))
        self.assertEqual(root.find('converter').attrib, {'name': 'transp2imas', 'version': '1.2.3'})
        items = {i.get('name'): i.get('value') for i in root.findall('transp_metadata/item')}
        self.assertEqual(items, {'Runid': '38530P94', 'TRANSP_DOI': 'doi'})        # bytes decoded, None skipped
        conv = {i.get('name'): i.get('value') for i in root.findall('coordinate_conventions/item')}
        self.assertEqual(conv, {'NLJCCW': 'False', 'NLBCCW': 'False', 'ip_sign': '-1.0', 'bphi_sign': '-1.0',
                                'q_sign': '1.0', 'psi_sign': '1.0', 'phi_sign': '-1.0'})

    def test_input_files_completeness_flags(self):
        files = self.build().find('input_files')
        self.assertEqual((files.get('complete'), files.get('discovery_truncated')), ('true', 'false'))
        files = self.build(inventory=Inventory(missing=['x'])).find('input_files')
        self.assertEqual((files.get('complete'), files.get('discovery_truncated')), ('false', 'false'))
        files = self.build(inventory=Inventory(truncated=True)).find('input_files')
        self.assertEqual((files.get('complete'), files.get('discovery_truncated')), ('false', 'true'))


if __name__ == '__main__':
    unittest.main(verbosity=1)
