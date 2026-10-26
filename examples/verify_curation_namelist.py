#!/usr/bin/python
"""Tests for the TRANSP namelist adapter (curation/namelist.py) and the declared timeline built on it.

Three stages are tested separately: source collection (statements, spans, references), syntax parsing (f90nml behind
the adapter) and state resolution (``declared_timeline``).  Synthetic fixtures cover every syntax and update item;
the bundled TR.DAT files and PT templates are mandatory (a missing file fails, it is not skipped).  The old scalar
reader is kept below as a frozen copy to prove semantic preservation on the bundled inputs.

No IMAS, no network.

    python examples/verify_curation_namelist.py
"""

import hashlib
import json
import os
from pathlib import Path
import re
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import f90nml

import fortran_literals
from input_provenance import strip_fortran_comment
from curation import namelist
from curation.inputs import declared_timeline, invalidate, literal

DATA = ROOT / 'examples' / 'data'
TR_DATS = [DATA / 'nstx' / '141711Y05' / '141711Y05TR.DAT', DATA / 'cmod' / '89846TR.DAT',
           DATA / 'iter' / '38530P94TR.DAT']
PT_TEMPLATES = [DATA / 'nstx' / '141711Y05' / 'y05_pt.nml', DATA / 'd3d' / '147634T16' / 't16_pt.nml',
                DATA / 'iter' / 'p94_pt.nml']
OLD_KEY = re.compile(r'^([A-Za-z][A-Za-z0-9_%]*)\s*(?:\(([^()]*)\))?\s*=\s*(.*)$')
OLD_INTEGER = re.compile(r'\s*([+-]?\d+)\s*\Z')
OLD_SLICE = re.compile(r'\s*([+-]?\d+)\s*:\s*([+-]?\d+)\s*\Z')


def old_scalar_reader(path):
    """Frozen copy of the pre-f90nml scalar reader: list of per-state ``(start, end, values)``."""
    def normalize(name):
        return re.sub(r'\s+', '', name).upper()

    def old_invalidate(values, base, first=None, last=None):
        for key in list(values):
            name, _, rest = key.partition('(')
            if name != base:
                continue
            if rest and first is not None and not first <= int(rest.rstrip(')')) <= last:
                continue
            del values[key]

    states, values, start = [], {}, None
    for line in Path(path).read_bytes().decode('utf-8', errors='replace').splitlines():
        clean = strip_fortran_comment(line).strip()
        if clean.startswith('&'):
            clean = re.sub(r'^&\w+\s*', '', clean)
        if not clean or clean == '/':
            continue
        update = re.match(r'^~\s*update_time\s*=\s*(.*)$', clean, re.I)
        if update:
            time = literal(update[1])
            if isinstance(time, bool) or not isinstance(time, (int, float)) or (start is not None and time <= start):
                continue
            states.append((start, time, dict(values)))
            start = time
            continue
        match = OLD_KEY.match(clean)
        if not match:
            continue
        base, index, raw = match.groups()
        base = normalize(base)
        value = literal(raw)
        single = OLD_INTEGER.match(index) if index is not None else None
        name = f'{base}({int(single[1])})' if single else base
        if index is not None and not single:
            span = OLD_SLICE.match(index)
            old_invalidate(values, base, *((int(span[1]), int(span[2])) if span else (None, None)))
        elif value is None:
            values.pop(name, None)
            if not single:
                old_invalidate(values, base)
        else:
            if not single:
                old_invalidate(values, base)
            values[name] = value
    states.append((start, None, dict(values)))
    return states


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class Scratch(unittest.TestCase):
    """A temporary directory with helpers to write inputs and read the declared timeline."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.dir = Path(self._tmp.name)

    def write(self, name, body):
        path = self.dir / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(body)
        return path

    def timeline(self, body, **kwargs):
        timeline = declared_timeline(self.write('TR.DAT', body), **kwargs)
        json.dumps(timeline, allow_nan=False)              # every timeline must be JSON-serializable
        return timeline

    def final(self, body):
        return self.timeline(body)['states'][-1]['values']


class SyntaxTests(Scratch):
    """Stage 2: what the adapter reads from one statement."""

    def test_multidimensional_arrays_follow_fortran_column_major_order(self):
        values = self.final('A(1:2,1:2)=1,2,3,4\n')
        self.assertEqual(values, {'A(1,1)': 1, 'A(2,1)': 2, 'A(1,2)': 3, 'A(2,2)': 4})
        nested = f90nml.reads('&g\n a(1:2,1:2)=1,2,3,4\n/')['g']['a']      # f90nml agrees: nested[j][i]
        self.assertEqual([[nested[j][i] for j in range(2)] for i in range(2)], [[1, 3], [2, 4]])
        three = self.final('B(1:2,1:2,1:2)=1,2,3,4,5,6,7,8\n')
        self.assertEqual((three['B(2,1,1)'], three['B(1,2,1)'], three['B(1,1,2)'], three['B(2,2,2)']), (2, 3, 5, 8))

    def test_a_section_with_a_fixed_subscript(self):
        self.assertEqual(self.final('A(2,1:2)=5,6\n'), {'A(2,1)': 5, 'A(2,2)': 6})

    def test_repeat_notation(self):
        self.assertEqual(self.final('A=3*1.0\n'), {'A(1)': 1, 'A(2)': 1, 'A(3)': 1})
        self.assertEqual(self.final('A=2*3.,20.\n'), {'A(1)': 3, 'A(2)': 3, 'A(3)': 20})
        self.assertEqual(self.final("S=2*'x'\n"), {'S(1)': 'x', 'S(2)': 'x'})

    def test_slices_and_element_continuation(self):
        self.assertEqual(self.final('A(2:3)=7,8\n'), {'A(2)': 7, 'A(3)': 8})
        self.assertEqual(self.final('A(3)=1,2\n'), {'A(3)': 1, 'A(4)': 2})          # NAME(k)=v1,v2 fills k, k+1
        self.assertEqual(self.final('A(1:3)=1,2\n'), {'A(1)': 1, 'A(2)': 2})        # fewer values leave the rest

    def test_sparse_one_based_arrays_keep_their_actual_indices(self):
        timeline = self.timeline('A(3)=1\n')
        self.assertEqual(timeline['states'][0]['values'], {'A(3)': 1})
        self.assertNotIn('A(1)', timeline['states'][0]['values'])
        record = timeline['assignments'][0]
        self.assertEqual((record['indices'], record['keys']), (['3'], ['A(3)']))
        self.assertEqual(self.final('A(0)=1\nA(-2:-1)=5,6\n'), {'A(0)': 1, 'A(-2)': 5, 'A(-1)': 6})

    def test_whole_array_without_subscript_assumes_lower_bound_one_and_says_so(self):
        timeline = self.timeline('A=1,2\nA=3\n')
        self.assertEqual(timeline['states'][0]['values'], {'A': 3})                 # a bare scalar replaces the elements
        kinds = [u['kind'] for u in timeline['unresolved']]
        self.assertEqual(kinds, ['assumed_lower_bound'])
        self.assertEqual(self.timeline('A=1,2\n')['unresolved'][0]['name'], 'A')

    def test_null_zero_omitted_and_invalid_elements_are_distinct(self):
        setup = 'A(1)=1\nA(2)=2\nA(3)=3\n'
        self.assertEqual(self.final(setup + 'A(1:3)=9,,8\n'), {'A(1)': 9, 'A(2)': 2, 'A(3)': 8})      # null: unchanged
        self.assertEqual(self.final(setup + 'A(1:3)=9,0,8\n'), {'A(1)': 9, 'A(2)': 0, 'A(3)': 8})     # explicit zero
        self.assertEqual(self.final(setup + 'A(1)=9\n'), {'A(1)': 9, 'A(2)': 2, 'A(3)': 3})           # omitted: unchanged
        timeline = self.timeline(setup + 'A(1:3)=9,1+2,8\n')                                         # invalid: unknown
        self.assertEqual(timeline['states'][0]['values'], {'A(1)': 9, 'A(3)': 8})
        self.assertEqual(timeline['gaps'][0]['invalidated'], ['A(2)'])
        null = self.timeline(setup + 'A(1:3)=9,,8\n')['assignments'][-1]
        self.assertEqual(null['value'], [9, {'null': True}, 8])
        self.assertEqual(self.timeline(setup + 'A(1:3)=9,1+2,8\n')['assignments'][-1]['value'][1], {'unknown': '1+2'})
        self.assertEqual(self.final('A(1:2)=,\n'), {})                                                # null on nothing
        self.assertEqual(self.final('A=,,3\n'), {'A(3)': 3})                                          # leading nulls

    def test_a_null_is_not_zero_in_a_fresh_state(self):
        values = self.final('A(1:2)=,5\n')
        self.assertNotIn('A(1)', values)
        self.assertEqual(values, {'A(2)': 5})

    def test_derived_type_components(self):
        values = self.final('pt%b%c=1\npt%b%d(2)=3\nPT%B%E=.TRUE.\n')
        self.assertEqual(values, {'PT%B%C': 1, 'PT%B%D(2)': 3, 'PT%B%E': True})
        self.assertEqual(self.timeline('a(1)%b=1\n')['gaps'][0]['kind'], 'subscripted_component')

    def test_exponents_quoted_bang_and_doubled_quotes(self):
        values = self.final("A=1.5D+2\nB=2.5e-1\nC=-.5d0\nS='it''s ! not a comment' ! real comment\nD=\"say \"\"hi\"\"\"\n")
        self.assertEqual(values, {'A': 150.0, 'B': 0.25, 'C': -0.5, 'S': "it's ! not a comment", 'D': 'say "hi"'})
        self.assertEqual(self.final('A=1.0D-3,2d-3\n'), {'A(1)': .001, 'A(2)': .002})

    def test_dotted_logicals_and_expressions(self):
        values = self.final('A=.T\nB=.F\nC=.TRUE.\nD=T\n')
        self.assertEqual(values, {'A': True, 'B': False, 'C': True, 'D': True})
        timeline = self.timeline('E=1+2\nF=$unsupported\nG=1 +\n')
        self.assertEqual(timeline['states'][0]['values'], {})                       # f90nml would read 1+2 as 100.0
        self.assertEqual([g['kind'] for g in timeline['gaps']].count('unsupported_literal'), 2)
        self.assertEqual(f90nml.reads('&g\n a=1+2\n/')['g']['a'], 100.0)             # why f90nml alone is not trusted
        for text in ('A=nan', 'A=1e999', 'A=(1.0,2.0)', 'A=0x10', "A='open"):
            self.assertEqual(self.final(text + '\n'), {}, text)

    def test_several_assignments_per_line_in_a_namelist_but_not_in_a_tr_dat(self):
        context = namelist.Context(f90nml)
        document = namelist.parse_document(context, '&g a=1, b=2 c=3,4 /\n', 'namelist')
        found = [(i['parsed']['name'], [e['key'] for e in i['parsed']['elements']])
                 for i in document['items'] if i['type'] == 'assignment']
        self.assertEqual(found, [('A', ['A']), ('B', ['B']), ('C', ['C(1)', 'C(2)'])])
        timeline = self.timeline('X(1)=5\nX(2)=1, Y=2\n')
        self.assertEqual(timeline['states'][0]['values'], {'X(1)': 5})              # TR.DAT: the line is not interpreted
        self.assertEqual(timeline['gaps'][0]['kind'], 'multiple_assignments_per_line')
        self.assertEqual(timeline['gaps'][0]['names'], ['X(2)', 'Y'])

    def test_continuation_lines(self):
        timeline = self.timeline('A=1,2,\n  3, 4\nB=5\n')
        self.assertEqual(timeline['states'][0]['values'], {'A(1)': 1, 'A(2)': 2, 'A(3)': 3, 'A(4)': 4, 'B': 5})
        self.assertEqual((timeline['assignments'][0]['line'], timeline['assignments'][0]['end_line']), (1, 2))
        gap = self.timeline('A=1\n  3, 4\n')['gaps'][0]                            # no trailing comma: not a continuation
        self.assertEqual(gap['kind'], 'unsupported_statement')

    def test_too_many_values_for_a_section_or_a_huge_repeat_do_not_resolve(self):
        timeline = self.timeline('A(1)=7\nA(1:2)=1,2,3\n')
        self.assertEqual(timeline['states'][0]['values'], {})
        self.assertEqual(timeline['gaps'][0]['kind'], 'too_many_values')
        timeline = self.timeline('A=1000000000*0.0\n')
        self.assertEqual(timeline['states'][0]['values'], {})
        self.assertEqual(timeline['gaps'][0]['kind'], 'syntax_error')

    def test_multidimensional_element_continuation_is_unresolved(self):
        timeline = self.timeline('P(1,1)=1\nP(2,1)=2\nP(1,1)=0.,180.\nQ(1,1)=9\n')
        self.assertEqual(timeline['states'][0]['values'], {'Q(1,1)': 9})
        self.assertEqual(timeline['gaps'][0]['kind'], 'unbounded_continuation')
        self.assertEqual(timeline['gaps'][0]['invalidated'], ['P(1,1)', 'P(2,1)'])

    def test_parser_configuration_and_version_are_recorded(self):
        parser = self.timeline('A=1\n')['parser']
        self.assertEqual((parser['kind'], parser['version']), ('f90nml', f90nml.__version__))
        options = parser['config']['parser_options']
        self.assertEqual((options['comment_tokens'], options['default_start_index'], options['row_major'],
                          options['sparse_arrays']), ('!', 1, False, False))
        for key in ('comments', 'start_index', 'array_order', 'dtype', 'repeat_and_null', 'scalars'):
            self.assertIn(key, parser['config'])
        self.assertEqual(namelist.parser_record()['version'], f90nml.__version__)

    def test_assignment_records_have_order_block_span_and_raw_text(self):
        timeline = self.timeline('A=1\n  B (2) = 7,8 ! c\n~update_time=1\nC=.T\n')
        records = timeline['assignments']
        self.assertEqual([r['order'] for r in records], [1, 2, 3])
        self.assertEqual([r['block'] for r in records], [0, 0, 1])
        second = records[1]
        self.assertEqual((second['name'], second['indices'], second['line'], second['column'], second['raw']),
                         ('B', ['2'], 2, 3, '7,8'))
        self.assertEqual(second['file'], timeline['source']['path'])
        self.assertEqual(records[2]['value'], True)

    def test_assignment_records_are_capped_with_a_flag(self):
        timeline = self.timeline('A=1\nB=2\nC=3\n', max_assignments=2)
        self.assertEqual((len(timeline['assignments']), timeline['assignments_truncated']), (2, True))
        self.assertEqual(len(self.timeline('A=1\n')['assignments']), 1)

    def test_synthetic_namelist_maps_back_to_source_lines(self):
        context = namelist.Context(f90nml)
        text = '! head\nA=1\n\n  B=2,\n 3\n~update_time=1\nC=4\n'
        document = namelist.parse_document(context, text, 'transp_dat')
        synthetic = document['synthetic']
        self.assertEqual(synthetic['text'].splitlines(), ['&tr_block_0', 'A=1', 'B=2, 3', '/', '&tr_block_1', 'C=4', '/'])
        self.assertEqual(synthetic['line_map'][2], {'line': 2, 'end_line': 2, 'column': 1})
        self.assertEqual(synthetic['line_map'][3], {'line': 4, 'end_line': 5, 'column': 3})
        self.assertIsNone(synthetic['line_map'][1])
        self.assertEqual(document['syntax_check']['synthetic_namelist'], 'ok')
        self.assertEqual(list(f90nml.reads(synthetic['text'])), ['tr_block_0', 'tr_block_1'])


class UpdateTests(Scratch):
    """Stage 3: ordered updates over the update blocks."""

    def test_explicit_declarations_persist_and_partial_updates_keep_untouched_elements(self):
        timeline = self.timeline('A=1,2,3\nB=5\n~update_time=1\nA(2)=9\n~update_time=2\nA(3:4)=7,8\n')
        states = [s['values'] for s in timeline['states']]
        self.assertEqual(states[0], {'A(1)': 1, 'A(2)': 2, 'A(3)': 3, 'B': 5})
        self.assertEqual(states[1], {'A(1)': 1, 'A(2)': 9, 'A(3)': 3, 'B': 5})
        self.assertEqual(states[2], {'A(1)': 1, 'A(2)': 9, 'A(3)': 7, 'A(4)': 8, 'B': 5})
        self.assertEqual([(s['start_time'], s['end_time']) for s in timeline['states']],
                         [(None, 1), (1, 2), (2, None)])
        self.assertEqual(timeline['states'][2]['source_lines']['A(2)'], 4)
        self.assertEqual(timeline['states'][2]['source_lines']['A(3)'], 6)

    def test_a_malformed_update_does_not_reuse_superseded_values(self):
        values = self.final('A=1,2,3\n~update_time=1\nA(2:3)=1+2\n')
        self.assertEqual(values, {'A(1)': 1})
        self.assertEqual(self.final('A=1,2,3\nA(1:3:2)=5,6\n'), {})                      # strided: not supported
        self.assertEqual(self.final('A=1,2,3\nA(:)=5,6\n'), {})                          # bound-less range
        self.assertEqual(self.final('A=1,2,3\nA(N)=5\n'), {})                            # expression subscript
        self.assertEqual(self.final('A=1,2,3\nB=1\nA(2)=garbage\n'), {'A(1)': 1, 'A(3)': 3, 'B': 1})
        timeline = self.timeline('A=1,2,3\nA(:)=5,6\n')
        self.assertEqual(sorted(timeline['gaps'][0]['invalidated']), ['A(1)', 'A(2)', 'A(3)'])
        self.assertEqual({u['kind'] for u in timeline['unresolved']}, {'assumed_lower_bound', 'bounds_less_range'})

    def test_a_scalar_update_after_elements_invalidates_the_elements(self):
        timeline = self.timeline('A(1)=1\nA(2)=2\nA=7\n')
        self.assertEqual(timeline['states'][0]['values'], {'A': 7})
        self.assertEqual(timeline['assignments'][-1]['invalidated'], ['A(1)', 'A(2)'])

    def test_update_time_rules_are_unchanged(self):
        timeline = self.timeline('A=1\n~update_time=2\n~update_time=1\n')
        self.assertFalse(timeline['temporal_valid'])
        self.assertEqual(timeline['gaps'][0]['line'], 3)
        self.assertFalse(self.timeline('~update_time=bad\n')['temporal_valid'])
        self.assertTrue(self.timeline('~update_time=1\n~update_time=2\n')['temporal_valid'])

    def test_group_markers_are_structure_not_gaps(self):
        timeline = self.timeline('$TRDATA\nA=1\n$END\n&INNML B=2\n/\n&nml\nC=3\n&end\n')
        self.assertEqual(timeline['states'][0]['values'], {'A': 1, 'B': 2, 'C': 3})
        self.assertEqual(timeline['gaps'], [])
        self.assertEqual(len(timeline['files']['entries'][0]['groups']), 6)

    def test_nothing_is_defaulted(self):
        self.assertEqual(self.final('A(2)=1\n'), {'A(2)': 1})
        self.assertEqual(self.timeline('')['states'][0]['values'], {})
        self.assertEqual(declared_timeline()['states'], [])

    def test_the_public_helpers_keep_their_contract(self):
        self.assertEqual(literal('1.0D-3'), .001)
        self.assertIsNone(literal('.T'))                                                  # fortran_literals contract
        self.assertTrue(namelist.scalar_value('.T'))
        values, locations = {'A(1)': 1, 'A(2)': 2, 'A(1,2)': 3, 'B(1)': 4}, {'A(1)': 1}
        self.assertEqual(invalidate(values, locations, 'A', 2, 3), ['A(1,2)', 'A(2)'])
        self.assertEqual(values, {'A(1)': 1, 'B(1)': 4})
        self.assertEqual(invalidate(values, locations, 'A'), ['A(1)'])

    def test_source_files_are_never_modified(self):
        path = self.write('TR.DAT', "PT_TEMPLATE='pt.nml'\nA=1,2\n")
        self.write('pt.nml', '&g\n a=1\n/\n')
        before = {p.name: (sha(p), p.stat().st_mtime_ns) for p in self.dir.iterdir()}
        declared_timeline(path)
        after = {p.name: (sha(p), p.stat().st_mtime_ns) for p in self.dir.iterdir()}
        self.assertEqual(before, after)
        for bundled in TR_DATS + PT_TEMPLATES:
            self.assertTrue(bundled.is_file(), f'missing bundled input {bundled}')
            before = (sha(bundled), bundled.stat().st_mtime_ns)
            declared_timeline(bundled)
            self.assertEqual(before, (sha(bundled), bundled.stat().st_mtime_ns))

    def test_the_size_limit_is_kept(self):
        path = self.write('big.DAT', '')
        with path.open('wb') as stream:
            stream.truncate(10 * 1024 * 1024 + 1)
        with self.assertRaisesRegex(ValueError, '10 MiB'):
            declared_timeline(path)

    def test_every_gap_and_unresolved_record_is_json_serializable(self):
        body = ("include 'x.nml'\n#include y.nml\nPT_TEMPLATE='nope.nml'\nA=1+2\nB(:)=1\nC(1:3:2)=1\nD(1)%e=1\n"
                'F(1,1)=1,2\nG=1,2\nG(1:2)=1,2,3\nH=1, I=2\nstray words\nJ=1000000000*0\n~update_time=bad\n')
        timeline = self.timeline(body)
        for record in timeline['gaps'] + timeline['unresolved']:
            json.dumps(record, allow_nan=False)
        kinds = {g['kind'] for g in timeline['gaps']}
        self.assertTrue({'include_not_executed', 'missing', 'unsupported_literal', 'bounds_less_range', 'strided_section',
                         'subscripted_component', 'unbounded_continuation', 'too_many_values',
                         'multiple_assignments_per_line', 'unsupported_statement', 'syntax_error',
                         'update_time'} <= kinds, kinds)
        self.assertFalse(timeline['temporal_valid'])


class ReferenceTests(Scratch):
    """Stage 1: file references."""

    def entries(self, timeline):
        return {(Path(e['path']).name, e['status']) for e in timeline['files']['entries']}

    def test_references_resolve_relative_to_the_including_document_and_record_sha256(self):
        self.write('sub/x.nml', "&g\n a=1\n/\nINCLUDE_FILE='y.nml'\n")        # not a standard namelist line; recorded only
        self.write('sub/y.nml', '&h\n b=2\n/\n')
        self.write('y.nml', '&wrong\n b=3\n/\n')                                 # the including-document directory wins
        timeline = self.timeline("PT_TEMPLATE='sub/x.nml'\n")
        files = timeline['files']
        self.assertEqual(files['allowed_roots'], [str(self.dir.resolve())])
        self.assertEqual(files['max_depth'], 4)
        top, first = files['entries'][0], files['entries'][1]
        self.assertEqual((top['role'], top['sha256']), ('top_level', sha(self.dir / 'TR.DAT')))
        self.assertEqual((first['path'], first['sha256'], first['status']),
                         (str((self.dir / 'sub' / 'x.nml').resolve()), sha(self.dir / 'sub' / 'x.nml'), 'parsed'))
        self.assertEqual((first['key'], first['line'], first['referenced_by'], first['depth']),
                         ('PT_TEMPLATE', 1, str((self.dir / 'TR.DAT').resolve()), 1))
        self.assertEqual([Path(d['path']) for d in timeline['documents']],
                         [(self.dir / 'sub' / 'x.nml').resolve(), (self.dir / 'sub' / 'y.nml').resolve()])
        self.assertEqual(files['entries'][2]['referenced_by'], str((self.dir / 'sub' / 'x.nml').resolve()))
        self.write('sub/z.nml', '&h\n b=2\n/\n')
        self.write('sub/w.nml', "&g\n pt_next = 'z.nml'\n/\n")
        timeline = self.timeline("PT_TEMPLATE='sub/w.nml'\n")
        self.assertEqual(self.entries(timeline) - {('TR.DAT', 'parsed')}, {('w.nml', 'parsed'), ('z.nml', 'parsed')})
        self.assertEqual({d['path'] for d in timeline['documents']},
                         {str((self.dir / 'sub' / 'w.nml').resolve()), str((self.dir / 'sub' / 'z.nml').resolve())})

    def test_a_referenced_template_is_its_own_document_and_is_not_interpreted(self):
        self.write('pt.nml', '1\n&pt_controls\n pt_update_time = -1.\n pt_a%b = 6*0.0\n/\n&pt_controls\n x=1\n/\n')
        timeline = self.timeline("PT_TEMPLATE='pt.nml'\nA=1\n")
        self.assertEqual(timeline['states'][0]['values'], {'PT_TEMPLATE': 'pt.nml', 'A': 1})   # no pt_ keys in the state
        document = timeline['documents'][0]
        self.assertEqual(document['parser']['kind'], 'f90nml')
        self.assertEqual(document['syntax_check']['template_count'], {'declared': 1, 'groups_found': 2, 'matches': False})
        self.assertEqual([r['name'] for r in document['assignments']], ['PT_UPDATE_TIME', 'PT_A%B', 'X'])
        self.assertEqual(document['assignments'][0]['value'], -1.0)
        self.assertEqual(document['assignments'][1]['keys'], [f'PT_A%B({i})' for i in range(1, 7)])
        self.assertIn('not interpreted', document['semantics'])

    def test_a_missing_referenced_file_is_recorded(self):
        timeline = self.timeline("PT_TEMPLATE='gone.nml'\nA=1\n")
        entry = timeline['files']['entries'][1]
        self.assertEqual((entry['status'], entry['sha256']), ('missing', None))
        self.assertEqual(timeline['gaps'][0]['kind'], 'missing')
        self.assertEqual(timeline['unresolved'][0]['kind'], 'missing')
        self.assertEqual(timeline['states'][0]['values']['A'], 1)

    def test_a_reference_cycle_is_detected(self):
        self.write('a.nml', "&g\n pt_next='b.nml'\n/\n")
        self.write('b.nml', "&g\n pt_next='a.nml'\n/\n")
        timeline = self.timeline("PT_TEMPLATE='a.nml'\n")
        self.assertEqual([e['status'] for e in timeline['files']['entries']], ['parsed', 'parsed', 'parsed', 'cycle'])
        self.assertEqual(len(timeline['documents']), 2)
        timeline = self.timeline("PT_TEMPLATE='TR.DAT'\n")                       # the top-level file referencing itself
        self.assertEqual(timeline['files']['entries'][1]['status'], 'cycle')

    def test_the_depth_limit(self):
        for index in range(7):
            self.write(f'n{index}.nml', f"&g\n pt_next='n{index + 1}.nml'\n/\n")
        timeline = self.timeline("PT_TEMPLATE='n0.nml'\n")
        self.assertEqual([e['status'] for e in timeline['files']['entries']],
                         ['parsed', 'parsed', 'parsed', 'parsed', 'parsed', 'depth_limit'])
        self.assertEqual(timeline['files']['max_depth'], 4)
        timeline = self.timeline("PT_TEMPLATE='n0.nml'\n", max_depth=1)
        self.assertEqual([e['status'] for e in timeline['files']['entries']], ['parsed', 'parsed', 'depth_limit'])
        self.assertEqual({g['kind'] for g in timeline['gaps']}, {'depth_limit'})

    def test_paths_outside_the_allowed_roots_are_not_read(self):
        outside = self.write('outside/o.nml', '&g\n a=1\n/\n')
        inner = self.dir / 'run'
        inner.mkdir()
        top = inner / 'TR.DAT'
        top.write_text("PT_TEMPLATE='../outside/o.nml'\nX_NML='" + str(outside) + "'\n")
        timeline = declared_timeline(top)
        self.assertEqual([e['status'] for e in timeline['files']['entries'][1:]], ['outside_allowed_roots'] * 2)
        self.assertEqual(timeline['documents'], [])
        self.assertIsNone(timeline['files']['entries'][1]['sha256'])
        widened = declared_timeline(top, allowed_roots=[inner, self.dir / 'outside'])
        self.assertEqual([e['status'] for e in widened['files']['entries'][1:]], ['parsed', 'parsed'])
        self.assertEqual(widened['files']['allowed_roots'], [str(inner.resolve()), str((self.dir / 'outside').resolve())])
        link = inner / 'link.nml'
        os.symlink(outside, link)                                                # a symlink cannot leave the roots
        top.write_text("PT_TEMPLATE='link.nml'\n")
        self.assertEqual(declared_timeline(top)['files']['entries'][1]['status'], 'outside_allowed_roots')

    def test_include_directives_and_expressions_are_recorded_never_executed(self):
        marker = self.dir / 'marker'
        self.write('other.nml', "&g\n a=1\n/\n")
        timeline = self.timeline(f"include 'other.nml'\n#include other.nml\nPT_TEMPLATE='`touch {marker}`.nml'\n"
                                 "PT_TEMPLATE='a*.nml'\nY=1+2\n")
        statuses = [(e['key'], e['status']) for e in timeline['files']['entries'][1:]]
        self.assertEqual(statuses[:2], [('INCLUDE', 'include_not_executed')] * 2)
        self.assertEqual(statuses[2], ('PT_TEMPLATE', 'missing'))
        self.assertEqual(len(statuses), 3)                                       # a glob is not a reference
        self.assertEqual(timeline['documents'], [])                              # other.nml was not parsed
        self.assertFalse(marker.exists())
        self.assertEqual([u['kind'] for u in timeline['unresolved']].count('include_not_executed'), 2)
        self.assertEqual([g['kind'] for g in timeline['gaps']].count('include_not_executed'), 2)
        self.assertNotIn('Y', timeline['states'][0]['values'])

    def test_a_recognised_reference_that_is_not_a_namelist_is_not_followed(self):
        self.write('x.inp', 'text')
        timeline = self.timeline("TORIC_INPUT_FILE='x.inp'\n")
        self.assertEqual(timeline['files']['entries'][1]['status'], 'recognised_not_followed')
        self.assertEqual(timeline['documents'], [])


class FallbackTests(Scratch):
    def test_without_f90nml_the_limited_reader_is_reported_not_disguised(self):
        body = 'A=1\nB=1,2\nC=.T\nPT_TEMPLATE=\'x.nml\'\n'
        with mock.patch.dict(sys.modules, {'f90nml': None}):                       # ``import f90nml`` raises ImportError
            self.assertIsNone(namelist.parser_record())
            timeline = self.timeline(body)
        parser = timeline['parser']
        self.assertEqual((parser['kind'], parser['version']), ('scalar_fallback', None))
        self.assertTrue(any('arrays' in text for text in parser['limitations']))
        self.assertTrue(any('referenced files' in text for text in parser['limitations']))
        self.assertEqual(timeline['states'][0]['values'], {'A': 1, 'PT_TEMPLATE': 'x.nml'})   # the old reader: no arrays, no .T
        self.assertEqual(timeline['gaps'][0]['kind'], 'parser_fallback')
        self.assertEqual(timeline['unresolved'][0]['kind'], 'parser_fallback')
        self.assertEqual((timeline['assignments'], timeline['documents'], timeline['files']), ([], [], None))
        self.assertEqual(self.timeline(body)['parser']['kind'], 'f90nml')           # and the mock did not leak


class BundledInputTests(unittest.TestCase):
    """The bundled inputs are mandatory; a missing file fails."""

    def test_the_bundled_inputs_exist(self):
        for path in TR_DATS + PT_TEMPLATES:
            self.assertTrue(path.is_file(), f'missing bundled input {path}')

    def test_every_bundled_timeline_is_json_serializable_and_hashes_its_files(self):
        for path in TR_DATS:
            timeline = declared_timeline(path)
            json.dumps(timeline, allow_nan=False)
            self.assertEqual(timeline['source']['sha256'], sha(path))
            self.assertEqual(timeline['parser']['kind'], 'f90nml')
            self.assertTrue(timeline['assignments'])
            self.assertFalse(timeline['assignments_truncated'])
            for entry in timeline['files']['entries'][1:]:
                self.assertEqual(entry['status'], 'parsed', (path, entry))
                self.assertEqual(entry['sha256'], sha(entry['path']))
            for record in timeline['gaps'] + timeline['unresolved']:
                json.dumps(record, allow_nan=False)

    def test_each_bundled_tr_dat_references_its_pt_template(self):
        expected = {'141711Y05TR.DAT': 'y05_pt.nml', '38530P94TR.DAT': 'p94_pt.nml'}
        for path in TR_DATS:
            names = [Path(e['path']).name for e in declared_timeline(path)['files']['entries'][1:]]
            self.assertEqual(names, [expected[path.name]] if path.name in expected else [], path)
        t16 = declared_timeline(DATA / 'd3d' / '147634T16' / '147634T16TR.DAT')
        self.assertIn('t16_pt.nml', [Path(e['path']).name for e in t16['files']['entries']])

    def test_bundled_pt_templates_parse_as_namelist_documents(self):
        context = namelist.Context(f90nml)
        for path in PT_TEMPLATES:
            text = path.read_text()
            document = namelist.parse_document(context, text, 'namelist')
            check = document['syntax_check']
            self.assertEqual((check['document'], check['synthetic_namelist']), ('ok', 'ok'), path)
            self.assertEqual(check['document_groups'], ['pt_controls'])
            self.assertEqual(check['template_count'], {'declared': 1, 'groups_found': 1, 'matches': True})
            assignments = [i for i in document['items'] if i['type'] == 'assignment']
            self.assertTrue(all(i['parsed']['status'] == 'resolved' for i in assignments), path)
            self.assertFalse([i for i in document['items'] if i['type'] == 'unsupported'], path)
            self.assertEqual(f90nml.read(str(path))['pt_controls']['pt_update_time'], -1.0)
            names = {i['parsed']['name'] for i in assignments}
            self.assertIn('PT_UPDATE_TIME', names)
        y05 = declared_timeline(DATA / 'nstx' / '141711Y05' / '141711Y05TR.DAT')['documents'][0]
        by_name = {r['name']: r for r in y05['assignments']}
        self.assertIs(by_name['PT_CONFINEMENT%CDBM%ACTIVE']['value'], True)
        self.assertEqual(by_name['PT_CONFINEMENT%CDBM%XANOM']['value'], [0.0] * 6)
        self.assertEqual(by_name['PT_CONFINEMENT%CDBM%XANOM']['keys'][:2], ['PT_CONFINEMENT%CDBM%XANOM(1)',
                                                                          'PT_CONFINEMENT%CDBM%XANOM(2)'])
        self.assertEqual(by_name['PT_RESIDUAL%RES_TI']['value'], 5e-3)                # 5D-3
        self.assertEqual(by_name['PT_UPDATE_TIME']['value'], -1.0)
        self.assertEqual([r['order'] for r in y05['assignments']], list(range(1, len(y05['assignments']) + 1)))

    def test_d_exponents_and_indexed_keys_of_the_real_inputs(self):
        nstx = declared_timeline(TR_DATS[0])['states'][-1]['values']
        self.assertEqual(nstx['TBONA(1)'], 9.20500010252e-2)
        self.assertEqual(nstx['PT_TEMPLATE'], 'y05_pt.nml')
        self.assertEqual([nstx[f'XYBAPA({i})'] for i in (1, 12)], [0.0, 0.0])
        self.assertEqual(nstx['NLMDIF'], False)                                         # .F, previously a gap

    def test_semantic_preservation_against_the_old_scalar_reader(self):
        for path in TR_DATS:
            old = old_scalar_reader(path)
            new = declared_timeline(path)['states']
            self.assertEqual(len(old), len(new), path)
            added, differences = 0, []
            for (start, end, old_values), state in zip(old, new):
                self.assertEqual((start, end), (state['start_time'], state['end_time']), path)
                for key, value in old_values.items():
                    if key not in state['values'] or state['values'][key] != value or \
                            type(state['values'][key]) is not type(value):
                        differences.append((key, value, state['values'].get(key)))
                added += len(set(state['values']) - set(old_values))
            print(f'{path.name}: {sum(len(v) for _s, _e, v in old)} old keys, {len(differences)} differences, '
                  f'{added} keys resolved in addition', file=sys.stderr)
            # the one known difference: the old reader took ``X='a','b'`` as the single string "a','b" (a bug); the new
            # parser reads the list and the old key is correctly gone, replaced by X(1), X(2), ...
            misread = [d for d in differences if isinstance(d[1], str) and "','" in d[1]]
            self.assertEqual([d for d in differences if d not in misread], [], f'{path}: old keys changed or vanished')
            self.assertTrue(all(d[0] in ('FB_NAMES', 'PFC_NAMES') and d[2] is None for d in misread), misread)
            for _key, old_value, _new in misread:
                self.assertEqual(old_value.count("','") + 1, len(old_value.split("','")))
            if misread:
                values = new[-1]['values']
                self.assertEqual([values[f'FB_NAMES({i})'] for i in (1, 2, 3)], ['iver', 'ihor', 'ioh1'])
                self.assertEqual(values['PFC_NAMES(11)'], 'pf6')
                self.assertEqual(len(misread), 2 * len(new))                                # both names, every block
            self.assertGreater(added, 0, path)                                          # arrays, repeats and .T are new

    def test_every_gap_of_the_old_reader_is_resolved_or_explained(self):
        for path in TR_DATS:
            timeline = declared_timeline(path)
            for gap in timeline['gaps']:
                self.assertTrue(gap.get('kind') and gap.get('reason'), (path, gap))


if __name__ == '__main__':
    unittest.main(verbosity=1)
