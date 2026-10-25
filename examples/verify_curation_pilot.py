#!/usr/bin/python
"""Run the curation pilot on the bundled real cases and check what must hold, strictly (no skipping).

The two cases are NSTX 141711Y05 and DIII-D 147634T16 (one TRANSP run each).  The facts used here are exactly
what the files themselves state plus the file hashes; they are the *mechanics* check of the pilot path and
are not campaign evidence: no origin, uncertainty or independence is claimed for any measurement, so those
findings must stay ``not_assessed``.  A real pilot supplies its own facts file (see
``examples/curation/pilot-template.json``).

    python examples/verify_curation_pilot.py
"""

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
CASES = (('nstx/141711Y05', '141711Y05', 'NSTX', 141711), ('d3d/147634T16', '147634T16', 'D3D', 147634))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class PilotOnRealCasesTests(unittest.TestCase):
    def curate(self, folder, name, device, shot, tmp, mutate=None):
        directory = ROOT / 'examples' / 'data' / folder
        cdf, trdat = directory / (name + '.CDF'), directory / (name + 'TR.DAT')
        for path in (cdf, trdat):
            self.assertTrue(path.is_file(), f'required bundled fixture missing: {path}')
        facts = {'schema_version': 'transp2imas.curation.pilot.v1', 'run_identity': {'device': device, 'shot': shot},
                 'files': [{'role': 'trdat', 'path': str(trdat)}]}
        if mutate:
            mutate(facts)
        pilot = Path(tmp) / (name + '.pilot.json')
        pilot.write_text(json.dumps(facts))
        before = (sha(cdf), sha(trdat))
        done = subprocess.run([sys.executable, str(ROOT / 'curate.py'), '-i', str(cdf), '--trdat', str(trdat),
                               '--pilot', str(pilot), '--output-dir', str(Path(tmp) / name)],
                              capture_output=True, text=True)
        self.assertEqual((sha(cdf), sha(trdat)), before, 'curation modified a source file')
        document = None
        if (Path(tmp) / name / 'curation.json').exists():
            document = json.loads((Path(tmp) / name / 'curation.json').read_text())
        return done, document

    def test_each_case_curates_with_identity_and_hashes_verified_and_nothing_certified(self):
        for folder, name, device, shot in CASES:
            with self.subTest(name), tempfile.TemporaryDirectory() as tmp:
                done, document = self.curate(folder, name, device, shot, tmp)
                self.assertIn(done.returncode, (0, 1), done.stderr)
                found = {f['id']: f for f in document['findings']}
                self.assertEqual(found['RUN.IDENTITY']['status'], 'pass')
                self.assertEqual(found['INPUT.RESOLVED.FILES']['status'], 'pass')
                self.assertEqual(found['INPUT.RESOLVED.FILES']['evidence']['files'][0]['status'], 'verified')
                for role in ('electron_temperature', 'ion_temperature', 'rotation'):
                    self.assertEqual(found['EVIDENCE.' + role]['status'], 'not_assessed', role)   # no origin was supplied
                self.assertFalse(document['summary']['scientifically_certified'])
                self.assertEqual(found['EVENTS.INVENTORY']['status'], 'not_assessed')

    def test_a_wrong_declared_shot_is_a_failed_identity_check(self):
        folder, name, device, shot = CASES[0]
        with tempfile.TemporaryDirectory() as tmp:
            done, document = self.curate(folder, name, device, shot + 1, tmp)
            found = {f['id']: f for f in document['findings']}
            self.assertEqual(found['RUN.IDENTITY']['status'], 'fail')
            self.assertEqual(done.returncode, 1)                  # a fail finding reaches the default threshold

    def test_a_reference_for_another_shot_range_is_not_applied(self):
        folder, name, device, shot = CASES[1]
        def add(facts):
            facts['input_expectations'] = [{'name': 'EINJA(1)', 'expected': 1.0, 'reference': 'test', 'units': 'keV',
                                            'applicability': {'shot_min': shot + 1000, 'shot_max': shot + 2000}}]
        with tempfile.TemporaryDirectory() as tmp:
            _done, document = self.curate(folder, name, device, shot, tmp, add)
            found = {f['id']: f for f in document['findings']}
            self.assertEqual(found['HEATING.REFERENCE.0']['status'], 'not_applicable')


if __name__ == '__main__':
    unittest.main(verbosity=1)
