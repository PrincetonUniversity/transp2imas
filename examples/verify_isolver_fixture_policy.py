#!/usr/bin/python
"""Regress strict/local versus allow-missing/CI ISolver fixture policies."""

from __future__ import annotations

from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest import mock

import verify_isolver_mapping as mapping


class TestISolverFixturePolicy(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        base = Path(self.temp.name)
        self.present_cdf = base / '38530P94.CDF'
        self.present_trdat = base / '38530P94TR.DAT'
        self.present_cdf.write_bytes(b'placeholder')
        self.present_trdat.write_text('LEVGEO=12\n')
        self.missing_cdf = base / '142202U01.CDF'
        self.missing_trdat = base / '142202U01TR.DAT'
        self.cases = (
            ('ITER', '38530P94', self.present_cdf, self.present_trdat, 11, 119),
            ('NSTX', '142202U01', self.missing_cdf, self.missing_trdat, 9, 82),
        )

    def invoke(self, *flags: str, fake_verify: bool = True) -> tuple[int, str]:
        output = StringIO()
        args = ['verify_isolver_mapping.py', *flags]
        rows = ([{'coil': 'CS3U', 'max_rel_cc_closure_error': '0.0'}],
                [{'cdf_variable': 'PS_1'}])

        def verify_stub(cdf: Path, trdat: Path, *unused):
            mapping.verify(cdf.is_file(), f'missing ISolver fixture: {cdf}')
            mapping.verify(trdat.is_file(), f'missing ISolver namelist: {trdat}')
            return rows

        with mock.patch.object(mapping, 'CASES', self.cases), \
                mock.patch.object(sys, 'argv', args), redirect_stdout(output):
            if fake_verify:
                with mock.patch.object(mapping, 'verify_run', side_effect=verify_stub):
                    status = mapping.main()
            else:
                status = mapping.main()
        return status, output.getvalue()

    def test_allow_missing_checks_available_and_reports_skipped(self) -> None:
        report = Path(self.temp.name) / 'reports'
        report.mkdir()
        stale = report / '142202U01-pf-active.csv'
        stale.write_text('stale previous output')
        status, output = self.invoke('--allow-missing', '--report-dir', str(report))
        self.assertEqual(status, 0)
        self.assertIn('ITER 38530P94: PASS', output)
        self.assertIn('NSTX 142202U01: SKIP', output)
        self.assertIn('PARTIAL (1/2 real fixtures verified', output)
        self.assertTrue((report / '38530P94-pf-active.csv').is_file())
        self.assertFalse(stale.exists(), 'skipped runs must not retain stale reports')

    def test_strict_local_mode_requires_missing_cdf(self) -> None:
        with self.assertRaisesRegex(AssertionError, 'missing ISolver fixture'):
            self.invoke()

    def test_present_cdf_still_requires_its_trdat(self) -> None:
        self.cases = (('ITER', '38530P94', self.present_cdf,
                       Path(self.temp.name) / 'missingTR.DAT', 11, 119),)
        with self.assertRaisesRegex(AssertionError, 'missing ISolver namelist'):
            self.invoke('--allow-missing', fake_verify=False)

    def test_allow_missing_never_passes_zero_fixtures(self) -> None:
        self.cases = (self.cases[1],)
        with self.assertRaisesRegex(AssertionError, 'No ISolver CDF fixtures'):
            self.invoke('--allow-missing')


if __name__ == '__main__':
    unittest.main(verbosity=2)
