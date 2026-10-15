#!/usr/bin/python
"""Check the Release Please feat/fix contract and workflow wiring locally.

This validates repository configuration only. Release Please still uses real
merged commit messages, which cannot be determined from a source tarball.
"""
import json
from pathlib import Path
import re
import unittest

ROOT = Path(__file__).resolve().parents[1]


class TestReleaseConfiguration(unittest.TestCase):
    def test_feat_minor_and_fix_patch_are_enabled(self):
        config = json.loads((ROOT/'release-please-config.json').read_text())
        pkg = config['packages']['.']
        self.assertEqual(pkg['release-type'], 'simple')
        self.assertEqual(pkg['versioning-strategy'], 'default')
        self.assertIs(pkg['bump-patch-for-minor-pre-major'], False)
        self.assertIs(pkg['bump-minor-pre-major'], False)
        sections = {item['type']: item for item in config['changelog-sections']}
        self.assertIn('feat', sections)
        self.assertIn('fix', sections)
        self.assertFalse(sections['feat'].get('hidden', False))
        self.assertFalse(sections['fix'].get('hidden', False))
        self.assertEqual(pkg['version-file'], 'version.txt')

    def test_workflow_runs_on_main_and_has_write_permissions(self):
        text = (ROOT/'.github/workflows/release.yml').read_text()
        self.assertIn('googleapis/release-please-action@v4', text)
        self.assertRegex(text, r'(?m)^\s+- main\s*$')
        for permission in ('contents: write', 'pull-requests: write'):
            self.assertIn(permission, text)
        self.assertIn('config-file: release-please-config.json', text)
        self.assertIn('manifest-file: .release-please-manifest.json', text)

    def test_conventional_commit_title_examples(self):
        # This is only a local title sanity check, NOT a replacement for the
        # Release Please parser, which uses the merged GitHub commit history.
        pattern = re.compile(r'^(feat|fix|perf|refactor|docs|test|chore)(\([\w.-]+\))?!?: .+')
        for title in ('feat: add GLF23 theory transport',
                      'feat(theory): support MMM07 components',
                      'fix: preserve TIME3', 'fix(imas): correct DD4 coordinates',
                      'perf: avoid redundant CDF scan', 'feat!: restructure IDS output'):
            self.assertRegex(title, pattern)
        for title in ('Feature: new physics', 'Merge pull request #13',
                      'feat add physics', 'Feat: update'):
            self.assertNotRegex(title, pattern)


if __name__ == '__main__':
    raise SystemExit(0 if unittest.TextTestRunner(verbosity=1).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(TestReleaseConfiguration)
    ).wasSuccessful() else 1)
