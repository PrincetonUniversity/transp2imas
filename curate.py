#!/usr/bin/python
"""Read-only TRANSP data curation CLI. Does not import IMAS or run a converter."""

import argparse
from pathlib import Path
import sys

from curation import __version__
from curation.engine import CurationEngine
from curation.inputs import declared_timeline
from curation.manifest import load_manifest
from curation.reporting import write_report
from curation.source import digest_file, open_cdf


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--input', required=True, type=Path, help='TRANSP CDF (never modified)')
    parser.add_argument('--trdat', type=Path, help='explicitly associated TR.DAT; omission is reported as unknown')
    parser.add_argument('--evidence', type=Path, help='strict JSON evidence/uncertainty/configuration sidecar')
    parser.add_argument('--output-dir', type=Path, help='default: build/curation/<CDF-stem>')
    parser.add_argument('--fail-on', choices=('error', 'fail', 'warning'), default='fail',
                        help='exit 1 at this severity or worse; exit 0 does not imply scientific acceptance')
    parser.add_argument('--require-complete', action='store_true', help='also exit 1 when any finding is not_assessed')
    parser.add_argument('--version', action='version', version='curation ' + __version__)
    args = parser.parse_args(argv)
    try:
        manifest = load_manifest(args.evidence)
        timeline = declared_timeline(args.trdat)
        with open_cdf(args.input) as source:
            report = CurationEngine(source, manifest, timeline).run()
        # Record original sidecar identity as well as normalized effective-manifest identity.
        report.source['evidence_file'] = (None if args.evidence is None else
                                         {'path': str(args.evidence.resolve()), 'sha256': digest_file(args.evidence)})
        report.source['converter_version'] = (Path(__file__).parent / 'version.txt').read_text().strip()
        report.source['curation_module_sha256'] = {
            p.name: digest_file(p) for p in sorted((Path(__file__).parent / 'curation').glob('*.py'))}
        directory = args.output_dir or Path('build') / 'curation' / args.input.stem
        targets, document = write_report(report, directory,
                                         protected_paths=(args.input, args.trdat, args.evidence))
    except (OSError, ValueError, TypeError, RuntimeError, ImportError, KeyError, IndexError, AttributeError) as exc:
        # exit 2 is "could not complete"; an unexpected data error must not look like exit 1, "findings at the threshold"
        print(f'Curation could not complete: {type(exc).__name__}: {exc}', file=sys.stderr)
        return 2
    counts = document['summary']['counts']
    print('Curation screening: ' + ', '.join(f'{key}={value}' for key, value in counts.items()))
    print('Scientific certification: NOT performed')
    for path in targets.values():
        print(path)
    blocking = {'error'}
    if args.fail_on in ('fail', 'warning'):
        blocking.add('fail')
    if args.fail_on == 'warning':
        blocking.add('warning')
    if args.require_complete:
        blocking.add('not_assessed')
    return int(any(counts[status] for status in blocking))


if __name__ == '__main__':
    raise SystemExit(main())
