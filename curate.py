#!/usr/bin/python
"""Read-only TRANSP data curation CLI. Does not import IMAS or run a converter."""

import argparse
from importlib.metadata import PackageNotFoundError, version as package_version
from pathlib import Path
import sys

from curation import __version__
from curation.engine import CurationEngine
from curation.inputs import declared_timeline
from curation.manifest import load_manifest, referenced_paths
from curation.pilot import load_pilot
from curation.policy import Policy, load_policy
from curation.provenance import code_identity
from curation.rules import default_registry
from curation.reporting import write_report
from curation.source import digest_file, open_cdf


def validate_only(args, manifest, policy):
    """Report what validation established, and say plainly what it did not check (nothing here reads a source)."""
    migration = manifest['migration']
    print(f"Evidence valid: declared {migration['from']}, canonical {manifest['schema_version']}")
    for step in migration['steps']:
        print('  migrated: ' + step)
    print(f"Policy valid: {policy.name} (fingerprint {policy.fingerprint()[:12]})")
    print('Not checked without a source: coordinates, units, grid compatibility, variable names, covariance shapes.')
    return 0


def build_parser():
    """The command-line parser (also read by the GUI, so its form cannot drift from the options that exist)."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--input', type=Path, help='TRANSP CDF (never modified); required unless --validate-only')
    parser.add_argument('--trdat', type=Path, help='explicitly associated TR.DAT; omission is reported as unknown')
    parser.add_argument('--evidence', type=Path, help='strict JSON evidence/uncertainty/configuration sidecar')
    parser.add_argument('--pilot', type=Path, help='analyst-supplied pilot facts (replaces --evidence); '
                        'input hashes are computed, placeholders are refused')
    parser.add_argument('--policy', type=Path, help='curation policy (rules to run, screening parameters, applications); '
                        'default: every rule, sidecar thresholds')
    parser.add_argument('--plugin', action='append', default=[], metavar='MODULE',
                        help='trusted external rule module exposing register(registry); opt-in, recorded with its hash')
    parser.add_argument('--output-dir', type=Path, help='default: build/curation/<CDF-stem>')
    parser.add_argument('--fail-on', choices=('error', 'fail', 'warning'), default='fail',
                        help='exit 1 at this severity or worse; exit 0 does not imply scientific acceptance')
    parser.add_argument('--validate-only', action='store_true',
                        help='validate the evidence/pilot/policy files (structure, migration, semantics, referenced '
                        'resources) without a CDF; checks that need source coordinates or units are NOT run')
    parser.add_argument('--require-complete', action='store_true', help='also exit 1 when any finding is not_assessed')
    parser.add_argument('--version', action='version', version='curation ' + __version__)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.input is None and not args.validate_only:
        parser.error('-i/--input is required unless --validate-only')
    try:
        if args.pilot is not None and args.evidence is not None:
            parser.error('--pilot and --evidence are alternatives')
        manifest = load_pilot(args.pilot) if args.pilot is not None else load_manifest(args.evidence)
        timeline = declared_timeline(args.trdat)
        policy = load_policy(args.policy) if args.policy is not None else Policy.default()
        if args.validate_only:
            return validate_only(args, manifest, policy)
        registry = default_registry()
        for module in args.plugin:
            registry.load_plugin(module)
        with open_cdf(args.input) as source:
            report = CurationEngine(source, manifest, timeline, policy, registry).run()
        # Record original sidecar identity as well as normalized effective-manifest identity.
        given = args.pilot or args.evidence
        report.source['evidence_file'] = (None if given is None else
                                         {'path': str(given.resolve()), 'sha256': digest_file(given),
                                          'kind': 'pilot_facts' if args.pilot else 'evidence_sidecar'})
        try:
            installed_version = package_version('transp2imas')
        except PackageNotFoundError:
            installed_version = (Path(__file__).parent / 'version.txt').read_text().strip()
        report.source['converter_version'] = installed_version
        identity = code_identity()                          # hashed once per invocation, recursively
        report.source['code_identity'] = identity
        report.source['curation_module_sha256'] = identity['modules']
        if args.policy is not None:
            report.source['policy_file'] = {'path': str(args.policy.resolve()), 'sha256': digest_file(args.policy)}
        directory = args.output_dir or Path('build') / 'curation' / args.input.stem
        targets, document = write_report(report, directory,
                                         protected_paths=(args.input, args.trdat, args.evidence, args.pilot, args.policy, *referenced_paths(manifest)))
    except (OSError, ValueError, TypeError, RuntimeError, ImportError, KeyError, IndexError, AttributeError) as exc:
        # exit 2 is "could not complete"; an unexpected data error must not look like exit 1, "findings at the threshold"
        print(f'Curation could not complete: {type(exc).__name__}: {exc}', file=sys.stderr)
        return 2
    counts = document['summary']['counts']
    print('Curation screening: ' + ', '.join(f'{key}={value}' for key, value in counts.items()))
    print('Scientific certification: NOT performed')
    for path in targets.values():
        print(path)
    if report.has_execution_errors():
        # a rule raised: the report is written, but the screening did not complete
        failed = [e['rule'] for e in report.ledger if e['state'] == 'execution_error']
        print('Curation could not complete: rule execution error in ' + ', '.join(failed), file=sys.stderr)
        return 2
    if args.require_complete and report.not_executed():
        print('Not all rules executed: ' + ', '.join(report.not_executed()), file=sys.stderr)
        return 1
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
