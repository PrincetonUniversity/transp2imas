#!/usr/bin/python
"""Planning a curation run: the option schema, validation and the command line.

Like ``convert_plan``, the form is derived from the program's own parser (``curate.build_parser``): every option is either
exposed here or explicitly managed, and a test fails when ``curate.py`` gains an option nobody decided about.  Pure functions:
nothing is run and nothing is written.  A request is validated field by field and the command is an argument *list*.

Curation reads the CDF, the TR.DAT and the evidence/policy files and writes only into the job directory; it never imports
IMAS.  Two things are deliberately not offered: ``--plugin`` (it imports a Python module by name, i.e. runs arbitrary code;
use the command line for that) and a free output directory (the report always goes into the job directory, so it can be
found again and never overwrites an input).
"""

import json
import os
from pathlib import Path
import shlex

from gui import browse

CDF_SUFFIXES = ('.cdf',)
JSON_SUFFIXES = ('.json',)
TRDAT_SUFFIXES = ('.dat',)
FAIL_ON = ('fail', 'error', 'warning')
TRDAT_MODES = ('auto', 'none', 'path')
EVIDENCE_KINDS = ('none', 'evidence', 'pilot')
MANAGED = {'input': 'chosen with the CDF picker', 'output_dir': 'the report is always written into the job directory',
           'plugin': 'a plugin is a Python module imported by name (arbitrary code): use curate.py --plugin on the command line',
           'version': 'not a run option', 'pilot': 'chosen with the evidence kind (evidence sidecar or pilot facts)'}
EXPOSED = ('trdat', 'evidence', 'policy', 'fail_on', 'validate_only', 'require_complete')
FORM_FIELDS = ('trdat_mode', 'trdat', 'evidence_kind', 'evidence', 'policy', 'fail_on', 'require_complete', 'validate_only')
MAX_CLASSIFY_BYTES = 4 * 1024 * 1024


class PlanError(ValueError):
    """The request is not acceptable; ``errors`` is ``[{'field', 'message'}]`` for the form."""

    def __init__(self, errors):
        self.errors = errors
        super().__init__('; '.join(f"{e['field']}: {e['message']}" for e in errors))


def availability():
    """``(True, None)`` when ``curate.py`` can run in this interpreter, else ``(False, reason)`` (missing dependencies)."""
    try:
        import curate                                                    # noqa: F401
    except ImportError as exc:
        return False, f'curation is not available in this interpreter: {exc}'
    return True, None


def schema_names():
    """The identifiers inside the JSON files, by what they are: ``{'evidence': (...), 'pilot': str, 'policy': str}``."""
    from curation import ACCEPTED_MANIFEST_SCHEMAS
    from curation.pilot import PILOT_SCHEMA
    from curation.policy import POLICY_SCHEMA
    return {'evidence': tuple(ACCEPTED_MANIFEST_SCHEMAS), 'pilot': (PILOT_SCHEMA,), 'policy': (POLICY_SCHEMA,)}


def classify(document):
    """``'evidence'``, ``'pilot'``, ``'policy'`` or ``None`` for a parsed JSON document, by its ``schema_version``."""
    if not isinstance(document, dict):
        return None
    version = document.get('schema_version')
    for kind, names in schema_names().items():
        if version in names:
            return kind
    return None


def option_schema():
    """The curation options as a form description, from curate.py's own parser."""
    import curate
    options = []
    for action in curate.build_parser()._actions:
        if action.dest == 'help':
            continue
        kind = ('flag' if action.option_strings and action.nargs == 0 else 'choice' if action.choices else
                'list' if action.__class__.__name__ == '_AppendAction' else 'path' if action.type is Path else 'text')
        options.append({'dest': action.dest, 'flags': list(action.option_strings), 'help': action.help, 'kind': kind,
                        'choices': list(action.choices) if action.choices else None,
                        'default': None if action.default is None or isinstance(action.default, list) else action.default,
                        'managed_by_gui': MANAGED.get(action.dest), 'exposed': action.dest in EXPOSED})
    available, reason = availability()
    return {'options': options, 'fail_on': list(FAIL_ON), 'evidence_kinds': list(EVIDENCE_KINDS), 'trdat_modes': list(TRDAT_MODES),
            'available': available, 'unavailable_reason': reason}


def trdat_beside(cdf):
    """The ``<stem>TR.DAT`` next to ``cdf`` (any letter case), or ``None``."""
    cdf = Path(cdf)
    wanted = (cdf.stem + 'TR.DAT').upper()
    try:
        for entry in sorted(cdf.parent.iterdir()):
            if entry.name.upper() == wanted and entry.is_file():
                return entry
    except OSError:
        return None
    return None


def scan_json(roots, max_depth=5, limit=300):
    """Evidence sidecars, pilot-facts files and policy files below the roots (recognised by ``schema_version``)."""
    found = {'evidence': [], 'pilot': [], 'policy': []}
    count = 0
    for root in roots:
        root = Path(os.path.realpath(root))
        base = len(root.parts)
        for directory, subdirs, files in os.walk(root):
            if len(Path(directory).parts) - base >= max_depth:
                subdirs[:] = []
            subdirs[:] = [d for d in subdirs if not d.startswith('.') and d not in ('__pycache__', 'build', 'jobs')]
            for name in sorted(files):
                if not name.lower().endswith('.json'):
                    continue
                path = Path(directory) / name
                try:
                    if path.stat().st_size > MAX_CLASSIFY_BYTES:
                        continue
                    kind = classify(json.loads(path.read_text(encoding='utf-8')))
                except (OSError, UnicodeDecodeError, ValueError):
                    continue
                if kind:
                    found[kind].append({'path': str(path), 'size': path.stat().st_size, 'name': name})
                    count += 1
                    if count >= limit:
                        return found
    return found


def _resolve_json(config, value, field, expected, errors):
    """A JSON input below the input roots whose ``schema_version`` says it is a ``expected`` file; else an error."""
    try:
        path = browse.resolve(config.input_roots, str(value), JSON_SUFFIXES)
    except browse.AccessError as exc:
        errors.append({'field': field, 'message': str(exc)})
        return None
    try:
        document = json.loads(path.read_text(encoding='utf-8'))
    except (OSError, UnicodeDecodeError, ValueError) as exc:
        errors.append({'field': field, 'message': f'not readable JSON: {exc}'})
        return None
    found = classify(document)
    if found != expected:
        label = {'evidence': 'an evidence sidecar', 'pilot': 'pilot facts', 'policy': 'a curation policy'}
        errors.append({'field': field, 'message': f'{path.name} is {"not a curation file" if found is None else label[found]}, '
                       f'not {label[expected]} (its schema_version decides)'})
        return None
    return path


def validate(request, config, output_dir='<job directory>/curation'):
    """A checked plan: ``argv``, the display command, where the report will be, warnings and what was resolved."""
    errors, warnings = [], []
    request = request or {}
    options = dict(request.get('options') or {})
    unknown = sorted(set(options) - set(FORM_FIELDS))
    if unknown:
        errors.append({'field': 'options', 'message': f'unknown or unmanaged option(s): {", ".join(unknown)}'})
    available, reason = availability()
    if not available:
        errors.append({'field': 'options', 'message': reason})

    def flag(key):
        value = options.get(key, False)
        if value in (None, False, 0, ''):
            return False
        if value is True:
            return True
        errors.append({'field': key, 'message': 'expected true or false'})
        return False

    validate_only = flag('validate_only')
    require_complete = flag('require_complete')
    fail_on = options.get('fail_on') or 'fail'
    if fail_on not in FAIL_ON:
        errors.append({'field': 'fail_on', 'message': f'must be one of {", ".join(FAIL_ON)}'})
        fail_on = 'fail'

    cdf = None
    requested_cdf = str(request.get('cdf') or '')
    if requested_cdf or not validate_only:
        try:
            cdf = browse.resolve(config.input_roots, requested_cdf, CDF_SUFFIXES, max_bytes=None)    # TRANSP output is large
        except browse.AccessError as exc:
            errors.append({'field': 'cdf', 'message': str(exc)})

    trdat_mode = options.get('trdat_mode') or 'auto'
    if trdat_mode not in TRDAT_MODES:
        errors.append({'field': 'trdat_mode', 'message': f'must be one of {", ".join(TRDAT_MODES)}'})
        trdat_mode = 'auto'
    trdat = None
    if trdat_mode == 'path':
        try:
            trdat = browse.resolve(config.input_roots, str(options.get('trdat') or ''), TRDAT_SUFFIXES, max_bytes=None)
        except browse.AccessError as exc:
            errors.append({'field': 'trdat', 'message': str(exc)})
    elif trdat_mode == 'auto' and cdf is not None:
        beside = trdat_beside(cdf)
        if beside is not None:
            try:
                trdat = browse.resolve(config.input_roots, str(beside), TRDAT_SUFFIXES, max_bytes=None)
            except browse.AccessError:
                trdat = None

    kind = options.get('evidence_kind') or 'none'
    if kind not in EVIDENCE_KINDS:
        errors.append({'field': 'evidence_kind', 'message': f'must be one of {", ".join(EVIDENCE_KINDS)}'})
        kind = 'none'
    evidence = None
    if kind != 'none':
        if not options.get('evidence'):
            errors.append({'field': 'evidence', 'message': 'choose the file'})
        else:
            evidence = _resolve_json(config, options['evidence'], 'evidence', kind, errors)
    policy = None
    if options.get('policy'):
        policy = _resolve_json(config, options['policy'], 'policy', 'policy', errors)

    if require_complete and validate_only:
        warnings.append('"Also fail when a rule did not run" has no effect when only validating.')
    if validate_only:
        warnings.append('Only the files are validated (structure, migration, semantics, referenced resources). Checks that need '
                        'the CDF (coordinates, units, grids, variable names, covariance shapes) are not run.')
        if evidence is None and policy is None:
            warnings.append('No evidence, pilot or policy file is chosen, so only the built-in defaults are validated.')
    else:
        if evidence is None:
            warnings.append('No evidence sidecar or pilot facts: rules that need facts you did not supply are reported as '
                            'not_assessed, never as passed.')
        if trdat is None:
            warnings.append('No TR.DAT is associated: the declared timeline is reported as unknown.')
        if policy is None:
            warnings.append('No policy: every rule runs with the sidecar thresholds.')
    warnings.append('Curation is a screening, not a scientific certification, whatever the exit status.')
    if errors:
        raise PlanError(errors)

    argv = []
    if cdf is not None:
        argv += ['-i', str(cdf)]
    if trdat is not None:
        argv += ['--trdat', str(trdat)]
    if evidence is not None:
        argv += ['--pilot' if kind == 'pilot' else '--evidence', str(evidence)]
    if policy is not None:
        argv += ['--policy', str(policy)]
    argv += ['--fail-on', fail_on]
    if require_complete and not validate_only:
        argv.append('--require-complete')
    if validate_only:
        argv.append('--validate-only')
    else:
        argv += ['--output-dir', str(output_dir)]
    report = None if validate_only else str(Path(output_dir) / 'curation.json')
    return {'argv': argv, 'display': shlex.join(['python', 'curate.py', *argv]), 'warnings': warnings, 'report': report,
            'validate_only': validate_only, 'cdf': None if cdf is None else str(cdf),
            'trdat': None if trdat is None else str(trdat), 'trdat_mode': trdat_mode,
            'evidence': None if evidence is None else {'kind': kind, 'path': str(evidence)},
            'policy': None if policy is None else str(policy), 'fail_on': fail_on, 'require_complete': require_complete and not validate_only}
