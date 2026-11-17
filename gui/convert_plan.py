#!/usr/bin/python
"""Planning a conversion: the option schema, validation, the command line and the target-entry guard.

Pure functions: no IMAS import, no subprocess, no file is written.  The form is *derived from the converter's own
argument parser* (``utilities.build_parser``), so it cannot drift from ``convert.py``: every option is either exposed
here or explicitly managed by the GUI, and a test fails when the converter gains an option nobody decided about.

A request is validated field by field; nothing is passed to the converter that did not come through a validator below,
and the command is an argument *list* (never a shell string).  The target entry is computed the way the converter
computes it (database root / machine / DD major / shot / run) so that an existing entry is shown, and replacing or
extending one needs the exact entry path typed back (``confirm_overwrite``).
"""

from pathlib import Path
import re
import shlex

import conversion_plan
import utilities
from gui import browse

BACKENDS = ('HDF5', 'HDF', 'MDS', 'MDSPLUS')
DEFAULT_BACKEND = 'HDF5'                      # the converter's own default is MDS; the form always sends one explicitly
DD_VERSION = re.compile(r'^\d+\.\d+\.\d+$')
MACHINE = re.compile(r'^[A-Za-z0-9_.-]{1,32}$')
INPUT_SUFFIXES = ('.cdf', '.pln')
JSON_SUFFIXES = ('.json',)
CHOICES = {'ip_direction': ('auto', 'ccw', 'cw'), 'bt_direction': ('auto', 'ccw', 'cw'), 'moments': ('yes', 'no'),
           'ec_nml_launch': ('yes', 'no'), 'log': ('debug', 'info', 'warn', 'error', 'critical')}
# parser dest -> how the GUI treats it.  'managed' options are set by the GUI itself (never taken from the form).
MANAGED = {'input_cdf': 'chosen with the TRANSP input picker', 'input_dir': 'batch conversion (-i all -d) is not offered',
           'user': 'database location is chosen in the Database section', 'report_json': 'always written into the job directory'}
EXPOSED = ('backend', 'run', 'new_db', 'occur', 'machine', 'input_format', 'pln_time_selection',
           'ip_direction', 'bt_direction', 'dd_version',
           'pf_passive_geometry', 'pf_active_geometry', 'ids', 'lenient_mapping', 'moments', 'nzones', 'ec_nml_launch',
           'debug', 'debug_max', 'verbose', 'log')
DD_CHOICES = ('auto', '4.1.1', '3.39.0')


class PlanError(ValueError):
    """The request is not acceptable; ``errors`` is ``[{'field', 'message'}]`` for the form."""

    def __init__(self, errors):
        self.errors = errors
        super().__init__('; '.join(f"{e['field']}: {e['message']}" for e in errors))


class Config:
    """What the server allows: where inputs may come from, where entries and jobs may go, and what is enabled."""

    def __init__(self, input_roots, output_roots, imas_python=None, enable_run=False, allow_default_database=False,
                 home=None, imas_note=None):
        self.input_roots = [Path(r) for r in input_roots]
        self.output_roots = [Path(r) for r in output_roots]
        self.imas_python = imas_python
        self.enable_run = enable_run
        self.allow_default_database = allow_default_database
        self.home = Path(home) if home else Path.home()
        self.imas_note = imas_note            # why there is no usable IMAS interpreter, shown on the page

    @property
    def scratch_database(self):
        return self.output_roots[0] / 'imasdb'

    @property
    def default_database(self):
        return self.home / 'public' / 'imasdb'


def option_schema():
    """The converter's options as a form description, from its own parser; managed options are listed, not exposed."""
    parser = utilities.build_parser()
    options = []
    for action in parser._actions:
        if action.dest == 'help':
            continue
        if action.option_strings and action.nargs == 0 and action.const is True:
            kind = 'flag'
        elif action.__class__.__name__ == '_CountAction':
            kind = 'count'
        elif action.choices:
            kind = 'choice'
        elif action.type is int:
            kind = 'int'
        else:
            kind = 'text'
        options.append({'dest': action.dest, 'flags': list(action.option_strings), 'help': action.help, 'kind': kind,
                        'choices': list(action.choices) if action.choices else None, 'default': action.default,
                        'managed_by_gui': MANAGED.get(action.dest), 'exposed': action.dest in EXPOSED})
    return {'options': options, 'ids': list(conversion_plan.WRITER_ORDER), 'backends': list(BACKENDS),
            'dd_versions': list(DD_CHOICES), 'default_backend': DEFAULT_BACKEND,
            'database_modes': ['scratch', 'custom', 'default']}


def _int(value, field, errors, low=0, high=None):
    if value in (None, ''):
        return None
    if isinstance(value, bool) or not isinstance(value, int):
        errors.append({'field': field, 'message': 'expected a whole number'})
        return None
    if value < low or (high is not None and value > high):
        errors.append({'field': field, 'message': f'must be at least {low}' + (f' and at most {high}' if high else '')})
        return None
    return value


def _choice(value, field, errors, allowed):
    if value in (None, ''):
        return None
    if value not in allowed:
        errors.append({'field': field, 'message': f'must be one of {", ".join(allowed)}'})
        return None
    return value


def _flag(value, field, errors):
    if value in (None, False, 0, ''):
        return False
    if value is True:
        return True
    errors.append({'field': field, 'message': 'expected true or false'})
    return False


def identify(input_path, machine_override=None, input_format='netcdf'):
    """Return machine, shot, run and run id from CDF metadata or a native PLN file name."""
    try:
        run_name = str(input_path)
        if input_format == 'pln':
            if not Path(run_name).name.upper().endswith('XF.PLN'):
                raise ValueError('choose the run XF.PLN file')
            run_name = str(Path(run_name).with_name(Path(run_name).name[:-6]))
        runno, shot, rid = utilities.find_runid_shot_code(run_name)
    except ValueError:
        raise PlanError([{'field': 'cdf', 'message': f'the file name {Path(input_path).name!r} does not give a shot and run id'}])
    device = None
    if input_format == 'netcdf':
        from curation.source import open_cdf
        try:
            with open_cdf(input_path) as source:
                device = source.global_attribute('device')
                attribute = source.global_attribute('shot')
                if isinstance(attribute, int):
                    shot = attribute
        except (OSError, ValueError, TypeError, RuntimeError, ImportError) as exc:
            raise PlanError([{'field': 'cdf', 'message': f'cannot read the CDF: {exc}'}])
    machine = machine_override or device
    if not machine:
        raise PlanError([{'field': 'machine', 'message': 'the input file does not identify the machine; give a machine name'}])
    return machine, shot, runno, rid


def database_root(request_database, config, errors):
    """``(root path, mode, passes -u?)`` for the database section; custom paths must lie below an output root."""
    mode = (request_database or {}).get('mode') or 'scratch'
    if mode not in ('scratch', 'custom', 'default'):
        errors.append({'field': 'database.mode', 'message': 'must be scratch, custom or default'})
        return None, mode, False
    if mode == 'scratch':
        return config.scratch_database, mode, True
    if mode == 'default':
        if not config.allow_default_database:
            errors.append({'field': 'database.mode', 'message': 'writing to the default IMAS database is not enabled on '
                           'this server (start it with --allow-default-database)'})
            return None, mode, False
        return config.default_database, mode, False
    requested = (request_database or {}).get('path')
    if not requested or '\0' in str(requested):
        errors.append({'field': 'database.path', 'message': 'give the database directory'})
        return None, mode, False
    import os
    real = Path(os.path.realpath(Path(requested) if Path(requested).is_absolute() else config.output_roots[0] / requested))
    if not any(real == Path(os.path.realpath(r)) or Path(os.path.realpath(r)) in real.parents for r in config.output_roots):
        errors.append({'field': 'database.path', 'message': 'must be inside an allowed output root: '
                       + ', '.join(str(r) for r in config.output_roots)})
        return None, mode, False
    return real, mode, True


def entry_candidates(root, machine, shot, run, backend, dd_version):
    """Where the converter would put the entry (HDF5: one directory per run; DD ``auto`` may be 3 or 4)."""
    majors = [dd_version.split('.')[0]] if dd_version and dd_version != 'auto' else ['3', '4']
    if backend in ('MDS', 'MDSPLUS'):
        return [root / machine / m / '0' for m in majors]
    return [root / machine / m / str(shot) / str(run) for m in majors]


def next_free_run(root, machine, shot, run, backend, dd_version, limit=1000):
    """The first run number above ``run`` with no entry (HDF5; MDSplus entries cannot be checked, so ``None``)."""
    if backend in ('MDS', 'MDSPLUS'):
        return None
    for candidate in range(run + 1, run + 1 + limit):
        if not any(_nonempty(c) for c in entry_candidates(root, machine, shot, candidate, backend, dd_version)):
            return candidate
    return None


def _nonempty(path):
    return path.is_dir() and any(path.iterdir())


def validate(request, config, report_json='<job directory>/conversion-report.json'):
    """A checked plan for ``request``: argv, display command, target entry, warnings and the overwrite requirement."""
    errors, warnings = [], []
    options = dict((request or {}).get('options') or {})
    unknown = sorted(set(options) - set(EXPOSED))
    if unknown:
        errors.append({'field': 'options', 'message': f'unknown or unmanaged option(s): {", ".join(unknown)}'})
    input_format = _choice(options.get('input_format') or 'netcdf', 'input_format', errors, ('netcdf', 'pln')) or 'netcdf'
    pln_time_selection = _choice(options.get('pln_time_selection') or 'all', 'pln_time_selection', errors, ('last', 'all')) or 'all'
    cdf = None
    try:
        cdf = browse.resolve(config.input_roots, str((request or {}).get('cdf') or ''), INPUT_SUFFIXES, max_bytes=None)
        if input_format == 'pln' and not cdf.name.upper().endswith('XF.PLN'):
            raise browse.AccessError('native TRANSP input must be the run XF.PLN file')
        if input_format == 'pln':
            prefix = cdf.name[:-6]
            browse.resolve(config.input_roots, str(cdf.with_name(prefix + 'YF.PLN')), ('.pln',), max_bytes=None)
            tf_path = cdf.with_name(prefix + 'TF.PLN')
            if not tf_path.is_file():
                browse.resolve(config.input_roots, str(cdf.with_name(prefix + '.CDF')), ('.cdf',), max_bytes=None)
        if input_format == 'netcdf' and cdf.suffix.lower() != '.cdf':
            raise browse.AccessError('NetCDF input must be a .CDF file')
    except browse.AccessError as exc:
        errors.append({'field': 'cdf', 'message': str(exc)})

    backend = _choice(options.get('backend') or DEFAULT_BACKEND, 'backend', errors, BACKENDS) or DEFAULT_BACKEND
    run = _int(options.get('run'), 'run', errors, 1)
    occur = _int(options.get('occur'), 'occur', errors, 0) or 0
    nzones = _int(options.get('nzones'), 'nzones', errors, 0, 4096) or 0
    verbose = _int(options.get('verbose'), 'verbose', errors, 0, 3) or 0
    dd_version = options.get('dd_version') or 'auto'
    if dd_version != 'auto' and not (isinstance(dd_version, str) and DD_VERSION.match(dd_version)):
        errors.append({'field': 'dd_version', 'message': 'use auto or a version like 4.1.1'})
        dd_version = 'auto'
    machine_override = options.get('machine') or None
    if machine_override is not None and not (isinstance(machine_override, str) and MACHINE.match(machine_override)):
        errors.append({'field': 'machine', 'message': 'letters, digits, _ . - only (at most 32)'})
        machine_override = None
    ip_direction = _choice(options.get('ip_direction'), 'ip_direction', errors, CHOICES['ip_direction'])
    bt_direction = _choice(options.get('bt_direction'), 'bt_direction', errors, CHOICES['bt_direction'])
    moments = _choice(options.get('moments'), 'moments', errors, CHOICES['moments'])
    ec_nml = _choice(options.get('ec_nml_launch'), 'ec_nml_launch', errors, CHOICES['ec_nml_launch'])
    log = _choice(options.get('log'), 'log', errors, CHOICES['log'])
    lenient = _flag(options.get('lenient_mapping'), 'lenient_mapping', errors)
    debug = _flag(options.get('debug'), 'debug', errors)
    debug_max = _flag(options.get('debug_max'), 'debug_max', errors)
    new_db = _flag(options.get('new_db', True), 'new_db', errors)

    ids = options.get('ids')
    if ids is not None:
        if not isinstance(ids, list) or not ids or any(i not in conversion_plan.WRITER_ORDER for i in ids):
            errors.append({'field': 'ids', 'message': 'choose one or more of: ' + ', '.join(conversion_plan.WRITER_ORDER)})
            ids = None
        else:
            ids = [i for i in conversion_plan.WRITER_ORDER if i in set(ids)]      # dependency-safe order, no duplicates
            if len(ids) == len(conversion_plan.WRITER_ORDER):
                ids = None
    geometry = {}
    for key in ('pf_passive_geometry', 'pf_active_geometry'):
        if options.get(key):
            try:
                geometry[key] = str(browse.resolve(config.input_roots, str(options[key]), JSON_SUFFIXES))
            except browse.AccessError as exc:
                errors.append({'field': key, 'message': str(exc)})

    root, mode, passes_user = database_root((request or {}).get('database'), config, errors)
    machine = shot = default_run = rid = None
    target = None
    if cdf is not None and not any(e['field'] in ('cdf', 'machine') for e in errors):
        try:
            machine, shot, default_run, rid = identify(cdf, machine_override, input_format)
        except PlanError as exc:
            errors.extend(exc.errors)
    if root is not None and machine and shot is not None:
        final_run = run if run else default_run
        candidates = entry_candidates(root, machine, shot, final_run, backend, dd_version)
        existing = [c for c in candidates if _nonempty(c)]
        mds = backend in ('MDS', 'MDSPLUS')
        target = {'database_root': str(root), 'database_mode': mode, 'machine': machine, 'shot': shot, 'run': final_run,
                  'occurrence': occur, 'backend': backend, 'dd_version': dd_version,
                  'entry_paths': [str(c) for c in candidates], 'existing': [str(c) for c in existing],
                  'exists': bool(existing), 'exists_unknown': mds,
                  'primary': str((existing or candidates)[0]),
                  'next_free_run': next_free_run(root, machine, shot, final_run, backend, dd_version) if existing else None}
        if mds:
            warnings.append('The MDSplus backend stores entries in a shared per-database directory; the page cannot tell '
                            'whether this run already exists, so it is treated as possibly existing and needs confirmation.')
        if dd_version == 'auto':
            warnings.append('DD version is auto: it follows the IMAS environment of the server process, and both major '
                            'versions are checked for an existing entry.')
        if existing and not new_db:
            warnings.append('The entry exists and "create new entry" is off: the selected IDSs will be written into the '
                            'existing entry (replacing those IDSs).')
        if existing and new_db:
            warnings.append('The entry exists: "create new entry" replaces it.')
        if not existing and not new_db and not mds:
            errors.append({'field': 'new_db', 'message': 'the entry does not exist yet; turn on "create new entry"'})
    if lenient:
        warnings.append('Variables that cannot be mapped (bad shape, unreadable) will be skipped and recorded with their reason '
                        'in the conversion report instead of stopping the conversion. By default the conversion stops at the '
                        'first one (converter flag --lenient-mapping).')
    if ids:
        warnings.append('Only the selected IDSs are written, so a new entry will be partial.')
    if occur:
        warnings.append('Occurrence > 0 is written beside the main occurrence 0 of the same entry.')
    if errors:
        raise PlanError(errors)

    argv = ['-i', str(cdf), '--input-format', input_format, '-b', backend]
    if input_format == 'pln':
        argv += ['--pln-time-selection', pln_time_selection]
    if run:
        argv += ['-r', str(run)]
    if occur:
        argv += ['-o', str(occur)]
    if passes_user:
        argv += ['-u', str(root)]
    if machine_override:
        argv += ['-m', machine_override]
    if ip_direction and ip_direction != 'auto':
        argv += ['--ip-direction', ip_direction]
    if bt_direction and bt_direction != 'auto':
        argv += ['--bt-direction', bt_direction]
    argv += ['--dd-version', dd_version]
    for key, flag in (('pf_passive_geometry', '--pf-passive-geometry'), ('pf_active_geometry', '--pf-active-geometry')):
        if key in geometry:
            argv += [flag, geometry[key]]
    if ids:
        argv += ['--ids', ','.join(ids)]
    if lenient:
        argv.append('--lenient-mapping')
    argv += ['--report-json', str(report_json)]
    if moments == 'yes':
        argv += ['-mom', 'yes']
    if nzones:
        argv += ['-zones', str(nzones)]
    if ec_nml == 'yes':
        argv += ['-ecnml', 'yes']
    if debug:
        argv.append('--debug')
    if debug_max:
        argv.append('--debug-max')
    argv += ['-v'] * verbose
    if log:
        argv += ['-l', log]
    if new_db:
        argv.append('-new')

    needs = bool(target['exists'] or target['exists_unknown'])
    return {'argv': argv, 'display': shlex.join(['python', 'convert.py', *argv]), 'target': target, 'warnings': warnings,
            'needs_confirmation': needs, 'confirm_phrase': target['primary'] if needs else None,
            'selected_ids': ids or list(conversion_plan.WRITER_ORDER), 'cdf': str(cdf),
            'input_format': input_format, 'pln_time_selection': pln_time_selection}


def confirm(plan, supplied):
    """Whether the typed phrase is exactly the entry path the plan says would be changed (or none is needed)."""
    if not plan['needs_confirmation']:
        return True
    return isinstance(supplied, str) and supplied == plan['confirm_phrase']


RUN_CDF = re.compile(r'[A-Za-z]*\d+(?:[A-Za-z]\d{2})?\Z')


def is_run_cdf(name):
    """Whether a file name is a TRANSP run output, ``<runid>.CDF`` (shot digits, then letter and two digits: 141711Y05.CDF).

    Companion files of the same run are not converted and are not offered: ``<runid>PH.CDF`` (plasma-history output),
    plasma-state files (``<runid>_ps.cdf``, ``_fi_*.cdf``, ...).  A file with another name can still be chosen by typing its path.
    """
    stem, dot, suffix = name.rpartition('.')
    return bool(dot) and suffix.lower() == 'cdf' and bool(RUN_CDF.match(stem))


def scan_inputs(roots, max_depth=5, limit=500):
    """Run CDF files (``<runid>.CDF``) below the input roots: ``[{'path', 'size', 'mtime', 'trdat'}]``, newest first."""
    import os
    found = []
    for root in roots:
        root = Path(os.path.realpath(root))
        base = len(root.parts)
        for directory, subdirs, files in os.walk(root):
            if len(Path(directory).parts) - base >= max_depth:
                subdirs[:] = []
            subdirs[:] = [d for d in subdirs if not d.startswith('.') and d not in ('__pycache__', 'build')]
            for name in sorted(files):
                if is_run_cdf(name) or name.upper().endswith('XF.PLN'):
                    path = Path(directory) / name
                    try:
                        stat = path.stat()
                    except OSError:
                        continue
                    shot_tr = [f for f in files if f.upper().endswith('TR.DAT') and f.upper().startswith(path.stem.upper())]
                    found.append({'path': str(path), 'size': stat.st_size, 'mtime': stat.st_mtime,
                                  'trdat': bool(shot_tr), 'format': 'pln' if name.upper().endswith('XF.PLN') else 'netcdf'})
                    if len(found) >= limit:
                        return sorted(found, key=lambda i: -i['mtime'])
    return sorted(found, key=lambda i: -i['mtime'])


def scan_entries(database_roots, max_depth=6, limit=500):
    """HDF5 entries below the database roots (a directory holding ``master.h5``), with the IDSs they contain."""
    import os
    found = []
    for root in database_roots:
        root = Path(os.path.realpath(root))
        if not root.is_dir():
            continue
        base = len(root.parts)
        for directory, subdirs, files in os.walk(root):
            if len(Path(directory).parts) - base >= max_depth:
                subdirs[:] = []
            if 'master.h5' not in files:
                continue
            path = Path(directory)
            names = sorted({re.sub(r'_\d+$', '', f[:-3]) for f in files if f.endswith('.h5') and f != 'master.h5'})
            parts = path.relative_to(root).parts
            found.append({'path': str(path), 'database_root': str(root),
                          'machine': parts[0] if len(parts) >= 4 else None, 'dd_major': parts[1] if len(parts) >= 4 else None,
                          'shot': parts[2] if len(parts) >= 4 else None, 'run': parts[3] if len(parts) >= 4 else None,
                          'ids': names, 'bytes': sum((path / f).stat().st_size for f in files if f.endswith('.h5')),
                          'mtime': max((path / f).stat().st_mtime for f in files if f.endswith('.h5'))})
            subdirs[:] = []
            if len(found) >= limit:
                break
    return sorted(found, key=lambda i: -i['mtime'])
