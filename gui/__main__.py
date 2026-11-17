#!/usr/bin/python
"""``python -m gui``: serve the report viewer (read-only) and, with ``--enable-run``, the conversion interface."""

import argparse
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import webbrowser

from gui import __version__, convert_plan
from gui.jobs import clean_environment
from gui.server import make_server

REPO = Path(__file__).resolve().parent.parent


def resolve_python(path):
    """The interpreter named by ``--imas-python``: a file, or a virtual-environment directory (its ``bin/python``)."""
    candidate = Path(path).expanduser()
    if candidate.is_dir():
        for name in ('python', 'python3'):
            if (candidate / 'bin' / name).is_file():
                candidate = candidate / 'bin' / name
                break
        else:
            raise ValueError(f'{path} is a directory without bin/python; give the interpreter itself or a virtual environment')
    if not candidate.is_file():
        raise ValueError(f'{path} does not exist')
    if not os.access(candidate, os.X_OK):
        raise ValueError(f'{candidate} is not executable')
    return str(candidate)


def check_imas(python):
    """The interpreter if it can import IMAS-Python, else ``None`` (the page then explains why running is off)."""
    try:
        done = subprocess.run([python, '-c', 'import imas'], capture_output=True, timeout=120, env=clean_environment())
    except (OSError, subprocess.TimeoutExpired):
        return None
    return python if done.returncode == 0 else None


def main(argv=None):
    parser = argparse.ArgumentParser(prog='python -m gui', description=__doc__)
    parser.add_argument('--root', action='append', type=Path, default=[],
                        help='directory whose reports may be opened and whose CDF files may be converted (repeatable; '
                        'default: the current directory)')
    parser.add_argument('--input-root', action='append', type=Path, default=[],
                        help='further directory CDF files and geometry JSON may be read from (repeatable)')
    parser.add_argument('--output-root', action='append', type=Path, default=[],
                        help='where conversion jobs and the scratch IMAS database are written (repeatable; the first is '
                        'used; default: build/gui in the repository)')
    parser.add_argument('--database-root', action='append', type=Path, default=[],
                        help='extra IMAS database directory whose HDF5 entries may be listed and plotted, read-only')
    parser.add_argument('--enable-run', action='store_true',
                        help='allow starting conversions (default: read-only viewer; planning is always available)')
    parser.add_argument('--imas-python', help='interpreter with IMAS-Python (or the virtual-environment directory that has it), '
                        'used to convert and to read entries (default: this interpreter if it has IMAS)')
    parser.add_argument('--allow-default-database', action='store_true',
                        help='allow writing to, and plotting from, your default IMAS database (~/public/imasdb)')
    parser.add_argument('--host', default='127.0.0.1', help='interface to bind (default: loopback only)')
    parser.add_argument('--port', type=int, default=0, help='port (default: a free one)')
    parser.add_argument('--allow-remote', action='store_true',
                        help='allow a non-loopback --host; prefer an SSH tunnel to the loopback port')
    parser.add_argument('--open', action='store_true', help='open the URL in a browser on this machine')
    parser.add_argument('--verbose', action='store_true', help='log each request')
    parser.add_argument('--version', action='version', version='gui ' + __version__)
    args = parser.parse_args(argv)
    roots = args.root or [Path.cwd()]
    note = None
    if args.imas_python:
        try:
            python = resolve_python(args.imas_python)
        except ValueError as exc:
            print(f'--imas-python: {exc}', file=sys.stderr)
            return 2
    else:
        python = sys.executable if importlib.util.find_spec('imas') else None
        if python is None:
            note = 'this interpreter has no IMAS-Python and --imas-python was not given; converting and plotting are off'
    if python:
        print(f'Checking that {python} can import IMAS-Python (the first import can take several seconds)...', flush=True)
        checked = check_imas(python)
        if checked is None:
            note = f'{python} cannot import imas; converting and plotting are off'
        python = checked
    output_roots = args.output_root or [REPO / 'build' / 'gui']
    config = convert_plan.Config(list(roots) + list(args.input_root), output_roots, python, args.enable_run,
                                 args.allow_default_database, imas_note=note)
    try:
        server = make_server(roots, args.host, args.port, allow_remote=args.allow_remote, verbose=args.verbose,
                             config=config, database_roots=args.database_root, repo_root=REPO)
    except (ValueError, OSError) as exc:
        print(f'Cannot start the viewer: {exc}', file=sys.stderr)
        return 2
    port = server.server_address[1]
    url = f'http://{args.host}:{port}/?token={server.token}'
    mode = 'conversion enabled' if args.enable_run and python else 'conversion off (no conversion can be started)'
    if python and not args.enable_run:
        mode += '; plotting is on, add --enable-run to allow conversions'
    print(f'Report viewer, curation and conversion interface: {mode}.', flush=True)
    print('Curation can be run from the page: it reads its inputs and writes only a report into the job directory below '
          f'{output_roots[0]}.', flush=True)
    if note:
        print('Note: ' + note, flush=True)
    elif args.enable_run and not python:
        print('Note: --enable-run needs an interpreter with IMAS-Python (--imas-python); conversion stays off', flush=True)
    print('Roots: ' + ', '.join(str(r) for r in server.roots), flush=True)
    print('Open: ' + url, flush=True)
    if args.host in ('127.0.0.1', 'localhost'):
        print(f'From another machine: ssh -L {port}:127.0.0.1:{port} <this host>, then open the same URL.', flush=True)
    if args.open:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print()
    finally:
        server.server_close()
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
