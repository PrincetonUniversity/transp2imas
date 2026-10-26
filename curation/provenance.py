#!/usr/bin/python
"""Identity of the code, artifacts and dependencies that produced a report, hashed once per invocation.

The hash covers every Python module of the ``curation`` package (subpackages included), the shared parser modules the
package imports from the repository root, the packaged schema/policy artifacts, and the versions of the dependencies
that affect results.  Results are cached for the life of the process, so a rule is never re-hashed per call.
"""

from functools import lru_cache
import hashlib
from importlib import metadata
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SHARED_MODULES = ('fortran_literals.py', 'input_provenance.py')
DEPENDENCIES = ('numpy', 'scipy', 'netCDF4', 'f90nml', 'jsonschema', 'referencing')


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@lru_cache(maxsize=None)
def code_identity():
    """``{'modules': {relative path: sha256}, 'artifacts': {...}, 'dependencies': {name: version or None}}``."""
    package = ROOT / 'curation'
    modules = {str(p.relative_to(ROOT)): _sha256(p) for p in sorted(package.rglob('*.py'))}
    for name in SHARED_MODULES:
        if (ROOT / name).is_file():
            modules[name] = _sha256(ROOT / name)
    artifacts = {str(p.relative_to(ROOT)): _sha256(p) for p in sorted(package.rglob('*.json'))}
    dependencies = {}
    for name in DEPENDENCIES:
        try:
            dependencies[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            dependencies[name] = None
    return {'modules': modules, 'artifacts': artifacts, 'dependencies': dependencies}
