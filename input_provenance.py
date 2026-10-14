#!/usr/bin/python
"""Lossless TRANSP/model namelist provenance and per-model XML, without IMAS imports.

The structured view is deliberately conservative: original file contents are
preserved verbatim (apart from forbidden XML-1.0 control characters); parsed
assignments are a searchable index, not a substitute for a Fortran parser.
"""

from dataclasses import dataclass
import hashlib
import os
from pathlib import Path
import re
import xml.etree.ElementTree as ET

MAX_INPUT_BYTES = 10 * 1024 * 1024
MAX_INPUT_FILES = 128
MAX_REFERENCE_DEPTH = 12

WAVE_MODELS = ('TORIC', 'FPP', 'CQL3D', 'GENRAY', 'TORBEAM', 'TORAY', 'LSC')
TRANSPORT_MODELS = ('PT_SOLVER', 'TGLF', 'GLF23', 'MMM')
MODEL_FAMILIES = {
    'waves': WAVE_MODELS,
    'core_transport': TRANSPORT_MODELS,
    'nbi': ('NUBEAM',),
}


@dataclass(frozen=True)
class InputRecord:
    role: str
    key: str
    reference: str
    path: Path
    depth: int
    status: str
    text: str
    sha256: str
    model: str = ''


def xml_safe_text(text):
    """Escape forbidden XML 1.0 code points, preserving all readable input."""
    return ''.join(ch if (ch in '\t\n\r' or 0x20 <= ord(ch) <= 0xD7FF
                          or 0xE000 <= ord(ch) <= 0xFFFD
                          or 0x10000 <= ord(ch) <= 0x10FFFF) else '\ufffd'
                   for ch in str(text))


def strip_fortran_comment(line):
    """Discard ! comments while respecting quote characters and doubled quotes."""
    quote = None
    out = []
    i = 0
    while i < len(line):
        ch = line[i]
        if quote and ch == quote and i + 1 < len(line) and line[i+1] == quote:
            out.extend((ch, ch))
            i += 2
            continue
        if ch in ("'", '"'):
            if quote is None:
                quote = ch
            elif quote == ch:
                quote = None
        if ch == '!' and quote is None:
            break
        out.append(ch)
        i += 1
    return ''.join(out)


def assignments(text):
    """Yield source-line, key, raw assignment value for simple Fortran inputs."""
    for line_no, line in enumerate(str(text).splitlines(), start=1):
        uncommented = strip_fortran_comment(line)
        match = re.match(r'^\s*([A-Za-z][A-Za-z0-9_%()]*?)\s*=\s*(.*)$', uncommented)
        if match:
            yield line_no, match.group(1), match.group(2).strip()


def parse_namelist(parent, text):
    """Append a non-evaluating XML index, including groups and source line nos."""
    parsed = ET.SubElement(parent, 'parsed_namelist')
    group = ET.SubElement(parsed, 'group', {'name': 'GLOBAL'})
    last = None
    for line_no, raw_line in enumerate(str(text).splitlines(), start=1):
        line = strip_fortran_comment(raw_line).strip()
        if not line:
            continue
        if line.startswith('&') or line.startswith('$'):
            group_name = line[1:].split()[0] if line[1:].split() else 'GLOBAL'
            if group_name.upper() in ('END',):
                group = ET.SubElement(parsed, 'group', {'name': 'GLOBAL'})
            else:
                group = ET.SubElement(parsed, 'group', {'name': group_name})
            last = None
            continue
        if line == '/' or line.startswith('/ '):
            group = ET.SubElement(parsed, 'group', {'name': 'GLOBAL'})
            last = None
            continue
        match = re.match(r'^([A-Za-z][A-Za-z0-9_%()]*)\s*=\s*(.*)$', line)
        if match:
            last = ET.SubElement(group, 'parameter', {
                'name': match.group(1), 'value': match.group(2),
                'source_line': str(line_no),
            })
        elif last is not None:
            # A continuation does not create a new assignment. Original text
            # remains the authority if the parser is insufficient.
            last.set('value', (last.get('value', '') + ' ' + line).strip())
    return parsed


def file_model(key, reference=''):
    """Best-evidence input owner from the reference name, not file contents."""
    key = str(key).upper()
    path = Path(str(reference)).name.upper()
    combined = key + ' ' + path
    for model, tags in (
        ('CQL3D', ('CQL3D', 'CQL_')),
        ('TORBEAM', ('TORBEAM',)),
        ('GENRAY', ('GENRAY',)),
        ('TORAY', ('TORAY',)),
        ('TORIC', ('TORIC',)),
        ('FPP', ('FPPRF', 'FPP_')),
        ('LSC', ('LSC_', 'LSC.')),
        ('TGLF', ('TGLF',)),
        ('GLF23', ('GLF23',)),
        ('MMM', ('MMM_', 'MMM.')),
        ('PT_SOLVER', ('PT_TEMPLATE', 'PT_SOLVER', 'PT_NML', 'PT_NAMELIST', 'PT.NML', '_PT.NML')),
        ('NUBEAM', ('NUBEAM',)),
    ):
        if any(tag in combined for tag in tags):
            return model
    return ''


def file_reference_key(key):
    key = key.upper()
    if 'TEMPLATE' in key:
        return True
    return bool(re.search(
        r'(?:^|_)(?:NML|NAMELIST|INFILE|INPUTFILE|INPUT_FILE|INPUT_NML|'
        r'INPUT_NAMELIST|CONFIG_FILE|PARFILE|PARAMETER_FILE)(?:_|$)', key
    ) or (key.startswith(('TORIC', 'CQL', 'GENRAY', 'TORBEAM', 'TORAY',
                           'LSC', 'FPP')) and key.endswith(('_FILE', '_INPUT'))))


def scalar_file_value(raw):
    """Accept a scalar filename, rejecting expression/list/numeric values."""
    s = str(raw).strip().rstrip(',').strip()
    if not s:
        return ''
    if s[0] in ("'", '"'):
        quote = s[0]
        end = s.find(quote, 1)
        if end < 1 or s[end+1:].strip().strip(','):
            return ''
        s = s[1:end]
    elif re.search(r'\s|,', s):
        return ''
    if not s or s.upper() in ('.TRUE.', '.FALSE.', 'T', 'F', 'NONE', 'NULL'):
        return ''
    if re.fullmatch(r'[+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+\-]?\d+)?', s):
        return ''
    if '*' in s or '?' in s or s.endswith(('/', '\\')):
        return ''
    return s


def references_from_text(text):
    """Yield (input key, filename) for explicit namelist/file references."""
    for _line, key, value in assignments(text):
        filename = scalar_file_value(value)
        if not filename:
            continue
        # A literal *.nml / *.namelist / *.inp reference is unambiguous
        # even when an older input uses a non-standard keyword name.
        named_namelist = bool(re.search(r'(?i)\.(?:nml|namelist|inp)$', filename))
        if file_reference_key(key) or named_namelist:
            yield key, filename
    # TRANSP/model namelists sometimes use #include or include directives.
    for line in str(text).splitlines():
        raw = strip_fortran_comment(line).strip()
        match = re.match(r'^#?\s*include\s+(?:[\'\"]([^\'\"]+)[\'\"]|(\S+))\s*$', raw, re.I)
        if match:
            filename = scalar_file_value(match.group(1) or match.group(2))
            if filename:
                yield 'INCLUDE', filename


def model_hints_from_text(text):
    """Identify configured model codes; never infer a run from comments alone."""
    active = {}
    for _line, key, value in assignments(text):
        active[key.upper()] = value.strip().strip("'\"").upper()
    models = set()
    for key in ('EC_MODEL', 'IC_MODEL', 'LH_MODEL'):
        value = active.get(key, '')
        for model in WAVE_MODELS:
            if model in value:
                models.add(model)
    # TRANSP legacy IC input can select TORIC without IC_MODEL; avoid
    # interpreting comments that merely mention TORIC as execution evidence.
    if 'NMDTORIC' in active or active.get('NTORIC_PSERVE', '') not in ('', '0'):
        models.add('TORIC')
    for key, value in active.items():
        if file_reference_key(key) and scalar_file_value(value):
            owner = file_model(key, value)
            if owner:
                models.add(owner)
    return models


def read_input(path):
    try:
        if not path.is_file():
            return '', '', 'missing'
        size = path.stat().st_size
        if size > MAX_INPUT_BYTES:
            return '', '', f'too_large:{size}'
        raw = path.read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        try:
            text = raw.decode('utf-8')
        except UnicodeDecodeError:
            text = raw.decode('ISO-8859-1')
        return xml_safe_text(text), digest, 'embedded'
    except OSError as exc:
        return '', '', f'error:{type(exc).__name__}'


def append_file_xml(parent, record):
    attrs = {
        'role': record.role, 'key': record.key,
        'reference': record.reference, 'resolved_path': str(record.path),
        'status': record.status,
    }
    if record.role == 'ADJACENT_UNREFERENCED':
        attrs['association'] = 'same_run_directory_only_not_explicitly_referenced'
    if record.model:
        attrs['model'] = record.model
    if record.sha256:
        attrs['sha256'] = record.sha256
    node = ET.SubElement(parent, 'file', attrs)
    if record.status == 'embedded':
        ET.SubElement(node, 'content', {'encoding': 'text'}).text = record.text
        parse_namelist(node, record.text)
    return node


class InputInventory:
    """Bounded, recursive inventory of one TRANSP run's text inputs."""

    def __init__(self, trdat, runid=''):
        self.trdat = Path(trdat).expanduser() if trdat else None
        self.runid = str(runid).strip()
        self.records = []
        self.truncated = False
        self._collect()

    def _resolve(self, filename, parent):
        candidate = Path(os.path.expandvars(filename)).expanduser()
        # The top-level TR.DAT path comes from convert.py, which may pass a
        # relative path containing parent directories (e.g. examples/data/...).
        # It is relative to the invocation working directory, NOT relative to
        # the TR.DAT's parent.  Only *referenced child* filenames are resolved
        # relative to their referring file.
        if parent is None or candidate.is_absolute():
            return candidate.resolve()
        roots = []
        if parent:
            roots.append(Path(parent).parent)
        if self.trdat:
            roots.append(self.trdat.resolve().parent)
        roots.append(Path.cwd())
        for root in roots:
            path = (root / candidate).resolve()
            if path.is_file():
                return path
        return (roots[0] / candidate).resolve()

    def _collect(self):
        if self.trdat is None:
            return
        # Each tuple: (role, key, filename, referring file, nesting depth,
        #              inherited model). Inherited ownership matters for
        # #include and generic nested-input keys inside model-specific files.
        pending = [('TRANSP_INPUT', 'TRANSP_NAMELIST', str(self.trdat), None, 0, '')]
        adjacent_pending = []
        # Discover adjunct files only when the directory is dedicated to this
        # run. They may not be referenced by TR.DAT: archive them, but label
        # their association as unverified instead of implying code execution.
        # This captures TORAY/GENRAY *.in, ECH namelists, PT_SOLVER/MMM inputs,
        # including orphan inputs when the run's TR.DAT was not supplied.
        if self.runid:
            directory = self.trdat.resolve().parent
            dedicated = directory.name.upper() == self.runid.upper()
            if dedicated:
                adjacent = sorted(path for path in directory.iterdir()
                                  if path.is_file() and path.suffix.lower()
                                  in ('.nml', '.namelist', '.inp', '.in'))
            else:
                adjacent = sorted(path for pattern in (f'{self.runid}*.nml',
                                       f'{self.runid}*.NML')
                                  for path in directory.glob(pattern))
            main_text, _, main_status = read_input(self.trdat)
            main_waves = (set(model_hints_from_text(main_text)) & set(WAVE_MODELS)
                          if main_status == 'embedded' else set())
            for path in adjacent:
                owner = file_model('AUTO_DISCOVERED', path.name)
                # Generic 'ech_inp.nml' in a dedicated run directory has one
                # possible owner only when TR.DAT explicitly selects one EC code.
                if not owner and path.name.lower().startswith('ech_') and len(main_waves) == 1:
                    owner = next(iter(main_waves))
                role = 'ADJACENT_UNREFERENCED' if dedicated else 'AUX_NAMELIST'
                adjacent_pending.append((role, 'AUTO_DISCOVERED', str(path),
                                         str(self.trdat), 1, owner))
        seen = set()
        i = 0
        while i < len(pending) and len(self.records) < MAX_INPUT_FILES:
            role, key, filename, parent, depth, inherited_model = pending[i]
            i += 1
            path = self._resolve(filename, parent)
            # One entry per physical file avoids cyclic inclusion or redundant
            # copies; its earliest parent reference provides the provenance.
            if path in seen:
                continue
            seen.add(path)
            text, digest, status = read_input(path)
            model = file_model(key, filename) or inherited_model
            self.records.append(InputRecord(role, key, filename, path, depth,
                                            status, text, digest, model))
            if status != 'embedded':
                if role == 'TRANSP_INPUT':
                    pending.extend(adjacent_pending)
                continue
            refs = list(references_from_text(text))
            if depth >= MAX_REFERENCE_DEPTH:
                if refs:
                    self.truncated = True
                continue
            for ref_key, ref_file in refs:
                if role == 'TRANSP_INPUT':
                    child_role = 'TEMPLATE' if 'TEMPLATE' in ref_key.upper() else 'REFERENCED_NAMELIST'
                else:
                    child_role = 'NESTED_NAMELIST'
                pending.append((child_role, ref_key, ref_file, str(path), depth+1,
                                model))
            if role == 'TRANSP_INPUT':
                # Explicit references take precedence over adjacent-only files.
                pending.extend(adjacent_pending)
        if i < len(pending):
            self.truncated = True

    @property
    def missing(self):
        return [item for item in self.records if item.status != 'embedded']

    def models(self, family=None):
        models = set()
        main = self.records[0] if self.records else None
        if main and main.role == 'TRANSP_INPUT' and main.status == 'embedded':
            models.update(model_hints_from_text(main.text))
        for record in self.records:
            if record.model:
                models.add(record.model)
            if record.status == 'embedded' and record.role != 'TRANSP_INPUT':
                models.update(m for m in model_hints_from_text(record.text)
                              if m in WAVE_MODELS or m in TRANSPORT_MODELS)
        if family:
            models.intersection_update(MODEL_FAMILIES.get(family, ()))
        return sorted(models)

    def scoped_xml(self, model):
        """Full TR.DAT plus model-specific referenced files (including missing)."""
        root = ET.Element('transp_model_inputs', {
            'format_version': '1', 'model': str(model),
            'note': 'source inputs only; execution of model not independently verified',
        })
        scoped_records = [record for record in self.records
                          if record.role == 'TRANSP_INPUT' or record.model == model]
        files = ET.SubElement(root, 'input_files', {
            'complete': str(not (any(r.status != 'embedded' for r in scoped_records)
                                 or self.truncated)).lower(),
            'discovery_truncated': str(self.truncated).lower(),
        })
        for record in scoped_records:
            append_file_xml(files, record)
        return ET.tostring(root, encoding='unicode')


def fill_model_libraries(code, inventory, family):
    """Populate code.library parameters when supported by the target DD."""
    models = inventory.models(family)
    if not models or not hasattr(code, 'library'):
        return []
    # Common IMAS code structure limits libraries to 10 for MDS+.
    models = models[:10]
    code.library.resize(len(models))
    for i, name in enumerate(models):
        lib = code.library[i]
        if hasattr(lib, 'name'):
            lib.name = name
        if hasattr(lib, 'description'):
            lib.description = ('Referenced TRANSP external-model inputs; '
                               'model execution not independently verified')
        if hasattr(lib, 'parameters'):
            lib.parameters = inventory.scoped_xml(name)
    return models
