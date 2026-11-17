#!/usr/bin/python
"""Provenance of a conversion: file records, ``code.parameters`` XML and the IDS properties text.

Everything here takes plain values (paths, a ``CocosSigns``, an ``InputInventory``) rather than the
converter object, so it can be tested without a database or a CDF.  The input files themselves are
discovered and embedded by ``input_provenance``; this module records what the conversion used and how
it was configured.

* ``source_file_records``: the short, portable list for ``ids_properties/provenance`` (the DD caps it
  at 10): name with a SHA-256 of the file, and the file's modification time,
* ``attach_source_provenance``: write those records into the IDS in the layout of its DD (``reference``
  from DD 3.42, ``sources`` before),
* ``build_code_parameters_xml``: the XML for ``code.parameters`` (run, converter, TRANSP metadata,
  coordinate conventions, embedded inputs),
* ``properties_comment``: the ``ids_properties/comment`` text.
"""

from datetime import datetime, timezone
import hashlib
from pathlib import Path
import xml.etree.ElementTree as ET

from input_provenance import append_file_xml

PROVENANCE_RECORD_CAP = 10        # the DD caps ids_properties/provenance/node/reference at 10
TRANSP_METADATA_KEYS = ('Runid', 'CDF_date', 'BUILD_date', 'TRANSP_version', 'TRANSP_DOI',
                        'TGLF_hash', 'GLF23_hash')


def file_timestamp(path):
    """UTC ISO-8601 modification time of a local file, or '' when it cannot be read."""
    try:
        return datetime.fromtimestamp(Path(path).stat().st_mtime, tz=timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
    except (OSError, ValueError, OverflowError):
        return ''


def file_sha256(path):
    """SHA-256 of a file, or '' when it cannot be read (the record then carries no digest)."""
    digest = hashlib.sha256()
    try:
        with open(path, 'rb') as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b''):
                digest.update(chunk)
    except OSError:
        return ''
    return digest.hexdigest()


def source_file_records(candidates, cap=PROVENANCE_RECORD_CAP):
    """``[(name, timestamp)]`` for ``(role, filename)`` candidates, unique by resolved path.

    The name is ``role:file:<basename>#sha256=<digest>``: no machine-specific absolute path (the exact
    resolved path stays in ``code.parameters``).  Empty file names are skipped, and at most ``cap``
    records are returned.
    """
    records, seen = [], set()
    for role, filename in candidates:
        if not filename:
            continue
        path = Path(filename).expanduser()
        try:
            resolved = path.resolve()
        except (OSError, RuntimeError):
            resolved = path
        if str(resolved) in seen:
            continue
        seen.add(str(resolved))
        digest = file_sha256(resolved)
        name = f'{role}:file:{Path(filename).name}'
        if digest:
            name += f'#sha256={digest}'
        records.append((name, file_timestamp(resolved)))
        if len(records) >= cap:
            break
    return records


def attach_source_provenance(ids, records):
    """Write ``records`` into ``ids_properties/provenance`` in the layout of the IDS's DD.

    DD >= 3.42 has the ``reference`` array of name/timestamp; earlier DDs have ``sources`` (names).
    Raises whatever the IDS raises: the caller decides that provenance is best-effort and records it.
    """
    provenance = ids.ids_properties.provenance
    provenance.node.resize(1)
    node = provenance.node[0]
    if hasattr(node, 'path'):
        node.path = ''
    if hasattr(node, 'reference'):
        node.reference.resize(len(records))
        for i, (name, timestamp) in enumerate(records):
            ref = node.reference[i]
            if hasattr(ref, 'name'):
                ref.name = name
            if timestamp and hasattr(ref, 'timestamp'):
                ref.timestamp = timestamp
    elif hasattr(node, 'sources'):
        node.sources = [name for name, _ in records]
    return ids


def properties_comment(runid, device, year, cdf_date, build_date, tglf, glf23, dd_version, signs):
    """``ids_properties/comment``: the run, build hashes, target DD and the sign policy."""
    comment = f'TRANSP results for run {runid} {device} year: {year} CDF date: {cdf_date}'
    return (comment + f' build: {build_date} tglf: {tglf} glf23: {glf23}'
            + f' DD={dd_version}; COCOS={signs.target_cocos};'
              f' NLJCCW={signs.nljccw}; NLBCCW={signs.nlbccw};'
              f' Ip_sign={signs.ip_sign:+.0f}; Bphi_sign={signs.bt_sign:+.0f};'
              f' q_sign={signs.q_sign:+.0f}; psi_sign={signs.psi_sign:+.0f}')


def build_code_parameters_xml(*, dd_version, signs, device, shot, imas_run, occurrence, source_cdf,
                              converter_version, metadata, inventory, source_format='netcdf', source_input=None):
    """XML for ``code.parameters``: complete TRANSP inputs verbatim (escaped) with hashes and CDF metadata.

    ``metadata`` maps TRANSP attribute names to values (bytes are decoded, ``None`` is skipped);
    ``inventory`` is the run's ``InputInventory``.
    """
    root = ET.Element('transp2imas_parameters', {
        'format_version': '1',
        'target_dd': str(dd_version),
        'target_cocos': str(signs.target_cocos),
    })
    run = ET.SubElement(root, 'run', {
        'device': str(device), 'shot': str(shot), 'imas_run': str(imas_run), 'occurrence': str(occurrence),
    })
    if source_format == 'netcdf':
        run.set('source_cdf', str(source_cdf))
    else:
        run.set('source_format', str(source_format))
        run.set('source_input', str(source_input or source_cdf))
    ET.SubElement(root, 'converter', {'name': 'transp2imas', 'version': converter_version})

    block = ET.SubElement(root, 'transp_metadata')
    for key in TRANSP_METADATA_KEYS:
        value = metadata.get(key)
        if value is None:
            continue
        if isinstance(value, bytes):
            value = value.decode('ISO-8859-1', errors='replace')
        ET.SubElement(block, 'item', {'name': key, 'value': str(value).strip()})

    conventions = ET.SubElement(root, 'coordinate_conventions')
    for key, value in (('NLJCCW', signs.nljccw), ('NLBCCW', signs.nlbccw),
                       ('ip_sign', signs.ip_sign), ('bphi_sign', signs.bt_sign),
                       ('q_sign', signs.q_sign), ('psi_sign', signs.psi_sign),
                       ('phi_sign', signs.phi_sign)):
        ET.SubElement(conventions, 'item', {'name': key, 'value': str(value)})

    inputs = ET.SubElement(root, 'input_files', {
        'complete': str(not (inventory.missing or inventory.truncated)).lower(),
        'discovery_truncated': str(inventory.truncated).lower(),
    })
    for record in inventory.records:
        append_file_xml(inputs, record)
    return ET.tostring(root, encoding='unicode')
