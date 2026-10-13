#!/usr/bin/python
"""Export an input-files XML preview; IMAS workflow wraps it in run metadata."""

import argparse
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from input_provenance import InputInventory, append_file_xml


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('trdat', type=Path, help='source TRANSP TR.DAT')
    parser.add_argument('--runid', help='CDF run identifier; defaults to TR.DAT prefix')
    parser.add_argument('--output-dir', type=Path, default=Path('build/input-xml'))
    args = parser.parse_args()
    runid = args.runid or args.trdat.name.removesuffix('TR.DAT')
    inv = InputInventory(args.trdat, runid)
    root = ET.Element('transp_source_inputs', {'runid': runid, 'format_version': '1'})
    files = ET.SubElement(root, 'input_files', {
        'complete': str(not (inv.missing or inv.truncated)).lower(),
        'discovery_truncated': str(inv.truncated).lower(),
    })
    for record in inv.records:
        append_file_xml(files, record)
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    (out / 'transp-inputs.xml').write_text(
        ET.tostring(root, encoding='unicode'), encoding='utf-8')
    for model in inv.models():
        (out / f'{model.lower()}-inputs.xml').write_text(
            inv.scoped_xml(model), encoding='utf-8')
    print('Input XML preview:', out.resolve())
    for record in inv.records:
        print(f'  {record.key}: {record.reference} [{record.status}]')
    print('Models:', ', '.join(inv.models()) or '(none identified)')
    if inv.missing or inv.truncated:
        print('WARNING: source input inventory is INCOMPLETE', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
