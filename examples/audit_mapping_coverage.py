#!/usr/bin/python
"""Audit TRANSP CDF variable coverage against the current transp2imas source tree.

The converter is intentionally extended group-by-group.  This tool provides a
repeatable inventory for that process by comparing every variable in a TRANSP
CDF with static references found in the converter, validation scripts, and a
small maintained registry of quantities that are deliberately deferred.

The result is a *source coverage* audit, not a proof that an IMAS mapping is
physically complete.  A variable can be referenced by the converter as a grid,
normalization, closure input, or feature-detection flag without being copied to
an IMAS field one-to-one.  Conversely, dynamically constructed variable names
can occasionally escape static discovery.  Use the report to choose the next
mapping group, then confirm semantics, units, grids, signs, and Data Dictionary
compatibility as described in docs/mapping-workflow.md.

Typical use::

    python examples/audit_mapping_coverage.py \
        examples/data/nstx/141711Y01/141711Y01.CDF

Write full machine-readable and Markdown inventories::

    python examples/audit_mapping_coverage.py run.CDF \
        --csv build/audit/run.csv \
        --markdown build/audit/run.md \
        --json build/audit/run.json
"""

from __future__ import annotations

import argparse
import ast
import csv
import fnmatch
import json
import re
import sys
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

from scipy.io import netcdf_file

# Dynamically generated NC/neutral channel names will not appear in the AST
# as exact source strings. Keep the static inventory consistent with the
# actual group-discovery logic, without implying a validated IMAS write.
REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))
from theory_transport_mapping import (
    discover_neoclassical_channels, TYPED_FIELDS, AUXILIARY_SPECS,
    GLF23_MODEL_SPEC, GLF23_AUXILIARY_SPECS,
)
from neutral_profiles_mapping import discover_neutral_channels


# CDF quantities that are intentionally used primarily to construct or
# normalize other mapped quantities.  Keeping these explicit avoids presenting
# every helper/grid variable as a direct one-to-one mapping.
DERIVED_INPUTS: dict[str, str] = {
    "TIME": "TRANSP time base used by actuator/hardware signals",
    "TIME3": "TRANSP profile time base",
    "X": "zone-centre normalized toroidal-flux coordinate",
    "XB": "zone-boundary normalized toroidal-flux coordinate",
    "DVOL": "zone volume used for cumulative/integrated quantities",
    "DAREA": "zone area used for cumulative/enclosed quantities",
    "SURF": "flux-surface area used to convert transport divergence to flux",
    "BZ": "vacuum toroidal-field normalization used in current conversions",
    "BZXR": "B*R normalization used by equilibrium/global quantities",
    "CONPLJB": "geometry factor used to convert toroidal current to <J.B>/B0",
    "RMAJM": "auxiliary radial grid used by selected TRANSP profiles",
    "RMJMP": "major-radius profile used by equilibrium geometry",
}

# Known quantities that have been reviewed but are deliberately not mapped to a
# current IDS field.  Add entries here only after the reason is understood and
# documented.  Pattern entries use shell-style wildcards so species/indexed
# families can be represented without enumerating each CDF variable.
DEFERRED_PATTERNS: dict[str, str] = {
    "T0RECYC": "neutral recycling source-state temperature; no audited IMAS destination yet",
    "T0GASFL": "gas source-state temperature; no audited IMAS destination yet",
    "OM0RECYC": "neutral recycling source-state rotation; no audited IMAS destination yet",
    "OM0GASFL": "gas source-state rotation; no audited IMAS destination yet",
    "E0INR": "neutral/source boundary-state quantity; semantics need a destination audit",
    "PEEDG": "electron energy crossing the plasma boundary, not a volumetric source",
    "PIEDG": "ion energy crossing the plasma boundary, not a volumetric source",
    "VPOLD_NC": "neoclassical poloidal velocity on the RMAJM grid; target/grid treatment deferred",
    "VTORD_NC": "neoclassical toroidal velocity on the RMAJM grid; target/grid treatment deferred",
    "GFLNC_I": "aggregate ion particle divergence lacks typed total-ion particle-flux destination",
    "GFLNC_X": "aggregate impurity particle divergence; do not invent an additional ion species",
    "GFLNC_TOK": "aggregate summed particle divergence; not an individual ion species",
    "QFLNC_X": "aggregate impurity heat divergence; not an additional ion species",
    "QFLNC_TOK": "aggregate heat divergence; not an additional ion species",
    "QFLNCC_X": "aggregate impurity NC-classical heat divergence; not an additional ion species",
    "QFLNCC_TOK": "aggregate NC-classical heat divergence; not an additional ion species",
    "PIC_NPHI": "aggregate IC n_phi power spectrum is not frequency-resolved for multi-wave mapping",
    "PICF??N??": "frequency/n_phi-resolved IC profile family; component-index semantics need a dedicated spectral mapping audit",
}

METADATA_MODULES = {
    "imas_dataset_fair.py",
    "imas_workflow.py",
}

# Files under examples that are intended to validate or inventory source data.
# inspect/plot examples are excluded because merely displaying a CDF variable
# should not make it look like part of converter validation coverage.
VALIDATION_PREFIXES = ("verify_",)


@dataclass(frozen=True)
class CoverageRow:
    """One CDF variable and its static transp2imas coverage classification."""

    name: str
    status: str
    dimensions: str
    shape: str
    units: str
    description: str
    source_files: str
    note: str


def _decode(value: object) -> str:
    """Return a printable NetCDF attribute value without byte-string noise."""

    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("utf-8", "replace").strip()
    if isinstance(value, (list, tuple)):
        return ", ".join(_decode(item) for item in value)
    return str(value).strip()


def _joined_string_pattern(node: ast.JoinedStr) -> str:
    """Convert an f-string AST node to a shell-style wildcard pattern."""

    parts: list[str] = []
    for value in node.values:
        if isinstance(value, ast.Constant) and isinstance(value.value, str):
            parts.append(value.value)
        else:
            parts.append("*")
    return "".join(parts)


def _string_patterns(path: Path) -> set[str]:
    """Extract exact/wildcard/regex-like string patterns from one Python file."""

    try:
        source = path.read_text(encoding="utf-8", errors="replace")
        tree = ast.parse(source, filename=str(path))
    except (OSError, SyntaxError):
        return set()

    patterns: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            text = node.value.strip()
            if text:
                patterns.add(text)
        elif isinstance(node, ast.JoinedStr):
            pattern = _joined_string_pattern(node).strip()
            if pattern:
                patterns.add(pattern)
    return patterns


def _looks_like_regex(pattern: str) -> bool:
    """Identify source strings that are plausibly regular expressions."""

    if len(pattern) > 160 or any(ch.isspace() for ch in pattern):
        return False
    return any(token in pattern for token in (r"\d", "(?:", "[", "^", "$", ".*", ".+"))


def _pattern_matches(pattern: str, variable: str) -> bool:
    """Return whether a source string plausibly names a CDF variable."""

    if pattern == variable:
        return True

    # TRANSP helper calls commonly use '?' to represent indexed/species names;
    # f-strings are converted to '*' by _joined_string_pattern().  Only treat a
    # wildcard string as a CDF-name pattern when it contains a meaningful
    # uppercase static token.  This prevents generic f-strings such as f"{x}"
    # (which reduce to "*") from matching every CDF variable.
    if "?" in pattern or "*" in pattern:
        if not re.search(r"[A-Z]{2,}", pattern) or any(ch.isspace() for ch in pattern):
            return False
        shell_pattern = pattern.replace("?", "*")
        if fnmatch.fnmatchcase(variable, shell_pattern):
            return True

    if _looks_like_regex(pattern):
        try:
            return re.fullmatch(pattern, variable) is not None
        except re.error:
            return False
    return False


def _references_for_files(files: Iterable[Path], variable_names: list[str]) -> dict[str, set[str]]:
    """Map CDF variable names to source files that statically reference them."""

    refs: dict[str, set[str]] = defaultdict(set)
    for path in files:
        patterns = _string_patterns(path)
        if not patterns:
            continue
        for variable in variable_names:
            if any(_pattern_matches(pattern, variable) for pattern in patterns):
                refs[variable].add(path.name)
    return refs


def _dynamic_source_references(variable_names: list[str]) -> dict[str, set[str]]:
    """Recognize variable families explicitly enumerated by converter logic.

    This tracks *source consumption* only, not a one-to-one typed IMAS
    mapping. Aggregate particle/impurity channels intentionally remain
    deferred (the writer does not fabricate aggregate ion species).
    """
    references: dict[str, set[str]] = defaultdict(set)
    channels = discover_neoclassical_channels(variable_names)
    for kind, found in channels.items():
        for token, name in found.items():
            if token in ("X", "TOK") or (kind == "particles" and token == "I"):
                continue
            references[name].add("imas_core_transport.py")
    for population in discover_neutral_channels(variable_names):
        for field in ("density", "temperature", "rotation"):
            source = population[field]
            if source:
                references[source].add("imas_core_profiles.py")
    theory_names = (set(TYPED_FIELDS) | set(AUXILIARY_SPECS)
                    | set(GLF23_MODEL_SPEC['coefficients'])
                    | set(GLF23_AUXILIARY_SPECS))
    for name in theory_names & set(variable_names):
        references[name].add('imas_core_transport.py')
    return references


def _deferred_reason(name: str) -> str:
    """Return the maintained deferral reason for a CDF variable, if any."""

    for pattern, reason in DEFERRED_PATTERNS.items():
        if fnmatch.fnmatchcase(name, pattern):
            return reason
    return ""


def _family(name: str) -> str:
    """Return a compact leading-name family for unclassified summary counts."""

    match = re.match(r"[A-Z]+", name.upper())
    return match.group(0) if match else name


def _classify(
    name: str,
    converter_refs: set[str],
    metadata_refs: set[str],
    validation_refs: set[str],
) -> tuple[str, str]:
    """Assign one stable coverage status and explanatory note."""

    deferred = _deferred_reason(name)
    non_metadata_converter_refs = converter_refs - METADATA_MODULES

    if non_metadata_converter_refs:
        if name in DERIVED_INPUTS:
            return "support/derived-input", DERIVED_INPUTS[name]
        return "mapped/consumed", "referenced by converter source"

    if metadata_refs:
        return "metadata/provenance", "referenced only by metadata/provenance writers"

    if deferred:
        return "deferred", deferred

    if validation_refs:
        return "validation-only", "referenced by source-side validation/inventory code only"

    return "unclassified", "no static converter, validation, or maintained deferral reference found"


def _build_rows(cdf: Path, repo_root: Path) -> list[CoverageRow]:
    """Read the CDF and build the full coverage table."""

    with netcdf_file(str(cdf), "r", mmap=False) as dataset:
        variable_names = sorted(dataset.variables)

        root_python = sorted(path for path in repo_root.glob("*.py") if path.is_file())
        metadata_python = [path for path in root_python if path.name in METADATA_MODULES]
        # Source-only contract helper contains generic regex fragments such as
        # "[A-Z0-9]+". Those are NOT CDF-name mappings and must not let the
        # AST matcher classify every variable as consumed. The actual
        # deliberately supported families are enumerated by dynamic_map.
        converter_python = [path for path in root_python
                            if path.name not in ('theory_transport_mapping.py', 'neutral_profiles_mapping.py')]
        validation_python = sorted(
            path
            for path in (repo_root / "examples").glob("*.py")
            if path.is_file() and path.name.startswith(VALIDATION_PREFIXES)
        )

        converter_map = _references_for_files(converter_python, variable_names)
        metadata_map = _references_for_files(metadata_python, variable_names)
        validation_map = _references_for_files(validation_python, variable_names)
        dynamic_map = _dynamic_source_references(variable_names)

        rows: list[CoverageRow] = []
        for name in variable_names:
            variable = dataset.variables[name]
            converter_refs = converter_map.get(name, set()) | dynamic_map.get(name, set())
            metadata_refs = metadata_map.get(name, set())
            validation_refs = validation_map.get(name, set())
            status, note = _classify(name, converter_refs, metadata_refs, validation_refs)
            source_files = sorted(converter_refs | metadata_refs | validation_refs)

            rows.append(
                CoverageRow(
                    name=name,
                    status=status,
                    dimensions=",".join(variable.dimensions),
                    shape="x".join(str(value) for value in variable.shape),
                    units=_decode(getattr(variable, "units", "")),
                    description=_decode(
                        getattr(variable, "long_name", getattr(variable, "description", ""))
                    ),
                    source_files=",".join(source_files),
                    note=note,
                )
            )
    return rows


def _summary(rows: list[CoverageRow]) -> Counter[str]:
    return Counter(row.status for row in rows)


def _print_text(cdf: Path, repo_root: Path, rows: list[CoverageRow], max_list: int) -> None:
    """Print a concise terminal report and the highest-value unclassified families."""

    counts = _summary(rows)
    total = len(rows)
    covered = sum(
        counts[name]
        for name in ("mapped/consumed", "support/derived-input", "metadata/provenance")
    )

    print("TRANSP mapping coverage audit")
    print(f"CDF:       {cdf}")
    print(f"repo root: {repo_root}")
    print(f"variables: {total}")
    print()
    for status in (
        "mapped/consumed",
        "support/derived-input",
        "metadata/provenance",
        "validation-only",
        "deferred",
        "unclassified",
    ):
        print(f"  {status:23s} {counts[status]:5d}")
    percent = (100.0 * covered / total) if total else 0.0
    print(f"  {'converter-referenced total':23s} {covered:5d}  ({percent:5.1f}%)")
    print()

    unclassified = [row for row in rows if row.status == "unclassified"]
    family_counts = Counter(_family(row.name) for row in unclassified)
    if family_counts:
        print("Largest unclassified variable families")
        for family, count in family_counts.most_common(20):
            print(f"  {family:18s} {count:5d}")
        print()

    deferred = [row for row in rows if row.status == "deferred"]
    if deferred:
        print("Known deferred variables")
        for row in deferred[:max_list]:
            print(f"  {row.name:18s} {row.note}")
        if len(deferred) > max_list:
            print(f"  ... {len(deferred) - max_list} more")
        print()

    if unclassified:
        print(f"First {min(max_list, len(unclassified))} unclassified variables")
        for row in unclassified[:max_list]:
            details = "; ".join(value for value in (row.units, row.description) if value)
            print(f"  {row.name:18s} ({row.dimensions}) {details}")
        if len(unclassified) > max_list:
            print(f"  ... {len(unclassified) - max_list} more")


def _write_csv(path: Path, rows: list[CoverageRow]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(asdict(rows[0]).keys()) if rows else [])
        if rows:
            writer.writeheader()
            writer.writerows(asdict(row) for row in rows)


def _markdown_escape(text: str) -> str:
    return text.replace("|", r"\|").replace("\n", " ")


def _write_markdown(path: Path, cdf: Path, repo_root: Path, rows: list[CoverageRow]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    counts = _summary(rows)
    total = len(rows)
    covered = sum(
        counts[name]
        for name in ("mapped/consumed", "support/derived-input", "metadata/provenance")
    )
    percent = (100.0 * covered / total) if total else 0.0
    family_counts = Counter(_family(row.name) for row in rows if row.status == "unclassified")

    lines = [
        "# TRANSP mapping coverage audit",
        "",
        f"- CDF: `{cdf}`",
        f"- Repository root: `{repo_root}`",
        f"- CDF variables: **{total}**",
        f"- Converter-referenced variables: **{covered} ({percent:.1f}%)**",
        "",
        "> This is a static source-coverage report, not proof of semantic completeness. ",
        "> Every new mapping still requires a physics/units/grid/sign/DD audit.",
        "",
        "## Summary",
        "",
        "| Status | Count |",
        "|---|---:|",
    ]
    for status in (
        "mapped/consumed",
        "support/derived-input",
        "metadata/provenance",
        "validation-only",
        "deferred",
        "unclassified",
    ):
        lines.append(f"| `{status}` | {counts[status]} |")

    lines.extend(["", "## Largest unclassified families", "", "| Family | Count |", "|---|---:|"])
    if family_counts:
        for family, count in family_counts.most_common(30):
            lines.append(f"| `{family}` | {count} |")
    else:
        lines.append("| _none_ | 0 |")

    lines.extend(
        [
            "",
            "## Full inventory",
            "",
            "| Variable | Status | Dimensions | Shape | Units | Description | Source files | Note |",
            "|---|---|---|---|---|---|---|---|",
        ]
    )
    for row in rows:
        values = [
            f"`{row.name}`",
            f"`{row.status}`",
            row.dimensions,
            row.shape,
            row.units,
            row.description,
            row.source_files,
            row.note,
        ]
        lines.append("| " + " | ".join(_markdown_escape(value) for value in values) + " |")

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_json(path: Path, cdf: Path, repo_root: Path, rows: list[CoverageRow]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    counts = _summary(rows)
    payload = {
        "cdf": str(cdf),
        "repo_root": str(repo_root),
        "summary": dict(sorted(counts.items())),
        "variables": [asdict(row) for row in rows],
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cdf", type=Path, help="TRANSP CDF file to audit")
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="transp2imas source root (default: inferred from this script)",
    )
    parser.add_argument("--csv", type=Path, help="write the full inventory as CSV")
    parser.add_argument("--markdown", type=Path, help="write the full inventory as Markdown")
    parser.add_argument("--json", type=Path, help="write the full inventory as JSON")
    parser.add_argument(
        "--max-list",
        type=int,
        default=40,
        help="maximum deferred/unclassified rows printed to the terminal (default: %(default)s)",
    )
    parser.add_argument(
        "--fail-on-unclassified",
        action="store_true",
        help="exit with status 2 if any CDF variables remain unclassified",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    cdf = args.cdf.resolve()
    repo_root = args.repo_root.resolve()

    if not cdf.is_file():
        raise SystemExit(f"CDF not found: {cdf}")
    if not (repo_root / "nc2imas.py").is_file():
        raise SystemExit(f"not a transp2imas source root: {repo_root}")

    rows = _build_rows(cdf, repo_root)
    _print_text(cdf, repo_root, rows, max(0, args.max_list))

    if args.csv:
        _write_csv(args.csv, rows)
        print(f"wrote CSV:      {args.csv}")
    if args.markdown:
        _write_markdown(args.markdown, cdf, repo_root, rows)
        print(f"wrote Markdown: {args.markdown}")
    if args.json:
        _write_json(args.json, cdf, repo_root, rows)
        print(f"wrote JSON:     {args.json}")

    if args.fail_on_unclassified and any(row.status == "unclassified" for row in rows):
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
