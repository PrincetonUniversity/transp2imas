# Schemas, migration and report changes

## Pipeline

```text
strict JSON decoding
  -> declared-version structural validation (JSON Schema, Draft 2020-12, offline)
  -> explicit migration and canonical normalization (to inputs.v3)
  -> resource/reference validation (files exist, hashes match)
  -> semantic validation (ordering, graph, units, numerical constraints)
  -> scientific rules and application policy
```

Each stage is separate. A schema-valid sidecar is **not** authenticated evidence, and a semantic error is never
labelled a schema failure (the message has no `[schema ...]` prefix).

* **Decoding** (`curation.schema.loads_strict`) rejects duplicate keys, `NaN`/`Infinity`, and numeric overflow such as
  `1e999` (Python would otherwise decode it to infinity). Documents built in code are scanned for nonfinite scalars
  before validation.
* **Structure** comes from the schemas in `curation/schemas/`: `evidence.v1`, `evidence.v2`, `evidence.v3`, `pilot.v1`,
  `policy.v1`, `report.v3`, with shared definitions in `common`. They are resolved through a local `referencing`
  registry; nothing is fetched. Unknown scientific keys are rejected. A deliberate extension point exists:
  `extensions` maps namespaced names (`vendor.rule`) to objects.
* **Migration** is by the module, never by adding fields to an old document: a v1 document may not contain v2 or v3
  fields, and a v2 document may not contain v3 fields or a migration record. The canonical form is always `inputs.v3`
  with a `migration` record (`from`, `raw_sha256`, `steps`). Normalization is idempotent and never mutates the
  caller's object; validating a canonical document returns it unchanged.
* **Resources**: relative paths in `resolved_inputs.inputs[].path` and `uncertainty.covariances[].path` resolve against
  the directory of the sidecar (not the working directory). The declared `path` is kept and `resolved_path` is added.
  A covariance resource must exist and match its SHA-256, otherwise loading fails before any assessment. Referenced
  resources are protected from being overwritten by the report, like the top-level evidence file.
* **Error text** names the schema, the JSON path, the offending value and its type, and the failed condition, for
  example: `[schema evidence.v3] events[0].start: 'soon' is not of type 'number' (value 'soon' of type str; failed 'type': "number")`.

### Validating without a CDF

```bash
python curate.py --validate-only --evidence my-evidence.json [--policy my-policy.json]
```

This runs decoding, structure, migration, resource and semantic validation and prints the declared and canonical
schema and the migration steps. It states what it did **not** check: coordinates, units, grid compatibility,
variable names and covariance shapes need a source and are only checked in a real run.

Editors can use the schemas directly: the example sidecars carry a `"$schema"` pointing at
`curation/schemas/*.schema.json`; that key is an editor hint and is dropped from the canonical form.

### Install

```bash
pip install -r requirements-curation.txt     # numpy, scipy, f90nml>=1.4, jsonschema>=4.18, referencing; no IMAS
```

## Versions

| Document | Current | Still accepted | What changed |
|---|---|---|---|
| Evidence sidecar | `inputs.v3` | `v1`, `v2` | v3 adds `sample_support`, `event_detectors`, `uncertainty`, `extensions`, the `migration` record, the `unclassified_transient` event type and event fields `origin`, `reference`, `diagnostic`, `detector`, `timing_uncertainty`, `padding`, `review_state`, `policy`, `applications`; evidence records and sources gain `selector` |
| Report | `report.v3` | (v2 is no longer written) | adds `execution_ledger`, `policy`, `registry`, `fingerprints`, `detected_events`, `uncertainty`, `summary.rule_coverage`, `summary.excused_findings`, `policy_exclusions`; each finding carries `rule` and `rule_version` |
| Policy | `policy.v1` | | new |
| Pilot facts | `pilot.v1` | | now validated against `pilot.v1`; produces a v3 manifest |

## Intentional behaviour differences

These are changes, not extractions; each has a regression test in `examples/verify_curation_gaps.py` (P0.x).

1. **v1 sidecars work through the CLI** (P0.1). Previously the engine revalidated the migrated document under the v1
   rule and exited 2. The raw document is validated once, as v1; the canonical document is v3.
2. **Resolved controls are compared over the whole support** (P0.2). A runtime change inside a declared interval is a
   mismatch on that part; gaps and unknown initial states are `not_assessed`; changes outside the assessed support
   (assessment window intersected with the native clocks) do not matter. `evidence.intervals` now lists atomic
   half-open pieces rather than one entry per declared state.
3. **Evidence time applicability** (P0.3). Generic evidence is covered only where its `time_min`/`time_max` meets the
   native support of its quantity (profiles on TIME3, heating on TIME). `evidence.coverage` records `full`,
   `partial` (with `unassessed_intervals`), `disjoint` or `unknown`; disjoint and unknown are `not_assessed`, partial
   downgrades a pass to a warning. A declared run identity that contradicts the CDF is not trusted, and the CDF value
   is not substituted: the identity item becomes unknown, so no reference is authorized by it.
4. **Lineage across aliases** (P0.4). Shared resources are detected throughout both ancestries (SHA-256 or source string
   on any pair of ancestors). With declared `selector`s the outcome is `same_resource` (warning), a shared container
   with different selectors (`pass` with an explicit "not proof of independence" message), or
   `same_resource_selector_unknown` (warning). Usage is scoped in time; ancestry is kept.
5. **Events no longer rewrite findings** (P0.5). Previously a declared event turned the negative-sample finding into a
   `pass`. Now the raw finding is unchanged and the exclusion is recorded in `policy_exclusions` with counts, masks
   (`*__excluded`) and origins. Source-integrity rules (`GRID.*`, grid mismatch, units, read errors, run identity) can
   never be touched by an event, including `rules: ["*"]`. Existing sidecar `events[].rules` migrate to
   `policy: {"kind": "legacy_exclusion"}` and origin `analyst_declaration`; to keep a v1/v2 expectation that an
   event makes a rule pass you must read `summary.excused_findings` instead.
6. **Unresolved array updates invalidate prior elements** (P0.6). See [namelist-resolution.md](namelist-resolution.md).
7. **Provenance keys**: `source.curation_module_sha256` now keys by repository-relative path
   (`curation/engine.py`), includes subpackages, and `source.code_identity` is added.
8. **Exit codes**: a rule that raises makes the CLI exit 2 after writing the report; `--require-complete` also exits 1
   when any rule did not execute.

## Migration instructions

* Existing v1/v2 sidecars need no change. Re-run to get a v3 report; read `effective_manifest.migration` for what was
  added.
* If you consumed `summary.applications` after relying on an event to make a finding pass, also read
  `summary.excused_findings` and `policy_exclusions`.
* If you parsed `source.curation_module_sha256` by file name, use the relative path keys.
* To pin a policy, keep the policy file next to the sidecar and compare `fingerprints.assessment_sha256`.
