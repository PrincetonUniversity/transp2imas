# Curation architecture: rules, policy and provenance

Curation is a standalone, read-only package (`curation/`). It never imports IMAS, never writes a source file, and
never claims scientific certification. This page describes how a run is organised and where each kind of decision
lives. For the behaviour changes that came with it, see [schema-migration.md](schema-migration.md).

## Layers

```text
curate.py                CLI: reads inputs, builds policy/registry, writes the report, chooses the exit code
curation/engine.py       orchestration only: plan -> execute -> ledger -> fingerprints
curation/registry.py     Rule, RuleRegistry: explicit registration, dependencies, deterministic plan, plugins
curation/policy.py       Policy: which rules run, screening parameters, applications; mandatory-rule patterns
curation/context.py      CurationContext: the read-only inputs every rule sees
curation/results.py      Finding, Report (+ledger, coverage, fingerprints), Recorder: the write side of rules
curation/rules/*.py      the checks (coordinates, profiles, gradients, evidence, composition, heating,
                         model_consistency, numerics, events)
curation/schema.py       offline JSON Schema validation and strict JSON decoding
curation/schemas/*.json  the schemas (evidence v1/v2/v3, pilot v1, policy v1, report v3)
curation/manifest.py     the evidence pipeline: raw validation -> migration -> semantics -> resources
curation/inputs.py       the declared-input timeline (see namelist-resolution.md)
curation/events.py       event membership and inventory;  curation/detectors.py: the optional candidate detector
curation/uncertainty.py  covariance propagation and sampling (offline, bounded)
curation/provenance.py   hashes of the code, artifacts and dependencies that produced a report
```

## Observation, policy and eligibility are three different things

* A **raw finding** is an observation about the source or the evidence: a value is negative, a sample is missing,
  a clock is reversed, a declared run identity disagrees with the file. Findings carry the rule id and version that
  produced them. Nothing downstream rewrites a finding.
* A **policy** decides what is asked (which rules run, which screening parameters, which applications are reported)
  and which findings *block* which application. Policy never changes a finding's status or message.
* **Eligibility** is a separate result: which evaluated samples an application may use. A declared event can make
  samples ineligible under its declared policy (`report.policy_exclusions`, `*__excluded` masks); the raw finding
  and its `__evaluated`/`__flagged` masks stay exactly as observed.

`Report.application_status` combines them: a finding blocks an application unless a recorded exclusion covers every
one of its flagged samples *for that application*; the finding is then listed under `summary.excused_findings`.

## Rules

A rule is a function `rule(ctx, rec)` plus a `Rule` record: a stable id, a version, a family, prerequisites
(`requires`), whether it is `mandatory`, and the screening parameters it reads. The built-in set is registered
explicitly in `curation/rules/__init__.py`:

| Rule id | Findings | Notes |
|---|---|---|
| `coordinates.grid` | `GRID.*` | mandatory |
| `identity.run` | `RUN.IDENTITY` | mandatory |
| `profiles.screen` | `PROFILE.*`, `FIT.*.WEAK_GRADIENT`, `FIT.*.GRADIENT_UNCERTAINTY` | requires `coordinates.grid` |
| `uncertainty.gradient_covariance` | `FIT.*.GRADIENT_COVARIANCE` | optional; needs a supplied covariance resource |
| `evidence.provenance` | `EVIDENCE.*`, `UNCERTAINTY.*`, `INDEPENDENCE.*` | requires `identity.run` |
| `evidence.rotation` | `ROTATION.AVAILABILITY` | |
| `composition.quasineutrality` | `COMPOSITION.QUASINEUTRALITY` | |
| `model.declared_inputs` | `INPUT.TIMELINE`, `MODEL.*`, `HEATING.FRACTIONS` | |
| `heating.reference` | `HEATING.REFERENCE.*` | requires `model.declared_inputs` |
| `inputs.resolved` | `INPUT.RESOLVED.*` | requires `model.declared_inputs` |
| `numerics.runtime` | `NUMERICS.*` | |
| `events.detect` | `EVENTS.DETECT.*` | optional; needs `event_detectors` configuration |
| `events.inventory` | `EVENTS.INVENTORY` | runs after every other rule |

Planning is deterministic (registration order, prerequisites first, `after_all` rules last). Duplicate ids, unknown
requested rules, unknown prerequisites and cycles are errors, never silent omissions. The extraction from the
former single-class engine kept every finding id; `examples/verify_curation_rules.py` and the finding-by-finding
comparison recorded in `docs/curation/schema-migration.md` protect that.

## The execution ledger

Every rule has a ledger entry whose `state` is one of

| State | Meaning |
|---|---|
| `executed` | the rule ran to completion (its findings may still be fail, warning or not assessed) |
| `policy_disabled` | the policy turned it off; the checks it would have made are **unassessed** |
| `missing_prerequisite` | a prerequisite did not execute, so the rule was not run |
| `execution_error` | the rule raised; the exception type and message are recorded and an `error` finding `RULE.<id>.EXECUTION` is added |

Two consequences follow. An exception is an execution error, never a convenient `not_assessed`, and the CLI then
exits 2 after writing the report. And any rule that did not execute keeps every application out of
`screening_only_not_certified`: a check that did not run is not a check that passed. Required (source-integrity)
and optional rules are counted separately in `summary.rule_coverage`; a policy cannot disable a required rule, and
the attempt is kept in the ledger note.

## Policy documents

```json
{
  "schema_version": "transp2imas.curation.policy.v1",
  "name": "grids-and-profiles-only",
  "rules": {"enable_only": ["coordinates.grid", "profiles.screen"]},
  "parameters": {"weak_gradient_relative": 0.002},
  "applications": ["heat_transport"]
}
```

`rules.disable` maps rule ids to a required reason; `rules.enable_only` lists the optional rules to run (required
rules always run). `parameters` may override only the screening thresholds `quasineutrality_relative`,
`weak_gradient_relative`, `gradient_sigma_multiplier` and `newton_iterations_warning`; they override the sidecar's
values and both are recorded. Everything else is evidence, not policy. Run one with
`python curate.py -i run.CDF --policy policy.json`.

## Plugins

External rules are trusted code. They are registered only by the operator, with `--plugin module` on the command
line (a module exposing `register(registry)`), never by a name inside evidence JSON (a sidecar naming a plugin is a
schema error). The report records each plugin's module name, `__version__`, file SHA-256 and the rule ids it added.
A plugin may read `ctx.manifest['extensions']['vendor.name']`, the namespaced extension object a sidecar may carry.

## Fingerprints and provenance

`report.fingerprints` keeps two identities apart:

* `evidence_sha256`: the scientific content of the canonical sidecar (its migration record excluded), so the same
  content declared as v1, v2 or v3 has one evidence identity;
* `assessment_sha256`: the evidence plus the policy, the effective screening parameters, every rule id and version
  and the plugins. Changing a threshold or disabling a rule changes it although the evidence is unchanged.

`source.code_identity` (written by the CLI, hashed once per invocation) lists the SHA-256 of every Python module of
the package including subpackages, the shared parser modules, the packaged schema files, and the versions of the
dependencies that can change a result.

## What this structure does not do

It does not make any rule's physics more authoritative, add TRANSP defaults, or resolve model-execution semantics.
A rule's `pass` still means only that its stated predicate passed.
