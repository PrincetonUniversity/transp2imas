# TRANSP data curation

The initial curation framework screens a TRANSP CDF and explicitly associated inputs **without
changing them, running TRANSP, writing an IMAS entry, or declaring a run scientifically validated**.
It is independent of the converter's strict/lenient mapping policy. A converter defect and an
uncertain experimental input are different problems.

## Quick start

```bash
python3 examples/verify_curation.py
python3 curate.py -i examples/data/iter/38530P94.CDF \
  --trdat examples/data/iter/38530P94TR.DAT \
  --output-dir build/curation/38530P94
```

Add an analyst-maintained evidence sidecar as it becomes available:

```bash
cp examples/curation/evidence-template.json my-evidence.json
# Edit the copy with actual diagnostic/configuration references and assumptions.
python3 curate.py -i /path/run.CDF --trdat /path/runTR.DAT \
  --evidence my-evidence.json --output-dir build/curation/run
```

There is no automatic search for a TR.DAT or evidence file: association must be explicit. Do not
reuse a template's numerical assumptions for an unrelated discharge. No external reference is
retrieved or authenticated by this first version.

A complete synthetic demonstration needs no experimental data or IMAS installation:

```bash
python3 examples/make_curation_demo.py --output-dir build/curation-demo
python3 curate.py -i build/curation-demo/synthetic.CDF \
  --trdat build/curation-demo/syntheticTR.DAT \
  --evidence build/curation-demo/evidence.json \
  --output-dir build/curation-demo/report
```

## Outputs

| File | Contents |
|---|---|
| `curation.json` | Source/evidence identities, native input timeline, thresholds, all findings, application-specific preliminary statuses, diagnostic inventory |
| `curation.md` | Readable summary, evidence gaps, numerical findings and recommendations |
| `curation.npz` | Exact evaluated/flagged masks and the associated native time/radius coordinates |

Each mask has a pair of keys ending in `__evaluated` and `__flagged`. A point outside the evaluated
mask is **unassessed**, not a pass. JSON includes index examples but explicitly records when they
are truncated; the NPZ preserves the complete masks. Load it with `numpy.load(..., allow_pickle=False)`.

Input hashes, effective-manifest hash, evaluator version, and curation-module hashes make the
report traceable. Outputs are separately and atomically replaced; the three-file report bundle is
not a transactional database update. Existing reports in the selected output directory are replaced.
Original CDF/TR.DAT/evidence files are protected against being used as output destinations.

## What is implemented

- Native `TIME`/`TIME3` and `X` screening, without equal-length grid substitution.
- Fill/missing-value preservation, physical profile screening, selectable time/radius windows.
- Small-gradient coefficient-identifiability warnings and optional **supplied derivative** uncertainty.
- Origin, uncertainty and independence declarations; missing rotation is not replaced by zero.
- Quasineutrality screening only for a supplied, explicitly complete charged-species inventory.
- A conservative declared-input timeline; scalar beam-current-fraction checks; comparisons to explicit
  campaign/configuration reference values; density/current coupling and rotation-control review flags.
- Native iteration/timing/residual/upwind-channel inventory and opt-in numerical thresholds.
- Strict JSON sidecar validation, read-only bounded-array caching, unit and CLI integration tests.

## What is not implemented

Raw diagnostic fitting; covariance estimation; consistent TRANSP rerun ensembles; a resolved
TRANSP-default/model execution state; species-completeness verification; device-geometry registry
retrieval; automatic correction; event segmentation; trained use-case acceptance policies; IMAS
curation IDS export; and production certification are **not** implemented.

`pass` means a particular screening predicate passed. It does not mean that the entire run is usable.
`not_assessed` is not an error bar of zero. Values declared `uncertainty: provided` are evidence
claims, not proof that uncertainties were propagated or calibrated.

## Runtime and exit codes

The engine needs NumPy; it does not import IMAS. The CLI prefers netCDF4. Without netCDF4 it can
read classic CDF1/CDF2 through SciPy; NetCDF4/HDF5 and CDF5 require netCDF4. The classic fallback
may buffer the file, so per-variable/cache limits are not a total process-memory bound.

- Exit `0`: execution completed with no finding at the chosen blocking severity. This is not certification.
- Exit `1`: blocking screening findings (`--fail-on fail`, by default), or unassessed findings when
  `--require-complete` is requested.
- Exit `2`: input/manifest/backend/output error prevented completion.

`--fail-on warning` is useful for review queues; `--fail-on error` records concerns without making
scientific findings fail a batch job. The report is still written for ordinary screening findings.

## Documents

- [Scientific design and curation ideas](design.md)
- [Evidence sidecar format](evidence-format.md)
- [Implemented rules and threshold semantics](rules.md)
- [Development roadmap and acceptance tests](roadmap.md)
- [Review of the t2i2 source snapshot](../review-t2i2.md)
