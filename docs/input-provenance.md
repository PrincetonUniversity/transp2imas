# TRANSP input XML and model-specific provenance

TRANSP2IMAS keeps the run inputs in the IMAS database, not just in the local
TRANSP source directory. This applies to the main TR.DAT and to local,
explicitly referenced model namelists available at conversion time.

## Where to find the XML

| Destination | Stored input content |
|---|---|
| `workflow.time_loop.component[0].parameters` (`TRANSP`) | Whole-run XML inventory: original TR.DAT, parsed assignments, referenced/nested namelists, SHA-256 hashes and missing-file markers |
| `workflow.time_loop.component[1:].parameters` | An individual, explicitly configured/referenced model's input subset (e.g. TORIC, CQL3D, GENRAY, TORBEAM, FPP, LSC, PT_SOLVER, TGLF, MMM) |
| `waves.code.parameters` | Whole-run XML inventory so RF results can be traced to the original TRANSP input |
| `waves.code.library[].parameters` | **Per-wave-model** XML: original TR.DAT (containing inline TORIC/RF controls) plus that solver's referenced input files |
| `core_transport.code.library[].parameters` | Referenced transport model/PT_SOLVER inputs where present |
| `nbi.code.library[].parameters` | Referenced NUBEAM inputs where present |
| Other IDS `code.parameters` | Whole-run XML when this DD/IDS has that field |

`code.name` remains `transp2imas`: the converter produced the IDS. External
model names are not substituted for the converter's identity. Model/library
entries describe **input references**, not proof that a given model executed.

## Discovery rules

- The original TR.DAT is preserved with complete text, with parsed groups and
  assignments as a convenient, non-authoritative XML index. The parser never
  evaluates Fortran expressions; source text is authoritative.
- Files named by `*TEMPLATE`, `*_NML`, `*_NAMELIST`, `*_INPUT_FILE`, or relevant
  model input filenames are resolved relative to the referring input, then to
  TR.DAT's directory, then the converter's current working directory. Literal
  Fortran `INCLUDE` / `#include` directives are also followed.
- References are followed recursively with cycle detection. Same-file repeats
  are embedded once. Explicit references take precedence over adjacent files.
  In per-run directories (e.g. `examples/data/d3d/147634G14/`), additional
  `*.nml`/`*.in`/`*.inp` files are preserved as
  `ADJACENT_UNREFERENCED` with an explicit unverified-association marker.
  In shared directories, only run-prefixed auxiliary namelists are discovered;
  unrelated files are not swept into the archive.
- Each embedded file has a SHA-256 hash. Missing, inaccessible, or too-large
  files remain in the XML with a `status` attribute, and generate a warning.
  `input_files/@complete` is `false` when any referenced file is missing or
  discovery was truncated. For safety, one input is capped at 10 MiB; recursion
  is capped at depth 12 and 128 distinct input files.
- The XML is saved *before* other writers, so a later RF/IDS conversion error
  does not prevent workflow input provenance from being written.

## Wave-code status and limitations

The system recognizes TORIC, FPP, CQL3D, GENRAY, TORBEAM, TORAY and LSC from
uncommented run controls and explicit file references, **without assuming**
that all codes were run. This preserves wave-model inputs when the source
files are present. CQL3D is not automatically inferred from a reference in a
comment. If a CQL3D input exists but is neither referenced in TR.DAT, nor
named by another referenced input, it cannot be reliably associated with the
run and is not silently archived.

CQL3D typically has its own namelist/control file. Ensure that file is
available at conversion time and is discoverable through an **actual, supported**
TRANSP namelist/template reference or a nested referenced input file. The
converter's synthetic tests use a `CQL3D_NML` key solely to test generic
reference discovery; this does **not** establish that `CQL3D_NML` is a valid
TRANSP input parameter, so do not add it to production TR.DAT without checking
TRANSP's namelist specification. If a referenced file is missing, restore it at
its original relative location before conversion and rerun.

## Testing and verification

```bash
make verify-inputs
make export-inputs TRDAT=examples/data/iter/38530P94TR.DAT CASE=38530P94
make check
make verify-dd
python3 convert.py -i examples/data/iter/38530P94.CDF \
    -b HDF5 -new --dd-version 4.1.1
```

The first two are source-only; `verify-dd` and the actual conversion require
IMAS-Python and the selected DD installed. After conversion, inspect
`workflow.time_loop.component[]` and `waves.code.library[]`; each model entry
has `name` and an XML `parameters` string. For the included ITER 38530P94
case, `PT_TEMPLATE=p94_pt.nml` is available and embedded. The existing C-Mod
TORIC and DIII-D GENRAY/TORBEAM examples retain inline RF controls in TR.DAT
but do not contain separate CQL3D files.

`make export-inputs` produces local preview XML files in `build/input-xml/`. These are inspection artifacts: their `input_files` contents match the IMAS inventory, but IMAS `workflow` wraps them with additional run/converter metadata.

## Newly supplied source examples

DIII-D `147634T13` contains an adjacent `toray.in`, `147634G14` contains
`genray_ech.in`, and `147634B13` contains `ech_inp.nml`. These are now embedded
with their relevant model names but remain explicitly marked as adjacent,
unreferenced inputs. NSTX `141711Y05` embeds `y05_pt.nml` via its actual
`PT_TEMPLATE` reference. DIII-D `147634T16` has `t16_pt.nml` and `mmm.nml`
but its `147634T16TR.DAT` was not supplied. The model files are preserved,
while the missing master input is visible as an incomplete-input warning.

Use `make verify-inputs` for all 10 examples or `make verify-all` to validate
input provenance together with CDF and physics checks.
