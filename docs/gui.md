# Workbench (GUI)

A local web page for TRANSP2IMAS: read curation and conversion reports, **run conversions**, and **plot the converted IMAS
data**. It is a thin layer over the existing scripts: the form is built from the converter's own option parser, the page
shows the exact command it will run, and `convert.py` stays the single implementation.

```bash
python3 -m gui --root build                     # read-only report viewer
python3 -m gui --root examples/data --enable-run --imas-python /path/to/imas/python    # also convert and plot
# Open: http://127.0.0.1:PORT/?token=...
```

On a remote or HPC host keep the default loopback binding and tunnel the port from your machine, then open the printed URL
(including the token): `ssh -L PORT:127.0.0.1:PORT user@host`.

| Option | Meaning |
|---|---|
| `--root DIR` (repeatable) | reports that may be opened, and CDF files that may be converted (default: current directory) |
| `--input-root DIR` | further directories CDF files and geometry JSON may come from |
| `--output-root DIR` | where jobs and the **scratch IMAS database** are written (default `build/gui`) |
| `--database-root DIR` | extra IMAS database whose HDF5 entries may be listed and plotted, read-only |
| `--enable-run` | allow starting conversions (default: read-only; planning always works) |
| `--imas-python PATH` | interpreter with IMAS-Python **or the virtual-environment directory that has it** (`~/venvs/imas/` works), used to convert and to read entries (default: this interpreter if it has IMAS) |
| `--allow-default-database` | allow writing to, and plotting from, your default IMAS database (`~/public/imasdb`) |
| `--host`, `--port`, `--allow-remote`, `--open`, `--verbose` | network and convenience options (loopback by default) |

The page runs in the curation-only environment; conversion and plotting run in the interpreter given to `--imas-python`.
If it cannot import IMAS the launcher and the page say so and conversion and plotting stay off.

**Running needs both** `--enable-run` and an interpreter with IMAS. Giving `--imas-python` alone turns plotting on but leaves
*Run conversion* disabled, and the page says so next to the button (the badge in the title bar reads *read-only* or
*conversion enabled*). A bad `--imas-python` path stops the launcher with the reason.

## Views

### Order of the views

Convert, Explore, Reports. The page opens on Convert (or on Reports when the URL has `?open=<report path>` or `?view=reports`);
opening Reports with exactly one report available opens it.

### Reports

Curation report: **Overview** (status counts, application summaries, fingerprints, rule coverage, migration), **Findings**
(filter by status, rule or text), **Rules** (the execution ledger and plugins), **Evidence** (lineage graph, time coverage,
propagated uncertainty), **Events** (timeline, origins, candidates, policy exclusions), **Masks** (evaluated, flagged and
excluded samples as a heat map) and **Inputs** (parser, followed files, unresolved constructs, code identity).
Conversion report: **Summary** (status, publication to the backend), **Writers**, **Mapping journal**. A report from an older
curation version is shown with what it has and a notice. A curation report is a screening record, not a certification.

### Convert

1. **What to convert**: a CDF found below the input roots (the TR.DAT beside it is used) or a typed path.
2. **How**: IMAS backend (HDF5 is the default the form always sends; the converter's own default is MDS), Data Dictionary
   version (`auto`, `4.1.1`, `3.39.0` or a custom `x.y.z`), machine, run number, occurrence, "create a new entry".
3. **Where to write**: the *scratch database* under the output root (default), a directory you choose **inside an output
   root**, or your default IMAS database (only with `--allow-default-database`).
4. **Which IDSs**: the 22 writers as checkboxes, with All, None, Core plasma and Hardware shortcuts; the converter runs
   them in its fixed dependency-safe order and a subset gives a partial entry (the page says so).
5. **Mapping and physics options**: *skip variables that cannot be mapped* (the converter's `--lenient-mapping`: such a variable
   is recorded with its reason in the conversion report instead of stopping the conversion; by default the conversion stops at
   the first one), plasma-current and toroidal-field directions, passive and active geometry JSON, equilibrium moments,
   poloidal points, EC launch from namelist.
6. **Diagnostics**: log level, verbosity, debug output.

Every option's label help comes from the converter's own parser; `tests` fail when `convert.py` gains an option the page has
not decided about. Beside the form the **plan** shows, live, the exact command (an argument list, shown with shell quoting
and copyable), the target entry (machine, shot, run, DD major, database, path), whether the entry exists, any warnings
(auto DD, skipped variables, partial entry, MDSplus), and the validation errors against the field that caused them.

**Overwrite guard.** The target is computed the way the converter computes it. If the entry already exists (or cannot be
checked, as for MDSplus), *Run* stays disabled until you confirm: type the **exact entry path** or press **Fill in the path**
(a deliberate click that puts the path in the box); the server checks it again and refuses anything else. Two ways forward are
offered side by side: confirm the overwrite, or press **Use run N instead**, which fills in the next run number with nothing
in it (HDF5 only; MDSplus entries cannot be checked).

What an overwrite does, as tested against the real converter: with *Create a new entry* on, the existing entry's contents are
**deleted** and only the IDSs converted now remain; with it off, the selected IDSs are written into the existing entry and
the other IDSs it already had stay. The page itself deletes nothing: the converter does the writing, and the page says which
of the two will happen before you confirm.

The form is never rebuilt while you type (only the plan beside it is), so digits and the caret stay where you put them; fields
with a problem are marked in place.

**Running.** One conversion runs at a time, in a job directory `<output root>/jobs/<id>/` holding `job.json` (command, plan,
outcome), `log.txt` (followed live) and the converter's `conversion-report.json`. A job can be cancelled (its process group
is stopped). A job left *running* by a server that died is shown as **lost**, never as succeeded. When it finishes, buttons
open the conversion report in the Reports view and the converted entry in Explore.

### Explore (visualization)

Choose a converted **HDF5** entry (found in the scratch database and any `--database-root`, or typed). Returning to the view
re-reads the entry ("Checking the entry for changes…") and drops the old plot only if the entry changed on disk, so an entry you
just overwrote or extended is shown as it is now. The page lists the IDSs
with their sizes and time ranges, offers curated **recipes** (80, grouped by quantity) for the IDSs that are present, and a **node
browser** for any numeric node of an IDS (read from its first time slice, so it is fast on large entries).  The recipes cover
electron and ion temperature, density and pressure (each ion species separately), safety factor, current components,
Zeff, toroidal rotation and angular momentum; the heat, particle and momentum diffusivities and fluxes of every transport
model and the energy, particle, momentum and current sources of every source; equilibrium profiles (q, shear, pressure,
F, FF', currents, shape, volume, |B|), global traces (Ip, q95, shape, axis positions, beta, li) and the maps of psi and of
B_R, B_phi and B_Z; the wall outline; PF coil currents; beam power and energy per unit; fast-ion density, pressure, current,
torque, collisional power and thermalization, and the fast-ion birth profiles of `distribution_sources`; radiation emissivity;
and the summary traces.

**Equilibrium maps** draw flux surfaces: a number box sets how many contour lines (default 10 for the flux map; for a flux
map whose axis and boundary flux are stored they are evenly spaced in flux from the axis to the boundary, otherwise evenly
inside the data range).  Check boxes hide the wall and the plasma boundary or switch the colouring off.  **The wall** is
drawn as an R,Z outline on its own and as an overlay on every equilibrium map.

**Fast-ion distribution function.**  When an entry has a `distributions` IDS the page offers a panel with an energy-pitch
map, an energy spectrum, a pitch distribution and the density of every zone drawn on the R,Z plane, with a choice of
distribution, an optional zone (click a zone on the R,Z picture to choose it) and optional energy and pitch windows, and the
time slider.  With no zone the plots sum over the plasma volume (so integrals are particle numbers); with a zone they show
the stored per-m^3 values of that zone.  The array is read as `(zone, energy, pitch)` with the zone index varying fastest;
this was verified on a real NUBEAM entry: the volume-weighted total equals `global_quantities.particles_n` and the row
averages equal `profiles_1d.density_fast` (`examples/verify_gui_plots.py`, `FastIonRealEntryTests`, which converts the local-only
`examples/data/nstx/120967E03_send` case (too large for the repository; skipped with that reason when absent, or set
`TRANSP2IMAS_DF_ENTRY` to an already converted NUBEAM entry).  The reader refuses an array whose size does not equal zones x energy x pitch.
Zones are drawn with straight chords between the stored nodes, so where the grid has few nodes (the innermost wedges, which
are triangles) a zone is smaller than its stored area (those wedges are 2/pi of it); the density values are not affected.

Plot options: a **time slider** with play/pause and step buttons (profile frames are fetched once and scrubbed locally;
maps load the slice on release), **overlay eight times** (one colour per quantity, fading in with time), a **log y axis**,
legend toggles, a hover readout, **compare with another entry** (overlaid, dashed; the page notes when the time bases
differ), and **export** as PNG, SVG or CSV. Wildcards expand to every source, model or unit and label each element.
A node the entry does not have is reported as not present, never drawn as zeros; long series are thinned for drawing and the
page says so; a nonfinite sample is a gap.

Plots show what is stored. They do not interpret it: `dim1`/`dim2` of the flux map follow the entry's grid type, and a
plot of converted data is not a validation of the conversion.  What the converter does not store cannot be drawn: there are
no X-point, strike-point or separatrix nodes (the plasma boundary outline and the magnetic axis are shown), ion heat
diffusivity is stored only as the total ion value, and PF coil geometry is not stored (their currents are plotted).

The CDF list offers only `<runid>.CDF` files (for example `141711Y05.CDF`). Companion files of the run, such as `<runid>PH.CDF` and
plasma-state CDFs (`<runid>_ps.cdf`, `_fi_*.cdf`), are not converted and are not listed; a file with another name can still be
chosen by typing its path.

## Security model

The page can now start a program, so the posture is stricter than a viewer's:

* it binds to `127.0.0.1` (another `--host` needs `--allow-remote`); every request needs the token and a `Host` header
  naming the listening address (against DNS rebinding);
* **running is off unless `--enable-run`**; without it every POST that starts or cancels a job is refused;
* the only POST endpoints are `convert/plan` (a dry run), `convert/start` and `jobs/<id>/cancel`. A POST needs the token in
  the **`X-GUI-Token` header** (a cookie or query token is refused, so another site cannot make the browser send it),
  `Content-Type: application/json`, a body of at most 64 KiB, and, if an `Origin` header is present, the viewer's own origin;
* the only program ever started is this repository's `convert.py`, from an **argument list** built by validators (never a
  shell string): unknown options are rejected, values are typed and range-checked, the CDF and geometry files must be
  below the input roots (symlinks that leave them are refused), a custom database must be inside an output root, and the
  default IMAS database needs `--allow-default-database`;
* plotting reads entries in a **read-only worker subprocess**; an entry must be an HDF5 entry directory below a readable
  database root; the worker never writes. One reader per entry at a time (HDF5 refuses concurrent opens);
* report files are read only as `.json`/`.npz` below the roots, with a **64 MiB cap because the viewer reads them into
  memory**. The cap does not apply to the TRANSP CDF you convert (it is only handed to the converter and routinely holds
  hundreds of MB; a 243 MB CDF plans normally); geometry JSON, which is read, keeps the cap. The page builds its DOM with `textContent` only and sends a
  strict Content-Security-Policy, so analyst-written text cannot run as script; no remote asset is loaded.

## Tests

| Command | What it covers |
|---|---|
| `make verify-gui` | `verify_gui.py` (53: view models, masks, path rules, server, static-page lint), `verify_gui_convert.py` (66: option schema vs the converter, planning and overwrite guard, jobs, POST API), `verify_gui_plots.py` (plot helpers, fast-ion layout and reductions, a synthetic `distributions` entry with hand-computed answers, a **real** conversion read back and compared with the source CDF and a direct IMAS read, and real overwrite/extend/free-run conversions on a private entry) |
| `make verify-gui-browser PLAYWRIGHT_PYTHON=...` | opt-in: Chromium drives a hostile report, a server with running off, then the whole workflow for real: typing digits one at a time, converting, the overwrite guard and free-run suggestion, a real overwrite, plotting (flux-surface contours, wall outline, ion profiles, and the fast-ion distribution function on a small synthetic entry) and export |

The real-conversion tests need IMAS-Python and are skipped with that reason otherwise. Guard tests were checked by breaking
the guard (confirmation, Origin check, header token) and confirming the suite fails.

## Layout of the code

```text
gui/reports.py, masks.py, browse.py    view models, mask reduction, the one place that decides which paths may be read
gui/convert_plan.py                    option schema, validation, command, target-entry guard (no IMAS, no subprocess)
gui/jobs.py                            running convert.py: log, cancel, persistent job records
gui/imas_worker.py, imas_client.py, recipes.py    reading converted entries for plotting (worker subprocess in the IMAS environment)
gui/server.py                          standard-library HTTP server and API
gui/static/                            core.js, charts.js, app.js, convert.js, explore.js, index.html, style.css (no build step, no assets)
```

No third-party dependency. The planning, job and view-model modules do not depend on the web layer.

## Not built yet

1. A curation runner (a form for the `curate.py` options with a live log), the counterpart of the Convert view.
2. Batch conversion (`-i all -d dir`), which the form deliberately does not offer.
3. Plotting MDSplus entries (HDF5 only), comparing maps, and overlaying a CDF variable on its converted counterpart.
4. A run history with fingerprints and one-click reruns beyond the recent-jobs list.
