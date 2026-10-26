# Namelist resolution (`curation/namelist.py`, `curation/inputs.py`)

Curation reads a TR.DAT through a TRANSP-specific adapter built on `f90nml` (supported: `f90nml>=1.4`; developed
against 1.5.0). No IMAS import, no network, no `eval`, no shell. Source files are only read, never written.

A TR.DAT is not a standard namelist: bare `NAME=value` assignments with no `&group`, `!` comments (also trailing, and
`!` inside quotes), `~update_time = <t>` lines that start a new time block, mixed case, tabs, D exponents, and lines such
as `PT_TEMPLATE='y05_pt.nml'` that reference a child file. A PT template is a standard namelist (`&pt_controls ... /`)
with derived-type names (`pt_confinement%cdbm%active`), repeat notation and arrays.

## Three stages

| Stage | Module | Does | Never does |
|---|---|---|---|
| 1. Source collection | `namelist.collect`, `references`, `load_references` | ordered statements with line/column spans, `~update_time` blocks, a synthetic standard namelist with a map back to source lines, file references | interpret values, execute includes |
| 2. Syntax parsing | `namelist.parse_statement` (f90nml) | literals, repeats, sections, nulls, derived-type paths, element placement with actual Fortran indices | resolve state, invent defaults |
| 3. State resolution | `inputs.declared_timeline` | ordered updates over the update blocks | read text, call f90nml |

Only stage 3 may speak of an "effective declared state", and even there it means *what the input explicitly declares*,
not what TRANSP ran with.

### Stage 1

* Comments are stripped respecting quotes and doubled quotes (`input_provenance.strip_fortran_comment`).
* `&name`, `$name`, `&end`, `$end` and `/` are group markers (structure, not gaps). A TR.DAT allows one assignment per
  line; a namelist document allows several per line (`a=1, b=2 c=3,4 /`). On a TR.DAT line with several assignments the
  line is not interpreted and every target on it is invalidated (the old reader never read such a line either).
* A non-assignment line continues the previous statement in a namelist document; in a TR.DAT only after a trailing comma.
  Otherwise it is an `unsupported_statement`.
* A bare integer before the first group of a namelist document is the PT-template group count; it is recorded and
  compared with the number of groups found (`syntax_check.template_count`), not interpreted.
* The synthetic namelist (`namelist.synthetic_namelist`) puts the assignments of each update block into one
  `&tr_block_N` group; `line_map[k]` gives `{line, end_line, column}` of synthetic line `k`. Statements whose value list
  cannot be tokenised (unbalanced quote, absurd repeat count) are left out so f90nml never receives them.
* References: only what `input_provenance.references_from_text` recognises. A key containing `TEMPLATE`
  (`PT_TEMPLATE`, `MMM_TEMPLATE`) or a value ending `.nml`/`.namelist` is followed and parsed as its own namelist
  document. Other recognised references (for example `TORIC_INPUT_FILE='x.inp'`) are listed with status
  `recognised_not_followed`. `include 'f'` / `#include f` are recorded (`include_not_executed`) and never executed; the
  bundled inputs contain none. A glob or expression value is not a reference.
* Resolution: relative to the *including* document's directory, symlinks resolved, then required to lie inside the
  allowed roots (default: the TR.DAT directory; `declared_timeline(path, allowed_roots=[...])`). Depth limit 4
  (`max_depth`), at most 32 files, 10 MiB per file, cycle detection. Each file entry records `path`, `sha256`, `size`,
  `role`, `depth`, `referenced_by`, `key`, `line` and a `status`: `parsed`, `missing`, `cycle`, `depth_limit`,
  `outside_allowed_roots`, `unreadable`, `file_limit`, `include_not_executed`, `recognised_not_followed`.
* A referenced PT template is its own document (`timeline['documents']`). Its `pt_update_time` and model semantics are
  preserved as assignments, not interpreted, and its values never enter the TR.DAT states.

### Stage 2

Every statement is passed to f90nml as its own one-assignment namelist, so ordering, spans and element-wise updates
stay under the adapter's control (f90nml merges repeated assignments and loses them). The adapter trusts an element only
when its raw token (`fortran_literals.parse_scalar`) and f90nml agree. f90nml alone is not trusted because it:
reads `1+2` as `100.0` (a Fortran exponent form), reads `$word` as a null, silently drops surplus values and strides,
drops all but the first value of `A(k)=v1,v2`, and with `strict_logical=False` accepts `Tx` as true.

Recorded parser configuration (`timeline['parser']`, `documents[i]['parser']`): kind `f90nml`, version, and
`config.parser_options` (`comment_tokens='!'`, `default_start_index=1`, `global_start_index=None`, `row_major=False`,
`sparse_arrays=False`, `strict_logical=False`) plus the adapter rules for comments, start index, array order, derived
types, repeats/nulls, scalars, statements and continuation.

| Item | Behaviour |
|---|---|
| scalars | `fortran_literals` (bool, number, quoted string with doubled quotes), plus the dotted logicals `.T .F .TRUE .FALSE`. The extension lives in `namelist.scalar_value`; `inputs.literal` stays the strict `fortran_literals` contract |
| D/E exponents | `9.205d-2`, `1.5E+2`, `-.5d0` |
| repeats | `3*1.0`, `2*3.,20.`, `2*'x'`; `n*` alone is n nulls |
| nulls | `A=1,,3`: the null leaves its element **unchanged** (standard namelist semantics). It is not zero, not unknown, not an omitted element in the record (`{"null": true}`) |
| sections | `A(2:3)=7,8`; fewer values than elements leave the rest; more is `too_many_values` |
| element continuation | `A(3)=1,2` assigns `A(3)`, `A(4)`; for one subscript f90nml receives the equivalent section |
| multidimensional | Fortran column-major (first subscript fastest): `A(1:2,1:2)=1,2,3,4` is `A(1,1)=1 A(2,1)=2 A(1,2)=3 A(2,2)=4`. Flat keys `A(i,j)` |
| several values from a multidimensional element start | `P(1,1)=0.,180.` needs the undeclared first extent: `unbounded_continuation`, unresolved |
| sparse arrays | `A(3)=1` stores `A(3)` only. No positional list, no `A(1)`, no `A(2)` |
| derived types | `a%b%c` is the flat key `A%B%C`; `a%b(2)` is `A%B(2)`; a subscript inside the path (`a(1)%b`) is `subscripted_component`, unsupported |
| repeated groups | statements apply in file order; f90nml merges repeated groups, `syntax_check.document_groups` lists the names |
| not supported | expression subscripts, strided sections (`A(1:5:2)`), bound-less sections (`A(:)`, `A(2:)`), complex numbers, `nan`/`inf`, `0x10`, unquoted strings, unterminated quotes |

An element whose literal is not a supported scalar is **unknown**: never 0, `False` or the previous value.

### Stage 3

* Blocks: `~update_time = t` closes the current state `[start, t)` and opens the next, exactly as before. A
  non-increasing or unparseable time records a gap, sets `temporal_valid` false and opens no block.
* A declaration persists into later blocks. A later partial update (`A(2)=9`, `A(3:4)=7,8`) keeps the earlier known
  elements it does not touch.
* An update that cannot be resolved invalidates what it may have changed (`inputs.invalidate`, P0.6), never reusing the
  superseded value: unknown element, only that element (and the rest of its section); whole-name statement
  (`A=garbage`) the whole array; section with known bounds, the section; bound-less, strided, expression-subscripted or
  multidimensional-continuation, the whole base array.
* A bare scalar (`A=7`) for a name that holds elements invalidates the elements (scalar or `A(1)` is ambiguous). A
  whole-array statement without a subscript (`A=1,2`) is `A(1)`, `A(2)` and supersedes a scalar `A`.
* No TRANSP defaults: omitted switches, omitted bounds and omitted array elements stay unknown. f90nml defaults
  (`default_start_index=1`) are parser settings, not TRANSP defaults. A whole-array assignment without a subscript
  assumes lower bound 1 and says so once per array and document (`unresolved`, kind `assumed_lower_bound`).
* No claim about additive versus replacing model semantics, preprocessing, namelist order inside TRANSP or runtime state.

## Timeline structure

Unchanged: `available`, `states[{start_time, end_time, values, source_lines}]`, `gaps`, `temporal_valid`,
`source{path, sha256}`, `semantics`. `values` is a flat dict of upper-case keys to bool/float/str; elements are `NAME(i)`
or `NAME(i,j)`; derived types `A%B%C`.

Added:

* `parser`: `{kind, version, config}` or `{kind: 'scalar_fallback', version: null, config: {}, limitations: [...]}`.
* `assignments`: ordered records `{order, block, file, line, end_line, column, name, indices, keys, value, values_truncated,
  raw, status, problems, assumption, invalidated}`. `value` is the scalar, or a list of per-element entries (value,
  `{"null": true}` or `{"unknown": raw}`); at most 256 elements per record. `status` is `resolved`, `partial` or
  `unresolved`. At most `max_assignments` (5000) records; `assignments_truncated` and `assignment_limit` say so.
* `files`: `{allowed_roots, max_depth, max_files, entries}` (the include chain with sha256). The top-level entry carries
  `syntax_check` (f90nml reading of the synthetic namelist) and `groups`.
* `documents`: each followed namelist with its own `parser`, `sha256`, `syntax_check`, `groups`, `assignments`,
  `unresolved`.
* `unresolved`: `{kind, file, line, name, reason, ...}` for every construct not resolved.

Everything is JSON-serializable with `allow_nan=False`.

## What is deliberately unresolved

TRANSP defaults and default catalogues; declared array bounds; `include` and any preprocessor; expressions (in values,
subscripts or file names); wildcard references; model execution semantics (PT solver, additive versus replacing
updates, `pt_update_time`); anything that happens at run time. Use supplied resolved run settings
(`resolved_inputs`) for those.

## Fallback

If `import f90nml` fails, `declared_timeline` runs the old scalar reader (`inputs.scalar_reader`) and says so:
`parser.kind == 'scalar_fallback'`, `parser.limitations` lists what it cannot do, a `parser_fallback` entry appears in
`gaps` and `unresolved`, `assignments` is empty, `files` is `null` (references are neither followed nor hashed).
It is a limited reader, never the full parser. Install with `pip install 'f90nml>=1.4'`; curation does not need IMAS.

## Migration notes

* `gaps` entries now always carry `kind` and `reason`; statement gaps also carry `name` and `invalidated` as before.
  New kinds: `unsupported_literal`, `syntax_error`, `f90nml_rejected`, `parser_disagreement`, `bounds_less_range`,
  `strided_section`, `expression_subscript`, `empty_section`, `subscripted_component`, `too_many_values`,
  `unbounded_continuation`, `multiple_assignments_per_line`, `unsupported_statement`, `include_not_executed`,
  `update_time`, `parser_fallback`, and the reference failures `missing`, `cycle`, `depth_limit`,
  `outside_allowed_roots`, `unreadable`, `file_limit`. A gap from a missing referenced file has `line` of the reference.
* Gaps that disappear because the construct is now resolved: arrays (`A=1,2,3`), repeats, sections, element
  continuation, `.T`/`.F` logicals, derived-type names, `$TRDATA`/`$END` markers and a `&group name=value` header line.
  On the bundled inputs the NSTX and ITER TR.DAT files go from 40 and 100 gaps to 0; C-Mod from 19 to 2 (the two
  `PHicha(1,j)=0.,180.` multidimensional continuations).
* `invalidated` lists the keys removed by that statement; for an unknown element it is only that element (plus the
  untouched remainder of its section).
* Keys that appear in addition: every element of an array/repeat statement and every `.T`/`.F` logical. No scalar or
  `NAME(i)` key the old reader produced changed value; the one difference is a correction (`X='a','b'` was read as the
  single string `a','b`; it is now `X(1)`, `X(2)`).
* `literal()` is unchanged (strict `fortran_literals`); `declared_timeline` accepts the extra optional arguments
  `allowed_roots`, `max_depth`, `max_assignments`.
* The engine's `INPUT.TIMELINE` evidence still reports only `gaps` and `source`; `parser`, `files` and `unresolved`
  are in the timeline for any consumer that wants them (the report writer dumps the whole timeline).
