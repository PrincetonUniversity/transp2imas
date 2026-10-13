# Real ISolver regression cases: ITER 38530P94 / NSTX 142202U01

These two independently generated TRANSP runs have `LEVGEO=12` (ISolver) and
contain the coil and passive-conductor current channels that were missing from
the earlier interpretive/RF regression cases.

## File layout

```text
examples/data/iter/
    38530P94.CDF
    38530P94TR.DAT
    p94_pt.nml
examples/data/nstx/142202U01/
    142202U01.CDF
    142202U01TR.DAT
```

The `TR.DAT` files list `PFC_NAMES` and select `LEVGEO=12`. The ITER case
includes PT_SOLVER settings in `p94_pt.nml` for input provenance. The special
`NLISODIF` and `NLINDUCT` choices differ between the cases; this regression
does not infer a circuit topology or conductor geometry from those switches.

| Source | Time samples | `KK_*` coil currents | `VV_*` terminal voltages | `VS_*` source voltages | `PS_n` passive currents |
| --- | ---: | ---: | ---: | ---: | ---: |
| ITER 38530P94 | 5 | 11 | 11 | 11 | 119 |
| NSTX 142202U01 | 486 | 9 | 9 | 9 | 82 |

All of these saved channel arrays have native CDF dimension `TIME`, with
unit attributes `A` or `V`. The writer also supports the independent `TIME3`
grid for other cases. The regression checks exact values, native times, channel
identities and the `PFC_NAMES` identifiers declared in the input file.

## Distinct meanings of coil-current-related outputs

* `KK_<coil>`: actual simulated *circuit current*, in amperes. This is the
  default source for `pf_active.coil[].current` when no validated machine
  mapping is available.
* `CC_<coil>`: current in the coil from the equilibrium solution; the two
  supplied CDFs demonstrate that it stores **effective ampere-turns**. For
  every matched coil, `CC_<coil> / KK_<coil>` is constant to the precision of
  the CDF floats and approximately an **integer**.
* `KT_<coil>`: target/reference circuit current, not an actual measurement.
* `KD_<coil>`: desired current for the corresponding control circuit.
* `KS_<coil>`: JIT sigma/uncertainty for the circuit.
* `VV_<coil>` and `VS_<coil>`: respectively terminal and source voltages.
* `PS_n`: current in passive-structure channel `n`; the CDF `long_name`
  includes a label such as `<iw51a>` or `<bvfb>`. These labels are retained
  verbatim, but cannot establish geometry or current orientation independently.

**Critical unit distinction.** In ITER `CC_CS3U / KK_CS3U ≈ 548`; in NSTX
`CC_OH / KK_OH ≈ 963`. Filling IMAS `coil.current` directly with `CC_`
would be wrong by factors of hundreds. A manifest can bind `CC_` only with an
**explicit, independently verified** `cc_effective_turns` divisor; do not use
the empirical ratios alone as a substitute for an approved device definition.
The converter never applies the plasma `Ip` or magnetic-field `Bphi` signs to
these engineering circuit currents without a manifest-provided sign mapping.

`make verify-isolver` also writes reports listing all active coil identifiers,
`CC/KK` ratios, and the names of every passive loop channel:

```text
build/pf-isolver/38530P94-pf-active.csv
build/pf-isolver/38530P94-pf-passive.csv
build/pf-isolver/142202U01-pf-active.csv
build/pf-isolver/142202U01-pf-passive.csv
```

## Verification

```bash
make check                 # fast checks, no large ISolver files needed
make verify-pf-active      # synthetic geometry / sign / DD3/DD4 tests
make verify-pf-passive     # synthetic geometry / sign / DD3/DD4 tests
make verify-isolver        # both real CDFs; requires data at paths above
make verify-all            # existing NSTX physics regressions
make verify-waves          # existing RF physics regressions
make verify-dd             # IMAS-installed environment only
```

`make verify-isolver` is a **source-side** test. It uses mock IMAS data
structures to verify channel-to-node shape, time and values; it does not
perform an actual `DBEntry.put()` in DD3/DD4 or verify the machine geometry.
Before an engineering-quality `pf_active`/`pf_passive` IDS is produced, obtain
the approved ISolver device coils/passive-element geometry, signed turns,
source-to-coil circuit associations, and loop orientation. Then supply an
explicit machine JSON manifest to the converter and run the end-to-end
`make verify-dd` and actual conversion on an IMAS-equipped installation.

## Large-fixture policy (GitHub)

`142202U01.CDF` is **about 166 MB**, exceeding GitHub's normal **100 MB
single-file limit**. It is included in the downloadable full source fixture
archive and expected at the path above for local regression, but is ignored
by `.gitignore` to prevent accidental plain-Git commits. The lightweight
`make check` and public GitHub CI do not require it.

If repository maintainers wish to store that CDF in GitHub, establish the
project's Git LFS policy and ensure LFS storage/quota is available, then run
something like:

```bash
git lfs install
git lfs track "examples/data/nstx/142202U01/142202U01.CDF"
git add .gitattributes
git add -f examples/data/nstx/142202U01/142202U01.CDF
```

Confirm the staged object is an LFS pointer before committing. Alternatively,
keep the CDF on the lab filesystem or release-data storage and use the local
`verify-isolver` target. Do not commit it as an ordinary Git blob.
