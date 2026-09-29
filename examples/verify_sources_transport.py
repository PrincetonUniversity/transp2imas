#!/usr/bin/python
"""Numerical sanity checks for audited source/transport/radiation/edge mappings.

This script reads a TRANSP CDF directly; it does not require IMAS.  It checks
SI conversion and the reconstruction formulas used by the converter.  It does
not validate the target Data Dictionary schema -- use verify_dd_compat.py for
that in an environment where the imas Python module is installed.
"""
import re
import sys
from pathlib import Path
import numpy as np
from scipy.io import netcdf_file

DEFAULT_CDF = (Path(__file__).resolve().parent / 'data' / 'nstx' /
               '141711Y01' / '141711Y01.CDF')
path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_CDF
f = netcdf_file(path, 'r', mmap=False)


def get(k):
    return np.asarray(f.variables[k][:], dtype=np.float64) if k in f.variables else np.array([])


def edge_flux(divname, scale, dv, surf):
    div = get(divname) * scale
    if not div.size:
        return np.array([])
    inside = np.cumsum(div * dv, axis=1)
    return np.divide(inside, surf, out=np.zeros_like(inside), where=np.abs(surf) > 0)


def radiation_integral(name, dv):
    q = get(name) * 1e6  # W/cm3 -> W/m3
    if not q.size:
        return np.array([]), np.array([])
    inside = np.cumsum(q * dv, axis=1)
    return q, inside


dv = get('DVOL') * 1e-6
surf = get('SURF') * 1e-4
print('file:', path)
print('grid:', dv.shape, 'surface:', surf.shape)

print('\ncore_sources volume-integral checks')
for prof, glob in [('PBE', 'BPTE'), ('PBI', 'BPTI'), ('PBTH', 'BPTH')]:
    q = get(prof) * 1e6
    if q.size:
        integ = np.sum(q * dv, axis=1)
        if glob in f.variables:
            ref = get(glob)
            rel = np.nanmax(np.abs(integ - ref) / np.maximum(np.abs(ref), 1.0))
            print(f'  {prof} -> {glob}: max relative integral error {rel:.3e}')
        else:
            print(f'  {prof}: integrated, no {glob} reference')

print('\ncore_transport observed/model divergence -> flux checks')
for divname, scale in [
    ('EPTR_OBS', 1e6), ('EETR_OBS', 1e6), ('IETR_OBS', 1e6), ('AMTR_OBS', 1e6),
    ('EPTR_MOD', 1e6), ('EETR_MOD', 1e6), ('IETR_MOD', 1e6), ('IPTR_MOD', 1e6),
]:
    flux = edge_flux(divname, scale, dv, surf)
    if flux.size:
        print(f'  {divname}: finite={np.isfinite(flux).all()} '
              f'edge range=[{flux[:, -1].min():.6e},{flux[:, -1].max():.6e}]')

species_mod = sorted(name for name in f.variables
                     if re.fullmatch(r'PTR(?:H|D|T|HE3|HE4|3|4)_MOD', name))
for divname in species_mod:
    flux = edge_flux(divname, 1e6, dv, surf)
    print(f'  {divname}: finite={np.isfinite(flux).all()} '
          f'edge range=[{flux[:, -1].min():.6e},{flux[:, -1].max():.6e}]')
if 'IPTR_MOD' in f.variables:
    print('  IPTR_MOD note: numerical reconstruction shown above, but no aggregate '
          'total-ion particle-flux node exists in core_transport; it is not exported directly.')


print('\nedge_sources neutral-process decomposition checks')
dv_cm3 = get('DVOL')


def max_rel(a, b, floor=1.0):
    if not (a.size and b.size):
        return np.nan
    den = np.maximum(np.maximum(np.abs(a), np.abs(b)), floor)
    return np.nanmax(np.abs(a-b) / den)


def scaled_residual(lhs, rhs, *terms):
    """Residual normalized by the largest term, robust when balance~0."""
    scale = np.ones_like(lhs, dtype=np.float64)
    for arr in terms:
        if arr.size:
            scale = np.maximum(scale, np.abs(arr))
    return np.nanmax(np.abs(lhs-rhs) / scale)


edge_groups = [
    dict(origin='gas_flow', e='SEGF_D', cx='SFCXGF_D', si='SIGF_D_D',
         agg='SISGF_D', global_rate='GASD', pi='PIGF_D', pcx='PCXGF_D',
         pbal='BALE0_GF_D', pfl='PFLX0GF_D', tq='TQIGF_D', tqcx='TQCXGF_D',
         tqbal='TQBA0_GF_D', tqfl='TQ0FLGF_D'),
    dict(origin='recycling', e='SERC_D', cx='SFCXRC_D', si='SIRC_D_D',
         agg='SISRC_D', global_rate='RCYD', pi='PIRC_D', pcx='PCXRC_D',
         pbal='BALE0_RC_D', pfl='PFLX0RC_D', tq='TQIRC_D', tqcx='TQCXRC_D',
         tqbal='TQBA0_RC_D', tqfl='TQ0FLRC_D'),
]
for d in edge_groups:
    e, cx, si, agg = get(d['e']), get(d['cx']), get(d['si']), get(d['agg'])
    print(f"  {d['origin']}")
    if si.size and e.size and cx.size:
        print(f"    {d['si']} = {d['e']} + {d['cx']}: max rel "
              f"{max_rel(si, e+cx):.3e}")
    if agg.size and si.size:
        print(f"    aggregate/species alias {d['agg']} vs {d['si']}: max rel "
              f"{max_rel(agg, si):.3e}")
    if si.size:
        integ = np.sum(si * dv_cm3, axis=1)
        print(f"    combined ion source integral range=[{integ.min():.6e},{integ.max():.6e}] s^-1")
        ref = get(d['global_rate'])
        if ref.size:
            print(f"    {d['si']} -> {d['global_rate']}: max relative integral error "
                  f"{max_rel(integ, ref):.3e}")

    pi, pcx, pbal, pfl = get(d['pi']), get(d['pcx']), get(d['pbal']), get(d['pfl'])
    if pi.size and pcx.size:
        net = pi - pcx
        print(f"    mapped net plasma ion energy source {d['pi']}-{d['pcx']}: "
              f"range=[{net.min():.6e},{net.max():.6e}] W/cm3")
        if pbal.size and pfl.size:
            rhs = pi - pcx + pfl
            print(f"    neutral power balance {d['pbal']} ~= {d['pi']}-{d['pcx']}+{d['pfl']}: "
                  f"scaled residual {scaled_residual(pbal, rhs, pi, pcx, pfl):.3e}")

    tq, tqcx, tqbal, tqfl = get(d['tq']), get(d['tqcx']), get(d['tqbal']), get(d['tqfl'])
    if tq.size and tqcx.size:
        net = tq - tqcx
        print(f"    mapped net plasma torque {d['tq']}-{d['tqcx']}: "
              f"range=[{net.min():.6e},{net.max():.6e}] N m/cm3")
        if tqbal.size and tqfl.size:
            rhs = tq - tqcx + tqfl
            print(f"    neutral torque balance {d['tqbal']} ~= {d['tq']}-{d['tqcx']}+{d['tqfl']}: "
                  f"scaled residual {scaled_residual(tqbal, rhs, tq, tqcx, tqfl):.3e}")

# Recombination: species-resolved neutral birth and recapture follow the same
# decomposition used for gas/recycling in the supplied CDF.
s0, s0d = get('S0RECO'), get('S0V0D')
if s0.size and s0d.size:
    print(f'  recombination S0RECO vs S0V0D: max rel {max_rel(s0, s0d):.3e}')
si, se, sfcx = get('SIV0D'), get('SERECO'), get('SFCXRECO')
if si.size and se.size and sfcx.size:
    print(f'  recombination SIV0D = SERECO + SFCXRECO: max rel '
          f'{max_rel(si, se+sfcx):.3e}')
for name, unit in [('S0RECO','s^-1'), ('P0RECO','W'), ('TQ0RECO','N m')]:
    a = get(name)
    if a.size:
        integ = np.sum(a * dv_cm3, axis=1)
        print(f'  {name} recombination-neutral source integral range='
              f'[{integ.min():.6e},{integ.max():.6e}] {unit}')

for name in ['N0SGF_D', 'N0SRC_D', 'N0GF_D_D', 'N0RC_D_D']:
    if name in f.variables:
        units = getattr(f.variables[name], 'units', b'')
        if isinstance(units, bytes):
            units = units.decode('latin1', 'replace').strip()
        print(f'  {name}: units={units!r}; neutral state/profile, NOT an edge_sources rate')

print('\nradiation emissivity -> cumulative/global power checks')
rad_names = ['PRAD_BR', 'PRAD_CY', 'PRAD_LI']
rad_integrals = {}
for name in rad_names:
    q, inside = radiation_integral(name, dv)
    if q.size:
        rad_integrals[name] = inside[:, -1]
        nonnegative = bool(np.all(q >= -1e-12 * np.maximum(1.0, np.nanmax(np.abs(q)))))
        print(f'  {name}: finite={np.isfinite(inside).all()} nonnegative~={nonnegative} '
              f'power range=[{inside[:, -1].min():.6e},{inside[:, -1].max():.6e}] W')

if 'PRAD' in f.variables:
    q, inside = radiation_integral('PRAD', dv)
    if q.size:
        print(f'  PRAD (used total): power range=[{inside[:, -1].min():.6e},'
              f'{inside[:, -1].max():.6e}] W')

line_names = sorted(name for name in f.variables
                    if re.fullmatch(r'PRL[XS]_[A-Z0-9]+(?:_[0-9]+)?', name))
if line_names:
    print('  species-resolved line-radiation candidates:', ', '.join(line_names))
    for name in line_names:
        q, inside = radiation_integral(name, dv)
        if q.size:
            print(f'    {name}: edge power range=[{inside[:, -1].min():.6e},'
                  f'{inside[:, -1].max():.6e}] W')

f.close()
