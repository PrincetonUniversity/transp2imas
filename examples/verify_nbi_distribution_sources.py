#!/usr/bin/python
"""Numerical checks for TRANSP NBI distribution_sources/core_sources mapping.

No IMAS installation is required.  The script validates the real CDF identities
used by v8.6.1: BDEP component closure, volume-integrated particle/kinetic power,
shinethrough bookkeeping, and the beam-driven-current reconstruction.
"""

import re
import sys
from pathlib import Path

import numpy as np
from scipy import constants
from scipy.io import netcdf_file

DEFAULT_CDF = (Path(__file__).resolve().parent / 'data' / 'nstx' /
               '141711Y01' / '141711Y01.CDF')
path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_CDF
f = netcdf_file(path, 'r', mmap=False)


def get(name):
    return np.asarray(f.variables[name][:], dtype=np.float64) if name in f.variables else np.array([])


def nrms(a, b):
    den = max(float(np.sqrt(np.nanmean(np.asarray(a) ** 2))), 1.0)
    return float(np.sqrt(np.nanmean((np.asarray(a) - np.asarray(b)) ** 2)) / den)


dvol_m3 = get('DVOL') * 1.0e-6
beams = sorted(int(re.match(r'BDEP(\d{2})_TOT', name).group(1))
               for name in f.variables if re.fullmatch(r'BDEP\d{2}_TOT', name))

print('file:', path)
print('configured BDEP beam arrays:', len(beams), beams)
active = []
for ib in beams:
    s = f'{ib:02d}'
    total = get(f'BDEP{s}_TOT')
    if not total.size:
        continue
    is_active = bool(np.any(np.isfinite(total) & (total != 0.0)))
    if is_active:
        active.append(ib)

    components = [get(f'BDEP{s}_E{k}') for k in (1, 2, 3)]
    if all(a.size for a in components):
        closure = nrms(total, sum(components))
    else:
        closure = np.nan

    if not is_active:
        print(f'  beam {ib:02d}: inactive/all-zero BDEP; distribution_sources skips it')
        continue

    particles_m3s = total * 1.0e6
    particles_global = np.sum(particles_m3s * dvol_m3, axis=1)

    power_density = np.zeros_like(total, dtype=np.float64)
    shine_power = np.zeros(total.shape[0], dtype=np.float64)
    shine_particles = np.zeros(total.shape[0], dtype=np.float64)
    have_energy = True
    for k in (1, 2, 3):
        dep = get(f'BDEP{s}_E{k}')
        einj = get(f'EINJ{s}_E{k}')
        if not (dep.size and einj.size):
            have_energy = False
            continue
        power_density += dep * 1.0e6 * einj[:, None] * constants.e
        bpsh = get(f'BPSH{s}_E{k}')
        if bpsh.size:
            shine_power += bpsh
            denom = einj * constants.e
            shine_particles += np.divide(
                bpsh, denom, out=np.zeros_like(bpsh),
                where=np.isfinite(denom) & (denom > 0.0))

    power_global = np.sum(power_density * dvol_m3, axis=1) if have_energy else np.array([])
    pinj = get(f'PINJ{s}')
    active_time = pinj > 1.0e4 if pinj.size else np.zeros(total.shape[0], dtype=bool)
    ratio = np.array([])
    if power_global.size and pinj.size and np.any(active_time):
        ratio = np.divide(power_global + shine_power, pinj,
                          out=np.full_like(pinj, np.nan), where=active_time)

    line = (f'  beam {ib:02d}: BDEP parts nRMS={closure:.3e}; '
            f'particles=[{particles_global.min():.3e},{particles_global.max():.3e}] s^-1')
    if power_global.size:
        line += f'; birth power=[{power_global.min():.3e},{power_global.max():.3e}] W'
    if ratio.size:
        line += (f'; median (birth+shine)/PINJ={np.nanmedian(ratio):.6f}'
                 f' (p10-p90={np.nanpercentile(ratio,10):.4f}-{np.nanpercentile(ratio,90):.4f})')
    print(line)

print('active distribution_sources beams:', len(active), active)

curb = get('CURB')
con = get('CONPLJB')
bz = get('BZ')
bdc_names = sorted(name for name in f.variables if re.fullmatch(r'BDC\d{2}', name))
if curb.size and bdc_names:
    bdc_sum = np.zeros_like(curb)
    for name in bdc_names:
        bdc_sum += get(name)
    print(f'CURB vs sum({len(bdc_names)} BDCnn): nRMS={nrms(curb, bdc_sum):.6e}; '
          f'corr={np.corrcoef(curb.ravel(), bdc_sum.ravel())[0,1]:.9f}')

if curb.size and con.size and bz.size:
    # 174783Q32 has NLJCCW=true, so ip_sign=+1. This is the exact reconstruction
    # written to core_sources/source[nbi]/profiles_1d/j_parallel; another run
    # must apply its ip_sign.
    jparallel = curb * 1.0e4 * np.divide(
        con, bz[:, None], out=np.zeros_like(curb), where=np.abs(bz[:, None]) > 0.0)
    print('NBI <J.B>/B0 reconstruction: finite=', bool(np.isfinite(jparallel).all()),
          'range=[%.6e,%.6e] A/m2' % (jparallel.min(), jparallel.max()))

# Per-unit total NBI torque is carried into distribution_sources.  The sum of
# TQTOTnn should reproduce TRANSP's total input torque profile TQIN for this
# beam-only run to within the small bookkeeping/numerical residual.
tq_names = [f'TQTOT{ib:02d}' for ib in beams if f'TQTOT{ib:02d}' in f.variables]
if tq_names:
    tq_sum = np.zeros_like(get(tq_names[0]))
    for name in tq_names:
        tq_sum += get(name)
    tqin = get('TQIN')
    if tqin.size and tqin.shape == tq_sum.shape:
        tq_scale = max(float(np.sqrt(np.nanmean(tqin ** 2))), np.finfo(float).tiny)
        tq_rel = float(np.sqrt(np.nanmean((tqin - tq_sum) ** 2)) / tq_scale)
        print(f'TQIN vs sum({len(tq_names)} TQTOTnn): rel-RMS={tq_rel:.6e}; '
              f'corr={np.corrcoef(tqin.ravel(), tq_sum.ravel())[0,1]:.9f}')

    # Convert N*m/cm^3 -> N/m^2 and integrate with m^3 shells, exactly as in
    # imas_distribution_sources. Report the global torque range as a sanity check.
    tq_global = np.sum(tq_sum * 1.0e6 * dvol_m3, axis=1)
    print('sum(TQTOTnn) integrated torque: range=[%.6e,%.6e] N m; finite=%s' %
          (tq_global.min(), tq_global.max(), bool(np.isfinite(tq_global).all())))

# PPHI is the time/radius total angular-momentum density.  Its TRANSP unit
# Nt*m*s/cm^3 converts by 1e6 directly to kg*m^-1*s^-1, the IMAS
# core_profiles momentum density unit.  AMOM is only a zero-dimensional
# placeholder in 174783Q32 and is deliberately not used.
pphi = get('PPHI')
if pphi.size:
    pphi_si = pphi * 1.0e6
    print('PPHI -> core_profiles momentum_phi: shape=%s finite=%s range=[%.6e,%.6e] kg m^-1 s^-1' %
          (pphi.shape, bool(np.isfinite(pphi_si).all()), pphi_si.min(), pphi_si.max()))
    amoi = get('AMOI')
    omega = get('OMEGA')
    if amoi.shape == pphi.shape and omega.shape == pphi.shape:
        pred = amoi * omega
        scale = max(float(np.sqrt(np.nanmean(pphi ** 2))), np.finfo(float).tiny)
        rel = float(np.sqrt(np.nanmean((pphi - pred) ** 2)) / scale)
        corr = float(np.corrcoef(pphi.ravel(), pred.ravel())[0, 1])
        print(f'PPHI = AMOI*OMEGA identity: rel-RMS={rel:.6e}; corr={corr:.12f}')

f.close()
