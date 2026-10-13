#!/usr/bin/python
"""Plot selected equilibrium IDS quantities after converting a TRANSP CDF.

Examples (run from repository root):
  python3 examples/plot_ids.py 38530P94 ITER HDF5 psirz --time-index 0
  python3 examples/plot_ids.py 38530P94 ITER HDF5 q --time-index 0 --save build/iter_q.png
  python3 examples/plot_ids.py 142202U01 NSTX HDF5 bfield --time-index 1

Legacy positional arguments (<runid> <database> [backend] [quantity]) are
accepted. This reads an existing IMAS database; it does not convert a CDF.
"""

import argparse
import os
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from utilities import find_runid_shot_code


FIELDS = ('psirz', 'bfield', 'q', 'psi1d', 'pressure', 'j_tor')


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('runid', help='TRANSP run identifier or CDF filename (e.g. 38530P94)')
    parser.add_argument('database', help='IMAS machine/database name (e.g. ITER, NSTX, D3D)')
    parser.add_argument('backend', nargs='?', default='HDF5', choices=('HDF5', 'HDF', 'MDSPLUS', 'MDS'))
    parser.add_argument('field', nargs='?', default='psirz', choices=FIELDS)
    parser.add_argument('--time-index', type=int, default=0,
                        help='zero-based equilibrium time_slice index (default: 0)')
    parser.add_argument('--run', type=int, default=None,
                        help='override the IMAS run number derived from the TRANSP suffix')
    parser.add_argument('--occurrence', type=int, default=0,
                        help='IMAS IDS occurrence (default: 0)')
    parser.add_argument('--user', default=os.getenv('USER'),
                        help='IMAS DBEntry user/path (default: $USER)')
    parser.add_argument('--dd-version', default=None,
                        help='IMAS DD version, e.g. 4.1.1 (default: IMAS runtime)')
    parser.add_argument('--save', type=Path, default=None,
                        help='write PNG/PDF instead of opening an interactive window')
    return parser.parse_args(argv)


def _time_slice(eq, index):
    count = len(eq.time_slice)
    if index < 0 or index >= count:
        raise IndexError(f'time index {index} outside valid range 0..{count-1}')
    return eq.time_slice[index]


def _data_2d(profile, field):
    if field == 'psirz':
        return [('Poloidal flux $\\psi$ (Wb)', np.asarray(profile.psi))]
    return [(label, np.asarray(getattr(profile, attr)))
            for label, attr in (('Radial field $B_R$ (T)', 'b_field_r'),
                                ('Toroidal field $B_\\phi$ (T)', 'b_field_tor'),
                                ('Vertical field $B_Z$ (T)', 'b_field_z'))]


def make_figure(eq, field, time_index):
    """Return a matplotlib Figure without showing it; suitable for unit tests."""
    from matplotlib import pyplot as plt

    ts = _time_slice(eq, time_index)
    time_value = None
    try:
        time_value = float(np.asarray(eq.time)[time_index])
    except (AttributeError, IndexError, ValueError, TypeError):
        pass
    time_label = (f'time slice {time_index}, t={time_value:.5g} s'
                  if time_value is not None else f'time slice {time_index}')

    if field in ('psirz', 'bfield'):
        if len(ts.profiles_2d) == 0:
            raise ValueError('equilibrium has no profiles_2d; choose q, pressure, or psi1d')
        profile = ts.profiles_2d[0]
        r = np.asarray(profile.grid.dim1, dtype=float)
        z = np.asarray(profile.grid.dim2, dtype=float)
        if r.ndim != 1 or z.ndim != 1 or r.size < 2 or z.size < 2:
            raise ValueError('equilibrium 2D grid has no usable R and Z axes')
        datasets = _data_2d(profile, field)
        fig, axes = plt.subplots(1, len(datasets), figsize=(6 * len(datasets), 5),
                                 constrained_layout=True)
        if len(datasets) == 1:
            axes = [axes]
        for ax, (label, values) in zip(axes, datasets):
            if values.shape == (r.size, z.size):
                arr = values.T
            elif values.shape == (z.size, r.size):
                arr = values
            else:
                raise ValueError(f'{label}: array shape {values.shape} incompatible '
                                 f'with grid R={r.size}, Z={z.size}')
            mesh = ax.pcolormesh(r, z, arr, shading='auto')
            fig.colorbar(mesh, ax=ax, label=label)
            if np.nanmax(arr) > np.nanmin(arr):
                ax.contour(r, z, arr, levels=12, colors='black', linewidths=0.4, alpha=0.5)
            try:
                boundary = ts.boundary.outline
                rr = np.asarray(boundary.r)
                zz = np.asarray(boundary.z)
                if rr.size > 2 and rr.shape == zz.shape:
                    ax.plot(rr, zz, 'w-', lw=1.2, label='LCFS')
            except (AttributeError, IndexError, ValueError):
                pass
            ax.set_xlabel('R (m)')
            ax.set_ylabel('Z (m)')
            ax.set_aspect('equal', adjustable='box')
            ax.set_title(label)
            ax.grid(alpha=0.15)
    else:
        field_name = {'q': 'q', 'psi1d': 'psi', 'pressure': 'pressure',
                      'j_tor': 'j_tor'}[field]
        values = np.asarray(getattr(ts.profiles_1d, field_name), dtype=float)
        if values.ndim != 1 or not values.size:
            raise ValueError(f'profiles_1d.{field_name} is not a populated 1D array')
        try:
            rho = np.asarray(ts.profiles_1d.rho_tor_norm, dtype=float)
            if rho.shape != values.shape:
                raise ValueError('rho_tor_norm grid length mismatch')
            xlabel = r'$\rho_{tor, norm}$'
        except (AttributeError, ValueError, TypeError):
            rho = np.arange(len(values))
            xlabel = 'Radial grid index'
        fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
        ax.plot(rho, values, lw=1.6)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(field_name)
        ax.grid(alpha=0.3)
    fig.suptitle(f'IMAS equilibrium — {time_label}')
    return fig


def main(argv=None):
    ns = parse_args(argv)
    inferred_run, shot, _code = find_runid_shot_code(ns.runid)
    run = inferred_run if ns.run is None else ns.run
    backend = ns.backend.upper()
    if ns.save:
        import matplotlib
        matplotlib.use('Agg')
    try:
        import imas
    except ImportError as exc:
        raise SystemExit('IMAS-Python is required to plot an existing IDS: '
                         'activate your IMAS environment first.') from exc
    backend_id = (imas.ids_defs.HDF5_BACKEND if backend in ('HDF5', 'HDF')
                  else imas.ids_defs.MDSPLUS_BACKEND)
    opts = {'dd_version': ns.dd_version} if ns.dd_version else {}
    entry = imas.DBEntry(backend_id, ns.database, shot, run, ns.user, **opts)
    print(f'Opening {ns.database} {shot}/{run}, backend={backend}, '
          f'occurrence={ns.occurrence}')
    try:
        entry.open()
        eq = entry.get('equilibrium', occurrence=ns.occurrence)
        fig = make_figure(eq, ns.field, ns.time_index)
    finally:
        entry.close()
    if ns.save:
        ns.save.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(ns.save, dpi=160)
        print('Saved:', ns.save.resolve())
    else:
        from matplotlib import pyplot as plt
        plt.show()
    from matplotlib import pyplot as plt
    plt.close(fig)
    return 0


if __name__ == '__main__':
    sys.exit(main())
