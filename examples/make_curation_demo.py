#!/usr/bin/python
"""Generate a tiny, explicitly synthetic curation demonstration (no IMAS needed)."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.io import netcdf_file

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from curation import MANIFEST_SCHEMA


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=Path('build/curation-demo'))
    parser.add_argument('--force', action='store_true', help='replace existing synthetic demonstration files')
    args = parser.parse_args(argv)
    directory = args.output_dir
    directory.mkdir(parents=True, exist_ok=True)
    paths = [directory / 'synthetic.CDF', directory / 'syntheticTR.DAT', directory / 'evidence.json',
             directory / 'covariance.npz', directory / 'policy.json']
    if not args.force and any(path.exists() for path in paths):
        parser.error('demo files already exist; choose a new directory or use --force')
    x = np.tile([.1, .5, .9], (2, 1))
    data = {
        'TIME': ([0., .25, .5, 1.], ('TIME',), 'SECONDS'),
        'TIME3': ([0., 1.], ('TIME3',), 'SECONDS'),
        'X': (x, ('TIME3', 'X'), ''),
        'NE': (np.full((2, 3), 1e13), ('TIME3', 'X'), 'N/CM**3'),
        'ND': (np.full((2, 3), .8e13), ('TIME3', 'X'), 'N/CM**3'),
        'NIMP_TEST': (np.full((2, 3), .1e13), ('TIME3', 'X'), 'N/CM**3'),
        'TE': ([[950, 750, 550], [1000, 1000, 1000]], ('TIME3', 'X'), 'EV'),
        'TI': (1500 - 700 * x, ('TIME3', 'X'), 'EV'),
        'OMEGA': (np.zeros_like(x), ('TIME3', 'X'), 'RAD/SEC'),
        'DTE_DX_SIGMA': (np.full_like(x, 600), ('TIME3', 'X'), 'EV'),
        'NPTRITER': ([0, 4, 65, 2], ('TIME',), ''),
        'RESL2TE': ([0, .005, .05, .001], ('TIME',), ''),
    }
    with netcdf_file(paths[0], 'w') as nc:
        nc.createDimension('TIME', 4); nc.createDimension('TIME3', 2); nc.createDimension('X', 3)
        nc.title = 'Synthetic curation demonstration; not an experimental TRANSP result'
        for name, (values, dimensions, units) in data.items():
            variable = nc.createVariable(name, 'd', dimensions)
            variable[:] = values
            variable.units = units
            variable.long_name = 'Synthetic example: ' + name
    paths[1].write_text('! SYNTHETIC EXAMPLE ONLY\nLPREDICT_NE=0\nNLMDIF=.TRUE.\n'
                        'EINJA(1)=95000\nFFULLA(1)=0.65\nFHALFA(1)=0.2\n'
                        '~update_time=0.5\nLPREDICT_NE=1\n')
    # EXPLICITLY SYNTHETIC covariance of the synthetic TE profile: independent 15 eV errors on each of the three
    # radii, one slice per TIME3 sample.  It is not a diagnostic covariance and calibrates nothing.
    covariance = np.stack([np.diag(np.full(3, 15.0 ** 2))] * 2)
    np.savez(paths[3], time=np.array([0., 1.]), covariance=covariance)
    evidence = {
        '$schema': '../../curation/schemas/evidence.v3.schema.json',
        'schema_version': MANIFEST_SCHEMA,
        'assessment_window': {'time_min': 0, 'time_max': 1, 'rho_min': .1, 'rho_max': .9},
        'evidence': {
            'rotation': {'origin': 'assumed', 'source': 'synthetic zero-rotation scenario', 'uncertainty': 'unknown',
                         'roles': ['input'], 'notes': 'Not a measurement; scenario sensitivity remains pending.'},
            'heating': {'origin': 'manual', 'source': 'synthetic campaign record', 'uncertainty': 'unknown', 'roles': ['input']}},
        'composition': {'complete': True, 'includes_fast_ions': True, 'source': 'synthetic inventory; no fast ions',
                        'species': [{'density': 'ND', 'charge': 1}, {'density': 'NIMP_TEST', 'charge': 2}]},
        'input_expectations': [{'name': 'EINJA(1)', 'expected': 80000, 'units': 'eV', 'atol': 1000,
                                'reference': 'SYNTHETIC reference, not an approved device configuration'}],
        'gradient_uncertainties': {'TE': {'sigma_variable': 'DTE_DX_SIGMA', 'coordinate': 'X',
                                         'method': 'numpy_gradient_edge_order_1', 'source': 'synthetic derivative ensemble'}},
        'thresholds': {'newton_iterations_warning': 40,
                       'residual_limits': {'RESL2TE': {'limit': .01, 'rationale': 'Synthetic demonstration threshold only'}}},
        'run_identity': {'device': 'SYNTHETIC', 'shot': 1},
        'sample_support': {'TIME3': {'averaging_interval': 0.1, 'alignment': 'center'}},
        'events': [{'id': 'synthetic_pellet', 'type': 'pellet', 'origin': 'actuator_command',
                    'start': 0.45, 'end': 0.55, 'timing_uncertainty': 0.02, 'rules': ['PROFILE.TE.*'],
                    'applications': ['heat_transport'], 'reference': 'SYNTHETIC actuator log, not a delivered pellet',
                    'note': 'A commanded event is not proof of a delivered one.'}],
        'uncertainty': {'covariances': [{
            'id': 'synthetic_te_covariance', 'quantity': 'TE', 'kind': 'covariance', 'path': 'covariance.npz',
            'sha256': hashlib.sha256(paths[3].read_bytes()).hexdigest(), 'units': 'EV', 'coordinate': 'X',
            'ordering': 'ascending_coordinate', 'estimator': 'numpy_gradient_edge_order_1',
            'provenance': 'SYNTHETIC diagonal covariance written by make_curation_demo.py; not diagnostic data',
            'included_sources': ['synthetic independent statistical error']}]},
        'notes': 'Every value and reference in this demonstration is synthetic; do not reuse as machine inputs.'}
    paths[2].write_text(json.dumps(evidence, indent=2) + '\n')
    paths[4].write_text(json.dumps({
        '$schema': '../../curation/schemas/policy.v1.schema.json', 'schema_version': 'transp2imas.curation.policy.v1',
        'name': 'synthetic-demo-subset', 'rules': {'disable': {'numerics.runtime': 'the demo has no solver runtime channels worth reviewing'}},
        'parameters': {'gradient_sigma_multiplier': 2.0}, 'applications': ['heat_transport']}, indent=2) + '\n')
    for path in paths:
        print(path)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
