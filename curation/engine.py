#!/usr/bin/python
"""Evidence-aware screening; findings are not a scientific certification.

Rules distinguish unavailable evidence, explicit concerns, internal consistency,
scientific applicability, and converter failures. No quantity is repaired here.
"""

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
import hashlib
import json
import re

import numpy as np

from curation import REPORT_SCHEMA, __version__
from curation.manifest import validate_manifest
from curation.inputs import normalize_name
from curation.source import SourceDataError, density_scale, unit_key

APPLICATIONS = ('heat_transport', 'momentum_transport', 'current_evolution', 'surrogate_training')
ROLE_APPS = {
    'electron_temperature': ('heat_transport', 'current_evolution', 'surrogate_training'),
    'ion_temperature': ('heat_transport', 'momentum_transport', 'surrogate_training'),
    'electron_density': APPLICATIONS, 'rotation': ('momentum_transport', 'heat_transport', 'surrogate_training'),
    'composition': APPLICATIONS, 'heating': APPLICATIONS,
}
STATUSES = ('pass', 'warning', 'fail', 'not_assessed', 'not_applicable', 'error')
PROFILES = {'NE': 'electron_density', 'TE': 'electron_temperature', 'TI': 'ion_temperature'}
SOLVER_CHANNELS = ('NPTRITER', 'PT_CPTIME', 'PT_WTIME', 'RESL2TE', 'RESL2TI', 'RESL2OMG',
                   'UPWIND_TE', 'UPWIND_TI', 'UPWIND_MO')


@dataclass
class Finding:
    id: str
    status: str
    message: str
    variables: list = field(default_factory=list)
    applications: list = field(default_factory=lambda: list(APPLICATIONS))
    evidence: dict = field(default_factory=dict)
    metrics: dict = field(default_factory=dict)
    scope: dict = field(default_factory=dict)
    recommendation: str = ''


@dataclass
class Report:
    source: dict
    manifest: dict
    timeline: dict
    findings: list = field(default_factory=list)
    masks: dict = field(default_factory=dict, repr=False)
    coordinates: dict = field(default_factory=dict, repr=False)
    channel_inventory: dict = field(default_factory=dict)
    read_counts: dict = field(default_factory=dict)

    def add(self, identifier, status, message, **kwargs):
        if status not in STATUSES or any(f.id == identifier for f in self.findings):
            raise ValueError(f'invalid or duplicate finding: {identifier}, {status}')
        finding = Finding(identifier, status, message, **kwargs)
        self.findings.append(finding)
        return finding

    def application_status(self, application):
        findings = [f for f in self.findings if application in f.applications]
        if any(f.status in ('fail', 'error') for f in findings):
            return 'blocked_by_screening_findings'
        if any(f.status == 'warning' for f in findings):
            return 'review_required'
        if any(f.status == 'not_assessed' for f in findings):
            return 'insufficient_evidence'
        return 'screening_only_not_certified'

    def as_dict(self):
        counts = {s: sum(f.status == s for f in self.findings) for s in STATUSES}
        return {'schema_version': REPORT_SCHEMA, 'curation_version': __version__,
                'generated_at_utc': datetime.now(timezone.utc).isoformat(),
                'source': self.source,
                'manifest_sha256': hashlib.sha256(json.dumps(self.manifest, sort_keys=True,
                                                             allow_nan=False).encode()).hexdigest(),
                'effective_manifest': self.manifest,
                'declared_input_timeline': self.timeline,
                'summary': {'counts': counts, 'scientifically_certified': False,
                            'scope': 'read-only source screening; no IDS backend validation',
                            'applications': {a: self.application_status(a) for a in APPLICATIONS}},
                'findings': [asdict(f) for f in self.findings],
                'channel_inventory': self.channel_inventory, 'source_read_counts': self.read_counts,
                'mask_storage': {'format': 'numpy_npz',
                                 'keys': sorted(set(self.masks) | set(self.coordinates)),
                                 'semantics': 'evaluated masks distinguish unassessed samples from unflagged samples'}}


class CurationEngine:
    """A standalone, injectable engine; input evidence remains analyst-declared."""

    def __init__(self, source, manifest, timeline=None):
        self.source = source
        self.manifest = validate_manifest(manifest)
        self.timeline = timeline or {'available': False, 'states': [], 'gaps': [], 'source': None}
        self.report = Report(source.provenance(), self.manifest, self.timeline)
        self.window = self.manifest['assessment_window']
        self.thresholds = self.manifest['thresholds']

    def _selection(self, name):
        """Select native samples; never infer a grid from equal array lengths."""
        data = self.source.array(name)
        meta = self.source.meta(name)
        dims = tuple(meta['dimensions'])
        mask = np.ones(data.shape, dtype=bool)
        for axis, dim in enumerate(dims):
            if dim not in ('TIME', 'TIME3'):
                continue
            if not self.source.has(dim):
                raise ValueError(f'{name}: missing native time coordinate {dim}')
            times = self.source.array(dim)
            if (times.ndim != 1 or len(times) != data.shape[axis] or
                    not np.all(np.isfinite(times)) or not np.all(np.diff(times) > 0)):
                raise ValueError(f'{name}: invalid native coordinate {dim}')
            if unit_key(self.source.meta(dim)['units']) not in ('S', 'SEC', 'SECS', 'SECOND', 'SECONDS'):
                raise ValueError(f'{dim}: native time units must be seconds')
            selected = np.ones(times.shape, dtype=bool)
            if 'time_min' in self.window:
                selected &= times >= self.window['time_min']
            if 'time_max' in self.window:
                selected &= times <= self.window['time_max']
            reshape = [1] * data.ndim
            reshape[axis] = len(times)
            mask &= selected.reshape(reshape)
            self.report.coordinates['coordinate_' + dim] = times
        # Save an explicit native radial coordinate even without a radial cut.
        for radial in ('X', 'XB'):
            if radial in dims and self.source.has(radial):
                radius = self.source.array(radial)
                if radius.shape == data.shape and tuple(self.source.meta(radial)['dimensions']) == dims:
                    self.report.coordinates['coordinate_' + radial] = radius
        if any(key.startswith('rho_') for key in self.window) and data.ndim >= 2:
            radial = next((d for d in dims if d in ('X', 'XB')), None)
            if radial is None or not self.source.has(radial):
                raise ValueError(f'{name}: requested radial window has no explicit X/XB coordinate')
            radius = self.source.array(radial)
            if radius.shape != data.shape or tuple(self.source.meta(radial)['dimensions']) != dims:
                raise ValueError(f'{name}: radial window requires the matching native X/XB grid')
            if not self._normalized_radius(radial):
                raise ValueError(f'{radial}: normalized-radius unit is not established')
            if 'rho_min' in self.window:
                mask &= radius >= self.window['rho_min']
            if 'rho_max' in self.window:
                mask &= radius <= self.window['rho_max']
            mask &= np.isfinite(radius)
            self.report.coordinates['coordinate_' + radial] = radius
        return mask

    def _normalized_radius(self, name):
        return unit_key(self.source.meta(name)['units']) in ('', '1', '-', 'DIMENSIONLESS')

    def _scope(self, name, identifier, evaluated, flagged):
        key = re.sub('[^A-Za-z0-9_]', '_', identifier)
        self.report.masks[key + '__evaluated'] = np.asarray(evaluated, dtype=bool)
        self.report.masks[key + '__flagged'] = np.asarray(flagged, dtype=bool)
        points = np.argwhere(flagged)
        return {'source_dimensions': self.source.meta(name)['dimensions'],
                'shape': list(evaluated.shape), 'evaluated_count': int(evaluated.sum()),
                'flagged_count': int(flagged.sum()), 'mask_prefix': key,
                'flagged_index_examples': points[:12].tolist(), 'index_examples_truncated': len(points) > 12,
                'assessment_window': self.window}

    def _samples(self, identifier, name, evaluated, flagged, message, *, severity='warning',
                 applications=APPLICATIONS, metrics=None, recommendation=''):
        n = int(evaluated.sum())
        status = 'not_assessed' if n == 0 else severity if np.any(flagged) else 'pass'
        return self.report.add(identifier, status, message if n else 'No assessable samples in the selected window.',
                               variables=[name], applications=list(applications), metrics=metrics or {},
                               scope=self._scope(name, identifier, evaluated, flagged), recommendation=recommendation)

    def _unavailable(self, identifier, message, names=(), apps=APPLICATIONS):
        return self.report.add(identifier, 'not_assessed', message, variables=list(names), applications=list(apps))

    def _source_problem(self, identifier, exc, names=(), apps=APPLICATIONS):
        # A backend/read failure is not merely absent scientific evidence.
        if isinstance(exc, (OSError, RuntimeError, TypeError, SourceDataError)):
            return self.report.add(identifier, 'error', str(exc), variables=list(names), applications=list(apps))
        return self._unavailable(identifier, str(exc), names, apps)

    def run(self):
        self.check_coordinates()
        self.check_profiles()
        self.check_evidence()
        self.check_rotation()
        self.check_composition()
        self.check_declared_inputs()
        self.check_numerics()
        self.report.read_counts = dict(self.source.read_counts)
        return self.report

    def check_coordinates(self):
        for name in ('TIME', 'TIME3'):
            identifier = 'GRID.' + name
            if not self.source.has(name):
                self._unavailable(identifier, 'Native time coordinate is absent.', [name])
                continue
            try:
                values = self.source.array(name)
                units = unit_key(self.source.meta(name)['units'])
                good = (values.ndim == 1 and values.size > 0 and np.all(np.isfinite(values))
                        and np.all(np.diff(values) > 0))
                self.report.add(identifier, 'pass' if good else 'fail',
                                'Native time coordinate must be finite and strictly increasing.', variables=[name],
                                metrics={'sample_count': int(values.size)},
                                evidence={'units': self.source.meta(name)['units']})
                if units not in ('S', 'SEC', 'SECS', 'SECOND', 'SECONDS'):
                    self._unavailable(identifier + '.UNITS', 'Time units are not recognized as seconds.', [name])
            except (ValueError, TypeError, OSError, RuntimeError) as exc:
                self._source_problem(identifier, exc, [name])
        if not self.source.has('X'):
            self._unavailable('GRID.X', 'No zone-centre normalized-radius grid.', ['X'])
            return
        try:
            x = self.source.array('X')
            if x.ndim != 2 or tuple(self.source.meta('X')['dimensions']) != ('TIME3', 'X'):
                raise ValueError('X is not a native (TIME3, X) grid')
            if not self._normalized_radius('X'):
                raise ValueError('X units do not establish a normalized radius')
            selected = self._selection('X')
            invalid = ~np.isfinite(x) | (x < -1e-6) | (x > 1 + 1e-6)
            invalid[:, 1:] |= np.diff(x, axis=1) <= 0
            self._samples('GRID.X', 'X', selected, selected & invalid,
                          'Normalized X must be finite, increasing and within [0, 1].', severity='fail')
            self.report.coordinates['coordinate_X'] = x
        except (ValueError, TypeError, OSError, RuntimeError) as exc:
            self._source_problem('GRID.X', exc, ['X'])

    def check_profiles(self):
        for name, role in PROFILES.items():
            apps = ROLE_APPS[role]
            if not self.source.has(name):
                self._unavailable('PROFILE.' + name, 'Profile is absent; no substitute is assumed.', [name], apps)
                continue
            try:
                data = self.source.array(name)
                if data.ndim != 2 or tuple(self.source.meta(name)['dimensions']) != ('TIME3', 'X'):
                    raise ValueError('Expected a native (TIME3, X) profile; no shape-based reinterpretation')
                if self.source.has('X') and data.shape != self.source.array('X').shape:
                    self.report.add('PROFILE.' + name + '.GRID_MISMATCH', 'fail',
                                    'Profile shape disagrees with the explicit native X grid.',
                                    variables=[name, 'X'], applications=list(apps))
                    continue
                selected = self._selection(name)
                self._samples('PROFILE.' + name + '.FINITE', name, selected, selected & ~np.isfinite(data),
                              'Missing, masked or nonfinite profile samples require review.', applications=apps)
                finite = selected & np.isfinite(data)
                self._samples('PROFILE.' + name + '.NEGATIVE', name, finite, finite & (data < 0),
                              'Negative density/temperature samples are not physically admissible.',
                              severity='fail', applications=apps)
                self._samples('PROFILE.' + name + '.ZERO', name, finite, finite & (data == 0),
                              'Zero values may be startup/off-plasma samples; do not use them as transport denominators.',
                              applications=apps)
                units = self.source.meta(name)['units']
                if name == 'NE':
                    density_scale(units)
                elif unit_key(units) not in ('EV', 'KEV'):
                    raise ValueError(f'{name}: unknown temperature units {units!r}')
                self.report.add('PROFILE.' + name + '.UNITS', 'pass', 'Native units are recognized.',
                                variables=[name], applications=list(apps), evidence={'units': units})
                if name in ('TE', 'TI'):
                    self.check_gradient(name)
            except (ValueError, TypeError, OSError, RuntimeError) as exc:
                self._source_problem('PROFILE.' + name + '.ASSESSMENT', exc, [name], apps)

    def check_gradient(self, name):
        apps = ('heat_transport', 'surrogate_training')
        ident = 'FIT.' + name
        if not self.source.has('X'):
            self._unavailable(ident + '.GRADIENT', 'No native radial grid.', [name, 'X'], apps)
            return
        data, x = self.source.array(name), self.source.array('X')
        if x.shape != data.shape or not self._normalized_radius('X'):
            self._unavailable(ident + '.GRADIENT', 'Profile and normalized radial grid are not aligned.', [name, 'X'], apps)
            return
        scale = 1000.0 if unit_key(self.source.meta(name)['units']) == 'KEV' else 1.0
        y = data * scale
        gradient = np.full(y.shape, np.nan)
        limit = np.full((len(y), 1), np.nan)
        for row in range(len(y)):
            if (y.shape[1] >= 2 and np.all(np.isfinite(y[row])) and np.all(np.isfinite(x[row]))
                    and np.all(np.diff(x[row]) > 0)):
                gradient[row] = np.gradient(y[row], x[row], edge_order=1)
                limit[row, 0] = (self.thresholds['weak_gradient_relative'] * np.max(np.abs(y[row]))
                                 / (x[row, -1] - x[row, 0]))
        evaluated = self._selection(name) & np.isfinite(gradient)
        weak = evaluated & (np.abs(gradient) <= limit)
        self._samples(ident + '.WEAK_GRADIENT', name, evaluated, weak,
                      'Small dT/dX makes a local diffusive interpretation ill-conditioned; this is not a bad-profile verdict.',
                      applications=apps, metrics={'relative_screening_threshold': self.thresholds['weak_gradient_relative']},
                      recommendation='Check flux and coefficient validity separately; assess gradient uncertainty.')
        spec = self.manifest['gradient_uncertainties'].get(name)
        if spec is None:
            self._unavailable(ident + '.GRADIENT_UNCERTAINTY',
                              'No supplied uncertainty for this derivative; none is inferred from pointwise errors.', [name], apps)
            return
        sigma_name = spec['sigma_variable']
        if not self.source.has(sigma_name):
            self._unavailable(ident + '.GRADIENT_UNCERTAINTY', 'Declared derivative uncertainty is absent.', [sigma_name], apps)
            return
        sigma = self.source.array(sigma_name)
        units = unit_key(self.source.meta(sigma_name)['units'])
        if (sigma.shape != y.shape or self.source.meta(sigma_name)['dimensions'] != self.source.meta(name)['dimensions']
                or units not in ('EV', 'KEV')):
            self._unavailable(ident + '.GRADIENT_UNCERTAINTY',
                              'Derivative uncertainty must share native dimensions and have eV/keV units per dimensionless X.',
                              [name, sigma_name], apps)
            return
        sigma = sigma * (1000.0 if units == 'KEV' else 1.0)
        bad_sigma = evaluated & (~np.isfinite(sigma) | (sigma < 0))
        self._samples(ident + '.SIGMA_VALIDITY', name, evaluated, bad_sigma,
                      'Derivative standard deviations must be finite and nonnegative.', severity='fail', applications=apps)
        usable = evaluated & np.isfinite(sigma) & (sigma >= 0)
        k = self.thresholds['gradient_sigma_multiplier']
        finding = self._samples(ident + '.GRADIENT_UNCERTAINTY', name, usable,
                               usable & (np.abs(gradient) <= k * sigma),
                               'The declared gradient uncertainty interval includes zero; coefficient sign/magnitude is not resolved.',
                               applications=apps, metrics={'sigma_multiplier': k})
        finding.evidence = dict(spec)

    def check_evidence(self):
        for role, apps in ROLE_APPS.items():
            record = self.manifest['evidence'].get(role, {})
            origin = record.get('origin', 'unknown')
            source = record.get('source', '')
            if origin in ('unknown', 'missing') or not source:
                status, message = 'not_assessed', 'Origin/reference is missing or unknown; automatic retrieval is not validation.'
            elif record.get('reliability') == 'unreliable':
                status, message = 'warning', 'Evidence is explicitly marked unreliable.'
            elif origin in ('assumed', 'inferred', 'predicted'):
                status, message = 'warning', 'This quantity is assumption/model dependent, not an independent measurement.'
            else:
                status, message = 'pass', 'An origin and reference were declared; their accuracy was not independently established.'
            self.report.add('EVIDENCE.' + role, status, message, applications=list(apps), evidence=record)
            uncertainty = record.get('uncertainty', 'unknown')
            self.report.add('UNCERTAINTY.' + role, 'pass' if uncertainty == 'provided' else 'not_assessed',
                            'Uncertainty is declared as provided; propagation/calibration remains unverified.' if uncertainty == 'provided'
                            else 'No usable uncertainty specification; unknown is not zero.', applications=list(apps))
        # Shared source strings are only a conservative exact-match warning, not a proof of independence.
        evidence = self.manifest['evidence']
        for role, record in evidence.items():
            if not record.get('independent_validation'):
                continue
            shared = [other for other, item in evidence.items()
                      if item.get('source') and item.get('source') == record.get('source')
                      and set(item.get('roles', ())) & {'input', 'calibration', 'initial_condition', 'boundary_condition'}]
            if shared:
                self.report.add('INDEPENDENCE.' + role, 'warning',
                                'The claimed independent validator shares a source used to constrain the run.',
                                applications=list(ROLE_APPS[role]), evidence={'shared_with': shared, 'source': record.get('source')})

    def check_rotation(self):
        names = [name for name in ('OMEGA', 'VTOR') if self.source.has(name)]
        apps = ('momentum_transport', 'heat_transport', 'surrogate_training')
        record = self.manifest['evidence'].get('rotation', {})
        if not names:
            self.report.add('ROTATION.AVAILABILITY', 'warning',
                            'No recognized rotation output; do not replace missing rotation by zero.',
                            applications=list(apps), recommendation='Provide measurement provenance or explicit rotation scenarios.')
        else:
            self.report.add('ROTATION.AVAILABILITY', 'pass',
                            'Rotation output exists; its presence does not establish measurement provenance.',
                            variables=names, applications=list(apps), evidence=record)

    def check_composition(self):
        ident = 'COMPOSITION.QUASINEUTRALITY'
        comp = self.manifest.get('composition')
        if not comp or comp.get('complete') is not True or comp.get('includes_fast_ions') is not True or not comp.get('species'):
            self._unavailable(ident, 'An explicit complete charged-species inventory, including fast-ion accounting, is required.')
            return
        ne_name = comp.get('electron_density', 'NE')
        names = [ne_name] + [s['density'] for s in comp['species']]
        names += [s['charge_variable'] for s in comp['species'] if 'charge_variable' in s]
        missing = [name for name in names if not self.source.has(name)]
        if missing:
            self._unavailable(ident, 'Declared composition variables are missing: ' + ', '.join(missing), missing)
            return
        try:
            ne = self.source.array(ne_name) * density_scale(self.source.meta(ne_name)['units'])
            dims = self.source.meta(ne_name)['dimensions']
            if ne.ndim != 2 or tuple(dims) != ('TIME3', 'X'):
                raise ValueError('Electron density must be on the native (TIME3, X) grid')
            charge_density = np.zeros_like(ne)
            valid = np.isfinite(ne) & (ne > 0)
            for species in comp['species']:
                name = species['density']
                density = self.source.array(name) * density_scale(self.source.meta(name)['units'])
                if density.shape != ne.shape or self.source.meta(name)['dimensions'] != dims:
                    raise ValueError(f'{name}: native composition grids do not match; no interpolation is assumed')
                charge = species.get('charge')
                if charge is None:
                    charge_name = species['charge_variable']
                    charge = self.source.array(charge_name)
                    if charge.shape != ne.shape or self.source.meta(charge_name)['dimensions'] != dims:
                        raise ValueError(f'{charge_name}: charge-state grid mismatch')
                    if unit_key(self.source.meta(charge_name)['units']) not in ('', '1', '-', 'DIMENSIONLESS'):
                        raise ValueError(f'{charge_name}: charge-state units are not dimensionless')
                valid &= np.isfinite(density) & (density >= 0) & np.isfinite(charge) & (charge > 0)
                charge_density += density * charge
            selected = self._selection(ne_name)
            self._samples(ident + '.INPUTS', ne_name, selected, selected & ~valid,
                          'Some composition samples are unavailable/nonphysical and cannot be tested for charge closure.')
            valid &= selected
            relative = np.full(ne.shape, np.nan)
            np.divide(np.abs(charge_density - ne), ne, out=relative, where=valid)
            tolerance = self.thresholds['quasineutrality_relative']
            metrics = {'relative_tolerance': tolerance, 'tolerance_kind': 'screening_not_statistical'}
            if valid.any():
                metrics['max_relative_residual'] = float(np.max(relative[valid]))
            finding = self._samples(ident, ne_name, valid, valid & (relative > tolerance),
                                   'Charge closure tested for the declared inventory; closure does not validate composition or Zeff.',
                                   metrics=metrics)
            finding.variables = names
            finding.evidence = {'source': comp.get('source', ''), 'inventory_is_analyst_declared': True}
        except (ValueError, TypeError, OSError, RuntimeError) as exc:
            self._source_problem(ident, exc, names)

    def _states_in_window(self):
        for state in self.timeline.get('states', []):
            start, end = state['start_time'], state['end_time']
            if end is not None and 'time_min' in self.window and end <= self.window['time_min']:
                continue
            if start is not None and 'time_max' in self.window and start > self.window['time_max']:
                continue
            yield state

    def check_declared_inputs(self):
        if not self.timeline.get('available'):
            self._unavailable('INPUT.TIMELINE', 'No TR.DAT supplied; controls and default values remain unknown.')
        else:
            gaps = self.timeline.get('gaps', [])
            self.report.add('INPUT.TIMELINE', 'warning' if gaps else 'pass',
                            'Only explicit scalar declarations were indexed; this is not the resolved TRANSP configuration.',
                            evidence={'parser_gaps': gaps, 'source': self.timeline.get('source')})
        temporal = self.timeline.get('temporal_valid', True)
        states = list(self._states_in_window()) if temporal else []
        coupled = []
        conflicts = []
        fractions = []
        for state in states:
            vals = state['values']
            density = [k for k in ('LPREDICT_NE', 'LPREDICT_NMAIN', 'LPREDICT_NIMP')
                       if isinstance(vals.get(k), (int, float)) and not isinstance(vals[k], bool) and vals[k] > 0]
            diffusion = [k for k in ('NLMDIF', 'NLISODIF') if vals.get(k) is True or (isinstance(vals.get(k), (int, float)) and vals[k] == 1)]
            if density and diffusion:
                coupled.append({'start_time': state['start_time'], 'end_time': state['end_time'],
                                'density_controls': density, 'diffusion_controls': diffusion})
            if vals.get('NLVWNC') in (True, 1) and vals.get('NLIVPO') in (True, 1):
                conflicts.append({'start_time': state['start_time'], 'end_time': state['end_time']})
            for name, value in vals.items():
                if re.fullmatch(r'FFULLA(?:\(\d+\))?', name):
                    other = name.replace('FFULLA', 'FHALFA')
                    half = vals.get(other)
                    if (isinstance(value, (int, float)) and not isinstance(value, bool)
                            and isinstance(half, (int, float)) and not isinstance(half, bool)):
                        fractions.append({'full_key': name, 'half_key': other, 'full': value, 'half': half,
                                          'start_time': state['start_time'], 'end_time': state['end_time'],
                                          'invalid': not (0 <= value <= 1 and 0 <= half <= 1 and value + half <= 1 + 1e-9)})
        if not temporal:
            # no declared state could be placed in time: say so instead of reporting a clean "nothing found"
            reason = ('Not evaluated: the TR.DAT update times are unusable, so no declared state could be placed '
                      'in time (see INPUT.TIMELINE).')
            for identifier in ('MODEL.DENSITY_CURRENT_COUPLING', 'MODEL.ROTATION_CONTROLS', 'HEATING.FRACTIONS'):
                self._unavailable(identifier, reason)
        else:
            self._report_declared_models(coupled, conflicts, fractions)
        self._report_heating_reference(states)

    def _report_declared_models(self, coupled, conflicts, fractions):
        self.report.add('MODEL.DENSITY_CURRENT_COUPLING', 'warning' if coupled else 'not_assessed',
                        'Density and magnetic diffusion are both declared active; inspect composition, resistivity, boundaries and update coupling.'
                        if coupled else 'No simultaneous activation established from explicit supported controls; consistency is not certified.',
                        evidence={'intervals': coupled}, recommendation='This combination is not prohibited; document the closure.')
        self.report.add('MODEL.ROTATION_CONTROLS', 'warning' if conflicts else 'not_assessed',
                        'Both radial-potential choices are declared; verify compatibility for the actual TRANSP version.'
                        if conflicts else 'No explicit double selection found; omitted/default controls are not resolved.',
                        evidence={'intervals': conflicts}, applications=list(ROLE_APPS['rotation']))
        self.report.add('HEATING.FRACTIONS', 'fail' if any(x['invalid'] for x in fractions) else 'pass' if fractions else 'not_assessed',
                        'FFULLA/FHALFA are current fractions; the implied third fraction must be nonnegative. No power-fraction equivalence is assumed.',
                        evidence={'declarations': fractions})

    def _report_heating_reference(self, states):
        expectations = self.manifest['input_expectations']
        if not expectations:
            self._unavailable('HEATING.REFERENCE', 'No campaign/configuration reference expectations were supplied.')
        for index, item in enumerate(expectations):
            key = normalize_name(item['name'])
            records = []
            missing = 0
            for state in states:
                value = state['values'].get(key)
                if not isinstance(value, (int, float)) or isinstance(value, bool):
                    missing += 1
                    continue
                matches = abs(value - item['expected']) <= item['atol'] + item['rtol'] * abs(item['expected'])
                records.append({'start_time': state['start_time'], 'end_time': state['end_time'],
                                'value': value, 'matches_reference': matches,
                                'source_line': state['source_lines'].get(key)})
            status = ('warning' if any(not r['matches_reference'] for r in records)
                      else 'not_assessed' if not records or missing else 'pass')
            self.report.add(f'HEATING.REFERENCE.{index}', status,
                            'Explicit native-unit declarations compared to the supplied reference; deviations are not automatically corrected.',
                            variables=[key], evidence={'expectation': item, 'intervals': records, 'unknown_intervals': missing})

    def check_numerics(self):
        present = []
        for name in SOLVER_CHANNELS:
            if not self.source.has(name):
                continue
            try:
                values = self.source.array(name)
                selected = self._selection(name)
                finite = selected & np.isfinite(values)
                entry = self.source.meta(name)
                entry.update({'selected_samples': int(selected.sum()), 'finite_samples': int(finite.sum()),
                              'nonzero_samples': int(np.count_nonzero(values[finite]))})
                if finite.any():
                    entry.update({'minimum': float(values[finite].min()), 'maximum': float(values[finite].max())})
                self.report.channel_inventory[name] = entry
                present.append(name)
                self._samples('NUMERICS.' + name + '.FINITE', name, selected, selected & ~np.isfinite(values),
                              'Runtime diagnostic contains unavailable/nonfinite samples; no convergence inference.')
                if name == 'NPTRITER':
                    invalid = (values < 0) | (np.abs(values - np.round(values)) > 1e-6)
                    self._samples('NUMERICS.ITERATION_DOMAIN', name, finite, finite & invalid,
                                  'Iteration counts should be nonnegative integers; zero is not proof of convergence.', severity='fail')
                    limit = self.thresholds.get('newton_iterations_warning')
                    if limit is not None:
                        self._samples('NUMERICS.ITERATION_LIMIT', name, finite, finite & (values > limit),
                                      'Iteration count exceeds the analyst screening threshold; this is not the solver stopping criterion.',
                                      metrics={'threshold': limit})
                elif name.startswith('RESL2'):
                    self._samples('NUMERICS.' + name + '.DOMAIN', name, finite, finite & (values < 0),
                                  'An L2 residual should not be negative.', severity='fail')
                    spec = self.thresholds['residual_limits'].get(name)
                    if spec:
                        finding = self._samples('NUMERICS.' + name + '.THRESHOLD', name, finite,
                                                finite & (values > spec['limit']),
                                                'Residual exceeds the supplied native-normalization threshold; check its rationale.',
                                                metrics={'threshold': spec['limit']})
                        finding.evidence = dict(spec)
            except (ValueError, TypeError, OSError, RuntimeError) as exc:
                self._source_problem('NUMERICS.' + name + '.READ', exc, [name])
        self.report.add('NUMERICS.AVAILABILITY', 'pass' if present else 'not_assessed',
                        'Runtime channels inventoried on their own coordinates; upwind adjustment is not a numerical heat flux.',
                        variables=present)
        self._unavailable('NUMERICS.CONVERGENCE',
                          'Per-equation stopping criteria, rejected steps and time-step/stabilization sensitivity are not certified by this first-stage screen.')
