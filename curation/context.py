#!/usr/bin/python
"""The read-only inputs every rule sees: source, evidence, declared timeline, events, thresholds and identity.

The context never initializes IMAS, never modifies a source array, and offers only derived read-only views (native
spans, normalized-radius test, declared states in the assessment window).  Rules write only through the recorder.
"""

import numpy as np

from curation.source import unit_key


class CurationContext:
    def __init__(self, source, manifest, timeline=None, parameters=None):
        self.source = source
        self.manifest = manifest
        self.timeline = timeline or {'available': False, 'states': [], 'gaps': [], 'source': None}
        self.window = manifest['assessment_window']
        # Screening thresholds: the sidecar's values, overridden by an explicit policy; both are recorded in the report.
        self.thresholds = {**manifest['thresholds'], **(parameters or {})}
        self.events = manifest['events']
        self.identity, self.identity_conflicts = self._resolve_identity()

    def _resolve_identity(self):
        """Device, shot and campaign used to decide applicability.

        The analyst's declaration is used only where it does not contradict the CDF.  A declaration that
        disagrees with the source's own attribute is not trusted, and neither is the source value taken in its
        place: the item becomes unknown, so no reference is authorized by a conflicting identity.
        """
        declared = dict(self.manifest['run_identity'])
        found = {'device': self.source.global_attribute('device'), 'shot': self.source.global_attribute('shot')}
        identity, conflicts = {}, []
        for key in ('device', 'shot', 'campaign'):
            have, want = found.get(key), declared.get(key)
            if have is not None and want is not None:
                same = (str(have).strip().lower() == str(want).strip().lower() if key == 'device'
                        else int(have) == int(want))
                if not same:
                    conflicts.append(key)
                    continue
            value = want if want is not None else have
            if value is not None:
                identity[key] = value
        return identity, conflicts

    def normalized_radius(self, name):
        return unit_key(self.source.meta(name)['units']) in ('', '1', '-', 'DIMENSIONLESS')


    def native_spans(self, clocks=('TIME', 'TIME3')):
        """``(first, last)`` of each usable native clock; an unusable or absent clock contributes nothing."""
        spans = []
        for dim in clocks:
            if self.source.has(dim):
                times = self.source.array(dim)
                if times.ndim == 1 and times.size and np.all(np.isfinite(times)):
                    spans.append((float(times.min()), float(times.max())))
        return spans

    def states_in_window(self):
        for state in self.timeline.get('states', []):
            start, end = state['start_time'], state['end_time']
            if end is not None and 'time_min' in self.window and end <= self.window['time_min']:
                continue
            if start is not None and 'time_max' in self.window and start > self.window['time_max']:
                continue
            yield state
