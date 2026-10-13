#!/usr/bin/python
"""Populate IMAS workflow/provenance metadata describing the TRANSP conversion."""

import numpy as np

from nc2imas import nc2imas


class imas_workflow:
    """Describe the TRANSP simulation and the source inputs that produced it."""

    def __init__(self, ncimas: nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger
        self.create_ids()

    def _attr(self, name):
        return self.im.get_ncattr_text(name, default='')

    def create_ids(self):
        self.logger.info('\nworkflow')
        try:
            ids = self.im.new_ids('workflow')
        except Exception as exc:
            self.logger.warning(f'workflow is unavailable in target DD: {exc}')
            return

        ids = self.im.fill_properties(ids)
        ids = self.im.fill_code_properties(ids)
        if hasattr(ids.ids_properties, 'homogeneous_time'):
            ids.ids_properties.homogeneous_time = 2

        if hasattr(ids, 'time_loop') and hasattr(ids.time_loop, 'component'):
            inventory = self.im.input_inventory()
            models = inventory.models()
            ids.time_loop.component.resize(1 + len(models))
            component = ids.time_loop.component[0]
            if hasattr(component, 'name'):
                component.name = 'TRANSP'
            if hasattr(component, 'description'):
                component.description = 'TRANSP integrated tokamak simulation producing the source CDF'
            version = self._attr('TRANSP_version')
            if version and hasattr(component, 'version'):
                component.version = version
            if hasattr(component, 'parameters'):
                component.parameters = self.im.build_code_parameters_xml()

            # Record referenced model input sets independently, without
            # claiming that every configured code was actually executed.
            for i, name in enumerate(models, start=1):
                external = ids.time_loop.component[i]
                if hasattr(external, 'name'):
                    external.name = name
                if hasattr(external, 'description'):
                    external.description = ('Referenced TRANSP model input set; '
                                            'execution not independently verified')
                if hasattr(external, 'parameters'):
                    external.parameters = inventory.scoped_xml(name)

            # Preserve model-code hashes exposed by TRANSP as workflow libraries.
            libraries = []
            for name, key, description in (
                ('TGLF', 'TGLF_hash', 'TGLF turbulent transport model referenced by TRANSP'),
                ('GLF23', 'GLF23_hash', 'GLF23 turbulent transport model referenced by TRANSP'),
            ):
                commit = self._attr(key)
                if commit:
                    libraries.append((name, commit, description))
            if libraries and hasattr(component, 'library'):
                component.library.resize(len(libraries))
                for i, (name, commit, description) in enumerate(libraries):
                    lib = component.library[i]
                    if hasattr(lib, 'name'):
                        lib.name = name
                    if hasattr(lib, 'description'):
                        lib.description = description
                    if hasattr(lib, 'commit'):
                        lib.commit = commit

        times = np.asarray(self.im.timesc, dtype=np.float64)
        if times.size and hasattr(ids, 'time_loop') and hasattr(ids.time_loop, 'time_end'):
            ids.time_loop.time_end = float(times[-1])

        # workflow_cycle is intentionally empty. A post-processed TRANSP CDF
        # does not retain the actual component execution schedule per time step.
        self.im.putdb_imas(ids)


if __name__ == '__main__':
    print('Please use: python convert.py')
