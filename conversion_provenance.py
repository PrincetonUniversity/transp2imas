#!/usr/bin/python
"""Provenance of one conversion: input inventory, source-file records, ``ids_properties`` and ``code``.

``ProvenanceRecorder`` owns what used to be three caches on the converter (the input inventory, the
source-file records and the ``code.parameters`` XML) and builds the provenance blocks every IDS carries.
It reads the converter's state at call time (file names, namelist, DD version, signs, run identity), so a
converter that changes its namelist calls :meth:`reset` and the next call rebuilds from the new inputs.
The text and XML rules themselves are in ``code_provenance`` and ``input_provenance``.

``source`` is the converter, or any object offering: ``ncfile``, ``nml_name``, ``device``, ``shot``,
``run``, ``occur``, ``dd_version``, ``signs``, ``logger``, ``mapping``, ``get_ncattr``, ``get_ncattr_text``
and ``get_nc_global_atts``.
"""

from datetime import datetime, timezone

from code_provenance import (TRANSP_METADATA_KEYS, attach_source_provenance, build_code_parameters_xml,
                             properties_comment, source_file_records)
from input_provenance import InputInventory, fill_model_libraries

REPOSITORY = 'https://github.com/PrincetonUniversity/transp2imas'


class ProvenanceRecorder:
    """Build and cache the provenance blocks of the conversion of ``source``."""

    def __init__(self, source, converter_version):
        self.source = source
        self.converter_version = converter_version
        self.reset()

    def reset(self):
        """Forget what was built from the previous inputs (a changed TR.DAT must not keep the old XML)."""
        self._code_parameters_xml = None
        self._file_records = None
        self._inventory = None

    def input_inventory(self):
        """Resolve and cache the run's TR.DAT and recursively referenced inputs.

        Missing or oversized files remain explicit provenance records; they are
        *not* silently counted as archived inputs.
        """
        source = self.source
        if self._inventory is None:
            self._inventory = InputInventory(source.nml_name, runid=source.get_ncattr_text('Runid', default=''))
            for item in self._inventory.missing:
                message = (f'Source input NOT embedded ({item.status}): '
                           f'{item.key} -> {item.path}')
                source.logger.warning(message)
                print(f'*** WARNING: {message}')
            if self._inventory.truncated:
                message = ('Source input discovery exceeded safety limits; '
                           'input inventory is incomplete')
                source.logger.warning(message)
                print(f'*** WARNING: {message}')
        return self._inventory

    def parameter_files(self):
        """``(role, key, reference, path)`` of every input in the inventory."""
        return [(item.role, item.key, item.reference, item.path)
                for item in self.input_inventory().records]

    def file_records(self):
        """Portable source-file records for IDS-level provenance (see ``code_provenance``).

        Kept small (the DD caps the reference list at 10); full file contents and recursively
        referenced templates stay preserved in ``code.parameters``.
        """
        source = self.source
        if self._file_records is not None:
            return list(self._file_records)
        if getattr(source, 'input_format', 'netcdf') == 'pln':
            prefix = source.ncfile[:-6]
            candidates = [('TRANSP XF.PLN', source.ncfile), ('TRANSP YF.PLN', prefix + 'YF.PLN')]
            pln_input = getattr(source, '_pln_input', None)
            if pln_input is not None and pln_input.tf_path.is_file():
                candidates.append(('TRANSP TF.PLN', str(pln_input.tf_path)))
            schema_cdf = getattr(pln_input, 'cdf_path', None)
            if schema_cdf:
                candidates.append(('TRANSP CDF schema', str(schema_cdf)))
        else:
            candidates = [('TRANSP CDF', source.ncfile)]
        if source.nml_name:
            candidates.append(('TRANSP TR.DAT', source.nml_name))
        # Include resolved model/template files when present, up to the DD cap.
        try:
            for role, key, filename, resolved in self.parameter_files():
                if role == 'TRDAT':
                    continue
                if resolved is not None:
                    candidates.append((f'{role}:{key}', str(resolved)))
        except Exception as exc:          # best-effort provenance: never abort, but never silent either
            source.mapping.skipped(f'could not collect model/template files: {type(exc).__name__}: {exc}',
                                   target='ids_properties/provenance')
        records = source_file_records(candidates)
        self._file_records = list(records)
        return list(records)

    def fill_source_provenance(self, ids):
        """Attach source CDF/TR.DAT provenance to the whole IDS (layout per DD: see ``code_provenance``)."""
        source = self.source
        try:
            attach_source_provenance(ids, self.file_records())
        except Exception as exc:          # best-effort provenance: never abort, but never silent either
            source.logger.warning(f'Could not attach IDS provenance: {exc}')
            source.mapping.skipped(f'could not attach IDS provenance: {type(exc).__name__}: {exc}',
                                   target='ids_properties/provenance')
        return ids

    def fill_properties(self, ids):
        source = self.source
        runid, device, year, cdf_date, build_date, tglf, glf = source.get_nc_global_atts()
        ids.ids_properties.comment = properties_comment(runid, device, year, cdf_date, build_date, tglf, glf,
                                                        source.dd_version, source.signs)
        ids.ids_properties.creation_date = datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
        if hasattr(ids.ids_properties, 'name'):
            ids.ids_properties.name = f'TRANSP {runid}'
        ids = self.fill_source_provenance(ids)

        # Do not fill the legacy top-level fields
        #   ids_properties.data_dictionary
        #   ids_properties.access_layer
        # They are not present in current IMAS Data Dictionary versions (including
        # DD 4.1.1).  In modern IMAS they live below ids_properties.version_put.
        # More importantly, IMAS-Python/DBEntry.put() fills version_put from the
        # Data Dictionary and backend actually used for the write.  Leaving this
        # metadata to DBEntry therefore avoids stale or inconsistent version data
        # and remains compatible with older DD versions as well.
        return ids

    def build_code_parameters_xml(self):
        """XML for IMAS code.parameters (see ``code_provenance.build_code_parameters_xml``); cached."""
        source = self.source
        if self._code_parameters_xml is not None:
            return self._code_parameters_xml
        metadata = {}
        for key in TRANSP_METADATA_KEYS:
            try:
                metadata[key] = source.get_ncattr(key)
            except (AttributeError, KeyError, RuntimeError):
                metadata[key] = None
        self._code_parameters_xml = build_code_parameters_xml(
            dd_version=source.dd_version, signs=source.signs, device=source.device, shot=source.shot,
            imas_run=source.run, occurrence=source.occur,
            source_cdf=source.ncfile if getattr(source, 'input_format', 'netcdf') == 'netcdf' else '',
            source_format=getattr(source, 'input_format', 'netcdf'), source_input=source.ncfile,
            converter_version=self.converter_version, metadata=metadata, inventory=self.input_inventory())
        return self._code_parameters_xml

    def fill_code_properties(self, ids):
        """Fill common code provenance, if this IDS contains a code node."""
        source = self.source
        if not hasattr(ids, 'code'):
            return ids
        ids.code.name = 'transp2imas'
        if hasattr(ids.code, 'description'):
            ids.code.description = ('TRANSP native PLN to IMAS converter'
                                    if getattr(source, 'input_format', 'netcdf') == 'pln'
                                    else 'TRANSP CDF to IMAS converter')
        if hasattr(ids.code, 'version'):
            ids.code.version = self.converter_version
        if hasattr(ids.code, 'repository'):
            ids.code.repository = REPOSITORY
        if hasattr(ids.code, 'parameters'):
            ids.code.parameters = self.build_code_parameters_xml()
        return ids

    def fill_external_model_inputs(self, ids, family):
        """Attach model-specific XML to the appropriate IDS code libraries."""
        if not hasattr(ids, 'code'):
            return []
        return fill_model_libraries(ids.code, self.input_inventory(), family)
