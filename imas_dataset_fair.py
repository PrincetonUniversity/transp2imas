#!/usr/bin/python
"""Populate dataset_fair metadata for a TRANSP-to-IMAS conversion when supported by the DD."""

from writer_base import IdsWriter


class imas_dataset_fair(IdsWriter):
    """Write dataset-level FAIR/provenance metadata without inventing PIDs.

    TRANSP CDF/TR.DAT inputs do not provide a persistent identifier, dataset
    license, rights holder, replacement chain, or validity interval for the
    resulting IMAS entry.  Those fields are therefore intentionally left empty.
    Whole-IDS source provenance is written through ids_properties/provenance by
    nc2imas.fill_properties().
    """

    def build(self):
        self.logger.info('\ndataset_fair')
        try:
            ids = self.dd.new_ids('dataset_fair')
        except ValueError as exc:          # imas.exception.IDSNameError: not in this DD
            self.ctx.mapping.unavailable_in_dd('dataset_fair', str(exc))
            self.logger.warning(f'dataset_fair is unavailable in target DD: {exc}')
            return []

        ids = self.fill.fill_properties(ids)
        # dataset_fair has no code node in DD 4.1.1, but call the common helper
        # for DD compatibility in case a supported schema provides one.
        ids = self.fill.fill_code_properties(ids)
        if hasattr(ids.ids_properties, 'homogeneous_time'):
            ids.ids_properties.homogeneous_time = 2

        # Do not use TRANSP_DOI here: it identifies TRANSP software/reference
        # material, not this converted dataset.  A dataset PID, rights holder,
        # license and validity interval must come from an authoritative catalog
        # or publication workflow and must not be fabricated by the converter.
        note = ('FAIR metadata generated from a TRANSP CDF conversion. Dataset PID, '
                'license, rights holder, replacement chain and validity interval '
                'are not present in the supplied TRANSP inputs and are intentionally blank.')
        if hasattr(ids.ids_properties, 'comment'):
            base = str(ids.ids_properties.comment or '').strip()
            ids.ids_properties.comment = (base + ' ' + note).strip()

        return [ids]


if __name__ == '__main__':
    print('Please use: python convert.py')
