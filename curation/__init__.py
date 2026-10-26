#!/usr/bin/python
"""Read-only TRANSP curation. It does not import IMAS and does not change inputs.

Schema identities (the documents are described by the offline JSON Schemas in ``curation/schemas``):

* evidence sidecars ``inputs.v1`` and ``inputs.v2`` are still accepted and migrated to the canonical ``inputs.v3``;
* reports are written as ``report.v3`` (v2 had no execution ledger, policy, registry or fingerprints).
"""

__version__ = '0.3.0'
REPORT_SCHEMA_V2 = 'transp2imas.curation.report.v2'
REPORT_SCHEMA = 'transp2imas.curation.report.v3'
MANIFEST_SCHEMA_V1 = 'transp2imas.curation.inputs.v1'
MANIFEST_SCHEMA_V2 = 'transp2imas.curation.inputs.v2'
MANIFEST_SCHEMA = 'transp2imas.curation.inputs.v3'
ACCEPTED_MANIFEST_SCHEMAS = (MANIFEST_SCHEMA_V1, MANIFEST_SCHEMA_V2, MANIFEST_SCHEMA)
