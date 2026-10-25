#!/usr/bin/python
"""Read-only scientific screening of TRANSP inputs and output; no IMAS import."""

__version__ = '0.2.0'
# Report schema: v2 adds per-rule excluded/unassessed accounting, the event inventory, the run identity, the
# resolved-input and evidence-dependency findings, and the manifest schema it was evaluated with.
REPORT_SCHEMA = 'transp2imas.curation.report.v2'
# Manifest schema: v2 is a superset of v1 (explicit evidence ids, parents and usage, applicability, overrides,
# resolved inputs, events, run identity).  A v1 document stays valid and is evaluated as a v2 document without them.
MANIFEST_SCHEMA = 'transp2imas.curation.inputs.v2'
MANIFEST_SCHEMA_V1 = 'transp2imas.curation.inputs.v1'
ACCEPTED_MANIFEST_SCHEMAS = (MANIFEST_SCHEMA_V1, MANIFEST_SCHEMA)
