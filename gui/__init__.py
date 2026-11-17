#!/usr/bin/python
"""A local web workbench for TRANSP2IMAS: read curation and conversion reports, run conversions, plot the results.

Default mode is a report viewer that can also run ``curate.py`` (read-only on its inputs; a report goes into the job directory)
but starts no conversion.  With ``--enable-run`` the Convert view can
start ``convert.py`` from a validated argument list into a scratch IMAS database, and the Explore view reads
converted HDF5 entries in a worker subprocess to plot them.  A curation report is a screening record, not a scientific
certification, and a passing conversion report is not a scientific validation of the converted data.

Structure: ``reports`` and ``masks`` turn documents into view models (pure functions, no web framework); ``browse`` finds
reports and decides which paths may be read; ``convert_plan`` and ``curate_plan`` turn a form into a checked command and a target-entry guard
(derived from the converter's own option parser); ``jobs`` runs the converter and keeps its log; ``imas_worker`` and
``imas_client`` read converted entries for plotting; ``recipes`` lists the curated plots; ``server`` is a
standard-library HTTP server for the static page in ``gui/static``.  Run it with ``python -m gui``.
"""

__version__ = '0.2.0'
