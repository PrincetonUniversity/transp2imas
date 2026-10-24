#!/usr/bin/python
"""Identity and schema-variant setters for IDS nodes whose layout depends on the Data Dictionary version.

Each function takes ``dd_at_least(major, minor, patch)``, the version predicate of the target
dictionary, instead of the converter, so the rules can be tested against any version.  ``nc2imas``
exposes them under their historical method names.
"""

import numpy as np


def set_species_name(dd_at_least, node, value):
    """Set the species/state text node for the selected DD version.

    DD 3.39 uses ``label``.  The species/state nodes were renamed to
    ``name`` in late DD3 (3.42) and remain ``name`` in DD4.
    """
    text = str(value)
    if dd_at_least(3, 42, 0):
        node.name = text
    else:
        node.label = text


def set_nuclear_charge(dd_at_least, node, value):
    """Set z_n using the DD-native scalar type.

    z_n changed from FLT_0D to INT_0D in DD4.
    """
    z = float(value)
    node.z_n = int(round(z)) if dd_at_least(4, 0, 0) else z


def set_ec_beam_identity(dd_at_least, beam, identifier, description=None):
    """Populate EC beam identity without mixing DD3/DD4 semantics.

    DD3.39: ``identifier`` is the short identifier and ``name`` is the
    descriptive beam name.  DD4: ``name`` is the short identifier and
    ``description`` carries descriptive text.
    """
    identifier = str(identifier)
    description = str(description if description is not None else identifier)
    if dd_at_least(4, 0, 0):
        beam.name = identifier
        beam.description = description
    else:
        beam.identifier = identifier
        beam.name = description


def set_nbi_unit_identity(dd_at_least, unit, identifier, description=None):
    """Populate NBI unit identity for DD3.39 or DD4."""
    identifier = str(identifier)
    description = str(description if description is not None else identifier)
    if dd_at_least(4, 0, 0):
        unit.name = identifier
        unit.description = description
    else:
        unit.identifier = identifier
        unit.name = description


def set_ec_mode(dd_at_least, beam, o_mode_fraction, times):
    """Write EC launch mode/polarization using the target DD schema.

    DD3 uses the integer ``mode`` node and, where available, the legacy
    beam-level ``o_mode_fraction`` array.  In DD4.1+, ``mode`` is removed
    and the fraction is stored under ``beam/polarization``.
    """
    frac = float(np.clip(float(o_mode_fraction), 0.0, 1.0))
    times = np.asarray(times, dtype=np.float64)
    vals = np.full(len(times), frac, dtype=np.float64)
    if dd_at_least(4, 1, 0):
        beam.polarization.time = times
        beam.polarization.o_mode_fraction = vals
    elif dd_at_least(4, 0, 0):
        # DD4.0 was transitional.  Prefer polarization when available;
        # fall back to the pre-4.1 beam-level fraction if present.
        if hasattr(beam, 'polarization') and hasattr(beam.polarization, 'o_mode_fraction'):
            if hasattr(beam.polarization, 'time'):
                beam.polarization.time = times
            beam.polarization.o_mode_fraction = vals
        elif hasattr(beam, 'o_mode_fraction'):
            beam.o_mode_fraction = vals
    else:
        # DD3 mode is the *main* excited wave mode.  For a mixed launch,
        # infer the dominant mode from the O-mode power fraction.
        beam.mode = 1 if frac >= 0.5 else -1
        if hasattr(beam, 'o_mode_fraction'):
            beam.o_mode_fraction = vals


class IdentityHelper:
    """The three setters a builder needs, bound to one DD predicate (``dd_at_least``)."""

    def __init__(self, dd_at_least):
        self.dd_at_least = dd_at_least

    def set_species_name(self, node, value):
        set_species_name(self.dd_at_least, node, value)

    def set_nuclear_charge(self, node, value):
        set_nuclear_charge(self.dd_at_least, node, value)
