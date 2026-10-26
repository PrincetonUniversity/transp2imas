#!/usr/bin/python
"""The built-in rule set, registered explicitly (no discovery).

Finding ids are stable across the extraction from the former single-class engine.  A rule is *mandatory* when it
screens source integrity (clocks, grids, units, identity); policy cannot disable those.
"""

from curation.registry import Rule, RuleRegistry
from curation.rules import (composition, coordinates, events, evidence, gradients, heating, model_consistency, numerics,
                            profiles)


def default_registry():
    registry = RuleRegistry()
    add = registry.register
    add(Rule('coordinates.grid', '1', 'coordinates', coordinates.check_coordinates,
             'Native time coordinates and the radial grid', mandatory=True))
    add(Rule('identity.run', '1', 'evidence', evidence.check_identity, 'Declared run identity versus the source',
             mandatory=True))
    add(Rule('profiles.screen', '1', 'profiles', profiles.check_profiles,
             'NE/TE/TI finite, negative, zero, unit and gradient screening', requires=('coordinates.grid',),
             parameters=('weak_gradient_relative', 'gradient_sigma_multiplier')))
    add(Rule('uncertainty.gradient_covariance', '1', 'gradients', gradients.check_gradient_covariance,
             'Gradient uncertainty propagated from a supplied fixed-grid covariance', requires=('coordinates.grid',),
             parameters=('gradient_sigma_multiplier',)))
    add(Rule('evidence.provenance', '1', 'evidence', evidence.check_evidence,
             'Origin, uncertainty declaration, time coverage and lineage of evidence', requires=('identity.run',)))
    add(Rule('evidence.rotation', '1', 'evidence', evidence.check_rotation, 'Rotation output availability'))
    add(Rule('composition.quasineutrality', '1', 'composition', composition.check_composition,
             'Charge closure of a declared complete inventory', parameters=('quasineutrality_relative',)))
    add(Rule('model.declared_inputs', '1', 'model_consistency', model_consistency.check_declared_inputs,
             'Declared TR.DAT timeline, model-control review flags and beam current fractions'))
    add(Rule('heating.reference', '1', 'heating', heating.check_heating_reference,
             'Declared heating configuration versus analyst references', requires=('model.declared_inputs',)))
    add(Rule('inputs.resolved', '1', 'model_consistency', model_consistency.check_resolved_inputs,
             'Resolved inputs versus declarations', requires=('model.declared_inputs',)))
    add(Rule('numerics.runtime', '1', 'numerics', numerics.check_numerics,
             'Runtime solver channel inventory and thresholds', parameters=('newton_iterations_warning',)))
    add(Rule('events.detect', '1', 'events', events.check_event_detectors,
             'Optional abrupt-drop candidate detection on configured diagnostics'))
    add(Rule('events.inventory', '1', 'events', events.check_events,
             'Declared events and the policy exclusions they produced', after_all=True))
    registry.validate()
    return registry
