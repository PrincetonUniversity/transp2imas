#!/usr/bin/python
"""Write physically defensible TRANSP active PF coil/supply data to IMAS.

The TRANSP CDF usually contains electrical waveforms but not an authoritative
coil winding model. CDF-only mode writes KK_ coil currents, VV_ coil voltages,
and VS_ supply voltages, while warning that coil geometry and polarity remain
unverified. An optional machine JSON manifest supplies verified winding signs,
coil elements, turns, resistances and DD-specific circuit connectivity.

KT_ is target/reference current: it must NOT overwrite supply.voltage and must
NOT be mistaken for actual supply.current. CC_ equilibrium-current channels
are opt-in via an authoritative manifest rather than silently substituting
for missing KK_ waveforms. Unknown/not saved channels are left empty, not
filled with zero waveforms or dummy -999 times.
"""

from __future__ import annotations

import hashlib
import xml.etree.ElementTree as ET

import numpy as np

from nc2imas import nc2imas
from pf_active_mapping import (
    all_time_coordinates, build_active_plan, populate_active_ids,
)
from writer_base import IdsWriter


class imas_pf_active(IdsWriter):
    """Create pf_active if actual CDF waveforms or machine geometry exist."""

    def __init__(self, ncimas: nc2imas, geometry_path: str | None = None):
        super().__init__(ncimas)
        self.geometry_path = geometry_path

    def build(self):
        """Write with correct per-signal time bases and optional sourced geometry."""
        self.logger.info("\npf_active")
        plan = build_active_plan(
            self.cdf.nc, geometry_path=self.geometry_path, machine=self.cdf.device,
        )
        if plan.ignored_targets:
            self.logger.warning(
                "pf_active: %d KT_* target-current channel(s) are NOT actual "
                "supply currents; pulse_schedule mapping is pending: %s",
                len(plan.ignored_targets), ", ".join(plan.ignored_targets),
            )
        if plan.ignored_controller:
            self.logger.info(
                "pf_active: ignoring %d KD_/KS_ control/uncertainty traces "
                "(not measured coil or supply current)", len(plan.ignored_controller),
            )
        if plan.ignored_equilibrium:
            self.logger.info(
                "pf_active: CC_* equilibrium-current candidates require "
                "explicit machine mapping to be used as coil currents: %s",
                ", ".join(plan.ignored_equilibrium),
            )
        if not plan.coils and not plan.supplies:
            self.logger.info("pf_active: no KK/VV/VS channels or sourced machine geometry; IDS not written")
            return []
        if plan.manifest_path is None:
            self.logger.warning(
                "pf_active: CDF-only PF waveforms; coil cross-sections, winding "
                "polarity, supply association and circuit connectivity are "
                "UNVERIFIED. Not suitable for magnetic-field reconstruction."
            )
        pfc = self.dd.new_ids("pf_active")
        self.fill.fill_properties(pfc)
        self.fill.fill_code_properties(pfc)
        dd_major = int(self.dd.dd_major_version)
        traces = all_time_coordinates(plan, dd_major)
        if traces:
            same = all(np.array_equal(traces[0], t) for t in traces[1:])
            pfc.ids_properties.homogeneous_time = 1 if same else 0
            if same:
                pfc.time = np.asarray(traces[0], dtype=np.float64)
            else:
                # Dynamic nodes each contain their own exact time arrays.
                # Generic pf_active.time serves only as a reference.
                generic = self.cdf.get_ncvar("TIME")
                pfc.time = (np.asarray(generic, dtype=np.float64)
                            if len(generic) else np.asarray(traces[0], dtype=np.float64))
        else:
            # A sourced device-only IDS is useful even when CDF currents are
            # absent; do not create any fake time points.
            pfc.ids_properties.homogeneous_time = 0

        count = populate_active_ids(pfc, plan, dd_major)
        if plan.manifest_text is None:
            pfc.ids_properties.comment += (
                "; pf_active: CDF waveforms only; no verified coil winding "
                "orientation, circuit connectivity or machine geometry"
            )
        else:
            digest = hashlib.sha256(plan.manifest_text.encode("utf-8")).hexdigest()
            pfc.ids_properties.comment += (
                f"; pf_active machine geometry={plan.manifest_path.name}; "
                f"sha256={digest}; source={plan.manifest_source}"
            )
            if hasattr(pfc, "code") and hasattr(pfc.code, "parameters"):
                try:
                    root = ET.fromstring(pfc.code.parameters)
                    source = ET.SubElement(root, "pf_active_geometry", {
                        "filename": plan.manifest_path.name,
                        "sha256": digest,
                        "source_reference": plan.manifest_source,
                    })
                    source.text = plan.manifest_text
                    pfc.code.parameters = ET.tostring(root, encoding="unicode")
                except (ET.ParseError, TypeError, ValueError) as exc:
                    self.logger.warning("pf_active: cannot embed machine provenance: %s", exc)

        self.logger.info(
            "pf_active: built %d coils, %d supplies, %d circuits, "
            "%d coil currents, %d coil voltages, %d supply currents, "
            "%d supply voltages",
            len(plan.coils), len(plan.supplies), count["circuits"],
            count["coil_current"], count["coil_voltage"],
            count["supply_current"], count["supply_voltage"],
        )
        return [pfc]


if __name__ == "__main__":
    print("Please use: python convert.py")
