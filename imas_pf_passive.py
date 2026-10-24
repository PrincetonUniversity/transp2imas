#!/usr/bin/python
"""Translate TRANSP PS_n currents and supplied machine geometry to pf_passive.

The TRANSP CDF contains current histories, not an authoritative passive-vessel
mechanical model. Therefore geometry, resistance, resistivity, area and the
signed number of turns are populated *only* from a validated, device-specific
JSON manifest supplied by the user (``--pf-passive-geometry``). Without it the
legacy current-only representation remains possible, with an explicit warning.

The writer supports DD3 and DD4 geometries; notably, DD4 requires the outline
polygon to repeat its initial vertex, whereas DD3 requires an open vertex list.
All loop time traces use their own ``loop.time`` regardless of CDF time base.
"""

import hashlib
import xml.etree.ElementTree as ET

import numpy as np

from nc2imas import nc2imas
from pf_passive_mapping import build_passive_plan, populate_passive_ids
from writer_base import IdsWriter


class imas_pf_passive(IdsWriter):
    """Create a passive-conductor IDS from available, validated source inputs."""

    def __init__(self, ncimas: nc2imas, geometry_path: str | None = None):
        super().__init__(ncimas)
        self.geometry_path = geometry_path

    def build(self):
        """Write pf_passive only when passive currents or sourced geometry exist."""
        self.logger.info("\npf_passive")
        plan = build_passive_plan(
            self.cdf.nc, geometry_path=self.geometry_path, machine=self.cdf.device
        )
        if not plan.loops:
            self.logger.info("pf_passive: no PS_n channels and no geometry; IDS not written")
            return []

        if plan.geometry_path is None:
            self.logger.warning(
                "pf_passive: %d current-only loops; no authoritative conductor "
                "geometry/resistance/turns supplied. IDS is INCOMPLETE and "
                "must not be used for electromagnetic modeling.",
                len(plan.loops),
            )
        else:
            self.logger.info(
                "pf_passive: %d geometry-defined loops loaded from %s",
                len(plan.loops), plan.geometry_path,
            )

        pfc = self.dd.new_ids("pf_passive")
        self.fill.fill_properties(pfc)
        self.fill.fill_code_properties(pfc)

        # IMAS pf_passive defines each current on loop[i].time, NOT simply on
        # the top-level pf_passive.time. Homogeneous_time can only be 1 if all
        # current-bearing loops have equal time arrays.
        current_times = [loop.channel.time for loop in plan.loops if loop.channel is not None]
        if current_times:
            shared_time = current_times[0]
            homogeneous = all(np.array_equal(shared_time, times) for times in current_times[1:])
            pfc.ids_properties.homogeneous_time = 1 if homogeneous else 0
            if homogeneous:
                pfc.time = np.array(shared_time, dtype=np.float64, copy=True)
            else:
                # Top-level time is informational in the heterogeneous case;
                # each dynamic quantity is explicitly coordinated by loop.time.
                scalar_time = self.cdf.get_ncvar("TIME")
                if len(scalar_time):
                    pfc.time = np.asarray(scalar_time, dtype=np.float64)
        else:
            # Geometry-only export is meaningful: it can provide the device
            # description even when a TRANSP run did not save passive currents.
            pfc.ids_properties.homogeneous_time = 1
            scalar_time = self.cdf.get_ncvar("TIME")
            if len(scalar_time):
                pfc.time = np.asarray(scalar_time, dtype=np.float64)

        count = populate_passive_ids(pfc, plan, int(self.dd.dd_major_version))

        if plan.geometry_json is not None:
            digest = hashlib.sha256(plan.geometry_json.encode("utf-8")).hexdigest()
            file_name = plan.geometry_path.name
            pfc.ids_properties.comment += (
                f"; pf_passive device geometry={file_name}; sha256={digest}; "
                f"source={plan.geometry_source}"
            )
            if hasattr(pfc, "code") and hasattr(pfc.code, "parameters"):
                try:
                    root = ET.fromstring(pfc.code.parameters)
                    source = ET.SubElement(root, "pf_passive_geometry", {
                        "filename": file_name,
                        "sha256": digest,
                        "source_reference": plan.geometry_source,
                    })
                    source.text = plan.geometry_json
                    pfc.code.parameters = ET.tostring(root, encoding="unicode")
                except (ET.ParseError, TypeError, ValueError) as exc:
                    self.logger.warning("pf_passive: cannot embed geometry provenance: %s", exc)
        else:
            pfc.ids_properties.comment += "; pf_passive: CDF currents only; no conductor geometry"

        self.logger.info(
            "pf_passive: built %d loops, %d current traces, %d defined geometries",
            len(plan.loops), count, sum(bool(loop.elements) for loop in plan.loops),
        )
        return [pfc]


if __name__ == "__main__":
    print("Please use: python convert.py")
