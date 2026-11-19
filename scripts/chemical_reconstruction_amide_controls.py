"""Frozen synthetic external-anchor H controls, not observed PDB reconstruction."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from scripts.chemical_reconstruction_ccd_controls import ROOT, _write
from scripts.chemical_reconstruction_naming_controls import check_pins as check_naming_pins
from scripts.chemical_reconstruction_naming_controls import raw_inventory
from scripts.chemical_reconstruction_peptide_controls import RTP, load_ccd

from dvbfixer.chemical_reconstruction.amide_authority import FFBONDED_SHA256, load_amide_authority
from dvbfixer.chemical_reconstruction.construction import propose
from dvbfixer.chemical_reconstruction.contracts import (
    ChemistryKind,
    Component,
    Connection,
    Constraint,
    ConstraintKind,
    LinkKind,
    Observation,
    Operation,
    Origin,
    Ownership,
    ReconstructionRequest,
    Unit,
    digest_record,
)
from dvbfixer.chemical_reconstruction.inventory_naming import transport_inventory
from dvbfixer.chemical_reconstruction.measurement import _dot, _sub, _unit
from dvbfixer.chemical_reconstruction.metric_geometry import MetricGeometryPlan
from dvbfixer.chemical_reconstruction.peptide_state import derive_local_peptide
from dvbfixer.chemical_reconstruction.proposals import (
    AuditPolicy,
    Placement,
    PlacementKind,
    Strategy,
)
from dvbfixer.domain.force_field_naming import ForceFieldTarget, NamingProfile
from dvbfixer.domain.nonprotein_reconstruction import ExternalLink
from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.rtp_parser import parse_rtp

PROTOCOL = ROOT / "docs/plans/chemical-reconstruction-amide-v1.md"
PROTOCOL_SHA256 = "f5faf25e32b9bb72370192b3060bf161f82981330c55cd126f8c8989b676b18e"
POLICY = AuditPolicy(.02, 5., .7, 0.)
FFBONDED = ROOT / "FF/amber19sb.ff/ffbonded.itp"


def synthetic_case(name: str, *, incoming_angle: float = 121.96):
    parameters = load_amide_authority(FFBONDED, expected_sha256=FFBONDED_SHA256)
    transport = transport_inventory(raw_inventory(name), target=ForceFieldTarget.AMBER19SB, profile=NamingProfile.STANDARD)
    identity = ComponentInstanceRef(1, "a", "82", "B", occurrence=2)
    authority = derive_local_peptide(load_ccd(name), load_ccd("ALA"), identity, transport.canonical)
    nitrogen, hydrogen, alpha, carbon = (ExactAtomRef(identity, n) for n in ("N", "H", "CA", "C"))
    incoming = ExactAtomRef(ComponentInstanceRef(1, "A", "81"), "C")
    points = {ExactAtomRef(identity, a.name): a.ideal_position for a in authority.state.graph.atoms if a.name != "H"}
    n = points[nitrogen]
    u = _unit(_sub(points[alpha], n))
    c = _unit(_sub(points[carbon], n))
    dot = _dot(c, u)
    v = _unit((c[0]-dot*u[0], c[1]-dot*u[1], c[2]-dot*u[2]))
    angle = math.radians(incoming_angle)
    # Synthetic anchor constructed without the absent H coordinate.
    incoming_point = tuple(n[i]+1.335*(math.cos(angle)*u[i]+math.sin(angle)*v[i]) for i in range(3))
    observations = tuple(Observation(ExactAtomRef(identity, a.name), a.element, points[ExactAtomRef(identity, a.name)], Unit.ANGSTROM)
                         for a in authority.state.graph.atoms if a.name != "H") + (
                             Observation(incoming, "C", incoming_point, Unit.ANGSTROM),)
    _, rtp = parse_rtp(RTP)
    types = {row[0]: row[1] for row in rtp[name].atoms}
    metrics = parameters.metrics(nitrogen, hydrogen, alpha, incoming, atom_types=(types["N"], types["H"], types["CA"]))
    constraints = (*authority.constraints(minimum_volume=.05, planar_tolerance_angstrom=.05),
                   Constraint(ConstraintKind.PLANARITY, (nitrogen, incoming, alpha, hydrogen), Unit.ANGSTROM, upper=.005))
    components = (Component(identity, ChemistryKind.PROTEIN, (authority.state,)), Component(incoming.component, ChemistryKind.UNKNOWN))
    source_digest = digest_record(("SYNTHETIC amide anchors, not original PDB bytes", name, incoming_angle, observations))
    request = ReconstructionRequest(source_digest, Operation.HYDROGENS, components, observations,
                                    frozenset({identity}), Ownership(frozenset({hydrogen})),
                                    (Connection(ExternalLink(incoming, nitrogen, 1.), LinkKind.COVALENT),), constraints,
                                    origin=Origin("synthetic-amide-metric-v1", parameters.digest,
                                                  (("inventory_transport", transport.digest), ("peptide_authority", authority.digest))))
    plan = MetricGeometryPlan(request.digest, Strategy.INTERNAL, POLICY,
                              (Placement(hydrogen, (nitrogen, incoming, alpha), PlacementKind.PLANAR),), metrics)
    return request, plan, {"parameters": parameters, "transport": transport, "authority": authority,
                           "synthetic_source": True, "incoming_angle_degrees": incoming_angle}


def run_controls(output: Path) -> dict:
    check_naming_pins()
    if hashlib.sha256(PROTOCOL.read_bytes()).hexdigest() != PROTOCOL_SHA256:
        raise ValueError("frozen amide protocol changed")
    load_amide_authority(FFBONDED, expected_sha256=FFBONDED_SHA256)
    output.mkdir(parents=True, exist_ok=False)
    records = []
    for name, angle in (*((name, 121.96) for name in ("ALA", "THR", "ILE", "LYS", "PHE")), ("ALA", 100.)):
        request, plan, evidence = synthetic_case(name, incoming_angle=angle)
        result = propose(request, plan)
        item = {"id": f"{name}-{angle:g}", "status": result.status, "reason": result.reason,
                "raw_passed": result.raw_audit.passed if result.raw_audit else False,
                "rounded_passed": result.rounded_audit.passed if result.rounded_audit else False}
        evidence.update(control=item, request=request, plan=plan, result=result)
        _write(output / f"{item['id']}.json", evidence)
        records.append(item)
    check_naming_pins()
    load_amide_authority(FFBONDED, expected_sha256=FFBONDED_SHA256)
    summary = {"schema": "dvbfixer.synthetic-amide-controls.v1", "protocol_sha256": PROTOCOL_SHA256,
               "declared": 6, "executed_entries": len(records), "controls": records,
               "observed_pdb_reconstruction": False, "h_or_complex_qualification": False, "production_enabled": False}
    _write(output / "summary.json", summary)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    print(json.dumps(run_controls(parser.parse_args().output), indent=2))
