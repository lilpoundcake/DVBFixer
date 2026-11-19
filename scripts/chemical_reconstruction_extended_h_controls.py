"""Frozen private planar-N/thiol controls; no implicit chemical-state fallback."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from scripts.chemical_reconstruction_ccd_controls import ROOT, _write
from scripts.chemical_reconstruction_hydrogen_controls import (
    FFBONDED,
    FIXED_POLICY,
    GENERATED_POLICY,
    _execute,
    backbone_types,
)
from scripts.chemical_reconstruction_joint_h_controls import check_pins as check_joint_pins
from scripts.chemical_reconstruction_naming_controls import raw_inventory
from scripts.chemical_reconstruction_peptide_controls import load_ccd
from scripts.chemical_reconstruction_source_controls import _bytes

from dvbfixer.chemical_reconstruction.amide_authority import FFBONDED_SHA256, load_amide_authority
from dvbfixer.chemical_reconstruction.ccd_authority import load_geometry_authority
from dvbfixer.chemical_reconstruction.contracts import (
    ChemistryKind,
    Component,
    Connection,
    LinkKind,
    Observation,
    Operation,
    Origin,
    Ownership,
    ReconstructionRequest,
    Unit,
    digest_record,
)
from dvbfixer.chemical_reconstruction.hydrogen_completion import (
    compile_complete_hydrogens,
    compile_extended_hydrogens,
    compile_joint_hydrogens,
)
from dvbfixer.chemical_reconstruction.inventory_naming import transport_inventory
from dvbfixer.chemical_reconstruction.measurement import _dot, _sub, _unit
from dvbfixer.chemical_reconstruction.observed_peptide import select_interior_target
from dvbfixer.chemical_reconstruction.peptide_state import PeptideStateError, derive_local_peptide
from dvbfixer.chemical_reconstruction.source_preparation import (
    audit_source_transport,
    prepare_source,
    prepare_source_probe,
)
from dvbfixer.domain.force_field_naming import ForceFieldTarget, NamingProfile
from dvbfixer.domain.nonprotein_reconstruction import ExternalLink
from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef

NAMES = ("ASN", "GLN", "ARG", "CYS", "TRP", "GLY", "PRO", "HIS")
PROTOCOL = ROOT / "docs/plans/chemical-reconstruction-extended-h-v1.md"
PROTOCOL_SHA256 = "7864f7e84b0667a89b97c8315a94a95800e0b88cb829b2c407c5ab0009ddca6f"
FIXTURES = ROOT / "tests/fixtures/chemical_reconstruction/ccd_h_groups_v1"
MANIFEST_SHA256 = "22c74384e25dbe735103b5d442d3c3c2cf6084492a68cf9b3ce887babc85b0a7"


def observed_path(name: str) -> str:
    return "tests/fixtures/regressions/1FR2.pdb" if name in {"CYS", "TRP"} else "tests/fixtures/protein_context/1UBQ.pdb"


def check_pins():
    check_joint_pins()
    for path, digest in ((PROTOCOL, PROTOCOL_SHA256), (FIXTURES / "manifest.json", MANIFEST_SHA256)):
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(f"frozen extended-H pin mismatch: {path}")
    for name, pin in json.loads((FIXTURES / "manifest.json").read_text())["components"].items():
        if hashlib.sha256((FIXTURES / f"{name}.cif").read_bytes()).hexdigest() != pin["sha256"]:
            raise ValueError(f"frozen extended-H CCD changed: {name}")


def chemical_authority(name, identity):
    pins = json.loads((FIXTURES / "manifest.json").read_text())["components"]
    parent = load_geometry_authority(FIXTURES / f"{name}.cif", name, expected_sha256=pins[name]["sha256"])
    try:
        raw = raw_inventory(name)
    except KeyError as exc:
        raise PeptideStateError(f"explicit-same-name-RTP-state-unavailable:{name}") from exc
    inventory = transport_inventory(raw, target=ForceFieldTarget.AMBER19SB, profile=NamingProfile.STANDARD)
    authority = derive_local_peptide(parent, load_ccd("ALA"), identity, inventory.canonical)
    return authority, inventory


def synthetic_base(name, evidence):
    identity = ComponentInstanceRef(1, "a", "82", "B", occurrence=2)
    authority, inventory = chemical_authority(name, identity)
    evidence.update(authority=authority, inventory_transport=inventory)
    points = {a.name: a.ideal_position for a in authority.state.graph.atoms if a.element not in {"H", "D"}}
    n = points["N"]
    u, c = _unit(_sub(points["CA"], n)), _unit(_sub(points["C"], n))
    projection = _dot(c, u)
    v = _unit(tuple(c[i]-projection*u[i] for i in range(3)))
    angle = math.radians(121.96)
    incoming = ExactAtomRef(ComponentInstanceRef(1, "A", "81"), "C")
    point = tuple(n[i]+1.335*(math.cos(angle)*u[i]+math.sin(angle)*v[i]) for i in range(3))
    observations = tuple(Observation(ExactAtomRef(identity, a.name), a.element, points[a.name], Unit.ANGSTROM)
                         for a in authority.state.graph.atoms if a.element not in {"H", "D"}) + (
                             Observation(incoming, "C", point, Unit.ANGSTROM),)
    pin = digest_record(("synthetic-extended-H-heavy-only-v1", name, observations))
    base = ReconstructionRequest(pin, Operation.HEAVY,
                                 (Component(identity, ChemistryKind.PROTEIN, (authority.state,)),
                                  Component(incoming.component, ChemistryKind.UNKNOWN)), observations,
                                 frozenset({identity}), Ownership(frozenset()),
                                 (Connection(ExternalLink(incoming, ExactAtomRef(identity, "N"), 1.), LinkKind.COVALENT),),
                                 authority.constraints(minimum_volume=.05, planar_tolerance_angstrom=.05),
                                 origin=Origin("synthetic-extended-H-heavy-only-v1", pin))
    return base, authority, incoming


def observed_base(name, directory, evidence):
    path = observed_path(name)
    original = (ROOT / path).read_bytes()
    _bytes(directory / "original.pdb", original)
    evidence.update(source=path, original_sha256=hashlib.sha256(original).hexdigest())
    prepared = prepare_source(original)
    _bytes(directory / "prepared.pdb", prepared.prepared_text.encode("ascii"))
    evidence["preparation"] = {"digest": prepared.digest, "original_sha256": prepared.original_sha256,
                               "prepared_sha256": prepared.prepared_sha256, "removed_metadata": prepared.removed_metadata}
    target = select_interior_target(original, name, chain="A", number=None)
    authority, inventory = chemical_authority(name, target)
    evidence.update(authority=authority, inventory_transport=inventory)
    probe = prepare_source_probe(prepared, authority, withheld=frozenset())
    evidence["source_transport_audit"] = audit_source_transport(prepared, probe)
    if not evidence["source_transport_audit"].passed:
        raise ValueError("extended-H-source-transport-failed")
    return probe.request, authority, probe.real_peptide_partners[0]


def compile_view(base, authority, incoming, evidence):
    parameters = load_amide_authority(FFBONDED, expected_sha256=FFBONDED_SHA256)
    evidence.update(base_request=base, parameters=parameters)
    kwargs = dict(backbone_types=backbone_types(authority.inventory.state_name), fixed_policy=FIXED_POLICY,
                  generated_policy=GENERATED_POLICY)
    try:
        old = compile_complete_hydrogens(base, authority, incoming, parameters, **kwargs)
    except ValueError as exc:
        evidence["old_compiler"] = {"status": "refused", "reason": str(exc)}
    else:
        evidence["old_compiler"] = {"status": "compiled", "view": old}
    view = compile_joint_hydrogens(compile_extended_hydrogens(base, authority, incoming, parameters, **kwargs))
    assert view.request.observations == base.observations
    return view


def run_controls(output: Path):
    check_pins()
    output.mkdir(parents=True, exist_ok=False)
    records = []
    for arm in ("synthetic", "observed"):
        for name in NAMES:
            directory = output / f"{arm}-{name}"
            directory.mkdir()
            item = {"id": directory.name, "arm": arm, "state": name, "proposal_entry_called": False}
            evidence = {"control": item}
            try:
                base, authority, incoming = synthetic_base(name, evidence) if arm == "synthetic" else observed_base(name, directory, evidence)
                view = compile_view(base, authority, incoming, evidence)
                item.update(target=authority.identity.sequence_number, source_observations=len(base.observations),
                            fractional_context=len(base.context_policy.fractional_atoms) if hasattr(base, "context_policy") else 0,
                            observations_unchanged=view.request.observations == base.observations)
                _execute(view, authority, item, evidence)
            except ValueError as exc:
                item.update(status="refused", reason=str(exc))
            _write(directory / "evidence.json", evidence)
            records.append(item)
    check_pins()
    summary = {"schema": "dvbfixer.extended-H-controls.v1", "protocol_sha256": PROTOCOL_SHA256,
               "manifest_sha256": MANIFEST_SHA256, "declared": 16, "executed_entries": len(records),
               "controls": records, "all20_or_native_H_qualified": False, "production_enabled": False}
    _write(output / "summary.json", summary)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    print(json.dumps(run_controls(parser.parse_args().output), indent=2))
