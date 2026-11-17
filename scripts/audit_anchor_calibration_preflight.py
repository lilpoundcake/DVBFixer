#!/usr/bin/env python3
"""Independent retained-particle readback; no Context, optimization or sampling."""
from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from importlib.metadata import version
from pathlib import Path

import numpy as np
from scripts.audit_retained_full_context_controls import atoms
from scripts.mac_reconstruction_comparison import digest, save
from scripts.retained_evidence import verify_bindings

from dvbfixer.model.diffusion.contract import DiffusionRequest


def key(identity):
    return tuple(identity[name] for name in ("chain", "residue_number", "insertion_code", "atom_name"))


def independent_windows(request, width):
    """Use observed target indices, independently of the worker's target-axis walk."""
    selected = set()
    for placement in request.sequence_placements:
        by_index = dict(zip(placement.observed_target_indices, placement.observed_residues, strict=True))
        for gap in request.gaps:
            if gap.chain != placement.chain or gap.left_anchor is None or gap.right_anchor is None:
                continue
            for first, direction in ((gap.target_interval.start - 1, -1), (gap.target_interval.stop, 1)):
                for offset in range(width):
                    index = first + direction * offset
                    if index not in by_index:
                        break
                    selected.add(by_index[index])
    return sorted(selected)


def audit_snapshot(case, snapshot, window):
    from openmm import NonbondedForce, XmlSerializer
    from openmm.unit import dalton

    request = DiffusionRequest.from_json((case / "request.json").read_text())
    proposed = atoms(case / "proposal.pdb")
    generated = {key(asdict(a)) for a in request.generated_atoms}
    residues = {(r.chain, r.residue_number, r.insertion_code) for r in independent_windows(request, window)}
    anchors = {a for a in proposed if a[:3] in residues and a not in generated}
    axis = json.loads((snapshot / "temporary-axis.json").read_text())["atoms"]
    identities = [key(a["identity"]) for a in axis]
    if len(set(identities)) != len(identities):
        raise ValueError("duplicate temporary identities")
    system = XmlSerializer.deserialize((snapshot / "temporary-system.xml").read_text())
    with np.load(snapshot / "initial-state.npz", allow_pickle=False) as data:
        positions = data["positions_nm"].copy()
        forces = data["full_force_kj_mol_nm"].copy()
    if positions.shape != (len(axis), 3) or forces.shape != positions.shape:
        raise ValueError("saved particle/force axes differ")
    if not np.isfinite(positions).all() or not np.isfinite(forces).all():
        raise ValueError("nonfinite initial coordinates or full forces")
    if system.getNumParticles() != len(axis):
        raise ValueError("serialized system/particle axis mismatch")
    base_path = snapshot / "base-system.xml"
    if base_path.is_file():
        base = XmlSerializer.deserialize(base_path.read_text())
        if base.getNumParticles() != len(axis) or base.getNumConstraints() or system.getNumConstraints():
            raise ValueError("base system particle/constraint mismatch")
        if any(base.getParticleMass(i).value_in_unit(dalton) <= 0 for i in range(len(axis))):
            raise ValueError("base system already contains frozen particles")
        if base.getNumForces() > system.getNumForces():
            raise ValueError("base system has extra force definitions")
        for index in range(base.getNumForces()):
            first, second = base.getForce(index), system.getForce(index)
            # Geometry conditioning changes only the full nonbonded force group.
            first.setForceGroup(0)
            second.setForceGroup(0)
            if XmlSerializer.serialize(first) != XmlSerializer.serialize(second):
                raise ValueError("base chemistry differs from captured full force field")
    indices = {identity: i for i, identity in enumerate(identities)}
    if not proposed.keys() <= indices.keys():
        raise ValueError("saved system lost full-context input atoms")
    for identity, atom in proposed.items():
        if not np.allclose(positions[indices[identity]]*10, atom["xyz"], rtol=0, atol=1.e-4):
            raise ValueError("saved initialization differs from proposal")
    actual_fixed = {identity for i, identity in enumerate(identities)
                    if system.getParticleMass(i).value_in_unit(dalton) == 0}
    if actual_fixed != proposed.keys() - generated - anchors:
        raise ValueError("independent fixed/mobile particle masks differ")
    for i, row in enumerate(axis):
        if row["mass_dalton"] != system.getParticleMass(i).value_in_unit(dalton):
            raise ValueError("reported/serialized particle masses differ")
        if row["temporary_addition"] != (identities[i] not in proposed):
            raise ValueError("temporary addition ownership mismatch")
    nonbonded = [f for f in system.getForces() if isinstance(f, NonbondedForce)]
    if len(nonbonded) != 1 or nonbonded[0].getNumParticles() != len(axis):
        raise ValueError("missing full nonbonded context participation")
    if nonbonded[0].getNonbondedMethod() != NonbondedForce.NoCutoff:
        raise ValueError("unexpected periodic/cutoff model")
    heavy_additions = [row for row in axis if row["temporary_addition"] and row["element"] not in {"H", "D"}]
    return {"independent_particle_readback_passed": True, "retained_heavy_atoms": len(proposed),
            "generated_mobile_heavy_atoms": len(generated), "observed_mobile_heavy_atoms": len(anchors),
            "fixed_heavy_atoms": len(actual_fixed), "temporary_heavy_additions": heavy_additions,
            "all_initial_forces_finite": True, "no_cutoff_nonbonded_particles": len(axis),
            "base_forcefield_matches_captured_system": base_path.is_file(),
            "system_and_initial_state_sha256": {p.name: digest(p) for p in (
                snapshot / "temporary-system.xml", snapshot / "temporary-axis.json", snapshot / "initial-state.npz",
                *([base_path] if base_path.is_file() else []))}}


def run(root, output):
    import openmm.app

    output.mkdir(parents=True, exist_ok=False)
    plan_path = root / "plan.json"
    plan = json.loads(plan_path.read_text())
    verify_bindings(plan["source_code_and_input_sha256"], (plan_path,))
    if digest(Path(plan["python"])) != plan["python_sha256"]:
        raise ValueError("declared interpreter changed")
    inventory = json.loads((root / "inventory.json").read_text())
    cases = []
    for row in inventory["cases"]:
        entry = row["entry"]
        item = {"entry": entry, "input_status": row["preparation_status"],
                "reason": row.get("reason"), "system_preflight_passed": False,
                "grid_execution_ready": False}
        snapshot = root / f"execution/{entry}-seed-7"
        result_path = snapshot / "result.json"
        if result_path.is_file():
            result = json.loads(result_path.read_text())
            item["worker_result_sha256"] = digest(result_path)
            item["system_preflight_passed"] = result["system_preflight_passed"]
            item["reason"] = result.get("reason")
            item["execution_receipt"] = json.loads((snapshot / "execution.json").read_text())
            if result["study_optimizer_calls"] != 0:
                raise ValueError("preflight unexpectedly executed study optimization")
            if result["system_preflight_passed"]:
                request = DiffusionRequest.from_json((root / f"cases/{entry}/request.json").read_text())
                for width in (0, 1, 2):
                    expected = [asdict(r) for r in independent_windows(request, width)]
                    if expected != row["windows"][str(width)]:
                        raise ValueError("independent sequence-window audit differs")
                try:
                    item["readback"] = audit_snapshot(root / "cases" / entry, snapshot, 2)
                except ValueError as exc:
                    item.update(system_preflight_passed=False, readback_error=str(exc))
        cases.append(item)
    data = Path(openmm.app.__file__).parent / "data"
    save(output / "receipt.json", {"cases": cases, "plan_sha256": digest(plan_path),
        "audit_code_sha256": digest(Path(__file__)), "independent_parser_sha256": digest(Path(__file__).with_name("audit_retained_full_context_controls.py")),
        "versions": {name: version(name) for name in ("openmm", "pdbfixer", "numpy")},
        "openmm_xml_resource_bindings": {str(p): digest(p) for p in sorted(data.rglob("*.xml"))},
        "study_optimizer_calls": 0, "grid_execution_ready": False,
        "qualification": False, "remaining_gate": "shared-initial-state adapter and full grid readback not implemented"})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    run(args.root.resolve(), args.output.resolve())
