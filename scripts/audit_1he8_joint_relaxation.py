#!/usr/bin/env python3
"""Re-audit retained 1HE8 mobility/outcomes without sampling or minimization."""
from __future__ import annotations

import argparse
import json
from collections import Counter
from dataclasses import asdict
from pathlib import Path

from scripts.mac_reconstruction_comparison import digest, save

from dvbfixer.ffutils.geometry import assert_all_l
from dvbfixer.model.diffusion.anchor_acceptance import _structure
from dvbfixer.model.diffusion.anchor_policy import (
    AnchorMovementLimits,
    AnchorPositionRestraint,
    atom_residue,
    check_anchor_movements,
)
from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest
from dvbfixer.model.diffusion.gap_portfolio import _coordinates
from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

PROJECT = Path(__file__).resolve().parents[1]
CASES = (
    {"name": "full-patchwork-direct-anchors", "root": "tmp/1HE8-patchwork-all-generated-anchors-v1/attempt-01",
     "declaration": "protocol.json", "input": "patchwork-input.pdb", "final": "patchwork-refined.pdb", "limits": "limits"},
    {"name": "richer-context-two-anchor-control",
     "root": "tmp/mac-reconstruction-2026-11-11/1he8-two-anchor-engineering-feasibility",
     "declaration": "declaration.json", "input": "input.pdb", "final": "refined.pdb", "limits": "engineering_limits"},
)


def read_restraints(declaration):
    return tuple(AnchorPositionRestraint(**{**r, "atom": AtomIdentity(**r["atom"]),
        "reference_angstrom": tuple(r["reference_angstrom"])}) for r in declaration["restraints"])


def inspect(case):
    root = PROJECT / case["root"]
    declaration = json.loads((root / case["declaration"]).read_text())
    request = DiffusionRequest.from_json((root / "request.json").read_text())
    source = (root / "source.pdb").read_text()
    input_text, final_text = (root / case["input"]).read_text(), (root / case["final"]).read_text()
    if digest(root / "source.pdb") != request.normalized_pdb.sha256 or declaration["source_sha256"] != request.normalized_pdb.sha256:
        raise ValueError("retained 1HE8 source binding differs")
    declared_input = declaration.get("input_sha256", declaration.get("patchwork_sha256"))
    if declared_input != digest(root / case["input"]):
        raise ValueError("retained input binding differs")
    restraints = read_restraints(declaration)
    generated = frozenset(request.generated_atoms)
    anchors = frozenset(r.atom for r in restraints)
    if generated & anchors or not anchors <= set(request.fixed_atoms):
        raise ValueError("generated and observed anchor ownership differs")
    coordinates, original_rows = _coordinates(input_text)
    final_coordinates, final_rows = _coordinates(final_text)
    if set(original_rows) != set(final_rows) or set(coordinates) != set(generated) | set(request.fixed_atoms):
        raise ValueError("retained atom inventory differs")
    if final_text != rewrite_generated_coordinates(input_text, {a: final_coordinates[a] for a in generated | anchors}):
        raise ValueError("nonmobile records or metadata changed")
    limits = AnchorMovementLimits(**declaration[case["limits"]])
    movements = check_anchor_movements(restraints, final_coordinates, limits)
    assess = assessment_adapter(request, source)
    initial, final = assess(input_text), assess(final_text, anchors)
    structure = _structure(final_text, "1HE8-retained-readback.pdb")
    assert_all_l(structure.pdb.topology, structure.pdb.positions)
    generated_residues = {atom_residue(a) for a in generated}
    anchor_residues = {atom_residue(a) for a in anchors}
    topology_atoms = {AtomIdentity(a.residue.chain.id, str(a.residue.id), a.residue.insertionCode.strip(), a.name): a
                      for a in structure.pdb.topology.atoms()}
    if not generated | anchors <= topology_atoms.keys():
        raise ValueError("OpenMM readback lost movable identities")
    issue_map = []
    for issue in final.issues:
        if not issue.actionable:
            continue
        residues = set(issue.residues)
        atoms = {a for a in topology_atoms if atom_residue(a) in residues}
        neighbors = set()
        for first, second in structure.pdb.topology.bonds():
            a = AtomIdentity(first.residue.chain.id, str(first.residue.id), first.residue.insertionCode.strip(), first.name)
            b = AtomIdentity(second.residue.chain.id, str(second.residue.id), second.residue.insertionCode.strip(), second.name)
            if a in atoms and b not in atoms:
                neighbors.add(b)
            if b in atoms and a not in atoms:
                neighbors.add(a)
        issue_map.append({"issue": asdict(issue), "residue_roles": [{"residue": asdict(r),
            "role": "generated" if r in generated_residues else "released-anchor" if r in anchor_residues else "fixed-observation"}
            for r in sorted(residues)], "external_bond_neighbors": [{"atom": asdict(a),
            "mobile": a in generated | anchors} for a in sorted(neighbors)]})
    historical = json.loads((root / "summary.json").read_text())
    diagnostics = historical["search_diagnostics"]
    return {"name": case["name"], "root": str(root), "generated_atoms_simultaneously_mobile": len(generated),
        "observed_anchor_atoms": len(anchors), "observed_anchor_residues": len(anchor_residues),
        "generated_regions": len(request.gaps), "limits": asdict(limits),
        "outside_movable_records_exact": True, "global_all_l": True,
        "initial_actionable_categories": dict(Counter(i.category for i in initial.issues if i.actionable)),
        "final_actionable_categories": dict(Counter(i.category for i in final.issues if i.actionable)),
        "final_admission_failures": list(final.admission_failures), "final_geometry_passed": final.geometry_passed,
        "anchor_movements": [asdict(m) for m in movements], "movement_violations": sum(not m.within_limits for m in movements),
        "accepted": final.geometry_passed and not final.admission_failures and all(m.within_limits for m in movements),
        "issue_map": issue_map, "historical_optimizer_diagnostics": diagnostics,
        "convergence_proven": False, "convergence_limit": "historical diagnostics contain phase caps/energies, not actual iteration counts, final gradient or termination reason",
        "temporary_hydrogenated_system_and_force_state_retained": False,
        "trajectory_no_passage_proven": False, "RMSD_to_crystal_used_for_acceptance": False,
        "new_sampling": False, "new_minimization": False,
        "input_sha256": digest(root / case["input"]), "final_sha256": digest(root / case["final"]),
        "request_sha256": digest(root / "request.json"), "source_sha256": digest(root / "source.pdb")}


def run(output):
    rows = [inspect(case) for case in CASES]
    output.mkdir(parents=True, exist_ok=False)
    save(output / "independent-readback.json", {"rows": rows, "auditor_sha256": digest(Path(__file__)),
        "no_new_sampling_or_optimization": True, "causal_explanation_complete": False})
    (output / "entry-script.py").write_bytes(Path(__file__).read_bytes())
    for row in rows:
        print(json.dumps({k: row[k] for k in ("name", "generated_atoms_simultaneously_mobile", "observed_anchor_atoms",
            "final_actionable_categories", "movement_violations", "accepted")}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.output.resolve())
