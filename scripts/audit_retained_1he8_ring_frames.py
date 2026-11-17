#!/usr/bin/env python3
"""Read the historical1HE8 TYR972 ring witness in saved frames, never repair it."""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from dataclasses import asdict
from pathlib import Path

from openmm import XmlSerializer
from scripts.continue_retained_1he8_segment_frames import particle_graph, saved_coordinates
from scripts.mac_reconstruction_comparison import digest, save
from scripts.run_mac_modeller_comparison import execute

from dvbfixer.model.diffusion.anchor_acceptance import _structure
from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest
from dvbfixer.model.diffusion.covalent_ring_probe import probe_covalent_ring
from dvbfixer.model.diffusion.gap_portfolio import _coordinates
from dvbfixer.model.diffusion.scope import _parse_coordinate_records

PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "tmp/mac-real-gap-continuation-2026-11-12"
HISTORICAL = PROJECT / "tmp/mac-reconstruction-2026-11-11"
CONTROL_NAMES = ("1he8-selected-raw", "1he8-context-seed7-refined",
                 "1he8-context-seed42-refined", "1he8-context-seed101-refined")


def focus_probe(points, bonds, ring, witness, ownership, *, partial_source=False):
    """Missing/undefined surfaces are not counted as no-piercing passes."""
    if any(len(edge) != 2 for edge in bonds):
        raise ValueError("Invalid topology edge")
    ring_edges = {frozenset((a, b)) for a, b in zip(ring, ring[1:] + ring[:1])}
    if not ring_edges <= bonds or frozenset(witness) not in bonds:
        raise ValueError("Declared ring/witness edges absent from explicit topology")
    required = {atom for edge in bonds for atom in edge}
    missing = required - points.keys()
    if missing and not partial_source:
        raise ValueError("Complete saved frame lacks topology coordinates")
    missing_ring = set(ring) - points.keys()
    available = frozenset(edge for edge in bonds if edge <= points.keys())
    base = {"topology_edges": len(bonds), "available_coordinate_edges": len(available),
            "missing_topology_coordinate_atoms": len(missing),
            "missing_ring_atoms": [asdict(atom) for atom in sorted(missing_ring)],
            "witness_atoms_present": all(atom in points for atom in witness),
            "partial_observed_source": partial_source,
            "no_continuous_path_or_causality_claim": True, "publication_gate": False}
    if missing_ring:
        return {**base, "status": "missing-observed-ring", "probe": None, "witness_piercings": 0}
    result = probe_covalent_ring(ring, points, available, ownership, planarity_max_angstrom=.10)
    status = "undefined-surface" if result.undefined_reason else "defined-surface"
    witness_hits = sum(frozenset(hit.bond) == frozenset(witness) for hit in result.piercings)
    return {**base, "status": status, "probe": asdict(result), "witness_piercings": witness_hits}


def summarize_frames(frames):
    counts = Counter(frame["screen"]["status"] for frame in frames)
    return {"frames": len(frames), "defined_surfaces": counts["defined-surface"],
            "undefined_surfaces": counts["undefined-surface"], "missing_observed_rings": counts["missing-observed-ring"],
            "frames_with_piercings": sum(bool(frame["screen"]["probe"] and frame["screen"]["probe"]["piercings"]) for frame in frames),
            "frames_with_witness_piercings": sum(frame["screen"]["witness_piercings"] > 0 for frame in frames),
            "undefined_or_missing_is_not_a_clean_pass": True}


def declare(output):
    source_root = ROOT / "1HE8-remaining-segment-frames-attempt01"
    receipt_path = source_root / "readback-receipt.json"
    receipt = json.loads(receipt_path.read_text())
    for name, expected in receipt["bindings"].items():
        if digest(Path(name)) != expected:
            raise ValueError(f"Previously audited frame evidence changed: {name}")
    if receipt["arms_checked"] != 12 or receipt["all_frames_recomputed"] != 467:
        raise ValueError("Twelve-arm saved-frame denominator changed")
    plans = [source_root / "plan.json", ROOT / "1HE8-accepted-endpoint-segment-frames-attempt01/plan.json"]
    runs = [row for path in plans for row in json.loads(path.read_text())["runs"]]
    if len(runs) != 12 or len({row["work"] for row in runs}) != 12:
        raise ValueError("Repeated/missing historical arm")
    witness_path = HISTORICAL / "1he8-972-1044-covalent-ring-investigation.json"
    history = json.loads(witness_path.read_text())["cases"]
    ring, witness = history[0]["ring"], history[0]["crossing_bond"]
    if any(row["ring"] != ring or row["crossing_bond"] != witness for row in history):
        raise ValueError("Historical ring witness changed between controls")
    controls_path = HISTORICAL / "small-ring-diagnostic-controls/declaration.json"
    controls = {row["name"]: row for row in json.loads(controls_path.read_text())["cases"]}
    files = [receipt_path, *plans, witness_path, controls_path, Path(__file__),
             PROJECT / "scripts/continue_retained_1he8_segment_frames.py",
             PROJECT / "scripts/run_mac_modeller_comparison.py",
             PROJECT / "src/dvbfixer/model/diffusion/covalent_ring_probe.py",
             PROJECT / "src/dvbfixer/model/diffusion/loop_surface_probe.py"]
    bindings = dict(receipt["bindings"])
    jobs, tasks = [], []
    for index, run in enumerate(runs):
        tasks.append({"kind": "saved-arm", **run})
        jobs.append({"case_id": f"ring-arm-{index:02d}", "case": run["work"], "seed": 0})
    for name in CONTROL_NAMES:
        control = controls[name]
        path = Path(control["path"])
        if digest(path) != control["sha256"]:
            raise ValueError("Historical positive/negative/undefined control changed")
        files.append(path)
        tasks.append({"kind": "historical-control", "name": name, "path": str(path)})
        jobs.append({"case_id": name, "case": str(path), "seed": 0})
    for index, job in enumerate(jobs):
        job["command"] = [sys.executable, "-m", "scripts.audit_retained_1he8_ring_frames", "worker", "--plan",
            str(output / "plan.json"), "--index", str(index), "--output",
            str(output / "execution" / f"{job['case_id']}-seed-0")]
    bindings.update({str(path): digest(path) for path in files})
    output.mkdir(parents=True, exist_ok=False)
    save(output / "plan.json", {"python": sys.executable, "python_sha256": digest(Path(sys.executable)),
        "jobs": jobs, "tasks": tasks, "count": len(jobs), "ring": ring, "witness_bond": witness,
        "source_code_and_input_sha256": bindings, "planarity_max_angstrom": .10,
        "expected_arm_frame_instances": 467, "historical_control_structures": 4,
        "scope": "One previously implicated covalent ring only, all12 retained trajectories plus four retrospective controls",
        "diagnostic_not_acceptance_rule": True, "new_restoration_or_optimizer_calls": 0,
        "no_threshold_tuning": True, "resource_guard": "pressure-v2 serial CPU executor,1800s/job,RSS8GiB,10GiB output cap",
        "limits": "Missing/nonplanar surfaces remain unassessed, no macrocycles/closure or unsaved-path proof; controls use endpoint PDB topology rather than unavailable historical force graphs"})
    print(json.dumps({"arms": len(runs), "frame_instances": 467, "controls": 4}), flush=True)


def worker(plan_path, index, output):
    if {p.name for p in output.iterdir()} != {"command.json", "output.log"}:
        raise ValueError("Ring worker needs a fresh executor reservation")
    plan = json.loads(plan_path.read_text())
    if not 0 <= index < len(plan["tasks"]):
        raise ValueError("Undeclared task")
    task = plan["tasks"][index]
    ring = tuple(AtomIdentity(**row) for row in plan["ring"])
    witness = tuple(AtomIdentity(**row) for row in plan["witness_bond"])
    frames = []
    if task["kind"] == "saved-arm":
        work = Path(task["work"])
        particles = json.loads((work / "temporary-particle-axis.json").read_text())
        system = XmlSerializer.deserialize((work / "temporary-system.xml").read_text())
        atoms, all_bonds, _ = particle_graph(particles, system)
        records, _, _ = _parse_coordinate_records((work / "input.pdb").read_text().splitlines())
        heavy = {record.identity for record in records if record.element not in {"H", "D"}}
        if not heavy <= set(atoms):
            raise ValueError("Temporary topology lost input heavy atoms")
        bonds = frozenset(frozenset(pair) for pair in all_bonds if set(pair) <= heavy)
        request = DiffusionRequest.from_json(Path(task["request"]).read_text())
        ownership = frozenset(request.generated_atoms)
        for stage, path in (("observed-source", Path(task["source"])), ("initial-proposal", work / "input.pdb"),
                            ("final", Path(task["final"]))):
            points, _ = _coordinates(path.read_text())
            frames.append({"stage": stage, "path": str(path), "screen": focus_probe(points, bonds, ring, witness, ownership,
                           partial_source=stage == "observed-source")})
        for name in task["frames"]:
            path = Path(name)
            points = saved_coordinates(path, atoms)
            frames.append({"stage": path.parent.name, "path": name,
                           "screen": focus_probe(points, bonds, ring, witness, ownership)})
        provenance = "Actual serialized temporary harmonic-bond graph restricted to original input heavy atoms; no proximity bonds"
        identity = task["work"]
    else:
        path = Path(task["path"])
        text = path.read_text()
        points, _ = _coordinates(text)
        structure = _structure(text, "historical-ring-control.pdb")
        ids = {atom.index: AtomIdentity(atom.residue.chain.id, atom.residue.id, atom.residue.insertionCode.strip(), atom.name)
               for atom in structure.pdb.topology.atoms()}
        bonds = frozenset(frozenset((ids[a.index], ids[b.index])) for a, b in structure.pdb.topology.bonds())
        frames.append({"stage": task["name"], "path": str(path),
                       "screen": focus_probe(points, bonds, ring, witness, frozenset())})
        provenance = "OpenMM canonical-residue endpoint PDB topology; not the unavailable historical temporary force graph"
        identity = task["name"]
    save(output / "result.json", {"identity": identity, "kind": task["kind"], "ring": plan["ring"],
        "witness_bond": plan["witness_bond"], "frames": frames, "summary": summarize_frames(frames),
        "topology_provenance": provenance, "new_restoration_or_optimizer_calls": 0,
        "scientific_acceptance_changed": False, "causal_timing_or_unsaved_trajectory_established": False})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run", "worker"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--plan", type=Path)
    parser.add_argument("--index", type=int)
    args = parser.parse_args()
    if args.mode == "declare":
        declare(args.output.resolve())
    elif args.mode == "run":
        execute(args.output.resolve())
    else:
        if args.plan is None or args.index is None:
            parser.error("worker requires --plan and --index")
        worker(args.plan.resolve(), args.index, args.output.resolve())
