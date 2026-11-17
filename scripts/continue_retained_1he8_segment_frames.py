#!/usr/bin/env python3
"""Finish the ten unreviewed saved1HE8 arms; never sample or optimize."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from openmm import HarmonicBondForce, XmlSerializer
from scripts.audit_retained_1he8_segment_frames import backbone_edges, frame_screen
from scripts.mac_reconstruction_comparison import digest, save
from scripts.run_mac_modeller_comparison import execute

from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest
from dvbfixer.model.diffusion.gap_portfolio import _coordinates

PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "tmp/mac-real-gap-continuation-2026-11-12"


def remaining_arms(rows, completed_works):
    by_work = {}
    for row in rows:
        work = str(Path(row["final_path"]).parent)
        if work in by_work or type(row["accepted"]) is not bool:
            raise ValueError("Repeated arm or invalid historical acceptance")
        by_work[work] = row
    if len(completed_works) != len(set(completed_works)) or not set(completed_works) <= by_work.keys():
        raise ValueError("Previously checked arms do not match retained history")
    if set(completed_works) != {work for work, row in by_work.items() if row["accepted"]}:
        raise ValueError("The existing narrow audit must cover exactly the two accepted arms")
    return [row for work, row in by_work.items() if work not in completed_works]


def declare(output):
    accepted_root = ROOT / "1HE8-accepted-endpoint-segment-frames-attempt01"
    accepted_plan_path, accepted_complete_path = [accepted_root / name for name in ("plan.json", "completed.json")]
    accepted_plan, accepted_complete = [json.loads(p.read_text()) for p in (accepted_plan_path, accepted_complete_path)]
    for name, expected in accepted_plan["bindings"].items():
        if digest(Path(name)) != expected:
            raise ValueError(f"Historical narrow-audit input changed: {name}")
    complete_works = [row["work"] for row in accepted_complete["runs"]]
    if complete_works != [row["work"] for row in accepted_plan["runs"]] or len(complete_works) != 2:
        raise ValueError("Historical narrow completion is incomplete")
    audits = [ROOT / name / "independent-readback.json" for name in
              ("1HE8-traced-mobility-independent-readback-attempt02", "1HE8-buffered-walls-independent-readback-attempt01")]
    rows = [row for path in audits for row in json.loads(path.read_text())["rows"]]
    selected = remaining_arms(rows, complete_works)
    if len(rows) != 12 or len(selected) != 10:
        raise ValueError("Historical12-arm/ten-remaining denominator changed")
    files = [*audits, accepted_plan_path, accepted_complete_path, Path(__file__),
             PROJECT / "scripts/audit_retained_1he8_segment_frames.py",
             PROJECT / "scripts/run_mac_modeller_comparison.py",
             PROJECT / "src/dvbfixer/model/diffusion/backbone_intersections.py"]
    jobs, runs = [], []
    for index, row in enumerate(selected):
        final = Path(row["final_path"])
        if digest(final) != row["final_sha256"]:
            raise ValueError("Previously audited final changed")
        work = final.parent
        execution_path = work / "execution.json"
        request_path = Path(json.loads(execution_path.read_text())["case"]) / "request.json"
        request = DiffusionRequest.from_json(request_path.read_text())
        source = request_path.parent / "source.pdb"
        if digest(source) != request.normalized_pdb.sha256:
            raise ValueError("Historical source/request custody changed")
        frames = []
        files += [execution_path, request_path, source, work / "input.pdb", final,
                  work / "temporary-particle-axis.json", work / "temporary-system.xml"]
        for phase in sorted(work.glob("minimize-*")):
            saved = sorted(phase.glob("callback-*.npz")) + [phase / "returned-state.npz"]
            metadata = [phase / name for name in ("declaration.json", "completed.json", "callbacks.jsonl")]
            if not all(path.is_file() for path in saved + metadata):
                raise ValueError("Missing retained optimizer evidence")
            files += saved + metadata
            frames += list(map(str, saved))
        if not frames or len(frames) != len(set(frames)):
            raise ValueError("Missing or repeated saved frames")
        runs.append({"work": str(work), "request": str(request_path), "source": str(source),
                     "final": str(final), "previous_joint_acceptance": row["accepted"], "frames": frames})
        jobs.append({"case_id": f"retained-arm-{index:02d}", "case": str(work), "seed": 0,
                     "command": [sys.executable, "-m", "scripts.continue_retained_1he8_segment_frames", "worker",
                         "--plan", str(output / "plan.json"), "--index", str(index), "--output",
                         str(output / "execution" / f"retained-arm-{index:02d}-seed-0")]})
    output.mkdir(parents=True, exist_ok=False)
    save(output / "plan.json", {"python": sys.executable, "python_sha256": digest(Path(sys.executable)),
        "jobs": jobs, "runs": runs, "count": len(jobs), "previously_checked_works": complete_works,
        "source_code_and_input_sha256": {str(p): digest(p) for p in files},
        "run_selector": "Exact complement of the previously audited two arms within immutable12-arm history",
        "new_sampler_or_optimizer_calls": 0, "threshold_angstrom": 0.75,
        "intersection_tolerance_angstrom": 1e-5, "fatal_gate": False,
        "scope": "Actual harmonic-bond backbone edges, source/proposal/final and all retained callback/returned frames",
        "resource_guard": "Existing pressure-v2 serial CPU executor:1800s/arm,RSS8GiB,10GiB output cap; seed0 is reservation only",
        "limits": "No ring/closure or unsaved-path qualification; no whole-context or acceptance change"})
    print(json.dumps({"runs": len(runs), "saved_frames": sum(len(r["frames"]) for r in runs)}), flush=True)


def particle_graph(particles, system):
    atoms = [AtomIdentity(p["chain"], p["residue_number"], p["insertion_code"], p["atom_name"]) for p in particles]
    if [p["index"] for p in particles] != list(range(len(atoms))) or len(set(atoms)) != len(atoms):
        raise ValueError("Particle identity/index axis changed")
    if system.getNumParticles() != len(atoms):
        raise ValueError("Particle count differs from serialized system")
    bonds = set()
    for force in system.getForces():
        if isinstance(force, HarmonicBondForce):
            for index in range(force.getNumBonds()):
                first, second, *_ = force.getBondParameters(index)
                bonds.add(tuple(sorted((atoms[first], atoms[second]))))
    if not bonds:
        raise ValueError("No explicit harmonic covalent graph")
    return atoms, bonds, backbone_edges(set(atoms), sorted(bonds))


def saved_coordinates(path, atoms):
    with np.load(path, allow_pickle=False) as archive:
        xyz = archive["positions_nm"] * 10
    if xyz.shape != (len(atoms), 3) or not np.isfinite(xyz).all():
        raise ValueError("Invalid saved coordinate frame")
    return dict(zip(atoms, xyz, strict=True))


def worker(plan_path, index, output):
    if {p.name for p in output.iterdir()} != {"command.json", "output.log"}:
        raise ValueError("Worker requires a fresh executor reservation")
    plan = json.loads(plan_path.read_text())
    if not 0 <= index < len(plan["runs"]):
        raise ValueError("Undeclared run index")
    job = plan["runs"][index]
    work = Path(job["work"])
    particles = json.loads((work / "temporary-particle-axis.json").read_text())
    system = XmlSerializer.deserialize((work / "temporary-system.xml").read_text())
    atoms, bonds, edges = particle_graph(particles, system)
    request = DiffusionRequest.from_json(Path(job["request"]).read_text())
    generated = frozenset(request.generated_atoms)
    if not generated <= set(atoms):
        raise ValueError("Serialized topology omits generated atom identities")
    frames = []
    for stage, path in (("observed-source", Path(job["source"])), ("initial-proposal", work / "input.pdb"),
                        ("final", Path(job["final"]))):
        coordinates, _ = _coordinates(path.read_text())
        frames.append({"stage": stage, "path": str(path), "screen": frame_screen(coordinates, edges, bonds, generated)})
    for name in job["frames"]:
        path = Path(name)
        metadata = json.loads((path.parent / "declaration.json").read_text())
        callback = None
        if path.name.startswith("callback-"):
            number = int(path.stem.split("-")[1])
            callbacks = [json.loads(line) for line in (path.parent / "callbacks.jsonl").read_text().splitlines()]
            matching = [row for row in callbacks if row["callback"] == number]
            if len(matching) != 1:
                raise ValueError("Saved callback lacks unique optimizer metadata")
            callback = matching[0]
        frames.append({"stage": path.parent.name, "path": name, "phase_metadata": metadata,
            "callback_metadata": callback, "screen": frame_screen(saved_coordinates(path, atoms), edges, bonds, generated)})
    result = {"work": str(work), "previous_joint_acceptance": job["previous_joint_acceptance"],
              "actual_covalent_backbone_segments": len(edges), "frames": frames,
              "frames_with_alarms": sum(bool(f["screen"]["alarms"]) for f in frames),
              "new_sampler_or_optimizer_calls": 0, "topological_defect_or_passage_established": False,
              "publication_gate": False, "unsaved_path_not_qualified": True}
    save(output / "result.json", result)


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
