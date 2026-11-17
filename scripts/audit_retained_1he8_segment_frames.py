#!/usr/bin/env python3
"""Read saved1HE8 frames against actual FF bonds; nonfatal diagnostic only."""
from __future__ import annotations

import argparse
import fcntl
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np
from openmm import HarmonicBondForce, XmlSerializer
from scripts.mac_reconstruction_comparison import digest, save

from dvbfixer.model.diffusion.backbone_intersections import BackboneSegment, screen_segments
from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest
from dvbfixer.model.diffusion.gap_portfolio import _coordinates

PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "tmp/mac-real-gap-continuation-2026-11-12"


def backbone_edges(axis, bonds):
    """Use supplied covalent bonds only, never residue numbering or proximity."""
    names = {"N", "CA", "C"}
    return tuple((a, b) for a, b in bonds if a in axis and b in axis
                 and a.atom_name in names and b.atom_name in names)


def frame_screen(coordinates, edges, bonds, generated):
    segments = tuple(BackboneSegment(a, b, tuple(coordinates[a]), tuple(coordinates[b]))
                     for a, b in edges if a in coordinates and b in coordinates)
    result = screen_segments(segments, generated, bonded_atoms=frozenset(frozenset(pair) for pair in bonds))
    row = asdict(result)
    row["segments"] = len(segments)
    for alarm in row["alarms"]:
        first, second = segments[alarm["first_index"]], segments[alarm["second_index"]]
        alarm["implicated_bonds"] = [[asdict(s.first), asdict(s.second)] for s in (first, second)]
    return row


def declare(output, *, accepted_only=False):
    output.mkdir(parents=True, exist_ok=False)
    bindings = {str(Path(__file__)): digest(Path(__file__))}
    runs = []
    for audit_name in ("1HE8-traced-mobility-independent-readback-attempt02",
                       "1HE8-buffered-walls-independent-readback-attempt01"):
        audit = ROOT / audit_name / "independent-readback.json"
        bindings[str(audit)] = digest(audit)
        for checked in json.loads(audit.read_text())["rows"]:
            if accepted_only and not checked["accepted"]:
                continue
            final = Path(checked["final_path"])
            if digest(final) != checked["final_sha256"]:
                raise ValueError("audited final structure changed")
            work = final.parent
            execution = json.loads((work / "execution.json").read_text())
            request_path = Path(execution["case"]) / "request.json"
            request = DiffusionRequest.from_json(request_path.read_text())
            # The historical private mobility worker explicitly reads source.pdb
            # beside request.json; normalized_pdb.path belongs to its old workspace.
            source = request_path.parent / "source.pdb"
            if digest(source) != request.normalized_pdb.sha256:
                raise ValueError("source/request binding changed")
            files = [work / "execution.json", work / "input.pdb", final, request_path, source,
                     work / "temporary-particle-axis.json", work / "temporary-system.xml"]
            frames = []
            for phase in sorted(work.glob("minimize-*")):
                phase_files = [phase / "declaration.json", phase / "completed.json", phase / "callbacks.jsonl"]
                saved = sorted(phase.glob("callback-*.npz")) + [phase / "returned-state.npz"]
                for p in phase_files + saved:
                    if not p.is_file():
                        raise ValueError("retained phase evidence incomplete")
                files += phase_files + saved
                frames += [str(p) for p in saved]
            bindings.update({str(p): digest(p) for p in files})
            runs.append({"work": str(work), "request": str(request_path), "source": str(source),
                         "final": str(final), "previous_joint_acceptance": checked["accepted"], "frames": frames})
            save(output / f"declared-run-{len(runs)-1:02d}.json", runs[-1])
            print(json.dumps({"retained_run_declared": len(runs), "frames": len(frames)}), flush=True)
    save(output / "plan.json", {"runs": runs, "bindings": bindings,
        "run_selector": "previously joint-accepted endpoints only; other10 remain unaudited" if accepted_only else "all12 retained arms",
        "scope": "N-CA/CA-C/actual C-N bond segments from serialized force-field graph; no inferred bridges",
        "threshold_angstrom": 0.75, "intersection_tolerance_angstrom": 1e-5,
        "all_retained_callback_and_returned_frames": True, "new_optimizer_or_model_calls": 0,
        "fatal_gate": False, "clean_sparse_frames_do_not_prove_no_passage": True,
        "timeout_seconds": 1800})


def run(output):
    plan = json.loads((output / "plan.json").read_text())
    for name, expected in plan["bindings"].items():
        if digest(Path(name)) != expected:
            raise ValueError(f"retained input changed: {name}")
    (output / "run-reservation").mkdir(exist_ok=False)
    # Blocking scheduler coordination, not progress polling.
    with (PROJECT / ".artifacts/mac-reconstruction-2026-11-11/ml.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        rows = []
        for index, job in enumerate(plan["runs"]):
            work = Path(job["work"])
            particles = json.loads((work / "temporary-particle-axis.json").read_text())
            atoms = [AtomIdentity(p["chain"], p["residue_number"], p["insertion_code"], p["atom_name"]) for p in particles]
            if [p["index"] for p in particles] != list(range(len(atoms))) or len(set(atoms)) != len(atoms):
                raise ValueError("temporary topology identity axis changed")
            system = XmlSerializer.deserialize((work / "temporary-system.xml").read_text())
            bonds = set()
            for force in system.getForces():
                if isinstance(force, HarmonicBondForce):
                    for i in range(force.getNumBonds()):
                        a, b, *_ = force.getBondParameters(i)
                        bonds.add(tuple(sorted((atoms[a], atoms[b]))))
            if system.getNumParticles() != len(atoms) or not bonds:
                raise ValueError("serialized force graph/particle axis mismatch")
            edges = backbone_edges(set(atoms), sorted(bonds))
            request = DiffusionRequest.from_json(Path(job["request"]).read_text())
            generated = frozenset(request.generated_atoms)
            frames = []
            for label, path in (("observed-source", Path(job["source"])), ("initial-proposal", work / "input.pdb"),
                                ("final", Path(job["final"]))):
                coordinates, _ = _coordinates(path.read_text())
                frames.append({"stage": label, "path": str(path), "screen": frame_screen(coordinates, edges, bonds, generated)})
            for name in job["frames"]:
                path = Path(name)
                xyz = np.load(path, allow_pickle=False)["positions_nm"]*10
                if xyz.shape != (len(atoms), 3) or not np.isfinite(xyz).all():
                    raise ValueError("retained optimizer frame invalid")
                metadata = json.loads((path.parent / "declaration.json").read_text())
                callback = None
                if path.name.startswith("callback-"):
                    number = int(path.stem.split("-")[1])
                    callbacks = [json.loads(line) for line in (path.parent / "callbacks.jsonl").read_text().splitlines()]
                    callback = next(row for row in callbacks if row["callback"] == number)
                frames.append({"stage": path.parent.name, "path": name, "phase_metadata": metadata,
                    "callback_metadata": callback, "screen": frame_screen(dict(zip(atoms, xyz, strict=True)), edges, bonds, generated)})
            result = {"work": str(work), "previous_joint_acceptance": job["previous_joint_acceptance"],
                "actual_covalent_backbone_segments": len(edges), "frames": frames,
                "frames_with_alarms": sum(bool(f["screen"]["alarms"]) for f in frames),
                "topological_defect_or_passage_established": False, "publication_gate": False,
                "unsaved_path_not_qualified": True}
            save(output / f"run-{index:02d}.json", result)
            rows.append({key: value for key, value in result.items() if key != "frames"} | {"frames": len(frames)})
        save(output / "completed.json", {"runs": rows, "new_optimizer_or_model_calls": 0,
            "fatal_policy_or_topology_qualification": False})
        print(json.dumps({"runs": len(rows), "frames": sum(row["frames"] for row in rows),
                          "runs_with_alarms": sum(bool(row["frames_with_alarms"]) for row in rows)}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--accepted-only", action="store_true")
    args = parser.parse_args()
    if args.mode == "declare":
        declare(args.output.resolve(), accepted_only=args.accepted_only)
    else:
        run(args.output.resolve())
