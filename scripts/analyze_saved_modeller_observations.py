#!/usr/bin/env python3
"""Localize observed drift to saved MODELLER files, without replay or vendor reads."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scripts.mac_reconstruction_comparison import digest, normalize_modeller, save

from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.gap_portfolio import _coordinates
from dvbfixer.model.diffusion.geometry import weighted_kabsch


def analyze(roots, output):
    rows, bindings = [], {}
    for root in roots:
        plan = json.loads((root / "plan.json").read_text())
        bindings[str(root / "plan.json")] = digest(root / "plan.json")
        for job in plan["jobs"]:
            case = Path(job["case"])
            request_path, source_path = case / "request.json", case / "source.pdb"
            request = DiffusionRequest.from_json(request_path.read_text())
            if digest(source_path) != request.normalized_pdb.sha256:
                raise ValueError("native observation reference changed")
            reference, _ = _coordinates(source_path.read_text())
            fixed = tuple(sorted(request.fixed_atoms))
            proposal = root / "execution" / f"{job['case_id']}-seed-{job['seed']}" / "proposal"
            engine = json.loads((proposal / "native-outputs.json").read_text())
            paths = [("saved-ini", proposal / "target.ini"), ("saved-IL", proposal / "target.IL00000001.pdb")]
            paths += [("reported-initial-output", proposal / engine["initial"][0]["name"]),
                      ("reported-loop-output", proposal / engine["loop"][0]["name"])]
            native_by_stage = {}
            for stage, path in paths:
                bindings[str(path)] = digest(path)
                try:
                    normalized, _ = normalize_modeller(request, path.read_text())
                    points, _ = _coordinates(normalized)
                except ValueError as exc:
                    rows.append({"case": job["case_id"], "seed": job["seed"], "stage": stage,
                                 "path": str(path), "transport_refusal": str(exc)})
                    continue
                native = np.array([points[a] for a in fixed])
                source = np.array([reference[a] for a in fixed])
                transform = weighted_kabsch(native, source)
                fitted = transform.apply(native)
                raw_distance = np.linalg.norm(native-source, axis=1)
                fit_distance = np.linalg.norm(fitted-source, axis=1)
                native_by_stage[stage] = native
                rows.append({"case": job["case_id"], "seed": job["seed"], "stage": stage, "path": str(path),
                    "observations_exact_in_source_frame": bool(np.array_equal(native, source)),
                    "source_frame_RMS_max_angstrom": [float(np.sqrt(np.mean(raw_distance**2))), float(raw_distance.max())],
                    "observed_global_fit_RMS_max_angstrom": [float(np.sqrt(np.mean(fit_distance**2))), float(fit_distance.max())],
                    "stage_metadata": "actual saved filename/engine output list, not an inferred unsaved optimizer phase"})
            if "reported-initial-output" in native_by_stage and "reported-loop-output" in native_by_stage:
                delta = np.linalg.norm(native_by_stage["reported-initial-output"]-native_by_stage["reported-loop-output"], axis=1)
                rows.append({"case": job["case_id"], "seed": job["seed"], "stage": "initial-to-loop-observations",
                    "coordinates_exactly_equal": bool(np.all(delta == 0)),
                    "RMS_max_angstrom": [float(np.sqrt(np.mean(delta**2))), float(delta.max())]})
            for p in (request_path, source_path, proposal / "native-outputs.json"):
                bindings[str(p)] = digest(p)
    output.mkdir(parents=True, exist_ok=False)
    save(output / "observation-stage-readback.json", {"rows": rows, "bindings": bindings,
        "analyst_sha256": digest(Path(__file__)), "new_model_or_optimizer_calls": 0,
        "native_initialization_routine_or_cause_not_established": True,
        "vendor_source_read_unavailable": True, "seed_confound_between_protocols_retained": True})
    print(json.dumps({"saved_stage_measurements": len(rows), "new_generation_calls": 0}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--roots", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    analyze([p.resolve() for p in args.roots], args.output.resolve())
