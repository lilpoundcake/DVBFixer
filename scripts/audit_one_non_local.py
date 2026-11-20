"""Read-only independent native-mask/readback for the private 1NON paired run."""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scripts.audit_retained_full_context_controls import atoms
from scripts.hybrid_gap_reconstruction import key
from scripts.mac_reconstruction_comparison import digest, save
from scripts.one_non_local_research import local_mask, movements
from scripts.retained_evidence import verify_bindings

from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter


def audit(root):
    plan = json.loads((root / "plan.json").read_text())
    verify_bindings(plan["source_code_and_input_sha256"], (root / "plan.json",))
    declaration = root / "refinement-declaration/plan.json"
    declared = json.loads(declaration.read_text())
    verify_bindings(declared["source_code_and_input_sha256"], (declaration,))
    if declared["parent_plan_sha256"] != digest(root / "plan.json"):
        raise ValueError("refinement declaration lineage changed")
    retry = root / "local-retry"
    if retry.is_dir():
        retry_plan = json.loads((retry / "plan.json").read_text())
        verify_bindings(retry_plan["source_code_and_input_sha256"], (retry / "plan.json",))
        if retry_plan["parent_plan_sha256"] != digest(root / "plan.json"):
            raise ValueError("local retry belongs to another run")
    request = DiffusionRequest.from_json((root / "request.json").read_text())
    if digest(root / "source.pdb") != request.normalized_pdb.sha256:
        raise ValueError("source checksum mismatch")
    source = atoms(root / "source.pdb")
    fixed, generated = set(source), set(map(key, request.generated_atoms))
    selected = set(local_mask(request, source))
    recorded = set(map(tuple, json.loads((root / "local-mask.json").read_text())["atoms"]))
    if selected != recorded:
        raise ValueError("source-only shell selection changed")
    results, initial = {}, {}
    for policy in ("fixed", "local"):
        allowed = selected if policy == "local" else set()
        for method in ("modeller", "rfd3"):
            name = f"{method}-{policy}"
            folder = (retry if policy == "local" and retry.is_dir() else root) / name
            result = json.loads((folder / "result.json").read_text())
            execution = json.loads((folder / "execution.json").read_text())
            if result["execution_status"] != "completed":
                results[name] = {"execution_status": result["execution_status"], "quality_status": "not_evaluated"}
                continue
            if execution["returncode"] or execution["guard_stop_reason"] is not None:
                raise ValueError("successful execution claim contradicts executor")
            before = json.loads((folder / "native-before.json").read_text())
            after = json.loads((folder / "native-after.json").read_text())
            native = json.loads((folder / "native.json").read_text())
            rows = before["atoms"]
            ids = [tuple(row["identity"]) for row in rows]
            x, y = np.array(before["coordinates"]), np.array(after["coordinates"])
            if len(set(ids)) != len(ids) or x.shape != y.shape or x.shape != (len(ids), 3) or not np.isfinite([x, y]).all():
                raise ValueError("native coordinate axis changed or is nonfinite")
            mapping = {a: i for i, a in enumerate(ids)}
            coords = {tuple(row["identity"]): row["xyz"] for row in native["atoms"]}
            if set(coords) != fixed | generated or not fixed | generated <= set(ids):
                raise ValueError("native requested atom axis changed")
            for row, a, p, q in zip(rows, ids, x, y, strict=True):
                expected = not row["virtual"] and (a in generated | allowed or row["hydrogen"])
                if row["mobile"] != expected or (not expected and np.linalg.norm(q-p) > 1e-7):
                    raise ValueError("native Cartesian mask/fixation failed")
                if a in coords and not np.array_equal(q, coords[a]):
                    raise ValueError("returned requested coordinates differ from native optimizer state")
                if a in fixed and not np.allclose(p, source[a]["xyz"], atol=1e-7, rtol=0):
                    raise ValueError("packing moved observations")
            initial[name] = {a: x[mapping[a]] for a in fixed | generated}
            candidate = atoms(folder / "candidate.pdb")
            actual_assessment = assessment_adapter(request, (root / "source.pdb").read_text())(
                (folder / "candidate.pdb").read_text(), frozenset(a for a in request.fixed_atoms if key(a) in allowed))
            if json.loads(json.dumps(asdict(actual_assessment))) != result["request_geometry"]:
                raise ValueError("common final-coordinate geometry assessment does not reproduce")
            if set(candidate) != fixed | generated:
                raise ValueError("serialized candidate inventory changed")
            for a, row in candidate.items():
                if row["xyz"] != tuple(float(f"{v:.3f}") for v in coords[a]):
                    raise ValueError("serialized coordinates differ from native minimizer return")
                if a in fixed:
                    old = source[a]["line"]
                    new = row["line"]
                    if (a not in allowed and old != new) or old[:30] != new[:30] or old[54:] != new[54:]:
                        raise ValueError("source metadata/outside-mask records changed")
            motion, outside = movements(source, candidate, allowed)
            if outside != 0 or motion != result["local_observed_movements"]:
                raise ValueError("independent movement audit does not reproduce")
            within = all(row["within_limits"] for row in motion)
            native_motion, native_outside = movements(source, {a: {"xyz": p} for a, p in coords.items()}, allowed)
            if native_motion != result["native_local_observed_movements"]:
                raise ValueError("native movement evidence does not reproduce")
            within = within and native_outside <= 1e-7 and all(row["within_limits"] for row in native_motion)
            if result["validation_passed"] and (not within or not result["geometry_passed"]
                    or not actual_assessment.geometry_passed or actual_assessment.admission_failures):
                raise ValueError("execution was incorrectly promoted to accepted structure")
            results[name] = {"execution_status": "completed", "quality_status": result["quality_status"],
                "native_outside_mask_max_A": max(math.dist(source[a]["xyz"], coords[a]) for a in fixed-allowed),
                "within_local_bounds": within, "candidate_sha256": digest(folder / "candidate.pdb")}
    pairs = {}
    for method in ("modeller", "rfd3"):
        a, b = f"{method}-fixed", f"{method}-local"
        if a in initial and b in initial:
            pairs[method] = all(np.array_equal(initial[a][k], initial[b][k]) for k in fixed | generated)
            if not pairs[method]:
                raise ValueError("fixed/local cells did not start from identical packed requested atoms")
    return {"plan_sha256": digest(root / "plan.json"), "cells": results,
            "paired_packed_initial_states_identical": pairs, "new_diffusion_calls": 0, "public_qualification": False}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    result = audit(args.root.resolve())
    save(args.root / "independent-readback.json", result)
    print(json.dumps(result))
