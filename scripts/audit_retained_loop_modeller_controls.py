#!/usr/bin/env python3
"""Recompute MODELLER transport/restoration without a complete reference."""
from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scripts.audit_mac_reconstruction_smoke import inventory
from scripts.audit_retained_loop_context_sampling import fixed_records
from scripts.mac_reconstruction_comparison import digest, normalize_modeller, save

from dvbfixer.ffutils.geometry import ChiralityError, assert_all_l
from dvbfixer.model.diffusion.anchor_acceptance import _assessment, _structure
from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest
from dvbfixer.model.diffusion.geometry import weighted_kabsch
from dvbfixer.model.diffusion.pdb_materialize import materialize_candidate_pdb
from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter


def all_l(text):
    structure = _structure(text, "audit.pdb")
    try:
        assert_all_l(structure.pdb.topology, structure.pdb.positions)
    except ChiralityError as exc:
        return False, str(exc)
    return True, None


def restore_generated_only(source_text, axis, aligned, native, request):
    points = dict(zip(axis, aligned, strict=True))
    generated = tuple(request.generated_atoms)
    return materialize_candidate_pdb(source_text, generated, [points[a] for a in generated],
        tuple(native[a][0] for a in generated), request, elements=tuple(native[a][1] for a in generated))


def audit(root, output, *, recover_restoration=False):
    plan = json.loads((root / "plan.json").read_text())
    source_snapshots = []
    for name, expected in plan["source_code_and_input_sha256"].items():
        if digest(Path(name)) != expected:
            snapshot = root / "entry-script-executed-snapshot.py"
            if Path(name).name != "run_retained_loop_modeller_controls.py" or not snapshot.is_file() or digest(snapshot) != expected:
                raise ValueError("MODELLER declaration code/input changed")
            source_snapshots.append({"declared_path": name, "executed_snapshot": str(snapshot), "sha256": expected})
    rows, bindings, recovered_payloads = [], {str(root / "plan.json"): digest(root / "plan.json")}, []
    for job in plan["jobs"]:
        work = root / "execution" / f"{job['case_id']}-seed-{job['seed']}"
        receipt_path = work / "execution.json"
        receipt = json.loads(receipt_path.read_text())
        if receipt["returncode"] or receipt["guard_stop_reason"] or not receipt["completed"]:
            rows.append({"case": job["case_id"], "technical_completion": False,
                         "returncode": receipt["returncode"], "resource_stop": receipt["guard_stop_reason"]})
            bindings[str(receipt_path)] = digest(receipt_path)
            continue
        case = Path(job["case"])
        request_path, source_path = case / "request.json", case / "source.pdb"
        request = DiffusionRequest.from_json(request_path.read_text())
        if request.seeds != (job["seed"],) or digest(source_path) != request.normalized_pdb.sha256:
            raise ValueError("MODELLER request/source mismatch")
        reported = json.loads((work / "result.json").read_text())
        if receipt["result"] != reported:
            raise ValueError("MODELLER receipt/result mismatch")
        raw_native_path = work / "proposal/native.pdb"
        normalized, _ = normalize_modeller(request, raw_native_path.read_text())
        if normalized != (work / "proposal/mapped-native.pdb").read_text():
            raise ValueError("native ordinal/name transport changed")
        if digest(case / "alignment.pir") != digest(work / "proposal/alignment.pir"):
            raise ValueError("native PIR differs from the preflighted declaration")
        source, native = inventory(source_path), inventory(work / "proposal/mapped-native.pdb")
        axis, fixed = tuple(sorted(native)), tuple(sorted(source))
        if set(axis) != set(request.fixed_atoms) | set(request.generated_atoms) or set(fixed) != set(request.fixed_atoms):
            raise ValueError("native/full-source atom ownership changed")
        transform = weighted_kabsch(np.array([native[a][2] for a in fixed]), np.array([source[a][2] for a in fixed]))
        aligned = transform.apply(np.array([native[a][2] for a in axis]))
        with np.load(work / "native-transport.npz", allow_pickle=False) as arrays:
            if (not np.array_equal(arrays["native"], np.array([native[a][2] for a in axis]))
                    or not np.allclose(arrays["aligned"], aligned, rtol=0, atol=1e-8)):
                raise ValueError("native observed-only alignment not reproduced")
        points = dict(zip(axis, aligned, strict=True))
        native_text = rewrite_generated_coordinates(normalized, points)
        source_text = source_path.read_text()
        if plan.get("restoration_policy") == "generated-only source copy":
            worker_restored = restore_generated_only(source_text, axis, aligned, native, request)
        else:
            worker_restored = materialize_candidate_pdb(source_text, axis, aligned, tuple(native[a][0] for a in axis), request,
                                                       elements=tuple(native[a][1] for a in axis))
        if native_text != (work / "aligned-native.pdb").read_text() or worker_restored != (work / "source-restored.pdb").read_text():
            raise ValueError("native or original worker materialization not reproduced")
        invalid_worker_restoration = fixed_records(worker_restored, set(request.fixed_atoms)) != fixed_records(source_text, set(request.fixed_atoms))
        if invalid_worker_restoration and not recover_restoration:
            raise ValueError("source observations changed after restoration")
        restored = restore_generated_only(source_text, axis, aligned, native, request)
        if fixed_records(restored, set(request.fixed_atoms)) != fixed_records(source_text, set(request.fixed_atoms)):
            raise ValueError("source observations changed after restoration")
        drift = np.array([np.linalg.norm(points[a]-source[a][2]) for a in fixed])
        native_assessment, _ = _assessment(request, _structure(source_text, "source.pdb"), _structure(native_text, "native.pdb"))
        restored_assessment = assessment_adapter(request, source_text)(restored)
        if json.loads(json.dumps(asdict(native_assessment))) != reported["native_assessment"]:
            raise ValueError("independent native assessment disagrees with worker")
        if not invalid_worker_restoration and json.loads(json.dumps(asdict(restored_assessment))) != reported["source_restored_assessment"]:
            raise ValueError("independent restored assessment disagrees with worker")
        native_l, native_error = all_l(native_text)
        restored_l, restored_error = all_l(restored)
        fixation_verified = False
        if plan.get("native_interface"):
            fixation_path = work / "native-fixation-evidence.json"
            fixation = json.loads(fixation_path.read_text())
            selections = sorted(work.glob("initial-selection-*.json"))
            if (fixation["seed"] != job["seed"] or fixation["selection_calls"] != len(selections)
                    or not selections or not fixation["no_posthoc_native_coordinate_clamp"]):
                raise ValueError("native fixation evidence incomplete")
            for path in selections:
                selected = tuple(AtomIdentity(**a) for a in json.loads(path.read_text())["atoms"])
                if len(set(selected)) != len(selected) or set(selected) != set(request.generated_atoms):
                    raise ValueError("native initial optimization mask changed")
                bindings[str(path)] = digest(path)
            bindings[str(fixation_path)] = digest(fixation_path)
            fixation_verified = True
        passed = restored_assessment.geometry_passed and not restored_assessment.admission_failures and restored_l
        if ((not invalid_worker_restoration and passed != reported["source_restored_geometry_and_L_passed"])
                or native_l != reported["native_global_all_l"]):
            raise ValueError("independent chirality/acceptance disagrees")
        corrected_path = output / job["case_id"] / "source-restored.pdb"
        recovered_payloads.append((corrected_path, restored))
        rows.append({"case": job["case_id"], "native_seed": job["seed"], "technical_completion": True,
            "native_geometry_passed": native_assessment.geometry_passed, "native_global_all_l": native_l,
            "native_chirality_failure": native_error, "restored_global_all_l": restored_l,
            "native_initial_generated_selection_verified": fixation_verified,
            "native_observed_coordinates_exact_before_fit": all(np.array_equal(native[a][2], source[a][2]) for a in fixed),
            "restored_chirality_failure": restored_error, "source_restored_passed": passed,
            "source_restored_assessment": asdict(restored_assessment), "source_fixed_records_exact": True,
            "source_restored_path": str(corrected_path),
            "original_worker_restoration_invalid": invalid_worker_restoration,
            "transport_only_recovery_no_generation": recover_restoration,
            "native_observed_drift_after_global_observed_fit_angstrom": {"rms": float(np.sqrt(np.mean(drift**2))), "max": float(drift.max())},
            "all_generated_gaps_requested_jointly": len(request.gaps), "source_context": "whole observed assembly",
            "not_context_or_seed_matched_to_Protenix": True, "no_crystal_reference_or_RMSD": True,
            "new_sampling_or_repair_in_audit": False, "manual_approval_or_publication": False})
        for p in (receipt_path, work / "result.json", request_path, source_path, raw_native_path,
                  work / "proposal/mapped-native.pdb", work / "proposal/alignment.pir", work / "aligned-native.pdb",
                  work / "source-restored.pdb", work / "native-transport.npz"):
            bindings[str(p)] = digest(p)
    output.mkdir(parents=True, exist_ok=False)
    for path, text in recovered_payloads:
        path.parent.mkdir(exist_ok=False)
        with path.open("x") as stream:
            stream.write(text)
        bindings[str(path)] = digest(path)
    save(output / "independent-readback.json", {"rows": rows, "bindings": bindings,
        "auditor_sha256": digest(Path(__file__)), "restoration_recovery": recover_restoration,
        "declaration_source_snapshots": source_snapshots,
        "original_worker_outputs_unchanged": True, "not_broad_method_superiority_or_final_evaluation": True})
    print(json.dumps({"cases": len(rows), "source_restored_passes": sum(row.get("source_restored_passed", False) for row in rows)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--recover-restoration", action="store_true")
    args = parser.parse_args()
    audit(args.root.resolve(), args.output.resolve(), recover_restoration=args.recover_restoration)
