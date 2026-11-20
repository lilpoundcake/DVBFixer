"""Read-only FrameDiPT native consumer/frame audit; no coordinate restoration."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scripts.audit_retained_full_context_controls import atoms
from scripts.mac_reconstruction_comparison import digest, save, target_axis
from scripts.retained_evidence import verify_bindings

from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.protein_atoms import ATOM37_ORDER
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter


def compare_frames(reference, actual, fixed):
    if reference.shape != actual.shape or reference.shape != (len(fixed), 7):
        raise ValueError("native frame axis changed")
    if not np.isfinite(actual).all():
        raise ValueError("nonfinite native frame")
    translation = float(np.linalg.norm(reference[fixed, 4:] - actual[fixed, 4:], axis=-1).max())
    # Quaternions q and -q denote the same rotation. No per-chain fitting.
    a, b = reference[fixed, :4], actual[fixed, :4]
    rotation = float(np.minimum(np.linalg.norm(a-b, axis=-1), np.linalg.norm(a+b, axis=-1)).max())
    if translation > 1e-4 or rotation > 1e-5:
        raise ValueError("native fixed frames moved before any copying")
    return {"translation_max_A": translation, "quaternion_max": rotation}


def audit(root):
    plan_path = root / "plan.json"
    plan = json.loads(plan_path.read_text())
    verify_bindings(plan["source_code_and_input_sha256"], (plan_path,))
    execution = json.loads((root / "execution.json").read_text())
    if execution["returncode"] or execution["guard_stop_reason"] is not None:
        raise ValueError("worker did not complete under its resource guard")
    request = DiffusionRequest.from_json((root / "request.json").read_text())
    if digest(root / "source.pdb") != request.normalized_pdb.sha256:
        raise ValueError("source digest mismatch")
    source = atoms(root / "source.pdb")
    axis = target_axis(request)
    order = [r for t in request.target_sequences for r in axis[t.chain]]
    identities = [(r.chain, r.residue_number, r.insertion_code) for r in order]
    gap_residues = {(r.chain, r.residue_number, r.insertion_code) for g in request.gaps for r in g.generated_residues}
    fixed = np.array([r not in gap_residues for r in identities])
    expected_pos = np.zeros((len(order), 37, 3))
    expected_mask = np.zeros((len(order), 37))
    index = {r: i for i, r in enumerate(identities)}
    for a, row in source.items():
        i, j = index[a[:3]], ATOM37_ORDER.index(a[3])
        expected_pos[i, j] = row["xyz"]
        expected_mask[i, j] = 1
    initial = dict(np.load(root / "native-initial.npz", allow_pickle=False))
    if (not np.array_equal(initial["fixed_mask"], fixed) or not np.all(initial["res_mask"] == 1)
            or not np.array_equal(initial["atom_positions"], expected_pos)
            or not np.array_equal(initial["atom_mask"], expected_mask)):
        raise ValueError("native consumer preparation dropped or invented source coordinates/masks")
    compare_frames(initial["rigids_0"], initial["rigids_t"], fixed)
    expected_seq, expected_chain = [], []
    offset = 0
    restypes = "ARNDCQEGHILKMFPSTWYV"
    for chain, target in enumerate(request.target_sequences):
        expected_seq.extend(range(offset, offset + len(target.sequence)))
        expected_chain.extend([chain] * len(target.sequence))
        offset += len(target.sequence) + 200
    sequence = "".join(t.sequence for t in request.target_sequences)
    if (initial["seq_idx"].tolist() != expected_seq or initial["chain_idx"].tolist() != expected_chain
            or initial["aatype"].tolist() != [restypes.index(a) for a in sequence]):
        raise ValueError("native known-sequence/chain/index axis changed")
    trace = [json.loads(line) for line in (root / "consumer-trace.jsonl").read_text().splitlines()]
    expected_calls = 1 if plan["mode"] == "preflight" else plan["steps"] + 1
    if len(trace) != expected_calls or any(r["call"] != i+1 or r["target_residues"] != len(order) for i, r in enumerate(trace)):
        raise ValueError("native model call trace is incomplete")
    if any(r["fixed_translation_max_A"] > 1e-4 or r["fixed_rotation_max_component"] > 1e-5 for r in trace):
        raise ValueError("native model saw moved observations")
    report = {"plan_sha256": digest(plan_path), "auditor_sha256": digest(Path(__file__)),
              "target_residues": len(order), "fixed_residues": int(fixed.sum()),
              "generated_residues": int((~fixed).sum()), "source_heavy_atoms": len(source),
              "native_model_calls": len(trace), "complete_sequence_and_residue_axis": True,
              "sidechain_conditioning": False, "strict_atomic_preservation": False,
              "public_qualification": False}
    if plan["mode"] == "sample":
        endpoint = dict(np.load(root / "native-endpoint.npz", allow_pickle=False))
        report["final_frames"] = compare_frames(initial["rigids_0"], endpoint["final_rigids"][0], fixed)
        policy = plan.get("frame_policy", "native")
        if policy == "reference-projection":
            projections = [json.loads(line) for line in (root / "projection-trace.jsonl").read_text().splitlines()]
            if len(projections) != plan["steps"]-1 or any(not p["generated_frames_exact"] or not p["fixed_tensor7_exact"] for p in projections):
                raise ValueError("fixed-frame projection trace is incomplete or changed generated frames")
            for i, row in enumerate(projections, 1):
                path = root / "projection-states" / f"{i:03d}.npz"
                if row["step"] != i or digest(path) != row["state_sha256"]:
                    raise ValueError("projection state binding changed")
                state = dict(np.load(path, allow_pickle=False))
                if (not np.array_equal(state["before"][..., ~fixed, :], state["after"][..., ~fixed, :])
                        or not np.array_equal(state["after"][..., fixed, :], initial["rigids_0"][None, fixed, :])):
                    raise ValueError("saved projection modified generated frames or failed to fix observations")
            report["fixed_frame_projection_steps"] = len(projections)
            report["unmodified_native_sampling"] = False
        else:
            report["unmodified_native_sampling"] = True
        bb = endpoint["backbone"]
        if bb.shape != expected_pos.shape or not np.isfinite(bb).all():
            raise ValueError("native backbone axis changed")
        boundary, internal = [], []
        cursor = 0
        for target in request.target_sequences:
            for i in range(cursor, cursor + len(target.sequence)-1):
                if fixed[i] and fixed[i+1]:
                    continue
                distance = float(np.linalg.norm(bb[i, 2] - bb[i+1, 0]))
                row = {"left": list(identities[i]), "right": list(identities[i+1]), "C_N_A": distance}
                (boundary if fixed[i] != fixed[i+1] else internal).append(row)
            cursor += len(target.sequence)
        report["native_backbone_peptides_touching_generation"] = {
            "junctions": boundary, "internal_min_A": min((r["C_N_A"] for r in internal), default=None),
            "internal_max_A": max((r["C_N_A"] for r in internal), default=None),
            "note": "distance inventory only; not full-atom chemical acceptance"}
        report["sampling_completed"] = True
        report["full_atom_completion"] = False
    else:
        report["sampling_completed"] = False
    return report


def audit_materialization(root):
    """Independently read native packing coordinates before source restoration."""
    plan_path = root / "plan.json"
    plan = json.loads(plan_path.read_text())
    if plan["protocol"] != "private-framedipt-backbone-plus-rosetta-sides-v1":
        raise ValueError("not a private FrameDiPT materialization")
    verify_bindings(plan["source_code_and_input_sha256"], (plan_path,))
    request = DiffusionRequest.from_json((root / "request.json").read_text())
    source, final = atoms(root / "source.pdb"), atoms(root / "candidate.pdb")
    if digest(root / "source.pdb") != request.normalized_pdb.sha256:
        raise ValueError("materialized source digest mismatch")
    fixed = {(a.chain, a.residue_number, a.insertion_code, a.atom_name) for a in request.fixed_atoms}
    generated = {(a.chain, a.residue_number, a.insertion_code, a.atom_name) for a in request.generated_atoms}
    native_rows = json.loads((root / "candidate-00/native.json").read_text())["atoms"]
    native = {tuple(r["identity"]): np.array(r["xyz"]) for r in native_rows}
    if len(native) != len(native_rows) or set(native) != fixed | generated or set(source) != fixed or set(final) != fixed | generated:
        raise ValueError("materialization ownership changed")
    if not np.isfinite(list(native.values())).all():
        raise ValueError("nonfinite native materialization")
    maximum = max(float(np.linalg.norm(native[a] - source[a]["xyz"])) for a in fixed)
    if maximum > 1e-7 or any(source[a]["line"] != final[a]["line"] for a in fixed):
        raise ValueError("materialization moved observations")
    if any(tuple(float(f"{v:8.3f}") for v in native[a]) != final[a]["xyz"] for a in generated):
        raise ValueError("generated coordinates differ from native packing result")
    neural = set(map(tuple, plan["frame_derived_backbone_CB"]))
    built = set(map(tuple, plan["rosetta_constructed_side_atoms"]))
    if neural & built or neural | built != generated:
        raise ValueError("FrameDiPT-derived versus constructed ownership mismatch")
    assessment = assessment_adapter(request, (root / "source.pdb").read_text())((root / "candidate.pdb").read_text())
    summary = json.loads((root / "summary.json").read_text())
    if assessment.geometry_passed != summary["geometry_passed"] or digest(root / "candidate.pdb") != summary["candidate_sha256"]:
        raise ValueError("materialized assessment/digest does not reproduce")
    return {"plan_sha256": digest(plan_path), "candidate_sha256": digest(root / "candidate.pdb"),
            "source_records_exact": True, "native_packing_fixed_max_A": maximum,
            "frame_derived_atoms": len(neural), "rosetta_constructed_atoms": len(built),
            "geometry_passed": assessment.geometry_passed, "native_diffusion_atomic_preservation": False,
            "public_qualification": False}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--materialized", action="store_true")
    args = parser.parse_args()
    result = audit_materialization(args.root.resolve()) if args.materialized else audit(args.root.resolve())
    save(args.root / "independent-readback.json", result)
    print(json.dumps(result))
