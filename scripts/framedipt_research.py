"""Private known-sequence FrameDiPT admission and CPU checkpoint controls.

No public capability, source repair, sequence design, cropping or all-atom
preservation claim. A residue-frame model does not condition on source sidechains.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scripts.audit_retained_full_context_controls import atoms
from scripts.hybrid_gap_reconstruction import launch_worker
from scripts.mac_reconstruction_comparison import digest, save, target_axis
from scripts.retained_evidence import freeze_code, verify_bindings

from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.protein_atoms import (
    _CANONICAL_ONE_TO_THREE,
    ATOM37_ORDER,
    CANONICAL_HEAVY_ATOMS,
)

REVISION = "74cf94cb7b4ed31d3a0d4230c5739d92cc061c50"
WEIGHTS_SHA256 = "b2d5768abdb0d03cfd40c4bd3fbc99cd0bda272b1f53d2b0d370eb1481790b51"
RESTYPES = "ARNDCQEGHILKMFPSTWYV"


def encode(request, source):
    """Source/request only: absent coordinates stay zero with masks false."""
    fixed = {(a.chain, a.residue_number, a.insertion_code, a.atom_name) for a in request.fixed_atoms}
    generated = {(a.chain, a.residue_number, a.insertion_code, a.atom_name) for a in request.generated_atoms}
    if set(source) != fixed or fixed & generated:
        raise ValueError("source differs from fixed ownership")
    if request.retained_explicit_links or request.heterogen_contexts:
        raise ValueError("initial FrameDiPT control refuses explicit links/nonprotein context")
    if any(g.left_anchor is None or g.right_anchor is None or not 3 <= len(g.generated_residues) <= 50 for g in request.gaps):
        raise ValueError("initial control supports only internal whole-residue gaps of 3..50")
    axis = target_axis(request)
    order = [r for target in request.target_sequences for r in axis[target.chain]]
    residue_to_index = {(r.chain, r.residue_number, r.insertion_code): i for i, r in enumerate(order)}
    if any(a[:3] not in residue_to_index for a in fixed | generated):
        raise ValueError("full source context must be represented on the target axis")
    n = len(order)
    positions, mask = np.zeros((n, 37, 3)), np.zeros((n, 37))
    aatype, seq_idx, chain_idx, fixed_mask = [], [], [], []
    offset = 0
    expected_fixed, expected_generated = set(), set()
    gap_residues = {(r.chain, r.residue_number, r.insertion_code) for gap in request.gaps for r in gap.generated_residues}
    rows = []
    for chain, target in enumerate(request.target_sequences):
        for ordinal, (r, letter) in enumerate(zip(axis[target.chain], target.sequence, strict=True)):
            rk = r.chain, r.residue_number, r.insertion_code
            resname = _CANONICAL_ONE_TO_THREE[letter]
            canonical = {(*rk, name) for name in CANONICAL_HEAVY_ATOMS[resname]}
            missing = rk in gap_residues
            if missing:
                if canonical - generated or any(a[:3] == rk for a in fixed):
                    raise ValueError("whole-gap atom ownership mismatch")
                expected_generated |= canonical
            else:
                actual = {a for a in fixed if a[:3] == rk}
                extras = actual - canonical
                if canonical - actual or extras - {(*rk, "OXT")} or (extras and ordinal != len(target.sequence) - 1):
                    raise ValueError("incomplete/noncanonical source residue or misplaced OXT")
                expected_fixed |= actual
                for a in actual:
                    row = source[a]
                    line = row["line"]
                    if (row["resname"] != resname or not line.startswith("ATOM  ")
                            or row["element"] in {"H", "D"} or not np.isfinite(row["xyz"]).all()):
                        raise ValueError("unsupported source chemistry/coordinates")
                    positions[residue_to_index[rk], ATOM37_ORDER.index(a[3])] = row["xyz"]
                    mask[residue_to_index[rk], ATOM37_ORDER.index(a[3])] = 1
            aatype.append(RESTYPES.index(letter))
            seq_idx.append(offset + ordinal)
            chain_idx.append(chain)
            fixed_mask.append(not missing)
            rows.append({**asdict(r), "resname": resname, "target_index": ordinal, "generated": missing})
        offset += len(target.sequence) + 200
    if expected_generated != generated or expected_fixed != fixed:
        raise ValueError("unaccounted fixed/generated atoms")
    arrays = {"atom_positions": positions, "atom_mask": mask, "aatype": np.array(aatype, dtype=np.int64),
              "seq_idx": np.array(seq_idx, dtype=np.int64), "chain_idx": np.array(chain_idx, dtype=np.int64),
              "fixed_mask": np.array(fixed_mask, dtype=np.float64), "res_mask": np.ones(n)}
    validate_transport(arrays)
    receipt = {"residues": rows, "chains": [{"chain": t.chain, "length": len(t.sequence)} for t in request.target_sequences],
               "target_residues": n, "fixed_residues": int(sum(fixed_mask)), "generated_residues": n-int(sum(fixed_mask)),
               "source_heavy_atoms": len(fixed), "generated_heavy_atoms_requested": len(generated),
               "frame_contributor_atoms": 3 * int(sum(fixed_mask)), "psi_contributor_oxygens": int(sum(fixed_mask)),
               "not_directly_conditioned_source_atoms": len(fixed) - 4 * int(sum(fixed_mask)),
               "sidechain_conditioning": False, "cropped": False, "reference_coordinates_read": False,
               "native_preservation_contract": "frames only; measure idealized backbone drift separately before copying"}
    return arrays, receipt


def validate_transport(arrays):
    n = len(arrays["aatype"])
    if arrays["atom_positions"].shape != (n, 37, 3) or arrays["atom_mask"].shape != (n, 37):
        raise ValueError("invalid atom37 transport shape")
    if not np.isfinite(arrays["atom_positions"]).all():
        raise ValueError("nonfinite atom37 transport")
    if not np.all(arrays["res_mask"] == 1) or not np.isin(arrays["fixed_mask"], [0, 1]).all():
        raise ValueError("missing targets must remain real, diffused residues")
    missing = arrays["fixed_mask"] == 0
    if np.any(arrays["atom_mask"][missing]) or np.any(arrays["atom_positions"][missing]):
        raise ValueError("unknown region contains coordinates or observation masks")
    if not np.all(arrays["atom_mask"][~missing, :3] == 1):
        raise ValueError("observed frames need N/CA/C")
    if not np.isin(arrays["aatype"], range(20)).all():
        raise ValueError("known sequence contains unknown amino acids")


def run(args):
    args.output.mkdir(parents=True, exist_ok=False)
    args.owned_output = True
    output = args.output.resolve()
    source_path, request_path = args.source.resolve(), args.request.resolve()
    upstream, weights = args.upstream.resolve(), args.weights.resolve()
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=upstream, text=True).strip()
    if revision != REVISION or digest(weights) != WEIGHTS_SHA256:
        raise ValueError("unrecognized FrameDiPT source revision/checkpoint")
    if subprocess.check_output(["git", "status", "--porcelain", "--untracked-files=no"], cwd=upstream, text=True).strip():
        raise ValueError("upstream source is modified")
    request = DiffusionRequest.from_json(request_path.read_text())
    if digest(source_path) != request.normalized_pdb.sha256:
        raise ValueError("source bytes differ from the frozen request")
    arrays, transport = encode(request, atoms(source_path))
    (output / "source.pdb").write_bytes(source_path.read_bytes())
    (output / "request.json").write_bytes(request_path.read_bytes())
    np.savez_compressed(output / "transport.npz", **arrays)
    save(output / "transport.json", transport)
    worker = Path("deploy/framedipt/research_worker.py").resolve()
    files = [source_path, request_path, weights, worker, Path(__file__).resolve(),
             Path("scripts/audit_framedipt_research.py").resolve(),
             Path(args.python).resolve(), Path("scripts/hybrid_gap_reconstruction.py").resolve(),
             Path("scripts/run_mac_reconstruction_guarded.py").resolve(),
             Path("scripts/mac_reconstruction_comparison.py").resolve(),
             Path("tests/test_framedipt_research.py").resolve(),
             *upstream.joinpath("framedipt").rglob("*.py"), *upstream.joinpath("openfold").rglob("*.py"),
             upstream / "experiments/utils.py", upstream / "config/base.yaml", upstream / "LICENSE.md",
             output / "source.pdb", output / "request.json", output / "transport.npz", output / "transport.json"]
    if args.paired_native is not None:
        files += [args.paired_native.resolve() / name for name in ("plan.json", "native-failure.json", "native-initial.npz")]
    plan = {"protocol": "private-framedipt-frame-inpainting-v1", "mode": args.mode, "device": "CPU",
            "frame_policy": args.frame_policy,
            "upstream": str(upstream), "upstream_revision": revision, "weights": str(weights),
            "runtime_python": str(args.python.absolute()),
            "paired_native_attempt": str(args.paired_native.resolve()) if args.paired_native else None,
            "weights_sha256": WEIGHTS_SHA256, "source_and_weights_license": "CC-BY-NC-SA-4.0",
            "source_code_and_input_sha256": {str(p): digest(p) for p in files},
            "seed": args.seed, "steps": 100, "min_t": .01, "noise_scale": .1,
            "native_inference_loop": "experiments.utils.inference_fn; native reverse or separately declared fixed-frame projection",
            "frame_translation_tolerance_A": 1e-4, "frame_rotation_tolerance": 1e-5,
            "native_atomic_preservation_claimed": False, "observed_coordinate_release": False,
            "sidechain_conditioning": False, "new_sequence_design": False,
            "sampling_requested": args.mode == "sample", "public_capability": False,
            "resources": {"timeout_seconds": args.timeout, "RSS_limit_GiB": 8, "threads": 1, "policy": "pressure-v2"}}
    save(output / "plan.json", plan)
    freeze_code(output / "plan.json")
    if args.mode == "sample":
        if args.preflight is None:
            raise ValueError("sample requires a retained matching preflight")
        preflight_plan = json.loads((args.preflight / "plan.json").read_text())
        preflight = json.loads((args.preflight / "preflight.json").read_text())
        independent = json.loads((args.preflight / "independent-readback.json").read_text())
        execution = json.loads((args.preflight / "execution.json").read_text())
        if (preflight_plan["mode"] != "preflight" or not preflight["consumer_passed"]
                or preflight_plan.get("frame_policy", "native") != args.frame_policy
                or execution["returncode"] or execution["guard_stop_reason"] is not None
                or independent["plan_sha256"] != digest(args.preflight / "plan.json")
                or independent["sampling_completed"] or not independent["complete_sequence_and_residue_axis"]
                or preflight_plan["source_code_and_input_sha256"][str(worker)] != digest(worker)
                or preflight_plan["source_code_and_input_sha256"][str(source_path)] != digest(source_path)
                or preflight_plan["source_code_and_input_sha256"][str(request_path)] != digest(request_path)):
            raise ValueError("preflight does not verify this complete source/request")
        verify_bindings(preflight_plan["source_code_and_input_sha256"], (args.preflight / "plan.json",))
        save(output / "preflight-binding.json", {"plan_sha256": digest(args.preflight / "plan.json"),
             "preflight_sha256": digest(args.preflight / "preflight.json")})
    env = {**os.environ, "PYTHONPATH": str(upstream), "PYTHONNOUSERSITE": "1", "PYTHONDONTWRITEBYTECODE": "1",
           "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
           "HF_HUB_OFFLINE": "1", "WANDB_MODE": "disabled", "PYTORCH_ENABLE_MPS_FALLBACK": "0"}
    # Preserve the venv entry path: resolving its interpreter symlink would
    # silently select the borrowed base prefix and lose the isolated overlay.
    command = [str(args.python.absolute()), str(worker), str(output)]
    save(output / "command.json", command)
    execution = launch_worker(command, output, env, args.timeout)
    verify_bindings(plan["source_code_and_input_sha256"])
    if args.paired_native is not None and (output / "native-initial.npz").is_file():
        previous = dict(np.load(args.paired_native / "native-initial.npz", allow_pickle=False))
        current = dict(np.load(output / "native-initial.npz", allow_pickle=False))
        matching = set(previous) == set(current) and all(np.array_equal(previous[k], current[k]) for k in previous)
        save(output / "paired-initial-state.json", {"native_plan_sha256": digest(args.paired_native / "plan.json"),
             "native_failure_sha256": digest(args.paired_native / "native-failure.json"),
             "all_initial_arrays_identical": matching, "resumed_native_run": False})
        if not matching:
            raise ValueError("paired projection arm did not reproduce the native initial state")
    print(json.dumps({"output": str(output), "returncode": execution["returncode"],
                      "guard_stop_reason": execution["guard_stop_reason"]}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("source", "request", "upstream", "weights", "python", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--mode", choices=("preflight", "sample"), default="preflight")
    parser.add_argument("--preflight", type=Path)
    parser.add_argument("--frame-policy", choices=("native", "reference-projection"), default="native")
    parser.add_argument("--paired-native", type=Path)
    parser.add_argument("--seed", type=int, default=23)
    parser.add_argument("--timeout", type=int, default=900, choices=range(1, 1801))
    args = parser.parse_args()
    try:
        run(args)
    except Exception as exc:
        import traceback
        if getattr(args, "owned_output", False):
            save(args.output / "failure.json", {"error": str(exc), "traceback": traceback.format_exc()})
        raise


if __name__ == "__main__":
    main()
