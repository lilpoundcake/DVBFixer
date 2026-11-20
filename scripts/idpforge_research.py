"""Private single N-terminal IDPForge control; no public profile or repair."""
from __future__ import annotations

import argparse
import json
import os
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scripts.audit_retained_full_context_controls import atoms
from scripts.hybrid_gap_reconstruction import launch_worker
from scripts.mac_reconstruction_comparison import digest, save, target_axis
from scripts.retained_evidence import freeze_code, verify_bindings

from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.protein_atoms import _CANONICAL_ONE_TO_THREE, CANONICAL_HEAVY_ATOMS

WEIGHTS = "a623756262138376672afcb17f2f8bd9e6780b265b57aa1b5fbd7e8b0f5d4769"


def admit(request, source):
    if len(request.target_sequences) != 1 or len(request.gaps) != 1:
        raise ValueError("initial control requires one chain and one N-terminal gap")
    if request.heterogen_contexts or request.retained_explicit_links:
        raise ValueError("initial control refuses nonprotein context and explicit links")
    target, gap = request.target_sequences[0], request.gaps[0]
    order = target_axis(request)[target.chain]
    keys = [(r.chain, r.residue_number, r.insertion_code) for r in order]
    missing = {(r.chain, r.residue_number, r.insertion_code) for r in gap.generated_residues}
    length = len(missing)
    if (not 1 <= length < len(order) or missing != set(keys[:length])
            or gap.left_anchor is not None or gap.right_anchor != order[length]):
        raise ValueError("requires an anchored whole-residue N-terminal gap")
    # Native process_pdb groups by integer resSeq alone. Do not silently pass
    # insertion codes, duplicated numbers or discontinuous numbering into it.
    numbers = [int(r.residue_number) for r in order]
    if any(r.insertion_code for r in order) or numbers != list(range(numbers[0], numbers[0] + len(order))):
        raise ValueError("native template parser requires contiguous numbers without insertion codes")
    fixed = {(a.chain, a.residue_number, a.insertion_code, a.atom_name) for a in request.fixed_atoms}
    generated = {(a.chain, a.residue_number, a.insertion_code, a.atom_name) for a in request.generated_atoms}
    if set(source) != fixed or fixed & generated:
        raise ValueError("source/fixed ownership mismatch")
    expected_fixed, expected_generated, rows = set(), set(), []
    for i, (residue, letter) in enumerate(zip(order, target.sequence, strict=True)):
        rk = keys[i]
        name = _CANONICAL_ONE_TO_THREE[letter]
        canonical = {(*rk, a) for a in CANONICAL_HEAVY_ATOMS[name]}
        if i < length:
            expected_generated |= canonical
        else:
            expected_fixed |= canonical
            for a in canonical:
                row = source.get(a)
                if (row is None or row["resname"] != name or not row["line"].startswith("ATOM  ")
                        or not np.isfinite(row["xyz"]).all() or row["line"][16].strip()
                        or float(row["line"][54:60]) != 1.):
                    raise ValueError("source must be complete canonical full-occupancy heavy atoms")
        rows.append({**asdict(residue), "resname": name, "generated": i < length})
    if fixed != expected_fixed or generated != expected_generated:
        raise ValueError("noncanonical atoms or incomplete whole-gap ownership")
    return {"sequence": target.sequence, "chain": target.chain, "gap_length": length,
            "residues": rows, "fixed_atoms": len(fixed), "generated_atoms": len(generated),
            "cropped": False, "reference_coordinates_read": False}


def run(args):
    args.output.mkdir(parents=True, exist_ok=False)
    args.owned_output = True
    root = args.output.resolve()
    source, req = args.source.resolve(), args.request.resolve()
    request = DiffusionRequest.from_json(req.read_text())
    if digest(source) != request.normalized_pdb.sha256 or digest(args.weights) != WEIGHTS:
        raise ValueError("unexpected source or user-provided checkpoint bytes")
    transport = admit(request, atoms(source))
    (root / "source.pdb").write_bytes(source.read_bytes())
    (root / "request.json").write_bytes(req.read_bytes())
    save(root / "transport.json", transport)
    runtime_src = args.runtime_src.resolve()
    upstream = args.upstream.resolve()
    worker = Path("deploy/idpforge/research_worker.py").resolve()
    bindings = [source, req, args.weights.resolve(), Path(__file__).resolve(), worker,
                Path("deploy/idpforge/prepare_runtime.py").resolve(),
                Path("scripts/audit_idpforge_research.py").resolve(),
                Path("scripts/hybrid_gap_reconstruction.py").resolve(),
                Path("scripts/mac_reconstruction_comparison.py").resolve(),
                Path("scripts/audit_retained_full_context_controls.py").resolve(),
                Path("scripts/run_mac_reconstruction_guarded.py").resolve(),
                Path("scripts/retained_evidence.py").resolve(),
                upstream / "mk_ldr_template.py", upstream / "configs/sample.yml", upstream / "LICENSE",
                runtime_src / "compatibility.json", *runtime_src.rglob("*.py"),
                root / "source.pdb", root / "request.json", root / "transport.json"]
    plan = {"protocol": "private-idpforge-nterminal-raw-v1", "mode": args.mode,
            "template_ss_policy": args.template_ss_policy,
            "weights": str(args.weights.resolve()), "weights_sha256": WEIGHTS,
            "weights_provenance": "user-provided official-archive contents; local digest, no independently retrieved publisher checksum",
            "upstream": str(upstream), "runtime_src": str(runtime_src), "seed": 42,
            "device": "CPU", "threads": 1, "training_steps": 200, "inference_steps": 40,
            "sampling_count": 1 if args.mode == "sample" else 0,
            "template_preparation": "unchanged mk_ldr_template.main, prepend known missing sequence, nsample=1; no truncation",
            "template_seed_floor_A": 6.46, "template_seed_ceiling_A": 9.12, "template_seed_skew": .3,
            "native_output_restoration": False, "relaxation": False, "automatic_repair": False,
            "potentials": False, "success_count_loop": False, "public_capability": False,
            "source_code_and_input_sha256": {str(p): digest(p) for p in bindings},
            "timeout_seconds": args.timeout, "resource_policy": "pressure-v2/RSS8GiB"}
    if args.mode == "sample":
        if args.preflight is None:
            raise ValueError("sample requires independently audited matching preflight")
        p = args.preflight.resolve()
        previous = json.loads((p / "plan.json").read_text())
        readback = json.loads((p / "independent-readback.json").read_text())
        if (previous["mode"] != "preflight" or not readback["consumer_passed"]
                or readback["plan_sha256"] != digest(p / "plan.json")
                or previous["seed"] != plan["seed"]
                or previous.get("template_ss_policy", "stock") != args.template_ss_policy
                or any(previous["source_code_and_input_sha256"].get(str(f)) != digest(f)
                       for f in (worker, source, req, args.weights.resolve()))):
            raise ValueError("preflight does not verify this worker/source/request/checkpoint")
        verify_bindings(previous["source_code_and_input_sha256"], (p / "plan.json",))
        plan["preflight"] = str(p)
        plan["source_code_and_input_sha256"].update({str(p / name): digest(p / name)
            for name in ("plan.json", "independent-readback.json", "native-initial.npz", "template.npz")})
    save(root / "plan.json", plan)
    freeze_code(root / "plan.json")
    env = {**os.environ, "PYTHONPATH": str(runtime_src), "PYTHONNOUSERSITE": "1", "PYTHONDONTWRITEBYTECODE": "1",
           "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
           "HF_HUB_OFFLINE": "1", "PYTORCH_ENABLE_MPS_FALLBACK": "0"}
    command = [str(args.python.absolute()), str(worker), str(root)]
    save(root / "command.json", command)
    execution = launch_worker(command, root, env, args.timeout)
    verify_bindings(plan["source_code_and_input_sha256"])
    print(json.dumps({"output": str(root), "returncode": execution["returncode"], "guard_stop_reason": execution["guard_stop_reason"]}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("source", "request", "weights", "upstream", "runtime-src", "python", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--mode", choices=("preflight", "sample"), default="preflight")
    parser.add_argument("--preflight", type=Path)
    parser.add_argument("--template-ss-policy", choices=("stock", "radian-consistent"), default="stock")
    parser.add_argument("--timeout", type=int, default=900, choices=range(1, 1801))
    arguments = parser.parse_args()
    try:
        run(arguments)
    except Exception as exc:
        import traceback
        if getattr(arguments, "owned_output", False):
            save(arguments.output / "failure.json", {"error": str(exc), "traceback": traceback.format_exc()})
        raise
