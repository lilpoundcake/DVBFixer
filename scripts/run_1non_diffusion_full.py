#!/usr/bin/env python3
"""Explicit private full-tetramer cc94 control; no MODELLER or local fallback."""
from __future__ import annotations

import argparse
import json
import os
import runpy
import subprocess
import sys
from dataclasses import asdict
from pathlib import Path

from scripts.mac_reconstruction_comparison import digest, save
from scripts.retained_evidence import freeze_code, verify_bindings

PROJECT = Path(__file__).resolve().parents[1]
SEED = 23


def planning(preparation, sequence_request, output):
    from dvbfixer.model.diffusion import cc94
    from dvbfixer.model.diffusion.contract import DiffusionRequest
    from dvbfixer.model.diffusion.joint_gap_research import (
        complete_context_admission,
        plan_complete_context,
    )
    from dvbfixer.model.diffusion.protein_context import prepare_cc94_request
    from dvbfixer.model.diffusion_cli import build_cli_diffusion_request

    output.mkdir(parents=True, exist_ok=False)
    original = DiffusionRequest.from_json(sequence_request.read_text())
    source = preparation / "protein-only.pdb"
    receipt = json.loads((preparation / "preparation.json").read_text())
    if digest(source) != receipt["protein_only_sha256"] or receipt["reattach_removed_context"]:
        raise ValueError("protein-only preparation changed")
    sequences = {t.chain: t.sequence for t in original.target_sequences}
    if set(sequences) != {"A", "B", "C", "D"} or {len(s) for s in sequences.values()} != {179}:
        raise ValueError("this declaration is only the retained 1NON tetramer")
    fresh = build_cli_diffusion_request(source, sequences, seeds=(SEED,), profile=cc94.MPS_PROFILE, no_terminal=True)
    # Public local contexts only validate the parent request contract. The actual
    # native view below is rebuilt from ALL observations and ALL four gaps.
    parent = prepare_cc94_request(fresh, source.read_bytes(), "auto")
    view = plan_complete_context(parent, source.read_bytes(), tuple(range(len(parent.gaps))), seed=SEED)
    admission = complete_context_admission(view)
    plan = view.request.protein_context
    if (len(parent.gaps) != 4 or any(len(g.generated_residues) != 18 for g in parent.gaps)
            or len(view.positions) != 716 or len(plan.spans) != 4
            or {(s.chain, s.start, s.stop) for s in plan.spans} != {(c, 0, 179) for c in sequences}
            or plan.preserved_only or {p.atom for p in plan.fixed_conditioning} != set(fresh.fixed_atoms)):
        raise ValueError("actual native conditioning is not the complete four-chain task")
    if fresh.fixed_atoms != original.fixed_atoms or fresh.generated_atoms != original.generated_atoms:
        raise ValueError("cleaning/fresh request changed the retained protein task identities")
    (output / "parent-request.json").write_text(parent.to_json())
    (output / "full-context-request.json").write_text(view.request.to_json())
    save(output / "engine-view.json", asdict(view))
    save(output / "data-admission.json", {**admission, "seed": SEED, "source_sha256": digest(source),
        "represented_chains": [{"chain": s.chain, "start": s.start, "stop": s.stop} for s in plan.spans],
        "engine_chain_count": len(view.engine_view.request.target_sequences),
        "gaps": [{"chain": g.chain, "length": len(g.generated_residues)} for g in fresh.gaps],
        "generated_atoms": len(fresh.generated_atoms), "sequence_request_sha256": digest(sequence_request),
        "actual_sampling_input_is_full_view_not_parent_local_views": True,
        "checkpoints_loaded": 0, "scientific_runs": 0, "protein_only": True,
        "source_protein": "deposited PyrR bifunctional protein; known unchanged sequence; coordinate restoration only"})


def native_arguments(plan, destination, prepare_only=False):
    argv = [str(PROJECT / "deploy/protpardelle-1c/joint_gap_research.py"),
        "--parent-request", plan["parent_request"], "--source", plan["source"],
        "--baseline", plan["source"],  # source only; never a MODELLER/gap proposal
        "--gap-indices", "0", "1", "2", "3", "--whole-structure", "--complete-observed-context",
        "--seed", str(SEED), "--checkpoint", plan["checkpoint"],
        "--refinement-python", plan["refinement_python"], "--output", str(destination),
        "--acknowledge-private-joint"]
    if prepare_only:
        argv.append("--prepare-only")
    if plan.get("fixed_scaffold_inpainting"):
        argv.append("--fixed-scaffold-inpainting")
    if plan.get("placement"):
        argv.append("--placement")
    return argv


def declare(preparation, planning_dir, output, *, inpainting=False, placement=False):
    from dvbfixer.model.diffusion import cc94

    output.mkdir(parents=True, exist_ok=False)
    admission = json.loads((planning_dir / "data-admission.json").read_text())
    if (not admission["supported"] or admission["required_tokens"] != 716
            or admission["engine_chain_count"] != 4 or admission["output_only_observed_atoms"]):
        raise ValueError("full-input planning did not pass")
    sampling = PROJECT / "dvbfixer/sampling/bin/python"
    refinement = PROJECT / "dvbfixer/refinement/bin/python"
    engine = PROJECT / ".artifacts/apple-cc94/protpardelle-source"
    checkpoint = PROJECT / ".artifacts/apple-cc94/model_params/weights" / cc94.CHECKPOINT_NAME
    if digest(checkpoint) != cc94.CHECKPOINT_SHA256:
        raise ValueError("checkpoint mismatch; no alternative model or download")
    paths = [Path(__file__), checkpoint, checkpoint.parent.parent / "configs/cc94.yaml", sampling, refinement]
    paths += list(planning_dir.glob("*"))+list(preparation.glob("*"))
    paths += list((PROJECT / "src/dvbfixer").rglob("*.py"))
    paths += list((PROJECT / "scripts").glob("*.py"))
    paths += list((PROJECT / "deploy/protpardelle-1c").glob("*.py"))
    paths += list((engine / "src").rglob("*.py"))
    plan = {"protocol": "private-1non-protein-only-complete-tetramer-cc94-v2-inpainting" if inpainting else "private-1non-protein-only-complete-tetramer-cc94-v1", "seed": SEED,
        "python": str(sampling), "refinement_python": str(refinement), "checkpoint": str(checkpoint),
        "engine_source": str(engine), "parent_request": str(planning_dir / "parent-request.json"),
        "source": str(preparation / "protein-only.pdb"),
        "source_code_and_input_sha256": {str(p): digest(p) for p in paths if p.is_file()},
        "model_calls": 1, "ode_updates": 500, "candidate_count": 1,
        "fixed_scaffold_inpainting": inpainting,
        "inpainting_policy": "private-cc94-fixed-scaffold-ve-inpainting-v1" if inpainting else None,
        "comparison": "same seed23/input/checkpoint/schedule as retained v1; new declared known-state intervention, not independent seeds or a silent retry",
        "output_mask": "canonical predictions plus fixed input-only OXT carried throughout the native process" if inpainting else "native canonical residue prediction mask",
        "input": "all four179-residue chains, all four18-residue gaps, all observed atoms conditioning; no preserved-only atoms",
        "native_fixed_gate": "audit raw native arrays after global rigid-frame fit BEFORE source copying; maximum displacement <=1e-4A",
        "native_float_tolerance_angstrom": 1.e-4, "generated_only_refinement": False,
        "trajectory_float_tolerance_angstrom": 1.e-4,
        "restarts": 0, "sampling_backend": "cc94 native MPS; fallback disabled",
        "resources": "pressure-v2; one ML process; wall1800s;RSS8GiB;MPS6GiB;10GiB evidence",
        "budget_scope": "one full-tetramer diffusion smoke; no cohort or automatic extension",
        "qualification": False, "production_changes": False}
    if placement:
        if not inpainting:
            raise ValueError('placement requires inpainting')
        plan.update(protocol='private-1non-cc94-residue-translation-v1', placement=True,
                    comparison='same seed23 raw fixed-scaffold control; assess raw generated-fixed, generated-generated clashes and eight junctions',
                    placement_rule='steps 250..499: nearest eligible fixed heavy atom within 3 A; same-chain adjacent rows excluded; translate all generated atoms of each residue by half the deficit, capped 0.35 A; consistent self-conditioning')
    save(output / "plan.json", plan)
    freeze_code(output / "plan.json")


def environment(plan):
    return {**os.environ, "PYTHONPATH": os.pathsep.join((str(PROJECT), str(PROJECT / "src"),
            str(Path(plan["engine_source"]) / "src"))),
        "PYTHONNOUSERSITE": "1", "PYTHONDONTWRITEBYTECODE": "1", "PYTORCH_ENABLE_MPS_FALLBACK": "0",
        "OPENMM_CPU_THREADS": "1", "OMP_NUM_THREADS": "1", "WANDB_MODE": "disabled"}


def load_plan(root):
    path = root / "plan.json"
    index = json.loads((root / "executed-code/index.json").read_text())
    if digest(path) != index["plan_sha256"]:
        raise ValueError("execution plan changed after archival")
    plan = json.loads(path.read_text())
    verify_bindings(plan["source_code_and_input_sha256"])
    return plan


def run(root, prepare_only=False):
    plan = load_plan(root)
    if prepare_only:
        command = [plan["python"], *native_arguments(plan, root / "native-data", True)]
        with (root / "native-data.log").open("xb") as log:
            subprocess.run(command, env=environment(plan), stdout=log, stderr=subprocess.STDOUT,
                           timeout=180, check=True)
        with (root / "native-reader.log").open("xb") as log:
            subprocess.run([plan["python"], "-m", "scripts.run_1non_diffusion_full", "reader", "--root", str(root)],
                           env=environment(plan), stdout=log, stderr=subprocess.STDOUT, timeout=120, check=True)
        return
    data = json.loads((root / "native-data/trace.json").read_text())
    if data["tokens"] != 716 or data["generated_gaps"] != 4 or data["status"] != "input-prepared":
        raise ValueError("pinned adapter data preparation differs from complete plan")
    reader = json.loads((root / "native-reader.json").read_text())
    if (reader["decoded_fixed_atoms"] != 5040 or reader["decoded_observed_residues"] != 644
            or reader["chain_count"] != 4 or reader["motif_sha256"] != digest(root / "native-data/motif.pdb")):
        raise ValueError("native reader did not verify the complete observed tetramer")
    command = [sys.executable, "-m", "scripts.run_mac_reconstruction_guarded", "--python", plan["python"],
        "--module", "scripts.run_1non_diffusion_full", "--output", str(root / "guarded"), "--timeout", "1800",
        "--resource-policy", "pressure-v2", "--", "worker", "--root", str(root)]
    save(root / "launch.json", {"command": command, "plan_sha256": digest(root / "plan.json")})
    result = subprocess.run(command, env=environment(plan), check=False)
    from scripts.audit_1non_diffusion_full import run as audit
    audit(root, root / "independent-readback")
    if result.returncode:
        raise RuntimeError("guarded diffusion did not complete; independent refusal report retained")


def native_reader(root):
    import numpy as np
    from protpardelle.common import residue_constants
    from protpardelle.data.pdb_io import read_pdb

    plan = load_plan(root)
    motif = root / "native-data/motif.pdb"
    protein, _, labels = read_pdb(motif)
    expected = {}
    for line in motif.read_text().splitlines():
        if line.startswith("ATOM  "):
            expected[(labels[line[21]], int(line[22:26]), line[12:16].strip())] = np.array(
                [float(line[i:i+8]) for i in (30, 38, 46)])
    decoded, errors = set(), []
    for row, (chain, number) in enumerate(zip(protein.chain_index, protein.residue_index, strict=True)):
        for name, column in residue_constants.atom_order.items():
            if protein.atom_mask[row, column]:
                key = (int(chain), int(number), name)
                if key not in expected:
                    raise ValueError("native reader invented an observation")
                decoded.add(key)
                errors.append(float(np.linalg.norm(protein.atom_positions[row, column]-expected[key])))
    if decoded != set(expected) or max(errors) > 1.e-5:
        raise ValueError("native reader changed the fixed atom mask/coordinates")
    save(root / "native-reader.json", {"decoded_fixed_atoms": len(decoded),
        "decoded_observed_residues": len(protein.residue_index), "chain_count": len(set(protein.chain_index)),
        "decode_max_error_angstrom": max(errors), "motif_sha256": digest(motif),
        "source_sha256": digest(Path(plan["source"])), "model_weights_loaded": False,
        "sampling_calls": 0, "native_reader_transport_only": True})


def worker(root):
    plan = load_plan(root)
    argv = native_arguments(plan, root / "native-sampling")
    sys.argv = argv
    runpy.run_path(argv[0], run_name="__main__")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("planning", "declare", "data", "reader", "run", "worker"))
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--preparation", type=Path)
    parser.add_argument("--sequence-request", type=Path)
    parser.add_argument("--planning", type=Path)
    parser.add_argument("--fixed-scaffold-inpainting", action="store_true")
    parser.add_argument("--placement", action="store_true")
    args = parser.parse_args()
    if args.mode == "planning":
        if args.preparation is None or args.sequence_request is None:
            parser.error("planning requires --preparation/--sequence-request")
        planning(args.preparation.resolve(), args.sequence_request.resolve(), args.root.resolve())
    elif args.mode == "declare":
        if args.preparation is None or args.planning is None:
            parser.error("declare requires --preparation/--planning")
        declare(args.preparation.resolve(), args.planning.resolve(), args.root.resolve(), inpainting=args.fixed_scaffold_inpainting, placement=args.placement)
    elif args.mode == "worker":
        worker(args.root.resolve())
    elif args.mode == "reader":
        native_reader(args.root.resolve())
    else:
        run(args.root.resolve(), prepare_only=args.mode == "data")
