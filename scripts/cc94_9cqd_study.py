"""Bounded research comparison of cc94 and cc89 ODE sampling on 9CQD M155–162.

This adapter reuses the maintained atom mapping/materialization driver without
changing public checkpoint locks, terminal admission, or production defaults.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import time
import traceback
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / ".artifacts/cc94-9cqd"
SEEDS = (7, 19, 42, 101, 202)
CELLS = {"N0": ("cc94", False), "N1": ("cc94", True), "B0": ("cc89", False), "B1": ("cc89", True)}
PINS = {
    "cc94": (
        "680096576f6bcec1da7e9c30a89c7278c04cf1ab16ed9a74c2cb99e37fa4db74",
        "a58c009ead22b6f79ef61d44e258b432b0fada1dfd8d44cb7a137d7459038a15",
        3100,
    ),
    "cc89": (
        "e9999ace79bf3044351cc982a624fc3459add3a4f91cbf97ff5b437b8942eb9d",
        "dfc9895b399ec4497bf6d646502168725f01dcbd55bc0b541fc3e3cc1f2f0483",
        415,
    ),
}
REQUEST_SHA = "3e28e7a3f061506b2eb6eebb9a98a921ddd56284aa92645dae4c92986d175021"
ENGINE_REVISION = "ee378400f25b801fa481028000f9060183d7fb4c"
ENGINE = ROOT / ".artifacts/9cqd-cc89-engine"
PYTHON = ROOT / ".artifacts/9cqd-linux-study/cc89-env/bin/python"
CORE = Path("/data/user/skripkins/micromamba/envs/dvbfixer/bin/python")
ARCHIVE = ROOT / ".artifacts/9cqd-protenix-handoff"


def digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n")


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def environment():
    env = dict(os.environ)
    env.update(
        PYTHONPATH=os.pathsep.join([str(ROOT), str(ROOT / "src"), str(ENGINE / "src")]),
        PATH=str(PYTHON.parent) + os.pathsep + str(CORE.parent) + os.pathsep + env.get("PATH", ""),
        CUDA_VISIBLE_DEVICES="0",
        OMP_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1",
        PROTPARDELLE_MODEL_PARAMS=str(WORK / "model_params"),
        PROTPARDELLE_OUTPUT_DIR=str(WORK / "runtime"),
        PYTHONDONTWRITEBYTECODE="1",
    )
    return env


def command(label, args, *, timeout=3600):
    log = WORK / "logs" / f"{label}.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    record = {
        "argv": [str(a) for a in args],
        "cwd": str(ROOT),
        "timeout_seconds": timeout,
        "started_unix": time.time(),
        "environment": {
            k: v
            for k, v in environment().items()
            if k
            in {
                "PYTHONPATH",
                "PATH",
                "CUDA_VISIBLE_DEVICES",
                "OMP_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "PROTPARDELLE_MODEL_PARAMS",
                "PROTPARDELLE_OUTPUT_DIR",
            }
        },
    }
    try:
        with log.open("w") as output:
            result = subprocess.run(
                record["argv"],
                cwd=ROOT,
                env=environment(),
                stdout=output,
                stderr=subprocess.STDOUT,
                timeout=timeout,
            )
        record["returncode"] = result.returncode
    except subprocess.TimeoutExpired:
        record.update(returncode=None, timed_out=True)
    record["finished_unix"] = time.time()
    write_json(WORK / "logs" / f"{label}.command.json", record)
    print(label, "exit", record["returncode"], log, flush=True)
    return record["returncode"] == 0


def assets(model):
    base = (
        WORK / "model_params"
        if model == "cc94"
        else ROOT / "test/protpardelle-1c-audit/model_params"
    )
    epoch = PINS[model][2]
    config, checkpoint = base / f"configs/{model}.yaml", base / f"weights/{model}_epoch{epoch}.pth"
    if (digest(config), digest(checkpoint)) != PINS[model][:2]:
        raise ValueError("checkpoint/config pin mismatch")
    return config, checkpoint


def inventory():
    import torch
    from protpardelle.core.models import Protpardelle, TrainingConfig, load_config

    rows = []
    for model in PINS:
        config, checkpoint = assets(model)
        network = Protpardelle(load_config(config, TrainingConfig), device="cpu")
        state = torch.load(checkpoint, map_location="cpu", weights_only=False)
        network.load_state_dict(state["model_state_dict"], strict=True)
        rows.append(
            {
                "model": model,
                "checkpoint_bytes": checkpoint.stat().st_size,
                "parameters": sum(p.numel() for p in network.parameters()),
                "task": network.task,
                "strict_state_dict_match": True,
                "position_encoding": network.config.model.struct_model.uvit.position_embedding_type,
                "sigma_data": float(network.sigma_data),
                "dummy_fill": network.config.data.dummy_fill_mode,
            }
        )
        del network, state
    tensor = torch.ones((4, 4), device="cuda")
    assert float((tensor @ tensor).sum().item()) == 64
    write_json(
        WORK / "model-inventory.json",
        {
            "models": rows,
            "python": sys.version,
            "torch": torch.__version__,
            "torch_cuda": torch.version.cuda,
            "gpu": torch.cuda.get_device_name(),
            "compute_capability": torch.cuda.get_device_capability(),
            "cuda_tensor_test": "passed",
        },
    )


def request_for(cell, seed, *, smoke=False):
    from dvbfixer.model.diffusion.contract import DiffusionRequest
    from dvbfixer.model.diffusion.protein_context import sampling_entities, validate_context_source

    if cell not in CELLS or seed not in SEEDS:
        raise ValueError("unregistered cell or seed")
    original = ARCHIVE / "inputs/original-gap0-seed7.json"
    if digest(original) != REQUEST_SHA:
        raise ValueError("archived selected-gap request changed")
    raw = json.loads(original.read_text())
    model, reinject = CELLS[cell]
    raw["seeds"] = [seed]
    raw["normalized_pdb"]["path"] = "source.pdb"
    options = {item["name"]: item["value"] for item in raw["backend_options"]}
    options["profile"] = f"9cqd-{model}-ode-{cell}-research-v1" + ("-hook-smoke" if smoke else "")
    options["cc89_chain_phase"] = "native" if model == "cc94" else "neutral-v1"
    options["cc89_position_features"] = "source-physical-v1"
    raw["backend_options"] = [{"name": k, "value": v} for k, v in options.items()]
    request = DiffusionRequest.from_json(json.dumps(raw))
    source = (ARCHIVE / "inputs/9CQD_split.pdb").read_bytes()
    if hashlib.sha256(source).hexdigest() != request.normalized_pdb.sha256:
        raise ValueError("source digest mismatch")
    validate_context_source(request, source)
    (gap,) = request.gaps
    target = next(t for t in request.target_sequences if t.chain == "M")
    assert target.sequence[gap.target_interval.start : gap.target_interval.stop] == "MGSKPNND"
    assert gap.left_anchor is None and gap.right_anchor.residue_number == "163"
    assert len(request.generated_atoms) == 58 and len(request.fixed_atoms) == 3441
    assert len(request.protein_context.fixed_conditioning) == 1399
    assert sum(len(e.residues) for e in sampling_entities(request)) == 189
    return request, source


def prepare():
    if (WORK / "protocol.json").exists():
        raise FileExistsError("protocol already frozen")
    hashes = {}
    for model in PINS:
        config, checkpoint = assets(model)
        hashes[model] = {"config_sha256": digest(config), "checkpoint_sha256": digest(checkpoint)}
    revision = subprocess.check_output(
        ["git", "-C", str(ENGINE), "rev-parse", "HEAD"], text=True
    ).strip()
    assert revision == ENGINE_REVISION
    if not command("model-inventory", [PYTHON, "-B", __file__, "inventory"]):
        raise RuntimeError("model inventory failed")
    from dvbfixer.model.diffusion.protein_context import (
        engine_context_view,
        source_context_position_features,
    )

    request, _ = request_for("N0", 7)
    view = engine_context_view(request)
    indices, chains = source_context_position_features(request)
    write_json(
        WORK / "context-axis.json",
        {
            "source_atom_map": [
                {"engine": asdict(a), "source": asdict(b)} for a, b in view.atom_to_source
            ],
            "residue_index": indices,
            "physical_chain_index": chains,
            "proxy_chain_count": len(view.request.target_sequences),
            "physical_chains": ["M", "W", "X"],
        },
    )
    for cell in CELLS:
        for seed in SEEDS:
            request, source = request_for(cell, seed)
            directory = WORK / "requests" / cell / f"seed-{seed}"
            directory.mkdir(parents=True, exist_ok=False)
            (directory / "request.json").write_text(request.to_json() + "\n")
            (directory / "source.pdb").write_bytes(source)
    sources = [
        Path(__file__),
        ROOT / "scripts/cc94_9cqd_metrics.py",
        ROOT / "deploy/protpardelle-1c/checkpoint_gap_smoke.py",
        ENGINE / "src/protpardelle/core/models.py",
        ROOT / "src/dvbfixer/model/diffusion/boundary_refinement.py",
        ROOT / "src/dvbfixer/model/diffusion/runner_refinement.py",
        ROOT / "src/dvbfixer/model/diffusion/validate.py",
    ]
    write_json(
        WORK / "protocol.json",
        {
            "schema_version": 1,
            "case": "9cqd-M155-162-one-anchor",
            "cells": {c: {"model": m, "reinjection": r} for c, (m, r) in CELLS.items()},
            "seeds": SEEDS,
            "steps": 500,
            "step_scale": 1.0,
            "s_churn": 0.0,
            "schedule": "ODE",
            "primary_samples": 20,
            "repeat_cells": ["N1", "B1"],
            "repeat_seed": 7,
            "additional_gap_filling": False,
            "generated_atom_count": 58,
            "fixed_atom_count": 3441,
            "context_tokens": 189,
            "conditioning_atoms": 1399,
            "sequence_fixed": True,
            "engine_revision": revision,
            "weights": hashes,
            "source_request_sha256": REQUEST_SHA,
            "metric_source_sha256": digest(ROOT / "scripts/cc94_9cqd_metrics.py"),
            "source_sha256": {str(p.relative_to(ROOT)): digest(p) for p in sources},
            "refinement": {
                "platform": "Reference",
                "max_iterations": 250,
                "restart_count": 8,
                "perturbation_angstrom": 0.75,
                "generated_only": True,
                "anchor_relaxation": False,
            },
            "limitations": [
                "cc94 differs in task, relative encoding, CA ghost filling and sigma_data, not just training mixture",
                "cc89 ODE controls are new; archived churn=200 arms are historical comparisons",
                "no native macOS hardware is available on this Linux host",
            ],
            "thresholds": {
                "severe_vdw_overlap_A": 0.9,
                "penetration_distance_A": 1.9,
                "segment_alarm_A": 0.5,
            },
        },
    )
    command("gpu-preflight", ["nvidia-smi", "-q"])
    command("python-packages", [PYTHON, "-m", "pip", "freeze"])
    print(WORK / "protocol.json")


def check_protocol():
    protocol = json.loads((WORK / "protocol.json").read_text())
    for path, expected in protocol["source_sha256"].items():
        if digest(ROOT / path) != expected:
            raise ValueError("frozen source changed: " + path)
    return protocol


def sample(cell, seed, attempt, *, smoke=False):
    import torch

    from dvbfixer.model.diffusion.benchmark import load_pdb_coordinates
    from dvbfixer.model.diffusion.contract import BackendProvenance, RunnerCandidate
    from dvbfixer.model.diffusion.protein_context import engine_context_view

    check_protocol()
    torch.cuda.init()
    model, reinject = CELLS[cell]
    request, source = request_for(cell, seed, smoke=smoke)
    directory = WORK / ("smoke" if smoke else "runs") / cell / f"seed-{seed}" / attempt
    directory.mkdir(parents=True, exist_ok=False)
    request_path = directory / "request.json"
    request_path.write_text(request.to_json() + "\n")
    (directory / "source.pdb").write_bytes(source)
    write_json(
        directory / "attempt.json",
        {
            "cell": cell,
            "seed": seed,
            "attempt": attempt,
            "smoke_only": smoke,
            "protocol_sha256": digest(WORK / "protocol.json"),
            "request_sha256": digest(request_path),
        },
    )
    module = load_module(
        "cc94_research_lowlevel", ROOT / "deploy/protpardelle-1c/checkpoint_gap_smoke.py"
    )
    config, checkpoint = assets(model)
    module._CONFIG_SHA256, module._CHECKPOINT_SHA256 = PINS[model][:2]
    profile = next(o.value for o in request.backend_options if o.name == "profile")
    module._sampler_profile = lambda *args, **kwargs: profile
    if model == "cc94":
        # Relative encoding retains physical chains; do not apply cc89's rotary
        # neutral-phase intervention. The helper below only admits this frozen case.
        def relative_position_policy(req, device, policy, phase, **kwargs):
            assert req == request and device == "cuda:0"
            assert policy == "source-physical-v1" and phase == "native"

        module._require_position_feature_policy = relative_position_policy

    def candidate(**kwargs):
        kwargs["candidate_id"] = kwargs["candidate_id"].replace("cc89", model)
        kwargs["score_provenance"] = kwargs["score_provenance"].replace("cc89", model)
        return RunnerCandidate(**kwargs)

    def provenance(**kwargs):
        kwargs["backend"] = kwargs["backend"].replace("cc89", model)
        return BackendProvenance(**kwargs)

    module.RunnerCandidate, module.BackendProvenance = candidate, provenance
    trace_builder = module.build_sampler_trace

    def build_trace(**kwargs):
        kwargs["profile"] = profile
        return trace_builder(**kwargs)

    module.build_sampler_trace = build_trace
    validation = module.validate_runner_result

    def persist(req, result, *, workspace, **kwargs):
        (directory / "raw-result.json").write_text(result.to_json() + "\n")
        validated = validation(req, result, workspace=workspace, **kwargs)
        write_json(directory / "raw-validation.json", [asdict(v) for v in validated])
        return validated

    module.validate_runner_result = persist
    mapping = dict(engine_context_view(request).atom_to_source)
    atom_mapper = module._model_atoms

    def capture(req, raw):
        axis, coords, names = atom_mapper(req, raw)
        np.savez_compressed(directory / "native-coordinates.npz", coordinates=coords)
        write_json(
            directory / "native-axis.json",
            {
                "atoms": [asdict(mapping.get(a, a)) for a in axis],
                "engine_atoms": [asdict(a) for a in axis],
                "residue_names": names,
            },
        )
        return axis, coords, names

    module._model_atoms = capture
    callbacks = []
    maker = module._make_per_step_reinjector

    def traced_maker(*args, **kwargs):
        callback, stats = maker(*args, **kwargs)

        def observe(coordinates, index, count, *context):
            projected = callback(coordinates, index, count, *context)
            callbacks.append(
                {
                    "index": index,
                    "count": count,
                    "device": str(coordinates.device),
                    "dtype": str(coordinates.dtype),
                    **dict(stats),
                }
            )
            return projected

        return observe, stats

    module._make_per_step_reinjector = traced_maker
    from protpardelle.core.models import Protpardelle

    original_sample = Protpardelle.sample

    def checked_sample(self, **kwargs):
        prediction = original_sample(self, **kwargs)
        sequence_matches = torch.equal(prediction["s"], kwargs["gt_aatype"])
        write_json(
            directory / "model-evidence.json",
            {
                "task": self.task,
                "position_encoding": self.config.model.struct_model.uvit.position_embedding_type,
                "parameters": sum(p.numel() for p in self.parameters()),
                "sequence_fixed": sequence_matches,
                "token_count": kwargs["gt_aatype"].shape[-1],
                "residue_index": kwargs["residue_index"].cpu().tolist(),
                "chain_index": kwargs["chain_index"].cpu().tolist(),
                "dummy_fill_mode": kwargs["dummy_fill_mode"],
                "step_scale": kwargs["step_scale"],
                "s_churn": kwargs["s_churn"],
                "uniform_steps": kwargs["uniform_steps"],
            },
        )
        if not sequence_matches:
            raise ValueError("sampler changed the fixed sequence")
        return prediction

    Protpardelle.sample = checked_sample
    try:
        (directory / "raw").mkdir()
        steps = 2 if smoke else 500
        result, _ = module._run_staged(
            request_path,
            directory / "raw",
            config,
            checkpoint,
            steps=steps,
            step_scale=1.0,
            s_churn=0.0,
            per_step_reinjection=reinject,
            device_name="cuda:0",
            mps_profile=False,
            chain_phase_policy="native" if model == "cc94" else "neutral-v1",
            position_feature_policy="source-physical-v1",
        )
        assert len(callbacks) == (steps if reinject else 0)
        assert [c["index"] for c in callbacks] == list(range(len(callbacks)))
        assert all(c["device"] == "cuda:0" and c["count"] == steps for c in callbacks)
        original = load_pdb_coordinates(directory, request.normalized_pdb)
        actual = load_pdb_coordinates(directory, result.candidates[0].coordinate_artifact)
        assert set(actual) == set(request.fixed_atoms) | set(request.generated_atoms)
        assert all(np.array_equal(actual[a], original[a]) for a in request.fixed_atoms)
        write_json(
            directory / "integrity.json",
            {
                "status": "passed",
                "fixed_atoms": 3441,
                "generated_atoms": 58,
                "callback_count": len(callbacks),
                "smoke_only": smoke,
            },
        )
    except BaseException as exc:
        write_json(
            directory / "failure.json",
            {"type": type(exc).__name__, "message": str(exc), "traceback": traceback.format_exc()},
        )
        raise
    finally:
        write_json(directory / "callback-evidence.json", callbacks)


def refine(cell, seed, attempt):
    from dvbfixer.model.diffusion.contract import DiffusionRequest, RunnerResult
    from dvbfixer.model.diffusion.runner_refinement import refine_runner_result
    from dvbfixer.model.diffusion.validate import validate_runner_result

    check_protocol()
    directory = WORK / "runs" / cell / f"seed-{seed}" / attempt
    request_path = directory / "request.json"
    request = DiffusionRequest.from_json(request_path.read_text())
    result = RunnerResult.from_json((directory / "raw-result.json").read_text())
    refined = refine_runner_result(
        request_path, result, directory / "refined", platform_name="Reference"
    )
    (directory / "refined-result.json").write_text(refined.to_json() + "\n")
    write_json(
        directory / "refined-validation.json",
        [asdict(v) for v in validate_runner_result(request, refined, workspace=directory)],
    )


def launch(cell, seed, attempt, *, smoke=False, refinement=False):
    args = [
        CORE if refinement else PYTHON,
        "-B",
        __file__,
        "refine" if refinement else "sample",
        "--cell",
        cell,
        "--seed",
        str(seed),
        "--attempt",
        attempt,
    ]
    if smoke:
        args.append("--smoke")
    return command(
        f"{cell}-seed{seed}-{attempt}-"
        + ("refine" if refinement else "smoke" if smoke else "sample"),
        args,
    )


def matrix(stage):
    check_protocol()
    if stage == "smoke":
        jobs = [(c, 7, "attempt-1") for c in CELLS]
    elif stage == "pilot":
        jobs = [(c, 7, "attempt-1") for c in CELLS]
    else:
        if not all(
            (WORK / "runs" / c / "seed-7" / "attempt-1" / "integrity.json").exists() for c in CELLS
        ):
            raise ValueError("pilot integrity not complete")
        jobs = [(c, seed, "attempt-1") for seed in SEEDS[1:] for c in CELLS]
        jobs += [(c, 7, "repeat-1") for c in ["N1", "B1"]]
    for cell, seed, attempt in jobs:
        base = WORK / ("smoke" if stage == "smoke" else "runs") / cell / f"seed-{seed}" / attempt
        if (base / "integrity.json").exists():
            continue
        if not launch(cell, seed, attempt, smoke=stage == "smoke"):
            raise RuntimeError("run failed; inspect " + str(base))


def refine_all():
    jobs = [
        p
        for p in (WORK / "runs").glob("*/*/*")
        if (p / "integrity.json").exists() and not (p / "refined-result.json").exists()
    ]

    def run(directory):
        c, seed, attempt = directory.parts[-3:]
        return launch(c, int(seed[5:]), attempt, refinement=True)

    with ThreadPoolExecutor(max_workers=22) as pool:
        outcomes = list(pool.map(run, jobs))
    write_json(WORK / "refinement-processes.json", outcomes)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "stage",
        choices=[
            "prepare",
            "inventory",
            "smoke",
            "pilot",
            "expand",
            "sample",
            "refine",
            "refine-all",
        ],
    )
    parser.add_argument("--cell", choices=CELLS, default="N1")
    parser.add_argument("--seed", type=int, choices=SEEDS, default=7)
    parser.add_argument("--attempt", default="attempt-1")
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    if args.stage == "prepare":
        prepare()
    elif args.stage == "inventory":
        inventory()
    elif args.stage in {"smoke", "pilot", "expand"}:
        matrix(args.stage)
    elif args.stage == "sample":
        sample(args.cell, args.seed, args.attempt, smoke=args.smoke)
    elif args.stage == "refine":
        refine(args.cell, args.seed, args.attempt)
    else:
        refine_all()


if __name__ == "__main__":
    main()
