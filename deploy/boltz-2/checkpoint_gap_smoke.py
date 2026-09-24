"""Run the patched Boltz-2 checkpoint with exact per-step reinjection."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import time
from dataclasses import asdict
from pathlib import Path

os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
os.environ.setdefault("CUEQ_DEFAULT_CONFIG", "1")
os.environ.setdefault("CUEQ_DISABLE_AOT_TUNING", "1")

import numpy as np
import torch
from atom_mapping import model_atom_axis, validate_feature_axis

from dvbfixer.model.diffusion.contract import (
    DIFFUSION_SCHEMA_VERSION,
    ArtifactReference,
    AtomIdentity,
    BackendProvenance,
    DiffusionRequest,
    DiffusionStatus,
    RunnerCandidate,
    RunnerDiagnostics,
    RunnerResourceMetrics,
    RunnerResult,
)
from dvbfixer.model.diffusion.validate import validate_runner_result

_BOLTZ_REVISION = "b1ebfc46ecf57f5414e0d1a6f9027bbb122c53bc"
_CHECKPOINT_SHA256 = "090e82ac8c92f5e943fa1b39e7410a44027bea7243c0bbb3caa67a77fc1428e1"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_fixed_coordinates(
    pdb_path: Path,
    fixed_atoms: tuple[AtomIdentity, ...],
) -> dict[AtomIdentity, np.ndarray]:
    requested = set(fixed_atoms)
    coordinates: dict[AtomIdentity, np.ndarray] = {}
    for line in pdb_path.read_text(encoding="ascii").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or line[16:17] not in {"", " ", "A"}:
            continue
        identity = AtomIdentity(
            line[21:22], line[22:26].strip(), line[26:27].strip(), line[12:16].strip()
        )
        if identity not in requested:
            continue
        if identity in coordinates:
            raise ValueError(f"duplicate fixed atom in normalized PDB: {identity}")
        coordinates[identity] = np.asarray(
            [float(line[30:38]), float(line[38:46]), float(line[46:54])],
            dtype=np.float64,
        )
    missing = requested - set(coordinates)
    if missing:
        raise ValueError(f"normalized PDB omits {len(missing)} fixed atoms")
    return coordinates


class FixedAtomCallback:
    """Synchronize each Boltz sample frame and restore observed atoms exactly."""

    def __init__(
        self,
        atom_axis: tuple[AtomIdentity, ...],
        padded_atom_count: int,
        fixed_coordinates: dict[AtomIdentity, np.ndarray],
    ) -> None:
        if padded_atom_count < len(atom_axis):
            raise ValueError("padded atom count is smaller than the mapped atom axis")
        if len(fixed_coordinates) < 3:
            raise ValueError("at least three fixed atoms are required")
        index_by_identity = {identity: index for index, identity in enumerate(atom_axis)}
        if set(fixed_coordinates) - set(index_by_identity):
            raise ValueError("fixed coordinate identity is absent from the Boltz axis")
        ordered = tuple(sorted(fixed_coordinates))
        target = np.asarray([fixed_coordinates[identity] for identity in ordered])
        if target.shape != (len(ordered), 3) or not np.isfinite(target).all():
            raise ValueError("fixed coordinates must be finite three-vectors")
        self._padded_atom_count = padded_atom_count
        self._indices = tuple(index_by_identity[identity] for identity in ordered)
        self._target = target
        self.steps: list[dict[str, float | int]] = []

    def __call__(
        self,
        coordinates: torch.Tensor,
        step_index: int,
        step_count: int,
    ) -> torch.Tensor:
        if step_count <= 0 or not 0 <= step_index < step_count:
            raise ValueError("invalid denoising step")
        if coordinates.ndim != 3 or coordinates.shape[-2:] != (self._padded_atom_count, 3):
            raise ValueError("coordinate tensor has an unexpected Boltz atom axis")
        if not coordinates.is_floating_point():
            raise TypeError("coordinate tensor must be floating point")

        indices = torch.as_tensor(self._indices, device=coordinates.device)
        target = torch.as_tensor(
            self._target, dtype=coordinates.dtype, device=coordinates.device
        )
        mobile = coordinates.index_select(-2, indices)
        mobile_center = mobile.mean(dim=-2)
        target_center = target.mean(dim=-2)
        mobile_centered = mobile - mobile_center[..., None, :]
        target_centered = target - target_center
        covariance = torch.einsum("...ki,kj->...ij", mobile_centered, target_centered)
        left, _singular_values, right_transpose = torch.linalg.svd(covariance)
        determinant = torch.linalg.det(
            right_transpose.transpose(-2, -1) @ left.transpose(-2, -1)
        )
        correction = torch.eye(
            3, dtype=coordinates.dtype, device=coordinates.device
        ).expand(*determinant.shape, 3, 3).clone()
        correction[..., 2, 2] = torch.where(
            determinant > 0.0,
            torch.ones_like(determinant),
            -torch.ones_like(determinant),
        )
        rotation = (
            right_transpose.transpose(-2, -1)
            @ correction
            @ left.transpose(-2, -1)
        )
        translation = target_center - torch.einsum(
            "...i,...ji->...j", mobile_center, rotation
        )
        projected = torch.einsum("...ni,...ji->...nj", coordinates, rotation)
        projected = projected + translation[..., None, :]
        projected[..., indices, :] = target
        errors = torch.linalg.vector_norm(
            projected.index_select(-2, indices) - target, dim=-1
        )
        self.steps.append(
            {
                "step_index": step_index,
                "step_count": step_count,
                "post_projection_max_error_angstrom": float(errors.max().item()),
            }
        )
        return projected


def _atom_name_field(atom_name: str, element: str) -> str:
    if len(atom_name) == 4 or (atom_name and atom_name[0].isdigit()) or len(element) == 2:
        return f"{atom_name:<4}"
    return f" {atom_name:<3}"


def _write_candidate(
    output_path: Path,
    coordinates: np.ndarray,
    atom_axis: tuple[AtomIdentity, ...],
    residue_names: tuple[str, ...],
    atomic_numbers: tuple[int, ...],
    request: DiffusionRequest,
) -> None:
    from rdkit import Chem

    if coordinates.shape != (len(atom_axis), 3) or not np.isfinite(coordinates).all():
        raise ValueError("prediction has invalid coordinate shape or non-finite values")
    included = set(request.fixed_atoms) | set(request.generated_atoms)
    periodic_table = Chem.GetPeriodicTable()
    lines: list[str] = []
    last: tuple[AtomIdentity, str] | None = None
    serial = 1
    for identity, residue_name, atomic_number, xyz in zip(
        atom_axis, residue_names, atomic_numbers, coordinates
    ):
        if identity not in included:
            continue
        try:
            residue_number = int(identity.residue_number)
        except ValueError as exc:
            raise ValueError("PDB smoke output requires integer residue numbers") from exc
        if not all(
            math.isfinite(float(value)) and len(f"{float(value):8.3f}") == 8 for value in xyz
        ):
            raise ValueError(
                "coordinate cannot be represented in PDB format: "
                f"{identity} -> {tuple(float(value) for value in xyz)}"
            )
        element = periodic_table.GetElementSymbol(atomic_number).upper()
        atom_field = _atom_name_field(identity.atom_name, element)
        lines.append(
            f"ATOM  {serial:5d} {atom_field} {residue_name:>3} "
            f"{identity.chain}{residue_number:4d}{identity.insertion_code or ' ':1}   "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
            f"  1.00  0.00          {element:>2}\n"
        )
        last = identity, residue_name
        serial += 1
    if last is None or serial - 1 != len(included):
        raise ValueError("candidate output did not contain every requested atom exactly once")
    identity, residue_name = last
    lines.append(
        f"TER   {serial:5d}      {residue_name:>3} {identity.chain}"
        f"{int(identity.residue_number):4d}{identity.insertion_code or ' ':1}\nEND\n"
    )
    output_path.write_text("".join(lines), encoding="ascii")


def _relative_artifact(workspace: Path, path: Path) -> ArtifactReference:
    try:
        relative = path.resolve().relative_to(workspace.resolve())
    except ValueError as exc:
        raise ValueError("candidate output must be inside the request workspace") from exc
    return ArtifactReference(relative.as_posix(), _sha256(path))


def run(
    request_path: Path,
    input_yaml: Path,
    output_dir: Path,
    cache_dir: Path,
    checkpoint: Path,
    *,
    recycling_steps: int,
    sampling_steps: int,
) -> None:
    from boltz.data.module.inferencev2 import (
        Boltz2InferenceDataModule,
        PredictionDataset,
        load_input,
    )
    from boltz.data.tokenize.boltz2 import Boltz2Tokenizer
    from boltz.data.types import Manifest
    from boltz.main import (
        Boltz2DiffusionParams,
        BoltzSteeringParams,
        MSAModuleArgs,
        PairformerArgsV2,
        process_inputs,
    )
    from boltz.model.models.boltz2 import Boltz2
    from pytorch_lightning import Trainer, seed_everything
    from rdkit import Chem

    if recycling_steps <= 0 or sampling_steps <= 0:
        raise ValueError("recycling_steps and sampling_steps must be positive")
    request_path = request_path.resolve()
    workspace = request_path.parent
    request = DiffusionRequest.from_json(request_path.read_text(encoding="utf-8"))
    if request.candidate_count != 1 or len(request.seeds) != 1:
        raise ValueError("Boltz checkpoint smoke requires one candidate and one seed")
    if _sha256(checkpoint) != _CHECKPOINT_SHA256:
        raise ValueError("Boltz checkpoint digest mismatch")
    source_path = workspace.joinpath(*Path(request.normalized_pdb.path).parts)
    if _sha256(source_path) != request.normalized_pdb.sha256:
        raise ValueError("normalized PDB digest mismatch")
    input_yaml = input_yaml.resolve()
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=False)

    seed = request.seeds[0]
    seed_everything(seed, workers=True)
    torch.use_deterministic_algorithms(True)
    torch.set_float32_matmul_precision("highest")
    torch.set_grad_enabled(False)
    Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)

    start = time.perf_counter()
    processed_root = output_dir / "boltz"
    process_inputs(
        data=[input_yaml],
        out_dir=processed_root,
        ccd_path=cache_dir / "ccd.pkl",
        mol_dir=cache_dir / "mols",
        msa_server_url="https://api.colabfold.com",
        msa_pairing_strategy="greedy",
        boltz2=True,
        preprocessing_threads=1,
    )
    processed = processed_root / "processed"
    manifest = Manifest.load(processed / "manifest.json")
    if len(manifest.records) != 1:
        raise RuntimeError("Boltz preprocessing did not produce exactly one record")
    record = manifest.records[0]
    chain_name_by_asym_id = {chain.chain_id: chain.chain_name for chain in record.chains}
    if len(chain_name_by_asym_id) != len(record.chains):
        raise ValueError("Boltz record has duplicate asym IDs")

    target_dir = processed / "structures"
    msa_dir = processed / "msa"
    constraints_dir = processed / "constraints"
    template_dir = processed / "templates"
    extra_mols_dir = processed / "mols"
    input_data = load_input(
        record,
        target_dir,
        msa_dir,
        constraints_dir,
        template_dir,
        extra_mols_dir,
    )
    tokenized = Boltz2Tokenizer().tokenize(input_data)
    atom_axis, residue_names, token_indices = model_atom_axis(
        tokenized,
        request,
        chain_name_by_asym_id=chain_name_by_asym_id,
    )
    preflight_dataset = PredictionDataset(
        manifest,
        target_dir,
        msa_dir,
        cache_dir / "mols",
        constraints_dir,
        template_dir,
        extra_mols_dir,
    )
    preflight_features = preflight_dataset[0]
    atomic_numbers = validate_feature_axis(atom_axis, token_indices, preflight_features)
    padded_atom_count = int(preflight_features["atom_pad_mask"].shape[0])
    fixed_coordinates = _read_fixed_coordinates(source_path, request.fixed_atoms)
    callback = FixedAtomCallback(atom_axis, padded_atom_count, fixed_coordinates)

    data_module = Boltz2InferenceDataModule(
        manifest=manifest,
        target_dir=target_dir,
        msa_dir=msa_dir,
        mol_dir=cache_dir / "mols",
        num_workers=0,
        constraints_dir=constraints_dir,
        template_dir=template_dir,
        extra_mols_dir=extra_mols_dir,
    )
    diffusion_params = Boltz2DiffusionParams()
    diffusion_params.step_scale = 1.5
    steering_args = BoltzSteeringParams()
    steering_args.fk_steering = False
    steering_args.physical_guidance_update = False
    steering_args.contact_guidance_update = False
    predict_args = {
        "recycling_steps": recycling_steps,
        "sampling_steps": sampling_steps,
        "diffusion_samples": 1,
        "max_parallel_samples": 1,
        "write_confidence_summary": True,
        "write_full_pae": False,
        "write_full_pde": False,
    }
    model_load_start = time.perf_counter()
    model_module = Boltz2.load_from_checkpoint(
        checkpoint,
        strict=True,
        predict_args=predict_args,
        map_location="cpu",
        diffusion_process_args=asdict(diffusion_params),
        ema=False,
        use_kernels=False,
        pairformer_args=asdict(PairformerArgsV2()),
        msa_args=asdict(MSAModuleArgs(use_paired_feature=True)),
        steering_args=asdict(steering_args),
    )
    model_load_seconds = time.perf_counter() - model_load_start
    model_module.dvbfixer_step_callback = callback
    model_module.eval()

    seed_everything(seed, workers=True)
    trainer = Trainer(
        default_root_dir=output_dir,
        accelerator="gpu",
        devices=1,
        precision="bf16-mixed",
        deterministic=True,
        logger=False,
        enable_checkpointing=False,
    )
    predictions = trainer.predict(
        model_module,
        datamodule=data_module,
        return_predictions=True,
    )
    if len(predictions) != 1 or predictions[0].get("exception"):
        raise RuntimeError("Boltz prediction failed or returned an unexpected batch count")
    prediction = predictions[0]
    coordinates = prediction["coords"]
    if coordinates.ndim != 3 or coordinates.shape[0] != 1:
        raise ValueError("Boltz checkpoint smoke expected exactly one coordinate sample")
    if len(callback.steps) != sampling_steps or any(
        step["step_count"] != sampling_steps for step in callback.steps
    ):
        raise RuntimeError("callback count does not match requested sampling steps")
    final = coordinates[0, : len(atom_axis)].detach().to(dtype=torch.float64, device="cpu").numpy()
    candidate_path = output_dir / "candidate.pdb"
    _write_candidate(
        candidate_path,
        final,
        atom_axis,
        residue_names,
        atomic_numbers,
        request,
    )

    candidate = RunnerCandidate(
        candidate_id=f"boltz-2-reinjection-seed-{seed}",
        seed=seed,
        coordinate_artifact=_relative_artifact(workspace, candidate_path),
        generated_atoms=request.generated_atoms,
        generated_residues=tuple(
            residue for gap in request.gaps for residue in gap.generated_residues
        ),
        raw_backend_score=None,
        score_provenance="boltz-2:unranked-checkpoint-smoke",
    )
    runner_result = RunnerResult(
        schema_version=DIFFUSION_SCHEMA_VERSION,
        status=DiffusionStatus.SUCCESS,
        candidates=(candidate,),
        runner_diagnostics=RunnerDiagnostics(exit_code=0, timed_out=False),
        backend_provenance=BackendProvenance(
            backend="boltz-2-hook-spike",
            runner_protocol_version=3,
            engine_repository="https://github.com/jwohlwend/boltz",
            engine_revision=_BOLTZ_REVISION,
            checkpoint_sha256=_CHECKPOINT_SHA256,
            device=torch.cuda.get_device_name(0),
            precision="bf16-mixed",
            framework="torch",
            framework_version=torch.__version__,
            cuda_version=str(torch.version.cuda or ""),
            deterministic_algorithms=torch.are_deterministic_algorithms_enabled(),
            deterministic_flags=("CUBLAS_WORKSPACE_CONFIG=:4096:8",),
        ),
        resource_metrics=RunnerResourceMetrics(
            wall_time_seconds=time.perf_counter() - start,
            model_load_seconds=model_load_seconds,
            peak_vram_bytes=torch.cuda.max_memory_allocated(),
        ),
    )
    validation = validate_runner_result(request, runner_result, workspace=workspace)[0]
    summary = {
        "candidate_sha256": candidate.coordinate_artifact.sha256,
        "seed": seed,
        "model_atom_count": len(atom_axis),
        "padded_atom_count": padded_atom_count,
        "written_atom_count": len(request.fixed_atoms) + len(request.generated_atoms),
        "recycling_steps": recycling_steps,
        "sampling_steps": sampling_steps,
        "callback_count": len(callback.steps),
        "callback_steps": callback.steps,
        "deterministic_algorithms": torch.are_deterministic_algorithms_enabled(),
        "validation": asdict(validation.summary),
        "ranking_metrics": [asdict(metric) for metric in validation.ranking_metrics],
        "model_load_seconds": model_load_seconds,
        "wall_time_seconds": runner_result.resource_metrics.wall_time_seconds,
        "peak_vram_bytes": runner_result.resource_metrics.peak_vram_bytes,
    }
    summary_path = output_dir / "summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(summary_path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("request", type=Path)
    parser.add_argument("input_yaml", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--recycling-steps", type=int, default=1)
    parser.add_argument("--sampling-steps", type=int, default=2)
    args = parser.parse_args()
    run(
        args.request,
        args.input_yaml,
        args.output_dir,
        args.cache,
        args.checkpoint,
        recycling_steps=args.recycling_steps,
        sampling_steps=args.sampling_steps,
    )


if __name__ == "__main__":
    main()
