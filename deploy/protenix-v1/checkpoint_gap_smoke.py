"""Run the patched Protenix checkpoint with exact per-step reinjection."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import time
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np
import torch
from reinjection import FixedAtomReinjector

from dvbfixer.model.diffusion.contract import (
    DIFFUSION_SCHEMA_VERSION,
    ArtifactReference,
    AtomIdentity,
    BackendProvenance,
    DiffusionRequest,
    DiffusionStatus,
    ResidueIdentity,
    RunnerCandidate,
    RunnerDiagnostics,
    RunnerResourceMetrics,
    RunnerResult,
)
from dvbfixer.model.diffusion.geometry import weighted_kabsch
from dvbfixer.model.diffusion.sampler import SamplingAblationMode
from dvbfixer.model.diffusion.validate import validate_runner_result

_ONE_TO_THREE = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
    "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
    "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
    "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _target_residue_map(request: DiffusionRequest) -> dict[int, ResidueIdentity]:
    if len(request.target_sequences) != 1 or len(request.sequence_placements) != 1:
        raise ValueError("checkpoint smoke requires exactly one target chain")
    target = request.target_sequences[0]
    placement = request.sequence_placements[0]
    mapping = dict(zip(placement.observed_target_indices, placement.observed_residues))
    for gap in request.gaps:
        if gap.chain != target.chain:
            raise ValueError("all gaps must belong to the single target chain")
        for index, residue in zip(
            range(gap.target_interval.start, gap.target_interval.stop),
            gap.generated_residues,
        ):
            if index in mapping:
                raise ValueError(f"target residue index {index} is mapped more than once")
            mapping[index] = residue
    expected = set(range(len(target.sequence)))
    if set(mapping) != expected:
        raise ValueError("request does not map every target sequence residue exactly once")
    return mapping


def _model_atom_axis(
    atom_array: Any,
    request: DiffusionRequest,
) -> tuple[tuple[AtomIdentity, ...], tuple[str, ...], tuple[str, ...]]:
    target = request.target_sequences[0]
    residue_map = _target_residue_map(request)
    chain_ids = {str(value) for value in atom_array.chain_id}
    if len(chain_ids) != 1:
        raise ValueError("featurized atom axis must contain exactly one chain")

    identities: list[AtomIdentity] = []
    residue_names: list[str] = []
    elements: list[str] = []
    for model_res_id, atom_name, residue_name, element in zip(
        atom_array.res_id,
        atom_array.atom_name,
        atom_array.res_name,
        atom_array.element,
    ):
        target_index = int(model_res_id) - 1
        if target_index not in residue_map:
            raise ValueError(f"model residue ordinal {model_res_id} is outside the target")
        expected_name = _ONE_TO_THREE[target.sequence[target_index]]
        if str(residue_name) != expected_name:
            raise ValueError(
                f"model residue {model_res_id} is {residue_name}, expected {expected_name}"
            )
        residue = residue_map[target_index]
        identities.append(
            AtomIdentity(
                residue.chain,
                residue.residue_number,
                residue.insertion_code,
                str(atom_name),
            )
        )
        residue_names.append(str(residue_name))
        elements.append(str(element).strip().upper())
    axis = tuple(identities)
    if len(set(axis)) != len(axis):
        raise ValueError("featurized atom axis maps to duplicate DVBFixer identities")
    required = set(request.fixed_atoms) | set(request.generated_atoms)
    missing = required - set(axis)
    if missing:
        raise ValueError(f"model atom axis omits {len(missing)} requested atoms")
    return axis, tuple(residue_names), tuple(elements)


def _read_fixed_coordinates(
    pdb_path: Path,
    fixed_atoms: tuple[AtomIdentity, ...],
) -> dict[AtomIdentity, np.ndarray]:
    requested = set(fixed_atoms)
    coordinates: dict[AtomIdentity, np.ndarray] = {}
    for line in pdb_path.read_text(encoding="ascii").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        altloc = line[16:17]
        if altloc not in {"", " ", "A"}:
            continue
        identity = AtomIdentity(
            line[21:22],
            line[22:26].strip(),
            line[26:27].strip(),
            line[12:16].strip(),
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


class _TracingCallback:
    def __init__(
        self,
        atom_axis: tuple[AtomIdentity, ...],
        fixed_coordinates: dict[AtomIdentity, np.ndarray],
    ) -> None:
        self._reinject = FixedAtomReinjector(atom_axis, fixed_coordinates)
        index_by_identity = {identity: index for index, identity in enumerate(atom_axis)}
        ordered = tuple(sorted(fixed_coordinates))
        self._indices = tuple(index_by_identity[identity] for identity in ordered)
        self._target = np.asarray([fixed_coordinates[identity] for identity in ordered])
        self.steps: list[dict[str, float | int]] = []

    def __call__(
        self, coordinates: torch.Tensor, step_index: int, step_count: int
    ) -> torch.Tensor:
        projected = self._reinject(coordinates, step_index, step_count)
        indices = torch.as_tensor(self._indices, device=projected.device)
        target = torch.as_tensor(
            self._target, dtype=projected.dtype, device=projected.device
        )
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
    if len(atom_name) == 4 or (atom_name and atom_name[0].isdigit()):
        return f"{atom_name:<4}"
    if len(element) == 2:
        return f"{atom_name:<4}"
    return f" {atom_name:<3}"


def _write_candidate(
    output_path: Path,
    coordinates: np.ndarray,
    atom_axis: tuple[AtomIdentity, ...],
    residue_names: tuple[str, ...],
    elements: tuple[str, ...],
    request: DiffusionRequest,
) -> None:
    if coordinates.shape != (len(atom_axis), 3) or not np.isfinite(coordinates).all():
        raise ValueError("prediction has invalid coordinate shape or non-finite values")
    included = set(request.fixed_atoms) | set(request.generated_atoms)
    lines: list[str] = []
    last: tuple[AtomIdentity, str] | None = None
    serial = 1
    for identity, residue_name, element, xyz in zip(
        atom_axis, residue_names, elements, coordinates
    ):
        if identity not in included:
            continue
        try:
            residue_number = int(identity.residue_number)
        except ValueError as exc:
            raise ValueError("PDB smoke output requires integer residue numbers") from exc
        if not all(
            math.isfinite(float(value)) and len(f"{float(value):8.3f}") == 8
            for value in xyz
        ):
            raise ValueError("coordinate cannot be represented in PDB format")
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


def _synchronize_frame(
    coordinates: np.ndarray,
    atom_axis: tuple[AtomIdentity, ...],
    fixed_coordinates: dict[AtomIdentity, np.ndarray],
) -> np.ndarray:
    index_by_identity = {identity: index for index, identity in enumerate(atom_axis)}
    ordered = tuple(sorted(fixed_coordinates))
    indices = [index_by_identity[identity] for identity in ordered]
    target = np.asarray([fixed_coordinates[identity] for identity in ordered])
    transform = weighted_kabsch(coordinates[indices], target)
    return transform.apply(coordinates)


def _relative_artifact(workspace: Path, path: Path) -> ArtifactReference:
    try:
        relative = path.resolve().relative_to(workspace.resolve())
    except ValueError as exc:
        raise ValueError("candidate output must be inside the request workspace") from exc
    return ArtifactReference(relative.as_posix(), _sha256(path))


def run(
    request_path: Path,
    input_json: Path,
    output_dir: Path,
    kalign: Path,
    *,
    cycles: int,
    steps: int,
    ablation_mode: SamplingAblationMode,
) -> None:
    os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
    from protenix.data.inference.infer_dataloader import get_inference_dataloader
    from protenix.utils.seed import seed_everything
    from runner.batch_inference import get_default_runner
    from runner.inference import update_inference_configs

    request_path = request_path.resolve()
    workspace = request_path.parent
    request = DiffusionRequest.from_json(request_path.read_text(encoding="utf-8"))
    if request.candidate_count != 1 or len(request.seeds) != 1:
        raise ValueError("checkpoint smoke requires exactly one candidate and seed")
    if cycles <= 0 or steps <= 0:
        raise ValueError("cycles and steps must be positive")
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=False)

    start = time.perf_counter()
    runner = get_default_runner(
        seeds=list(request.seeds),
        n_cycle=cycles,
        n_step=steps,
        n_sample=1,
        dtype="bf16",
        use_msa=False,
        trimul_kernel="torch",
        triatt_kernel="torch",
        enable_cache=False,
        enable_fusion=False,
        use_template=True,
        kalign_binary_path=str(kalign.resolve()),
    )
    runner.configs.deterministic = True
    runner.configs.deterministic_seed = True
    runner.configs.input_json_path = str(input_json.resolve())
    seed = request.seeds[0]
    seed_everything(seed=seed, deterministic=runner.configs.deterministic)
    batch = next(iter(get_inference_dataloader(configs=runner.configs)))
    data, atom_array, error_message = batch[0]
    if error_message:
        raise RuntimeError(error_message)
    atom_axis, residue_names, elements = _model_atom_axis(atom_array, request)
    source_path = workspace.joinpath(*Path(request.normalized_pdb.path).parts)
    if _sha256(source_path) != request.normalized_pdb.sha256:
        raise ValueError("normalized PDB digest mismatch")
    fixed_coordinates = _read_fixed_coordinates(source_path, request.fixed_atoms)
    callback = None
    if ablation_mode is SamplingAblationMode.REINJECTION:
        callback = _TracingCallback(atom_axis, fixed_coordinates)
        runner.model.dvbfixer_step_callback = callback

    runner.update_model_configs(
        update_inference_configs(runner.configs, data["N_token"].item())
    )
    prediction = runner.predict(data)
    coordinates = prediction["coordinate"]
    if coordinates.shape[0] != 1:
        raise ValueError("checkpoint smoke expected exactly one sample")
    callback_steps = [] if callback is None else callback.steps
    expected_callback_count = (
        0 if ablation_mode is SamplingAblationMode.TEMPLATE_ONLY else steps
    )
    if len(callback_steps) != expected_callback_count or any(
        step["step_count"] != steps for step in callback_steps
    ):
        raise RuntimeError("callback count does not match the requested ablation")

    final = coordinates[0].detach().to(dtype=torch.float64, device="cpu").numpy()
    final_frame_synchronization = (
        ablation_mode is SamplingAblationMode.TEMPLATE_ONLY
    )
    if final_frame_synchronization:
        final = _synchronize_frame(final, atom_axis, fixed_coordinates)
    output_path = output_dir / "candidate.pdb"
    _write_candidate(output_path, final, atom_axis, residue_names, elements, request)
    candidate = RunnerCandidate(
        candidate_id=f"protenix-v1-{ablation_mode.value}-seed-{seed}",
        seed=seed,
        coordinate_artifact=_relative_artifact(workspace, output_path),
        generated_atoms=request.generated_atoms,
        generated_residues=tuple(
            residue for gap in request.gaps for residue in gap.generated_residues
        ),
        raw_backend_score=None,
        score_provenance="protenix-v1:unranked-checkpoint-smoke",
    )
    runner_result = RunnerResult(
        schema_version=DIFFUSION_SCHEMA_VERSION,
        status=DiffusionStatus.SUCCESS,
        candidates=(candidate,),
        runner_diagnostics=RunnerDiagnostics(exit_code=0, timed_out=False),
        backend_provenance=BackendProvenance(
            backend="protenix-v1-hook-spike",
            runner_protocol_version=3,
            engine_repository="https://github.com/bytedance/Protenix",
            engine_revision="85767b811c40ed46e73a9b39519cf6bfca8701ba",
            checkpoint_sha256=(
                "2b7d5a8b30494514fc47fd2271a16260528cdba170ba09cc112fdecd8f85ec04"
            ),
            device=str(runner.device),
            precision=runner.configs.dtype,
            framework="torch",
            framework_version=torch.__version__,
            cuda_version=str(torch.version.cuda or ""),
            deterministic_algorithms=bool(runner.configs.deterministic),
        ),
        resource_metrics=RunnerResourceMetrics(wall_time_seconds=time.perf_counter() - start),
    )
    validation = validate_runner_result(request, runner_result, workspace=workspace)[0]
    summary = {
        "candidate_sha256": candidate.coordinate_artifact.sha256,
        "ablation_mode": ablation_mode.value,
        "model_atom_count": len(atom_axis),
        "written_atom_count": len(request.fixed_atoms) + len(request.generated_atoms),
        "cycles": cycles,
        "diffusion_steps": steps,
        "deterministic_algorithms": torch.are_deterministic_algorithms_enabled(),
        "callback_count": len(callback_steps),
        "callback_steps": callback_steps,
        "final_frame_synchronization": final_frame_synchronization,
        "validation": asdict(validation.summary),
        "ranking_metrics": [asdict(metric) for metric in validation.ranking_metrics],
        "wall_time_seconds": runner_result.resource_metrics.wall_time_seconds,
    }
    summary_path = output_dir / "summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(summary_path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("request", type=Path)
    parser.add_argument("input_json", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--kalign", required=True, type=Path)
    parser.add_argument("--cycles", default=1, type=int)
    parser.add_argument("--steps", default=2, type=int)
    parser.add_argument(
        "--ablation-mode",
        choices=(
            SamplingAblationMode.TEMPLATE_ONLY.value,
            SamplingAblationMode.REINJECTION.value,
        ),
        default=SamplingAblationMode.REINJECTION.value,
    )
    args = parser.parse_args()
    run(
        args.request,
        args.input_json,
        args.output_dir,
        args.kalign,
        cycles=args.cycles,
        steps=args.steps,
        ablation_mode=SamplingAblationMode(args.ablation_mode),
    )


if __name__ == "__main__":
    main()
