"""Run Protpardelle-1c cc89 on one materialized DVBFixer gap request."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
import tempfile
import time
from collections.abc import Callable
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np

from dvbfixer.model.diffusion.benchmark import candidate_quality, load_pdb_coordinates
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
from dvbfixer.model.diffusion.geometry import synchronize_and_reinject, weighted_kabsch
from dvbfixer.model.diffusion.validate import validate_runner_result

_REVISION = "ee378400f25b801fa481028000f9060183d7fb4c"
_CHECKPOINT_SHA256 = "dfc9895b399ec4497bf6d646502168725f01dcbd55bc0b541fc3e3cc1f2f0483"
_CONFIG_SHA256 = "e9999ace79bf3044351cc982a624fc3459add3a4f91cbf97ff5b437b8942eb9d"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _target_residue_map(request: DiffusionRequest) -> dict[int, ResidueIdentity]:
    if len(request.target_sequences) != 1 or len(request.sequence_placements) != 1:
        raise ValueError("cc89 benchmark requires exactly one target chain")
    target = request.target_sequences[0]
    placement = request.sequence_placements[0]
    if placement.chain != target.chain:
        raise ValueError("target sequence and placement chains differ")
    if placement.target_length != len(target.sequence):
        raise ValueError("sequence placement length differs from target sequence")
    mapping = dict(zip(placement.observed_target_indices, placement.observed_residues))
    for gap in request.gaps:
        for index, residue in zip(
            range(gap.target_interval.start, gap.target_interval.stop),
            gap.generated_residues,
        ):
            if index in mapping:
                raise ValueError(f"target index {index} is mapped more than once")
            mapping[index] = residue
    if set(mapping) != set(range(len(target.sequence))):
        raise ValueError("request does not map every target residue exactly once")
    return mapping


def _read_fixed_coordinates(
    path: Path,
    fixed_atoms: tuple[AtomIdentity, ...],
) -> dict[AtomIdentity, np.ndarray]:
    expected = set(fixed_atoms)
    coordinates: dict[AtomIdentity, np.ndarray] = {}
    for line in path.read_text(encoding="ascii").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or line[16:17] not in {"", " ", "A"}:
            continue
        identity = AtomIdentity(
            line[21:22], line[22:26].strip(), line[26:27].strip(), line[12:16].strip()
        )
        if identity in expected:
            if identity in coordinates:
                raise ValueError(f"duplicate fixed atom: {identity}")
            coordinates[identity] = np.asarray(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                dtype=np.float64,
            )
    missing = expected - set(coordinates)
    if missing:
        raise ValueError(f"normalized PDB omits {len(missing)} fixed atoms")
    return coordinates


def _write_motif_input(
    source: Path,
    destination: Path,
    request: DiffusionRequest,
) -> None:
    placement = request.sequence_placements[0]
    fixed_atom_set = set(request.fixed_atoms)
    target_index = {
        residue: index
        for index, residue in zip(placement.observed_target_indices, placement.observed_residues)
    }
    lines: list[str] = []
    serial = 1
    seen: set[AtomIdentity] = set()
    for line in source.read_text(encoding="ascii").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or line[16:17] not in {"", " ", "A"}:
            continue
        identity = AtomIdentity(
            line[21:22], line[22:26].strip(), line[26:27].strip(), line[12:16].strip()
        )
        residue = ResidueIdentity(identity.chain, identity.residue_number, identity.insertion_code)
        if identity not in fixed_atom_set or residue not in target_index:
            continue
        if identity in seen:
            raise ValueError(f"duplicate motif atom: {identity}")
        seen.add(identity)
        ordinal = target_index[residue] + 1
        record = "ATOM  " + f"{serial:5d}" + line[11:21] + "A" + f"{ordinal:4d}" + " " + line[27:]
        lines.append(record + "\n")
        serial += 1
    if seen != fixed_atom_set:
        raise ValueError(
            "motif input atom set mismatch: "
            f"missing={len(fixed_atom_set - seen)}, extra={len(seen - fixed_atom_set)}"
        )
    lines.append(f"TER   {serial:5d}\nEND\n")
    destination.write_text("".join(lines), encoding="ascii")


def _model_atoms(
    request: DiffusionRequest,
    coordinates: np.ndarray,
) -> tuple[tuple[AtomIdentity, ...], np.ndarray, tuple[str, ...]]:
    from protpardelle.common import residue_constants

    target = request.target_sequences[0]
    residue_map = _target_residue_map(request)
    if coordinates.shape != (len(target.sequence), 37, 3):
        raise ValueError(f"unexpected cc89 coordinate shape: {coordinates.shape}")
    identities: list[AtomIdentity] = []
    values: list[np.ndarray] = []
    residue_names: list[str] = []
    for index, amino_acid in enumerate(target.sequence):
        aatype = residue_constants.restype_order[amino_acid]
        residue = residue_map[index]
        residue_name = residue_constants.restype_1to3[amino_acid]
        for atom_index, present in enumerate(residue_constants.restype_atom37_mask[aatype]):
            if not present:
                continue
            identities.append(
                AtomIdentity(
                    residue.chain,
                    residue.residue_number,
                    residue.insertion_code,
                    residue_constants.atom_types[atom_index],
                )
            )
            values.append(coordinates[index, atom_index])
            residue_names.append(residue_name)
    return tuple(identities), np.asarray(values, dtype=np.float64), tuple(residue_names)


def _synchronize_and_restore_fixed_atoms(
    identities: tuple[AtomIdentity, ...],
    coordinates: np.ndarray,
    residue_names: tuple[str, ...],
    fixed_coordinates: dict[AtomIdentity, np.ndarray],
    fixed_residue_names: dict[AtomIdentity, str],
) -> tuple[tuple[AtomIdentity, ...], np.ndarray, tuple[str, ...]]:
    represented = set(identities)
    alignment_fixed = {
        identity: coordinate
        for identity, coordinate in fixed_coordinates.items()
        if identity in represented
    }
    synchronized = synchronize_and_reinject(identities, coordinates, alignment_fixed)
    source_only = tuple(sorted(set(fixed_coordinates) - represented))
    if not source_only:
        return identities, synchronized, residue_names
    if set(fixed_residue_names) != set(fixed_coordinates):
        raise ValueError("fixed residue names must cover every fixed atom")
    restored_coordinates = np.vstack(
        [synchronized, *(fixed_coordinates[identity] for identity in source_only)]
    )
    return (
        identities + source_only,
        restored_coordinates,
        residue_names + tuple(fixed_residue_names[identity] for identity in source_only),
    )


def _native_conditioning_fit(
    identities: tuple[AtomIdentity, ...],
    coordinates: np.ndarray,
    fixed_coordinates: dict[AtomIdentity, np.ndarray],
) -> dict[str, float | int]:
    index_by_identity = {identity: index for index, identity in enumerate(identities)}
    represented = tuple(sorted(set(fixed_coordinates) & set(index_by_identity)))
    if len(represented) < 3:
        raise ValueError("native conditioning fit requires at least three represented fixed atoms")
    mobile = np.vstack([coordinates[index_by_identity[identity]] for identity in represented])
    target = np.vstack([fixed_coordinates[identity] for identity in represented])
    aligned = weighted_kabsch(mobile, target).apply(mobile)
    displacements = np.linalg.norm(aligned - target, axis=1)
    return {
        "represented_fixed_atom_count": len(represented),
        "source_only_fixed_atom_count": len(fixed_coordinates) - len(represented),
        "rmsd_angstrom": float(np.sqrt(np.mean(np.square(displacements)))),
        "max_displacement_angstrom": float(np.max(displacements)),
    }


def _make_per_step_reinjector(
    request: DiffusionRequest,
    fixed_coordinates: dict[AtomIdentity, np.ndarray],
) -> tuple[Callable[..., Any], dict[str, float | int]]:
    import torch
    from protpardelle.common import residue_constants

    target = request.target_sequences[0]
    target_index_by_residue = {
        residue: index for index, residue in _target_residue_map(request).items()
    }
    represented: list[tuple[int, int, np.ndarray]] = []
    for identity in sorted(fixed_coordinates):
        residue = ResidueIdentity(
            identity.chain,
            identity.residue_number,
            identity.insertion_code,
        )
        target_index = target_index_by_residue.get(residue)
        atom_index = residue_constants.atom_order.get(identity.atom_name)
        if target_index is None or atom_index is None:
            continue
        aatype = residue_constants.restype_order[target.sequence[target_index]]
        if not residue_constants.restype_atom37_mask[aatype][atom_index]:
            continue
        represented.append((target_index, atom_index, fixed_coordinates[identity]))
    if len(represented) < 3:
        raise ValueError("per-step reinjection requires at least three represented fixed atoms")

    target_indices = tuple(item[0] for item in represented)
    atom_indices = tuple(item[1] for item in represented)
    stats: dict[str, float | int] = {
        "callback_count": 0,
        "represented_fixed_atom_count": len(represented),
        "maximum_post_projection_error_angstrom": 0.0,
    }

    def callback(
        state: Any,
        step_index: int,
        step_count: int,
        motif_indices: list[list[int]] | None,
        motif_coordinates: Any,
        motif_mask: Any,
    ) -> Any:
        if state.ndim != 4 or state.shape[0] != 1 or state.shape[-1] != 3:
            raise ValueError("per-step reinjection expects coordinate shape (1, L, A, 3)")
        if step_index != stats["callback_count"] or step_count <= step_index:
            raise ValueError("per-step callback sequence is inconsistent")
        if motif_indices is None or motif_coordinates is None or motif_mask is None:
            raise ValueError("per-step reinjection requires internal motif coordinates")
        if len(motif_indices) != state.shape[0]:
            raise ValueError("motif index batch does not match coordinate batch")
        with torch.no_grad():
            projected = state.clone()
            positions = motif_indices[0]
            coordinates = motif_coordinates[0].to(device=state.device, dtype=state.dtype)
            mask = motif_mask[0].to(device=state.device, dtype=torch.bool)
            if len(positions) != coordinates.shape[0] or mask.shape != coordinates.shape[:-1]:
                raise ValueError("motif callback context has inconsistent shapes")
            motif_row_by_target = {target_index: row for row, target_index in enumerate(positions)}
            try:
                motif_rows = tuple(motif_row_by_target[index] for index in target_indices)
            except KeyError as exc:
                raise ValueError("represented fixed residue is absent from motif placement") from exc
            if not bool(mask[motif_rows, atom_indices].all()):
                raise ValueError("represented fixed atom is absent from internal motif mask")
            authoritative = coordinates[motif_rows, atom_indices]
            projected[0, target_indices, atom_indices] = authoritative
            error = torch.linalg.vector_norm(
                projected[0, target_indices, atom_indices] - authoritative,
                dim=-1,
            ).max()
        stats["callback_count"] = int(stats["callback_count"]) + 1
        stats["maximum_post_projection_error_angstrom"] = max(
            float(stats["maximum_post_projection_error_angstrom"]),
            float(error.item()),
        )
        return projected

    return callback, stats


def _atom_name_field(atom_name: str, element: str) -> str:
    if len(atom_name) == 4 or (atom_name and atom_name[0].isdigit()) or len(element) == 2:
        return f"{atom_name:<4}"
    return f" {atom_name:<3}"


def _write_candidate(
    path: Path,
    identities: tuple[AtomIdentity, ...],
    coordinates: np.ndarray,
    residue_names: tuple[str, ...],
    request: DiffusionRequest,
) -> None:
    included = set(request.fixed_atoms) | set(request.generated_atoms)
    lines: list[str] = []
    seen: set[AtomIdentity] = set()
    serial = 1
    last: tuple[AtomIdentity, str] | None = None
    for identity, xyz, residue_name in zip(identities, coordinates, residue_names):
        if identity not in included:
            continue
        if identity in seen:
            raise ValueError(f"duplicate output atom: {identity}")
        seen.add(identity)
        try:
            residue_number = int(identity.residue_number)
        except ValueError as exc:
            raise ValueError("cc89 benchmark output requires integer residue numbers") from exc
        if not all(math.isfinite(float(value)) and len(f"{float(value):8.3f}") == 8 for value in xyz):
            raise ValueError(f"coordinate cannot be represented in PDB: {identity}")
        element = identity.atom_name[0]
        atom_field = _atom_name_field(identity.atom_name, element)
        lines.append(
            f"ATOM  {serial:5d} {atom_field} {residue_name:>3} "
            f"{identity.chain}{residue_number:4d}{identity.insertion_code or ' ':1}   "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
            f"  1.00  0.00          {element:>2}\n"
        )
        serial += 1
        last = identity, residue_name
    if seen != included or last is None:
        raise ValueError(f"candidate atom set mismatch: missing={len(included - seen)}")
    identity, residue_name = last
    lines.append(
        f"TER   {serial:5d}      {residue_name:>3} {identity.chain}"
        f"{int(identity.residue_number):4d}{identity.insertion_code or ' ':1}\nEND\n"
    )
    path.write_text("".join(lines), encoding="ascii")


def _artifact(workspace: Path, path: Path) -> ArtifactReference:
    return ArtifactReference(path.resolve().relative_to(workspace.resolve()).as_posix(), _sha256(path))


def _run_staged(
    request_path: Path,
    output_dir: Path,
    config_path: Path,
    checkpoint_path: Path,
    *,
    steps: int,
    step_scale: float,
    s_churn: float,
    per_step_reinjection: bool,
) -> Path:
    import torch
    from protpardelle.common import residue_constants
    from protpardelle.core.models import load_model
    from protpardelle.utils import seed_everything

    if steps <= 0:
        raise ValueError("steps must be positive")
    if _sha256(config_path) != _CONFIG_SHA256 or _sha256(checkpoint_path) != _CHECKPOINT_SHA256:
        raise ValueError("cc89 config or checkpoint digest mismatch")
    request_path = request_path.resolve()
    workspace = request_path.parent
    request = DiffusionRequest.from_json(request_path.read_text(encoding="utf-8"))
    if len(request.gaps) != 1 or request.candidate_count != 1 or len(request.seeds) != 1:
        raise ValueError("cc89 benchmark requires one internal gap and one candidate")
    gap = request.gaps[0]
    target = request.target_sequences[0]
    if gap.target_interval.start <= 0 or gap.target_interval.stop >= len(target.sequence):
        raise ValueError("cc89 benchmark requires a two-anchor internal gap")
    if len(target.sequence) > 512:
        raise ValueError("target exceeds the cc89 fixed-size limit")
    if any(amino_acid not in residue_constants.restype_order for amino_acid in target.sequence):
        raise ValueError("target sequence contains a non-canonical residue")

    source = workspace.joinpath(*Path(request.normalized_pdb.path).parts)
    if _sha256(source) != request.normalized_pdb.sha256:
        raise ValueError("normalized PDB digest mismatch")
    fixed_coordinates = _read_fixed_coordinates(source, request.fixed_atoms)
    residue_map = _target_residue_map(request)
    residue_name_by_residue = {
        residue_map[index]: residue_constants.restype_1to3[amino_acid]
        for index, amino_acid in enumerate(target.sequence)
    }
    fixed_residue_names = {
        atom: residue_name_by_residue[
            ResidueIdentity(atom.chain, atom.residue_number, atom.insertion_code)
        ]
        for atom in request.fixed_atoms
    }
    motif_path = output_dir / "motif.pdb"
    _write_motif_input(source, motif_path, request)

    start = time.perf_counter()
    torch.cuda.reset_peak_memory_stats()
    seed = request.seeds[0]
    seed_everything(seed)
    model = load_model(config_path, checkpoint_path, device="cuda")
    seq_mask, residue_index, chain_index = model.make_seq_mask_for_sampling(
        prot_lens_per_chain=torch.tensor([[len(target.sequence)]])
    )
    gt_aatype = torch.tensor(
        [[residue_constants.restype_order[amino_acid] for amino_acid in target.sequence]],
        dtype=torch.long,
        device=model.device,
    )
    fixed_positions = [
        *range(gap.target_interval.start),
        *range(gap.target_interval.stop, len(target.sequence)),
    ]
    placement = (
        f"A1-{gap.target_interval.start}/"
        f"{gap.target_interval.stop - gap.target_interval.start}/"
        f"A{gap.target_interval.stop + 1}-{len(target.sequence)}"
    )
    conditional = {
        "enabled": True,
        "discontiguous_motif_assignment": {
            "enabled": True,
            "strategy": "fixed",
            "fixed_motif_pos": fixed_positions,
        },
        "num_recurrence_steps": 1,
        "crop_conditional_guidance": {
            "enabled": True,
            "start": 0.0,
            "end": 2.0,
            "freq": 1,
            "freq_start": 0.0,
            "freq_end": 0.0,
            "strategy": "backbone-sidechain",
        },
        "reconstruction_guidance": {
            "enabled": False,
            "start": 0.0,
            "end": 2.0,
            "schedule": "custom",
            "max_scale": 10.0,
            "loss_weights": {"motif": 1.0},
        },
        "replacement_guidance": {"enabled": False, "start": 0.0, "end": 0.92},
    }
    callback_stats: dict[str, float | int] = {
        "callback_count": 0,
        "represented_fixed_atom_count": 0,
        "maximum_post_projection_error_angstrom": 0.0,
    }
    sample_options: dict[str, Any] = {}
    if per_step_reinjection:
        callback, callback_stats = _make_per_step_reinjector(request, fixed_coordinates)
        sample_options["step_callback"] = callback
    sampled = model.sample(
        seq_mask=seq_mask,
        residue_index=residue_index,
        chain_index=chain_index,
        gt_aatype=gt_aatype,
        num_steps=steps,
        step_scale=step_scale,
        s_churn=s_churn,
        motif_file_path=str(motif_path),
        motif_placements_full=[placement],
        dummy_fill_mode=model.config.data.dummy_fill_mode,
        partial_diffusion={"enabled": False, "pdb_file_path": None, "num_steps": 100},
        conditional_cfg=conditional,
        sidechain_mode=False,
        skip_mpnn_proportion=1.0,
        use_fullmpnn=False,
        use_fullmpnn_for_final=False,
        jump_steps=False,
        uniform_steps=True,
        **sample_options,
    )
    if per_step_reinjection and callback_stats["callback_count"] != steps:
        raise RuntimeError(
            "per-step callback count mismatch: "
            f"expected {steps}, observed {callback_stats['callback_count']}"
        )
    raw = sampled["x"][0].detach().to(dtype=torch.float64, device="cpu").numpy()
    identities, coordinates, residue_names = _model_atoms(request, raw)
    native_conditioning_fit = _native_conditioning_fit(
        identities,
        coordinates,
        fixed_coordinates,
    )
    identities, coordinates, residue_names = _synchronize_and_restore_fixed_atoms(
        identities,
        coordinates,
        residue_names,
        fixed_coordinates,
        fixed_residue_names,
    )
    candidate_path = output_dir / "candidate.pdb"
    _write_candidate(candidate_path, identities, coordinates, residue_names, request)

    candidate = RunnerCandidate(
        candidate_id=f"protpardelle-1c-cc89-seed-{seed}",
        seed=seed,
        coordinate_artifact=_artifact(workspace, candidate_path),
        generated_atoms=request.generated_atoms,
        generated_residues=gap.generated_residues,
        raw_backend_score=None,
        score_provenance="protpardelle-1c-cc89:unranked",
    )
    elapsed = time.perf_counter() - start
    result = RunnerResult(
        schema_version=DIFFUSION_SCHEMA_VERSION,
        status=DiffusionStatus.SUCCESS,
        candidates=(candidate,),
        runner_diagnostics=RunnerDiagnostics(exit_code=0, timed_out=False),
        backend_provenance=BackendProvenance(
            backend=(
                "protpardelle-1c-cc89-proxy-only+per-step-reinjection"
                if per_step_reinjection
                else "protpardelle-1c-cc89-proxy-only"
            ),
            runner_protocol_version=1,
            engine_repository="https://github.com/ProteinDesignLab/protpardelle-1c",
            engine_revision=_REVISION,
            checkpoint_sha256=_CHECKPOINT_SHA256,
            device=str(model.device),
            precision="float32",
            framework="torch",
            framework_version=torch.__version__,
            cuda_version=str(torch.version.cuda or ""),
            deterministic_algorithms=False,
            deterministic_flags=(
                f"seed={seed}",
                "fixed-sequence",
                "mpnn-disabled",
                f"per-step-reinjection={per_step_reinjection}",
            ),
        ),
        resource_metrics=RunnerResourceMetrics(
            wall_time_seconds=elapsed,
            peak_vram_bytes=torch.cuda.max_memory_allocated(),
        ),
    )
    validation = validate_runner_result(request, result, workspace=workspace)[0]
    reference_path = workspace / "reference.pdb"
    reference = load_pdb_coordinates(workspace, ArtifactReference("reference.pdb", _sha256(reference_path)))
    candidate_coordinates = load_pdb_coordinates(workspace, candidate.coordinate_artifact)
    quality = candidate_quality(
        candidate.candidate_id,
        reference,
        candidate_coordinates,
        generated_residues=set(gap.generated_residues),
        fixed_atoms=set(request.fixed_atoms),
        anchor_residues={gap.left_anchor, gap.right_anchor},
        closure_passed="junction-peptide-connectivity" not in validation.summary.hard_gate_failures,
        validation_passed=validation.summary.passed,
    )
    summary = {
        "candidate_sha256": candidate.coordinate_artifact.sha256,
        "config_sha256": _CONFIG_SHA256,
        "checkpoint_sha256": _CHECKPOINT_SHA256,
        "diffusion_steps": steps,
        "step_scale": step_scale,
        "s_churn": s_churn,
        "motif_placement": placement,
        "native_conditioning_fit": native_conditioning_fit,
        "per_step_reinjection": per_step_reinjection,
        "per_step_callback": callback_stats,
        "quality": {
            "gap_backbone_rmsd_angstrom": quality.gap_backbone_rmsd_angstrom,
            "gap_all_heavy_rmsd_angstrom": quality.gap_all_heavy_rmsd_angstrom,
            "fixed_heavy_rmsd_angstrom": quality.fixed_heavy_rmsd_angstrom,
        },
        "validation": asdict(validation.summary),
        "ranking_metrics": [asdict(metric) for metric in validation.ranking_metrics],
        "wall_time_seconds": elapsed,
        "peak_vram_bytes": result.resource_metrics.peak_vram_bytes,
    }
    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary_path


def run(
    request_path: Path,
    output_dir: Path,
    config_path: Path,
    checkpoint_path: Path,
    *,
    steps: int,
    step_scale: float,
    s_churn: float,
    per_step_reinjection: bool = False,
) -> None:
    request_path = request_path.resolve()
    workspace = request_path.parent
    output_dir = output_dir.resolve()
    try:
        output_dir.relative_to(workspace)
    except ValueError as exc:
        raise ValueError("cc89 benchmark output must be inside the request workspace") from exc
    if output_dir.exists():
        raise FileExistsError(f"output directory already exists: {output_dir}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent))
    try:
        _run_staged(
            request_path,
            staging,
            config_path,
            checkpoint_path,
            steps=steps,
            step_scale=step_scale,
            s_churn=s_churn,
            per_step_reinjection=per_step_reinjection,
        )
        staging.rename(output_dir)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    print(output_dir / "summary.json")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("request", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--steps", default=500, type=int)
    parser.add_argument("--step-scale", default=1.2, type=float)
    parser.add_argument("--s-churn", default=200.0, type=float)
    parser.add_argument("--per-step-reinjection", action="store_true")
    args = parser.parse_args()
    run(
        args.request,
        args.output_dir,
        args.config,
        args.checkpoint,
        steps=args.steps,
        step_scale=args.step_scale,
        s_churn=args.s_churn,
        per_step_reinjection=args.per_step_reinjection,
    )


if __name__ == "__main__":
    main()
