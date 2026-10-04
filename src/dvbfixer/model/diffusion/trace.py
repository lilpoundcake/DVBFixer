"""Helpers for writing digest-linked sampler trace artifacts."""

from __future__ import annotations

import hashlib
from pathlib import Path

from dvbfixer.model.diffusion.contract import (
    DIFFUSION_SCHEMA_VERSION,
    ArtifactReference,
    AtomIdentity,
    BackendProvenance,
    DiffusionContractError,
    DiffusionRequest,
    RunnerResourceMetrics,
    SamplerStepTrace,
    SamplerTrace,
    atom_identity_digest,
    fixed_mask_digest,
)
from dvbfixer.model.diffusion.sampler import SamplerCapabilities, SamplingAblationMode


def build_sampler_trace(
    *,
    profile: str,
    engine_repository: str,
    engine_revision: str,
    patch_identity: str,
    atom_order: tuple[AtomIdentity, ...],
    fixed_atoms: tuple[AtomIdentity, ...],
    represented_fixed_atoms: tuple[AtomIdentity, ...] | None = None,
    device: str,
    fallback_disabled: bool,
    denoising_update_count: int | None,
    projection_errors_angstrom: tuple[float, ...] = (),
    final_fixed_coordinate_restoration: bool,
    fixed_tolerance_angstrom: float = 0.01,
    resource_metrics: RunnerResourceMetrics | None = None,
    sampler_evidence_complete: bool = True,
) -> SamplerTrace:
    """Build the common no-refinement trace from observed sampler evidence."""
    reinjection = bool(projection_errors_angstrom)
    represented = fixed_atoms if represented_fixed_atoms is None else represented_fixed_atoms
    capabilities = SamplerCapabilities(
        mutable_state_each_step=reinjection and sampler_evidence_complete,
        identity_mapping_each_step=reinjection and sampler_evidence_complete,
        fixed_coordinate_overwrite_each_step=reinjection and sampler_evidence_complete,
        localized_boundary_refinement=sampler_evidence_complete,
    )
    identity_sha256 = atom_identity_digest(atom_order)
    steps = tuple(
        SamplerStepTrace(
            step_index=index,
            atom_identity_sha256=identity_sha256,
            post_projection_max_error_angstrom=error,
        )
        for index, error in enumerate(projection_errors_angstrom)
    )
    return SamplerTrace(
        schema_version=DIFFUSION_SCHEMA_VERSION,
        profile=profile,
        capabilities=capabilities,
        engine_repository=engine_repository,
        engine_revision=engine_revision,
        patch_identity=patch_identity,
        atom_order=atom_order,
        atom_order_sha256=identity_sha256,
        fixed_mask_sha256=fixed_mask_digest(atom_order, represented),
        ablation_mode=(
            SamplingAblationMode.REINJECTION
            if reinjection
            else SamplingAblationMode.TEMPLATE_ONLY
        ),
        fixed_atoms=fixed_atoms,
        represented_fixed_atoms=represented,
        fixed_tolerance_angstrom=fixed_tolerance_angstrom,
        steps=steps,
        callback_update_count=len(steps) if sampler_evidence_complete else None,
        denoising_update_count=denoising_update_count,
        max_projection_error_angstrom=(
            max(projection_errors_angstrom) if projection_errors_angstrom else None
        ),
        sampler_evidence_complete=sampler_evidence_complete,
        final_fixed_coordinate_restoration=final_fixed_coordinate_restoration,
        refinement_mode="none",
        refinement_parameters=(),
        device=device,
        fallback_disabled=fallback_disabled,
        resource_metrics=resource_metrics or RunnerResourceMetrics(),
        final_heavy_coordinate_operations=(
            ("final-fixed-coordinate-restoration",)
            if final_fixed_coordinate_restoration
            else ()
        ),
    )


def write_sampler_trace(
    workspace: Path,
    path: Path,
    trace: SamplerTrace,
) -> ArtifactReference:
    """Write one trace and return its workspace-relative artifact reference."""
    path.write_text(trace.to_json(), encoding="utf-8")
    relative = path.resolve().relative_to(workspace.resolve()).as_posix()
    return ArtifactReference(relative, hashlib.sha256(path.read_bytes()).hexdigest())


def validate_sampler_trace_context(
    trace: SamplerTrace,
    request: DiffusionRequest,
    provenance: BackendProvenance,
) -> None:
    """Match standalone trace evidence to its request and backend identity."""
    if trace.fixed_atoms != request.fixed_atoms:
        raise DiffusionContractError("sampler trace changed the fixed atom mask")
    expected_sampler_atoms = set(trace.represented_fixed_atoms) | set(
        request.generated_atoms
    )
    if set(trace.atom_order) != expected_sampler_atoms:
        raise DiffusionContractError(
            "sampler trace atom order must contain exactly the represented fixed "
            "and requested generated atoms"
        )
    expected_profile = next(
        (option.value for option in request.backend_options if option.name == "profile"),
        None,
    )
    if expected_profile is not None and trace.profile != expected_profile:
        raise DiffusionContractError("sampler trace changed the requested profile")
    if (
        trace.engine_repository != provenance.engine_repository
        or trace.engine_revision != provenance.engine_revision
    ):
        raise DiffusionContractError("sampler trace engine identity mismatch")
    if trace.device != provenance.device:
        raise DiffusionContractError("sampler trace device does not match backend provenance")
