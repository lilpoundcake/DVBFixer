"""Backend-neutral conformance records for constrained diffusion samplers."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum


class SamplingAblationMode(StrEnum):
    """Predeclared constrained-sampling ablations for later GPU benchmarks."""

    TEMPLATE_ONLY = "template-conditioning-only"
    REINJECTION = "template-conditioning-plus-reinjection"
    REINJECTION_BOUNDARY_REFINEMENT = "reinjection-plus-boundary-refinement"


@dataclass(frozen=True, slots=True)
class SamplerCapabilities:
    """Evidence-backed controls exposed by one pinned sampler adapter."""

    mutable_state_each_step: bool
    identity_mapping_each_step: bool
    fixed_coordinate_overwrite_each_step: bool
    localized_boundary_refinement: bool

    def supports_per_step_reinjection(self) -> bool:
        return (
            self.mutable_state_each_step
            and self.identity_mapping_each_step
            and self.fixed_coordinate_overwrite_each_step
        )

    def supports_mode(self, mode: SamplingAblationMode) -> bool:
        if mode is SamplingAblationMode.TEMPLATE_ONLY:
            return True
        if mode is SamplingAblationMode.REINJECTION:
            return self.supports_per_step_reinjection()
        if mode is SamplingAblationMode.REINJECTION_BOUNDARY_REFINEMENT:
            return (
                self.supports_per_step_reinjection()
                and self.localized_boundary_refinement
            )
        raise ValueError(f"unknown sampling ablation mode: {mode!r}")


@dataclass(frozen=True, slots=True)
class SamplerConformance:
    """Conformance decision that never infers unsupported denoising control."""

    requested_mode: SamplingAblationMode
    supported: bool
    missing_capabilities: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.supported and self.missing_capabilities:
            raise ValueError("supported sampler conformance cannot name missing capabilities")
        if not self.supported and not self.missing_capabilities:
            raise ValueError("unsupported sampler conformance must name missing capabilities")


def assess_sampler_conformance(
    capabilities: SamplerCapabilities,
    mode: SamplingAblationMode,
) -> SamplerConformance:
    """Require explicit mutable-state and identity control for reinjection claims."""
    missing: list[str] = []
    if mode is not SamplingAblationMode.TEMPLATE_ONLY:
        if not capabilities.mutable_state_each_step:
            missing.append("mutable-state-each-step")
        if not capabilities.identity_mapping_each_step:
            missing.append("identity-mapping-each-step")
        if not capabilities.fixed_coordinate_overwrite_each_step:
            missing.append("fixed-coordinate-overwrite-each-step")
    if (
        mode is SamplingAblationMode.REINJECTION_BOUNDARY_REFINEMENT
        and not capabilities.localized_boundary_refinement
    ):
        missing.append("localized-boundary-refinement")
    return SamplerConformance(
        requested_mode=mode,
        supported=not missing,
        missing_capabilities=tuple(missing),
    )
