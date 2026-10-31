"""Strict, separate completion trace contract; not a schema-4 SamplerTrace."""

from __future__ import annotations

from dataclasses import dataclass, fields
from typing import Any, ClassVar

from dvbfixer.model.diffusion import cc94
from dvbfixer.model.diffusion.contract import (
    ArtifactReference,
    DiffusionContractError,
    _JsonContract,
)


@dataclass(frozen=True, slots=True)
class CompletionSamplerTrace(_JsonContract):
    request: ArtifactReference
    source: ArtifactReference
    motif: ArtifactReference
    reader: ArtifactReference
    native: ArtifactReference
    raw: ArtifactReference
    refined: ArtifactReference
    engine_view: ArtifactReference
    preflight: ArtifactReference
    refinement: ArtifactReference
    progress: ArtifactReference
    adapter_source: ArtifactReference
    seed: int
    device: str
    checkpoint_sha256: str = cc94.CHECKPOINT_SHA256
    config_sha256: str = cc94.CONFIG_SHA256
    engine_revision: str = cc94.ENGINE_REVISION
    ode_updates: int = 500
    step_scale: float = 1.0
    s_churn: float = 0.0
    chain_phase: str = "native"
    fallback_disabled: bool = True
    refinement_mode: str = "generated-only"

    _schema_version: ClassVar[int] = 1

    def __post_init__(self) -> None:
        super(CompletionSamplerTrace, self).__post_init__()
        if type(self.schema_version) is not int:
            raise DiffusionContractError("completion trace schema must be an integer")
        if type(self.seed) is not int or not 0 <= self.seed < 2**32 or self.device not in {"mps", "cuda:0"}:
            raise DiffusionContractError("invalid completion trace seed/device")
        if (self.checkpoint_sha256 != cc94.CHECKPOINT_SHA256 or self.config_sha256 != cc94.CONFIG_SHA256
                or self.engine_revision != cc94.ENGINE_REVISION):
            raise DiffusionContractError("completion trace requires pinned cc94 identity")
        if (type(self.ode_updates) is not int or self.ode_updates != 500 or type(self.step_scale) not in {int, float}
                or type(self.s_churn) not in {int, float} or self.step_scale != 1.0
                or self.s_churn != 0.0 or self.chain_phase != "native" or self.fallback_disabled is not True
                or self.refinement_mode != "generated-only"):
            raise DiffusionContractError("completion trace changed frozen sampling/refinement policy")
        paths = [artifact.path for artifact in self.artifacts().values()]
        if len(set(paths)) != len(paths):
            raise DiffusionContractError("completion trace evidence paths must be unique")

    def artifacts(self) -> dict[str, ArtifactReference]:
        return {field.name: value for field in fields(self)
                if isinstance(value := getattr(self, field.name), ArtifactReference)}


def completion_conditioning(fixed_positions: tuple[int, ...]) -> dict[str, Any]:
    """Frozen native cc94 conditioning shared by execution and offline auditing."""
    return {
        "enabled": True,
        "discontiguous_motif_assignment": {"enabled": True, "strategy": "fixed", "fixed_motif_pos": list(fixed_positions)},
        "num_recurrence_steps": 1,
        "crop_conditional_guidance": {"enabled": True, "start": 0.0, "end": 2.0, "freq": 1,
                                      "freq_start": 0.0, "freq_end": 0.0, "strategy": "backbone-sidechain"},
        "reconstruction_guidance": {"enabled": False, "start": 0.0, "end": 2.0, "schedule": "custom",
                                    "max_scale": 10.0, "loss_weights": {"motif": 1.0}},
        "replacement_guidance": {"enabled": False, "start": 0.0, "end": 0.92},
    }
