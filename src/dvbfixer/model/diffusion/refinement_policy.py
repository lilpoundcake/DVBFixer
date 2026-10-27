"""Request-bound refinement policy; absence preserves historical semantics.

The version binds the approved movement limits, not mutable runner defaults.
Encoding is explicit so adding support never changes historical request hashes.
Public request builders must enable this only once staged publication is wired.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

from dvbfixer.model.diffusion.anchor_policy import AnchorMovementLimits
from dvbfixer.model.diffusion.contract import BackendOption, DiffusionContractError

POLICY_REVISION = "geometry-gated-anchors-v1"
_PREFIX = "geometry_refinement."
_VERSION = _PREFIX + "policy"
_ANCHORS = _PREFIX + "relax_anchors"
_BFACTORS = _PREFIX + "bfactor_weighting"
_KEYS = frozenset({_VERSION, _ANCHORS, _BFACTORS})


def add_refinement_arguments(group) -> None:
    group.add_argument("--no-relax-anchors", action="store_true",
        help="Disable stage-2 observed-anchor relaxation (diffusion only; enabled by default)")
    group.add_argument("--no-anchor-bfactor-weighting", action="store_true",
        help="Use uniform observed-anchor restraints instead of experimental B weighting (diffusion only)")


@dataclass(frozen=True)
class RefinementPolicy:
    """Explicit staged policy, with stage 2 and B weighting enabled by default."""

    relax_anchors: bool = True
    bfactor_weighting: bool = True

    def __post_init__(self) -> None:
        if type(self.relax_anchors) is not bool or type(self.bfactor_weighting) is not bool:
            raise DiffusionContractError("refinement policy switches must be booleans")

    @property
    def movement_limits(self) -> AnchorMovementLimits:
        # Literal version-bound values: a future change requires a new version.
        return AnchorMovementLimits(
            backbone_rmsd_angstrom_max=0.5,
            ca_displacement_angstrom_max=0.5,
            heavy_displacement_angstrom_max=1.25,
        )

    def backend_options(self) -> tuple[BackendOption, ...]:
        return (
            BackendOption(_VERSION, POLICY_REVISION),
            BackendOption(_ANCHORS, "true" if self.relax_anchors else "false"),
            BackendOption(_BFACTORS, "true" if self.bfactor_weighting else "false"),
        )


def read_refinement_policy(options: Iterable[BackendOption]) -> RefinementPolicy | None:
    """Decode strictly; None means legacy, never implicitly authorizes motion."""
    selected: dict[str, str] = {}
    for option in options:
        if not option.name.startswith(_PREFIX):
            continue
        if option.name not in _KEYS or option.name in selected:
            raise DiffusionContractError("unknown or duplicate refinement policy option")
        selected[option.name] = option.value
    if not selected:
        return None
    if set(selected) != _KEYS or selected[_VERSION] != POLICY_REVISION:
        raise DiffusionContractError("incomplete or unsupported refinement policy")
    if any(selected[key] not in {"true", "false"} for key in (_ANCHORS, _BFACTORS)):
        raise DiffusionContractError("refinement policy switches require true or false")
    return RefinementPolicy(selected[_ANCHORS] == "true", selected[_BFACTORS] == "true")
