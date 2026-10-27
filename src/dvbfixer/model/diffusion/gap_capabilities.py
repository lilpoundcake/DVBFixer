"""Profile-specific evidence boundaries, separate from research input capacity."""

from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType


@dataclass(frozen=True)
class GapCapability:
    stratum: str
    accepted: bool
    research_eligible: bool
    reason: str


ACCEPTED_INTERNAL_RANGES = MappingProxyType({
    "protpardelle-1c-cuda": ((5, 10),),
    "protpardelle-1c-mps-local-v1": ((5, 10),),
    "protenix-v1-cuda": ((3, 12),),
    "protpardelle-1c-mps": ((3, 12),),
})


def accepted_gap_ranges(profile: str) -> tuple[tuple[int, int], ...]:
    """Return the immutable evidence range used for extrapolation warnings.

    Unnamed/internal test adapters retain the initial baseline. They cannot
    inherit a future expansion accepted for a named hardware profile.
    """
    return ACCEPTED_INTERNAL_RANGES.get(profile, ((3, 12),))


def accepted_gap_range_label(profile: str) -> str:
    return ", ".join(f"{start}-{stop}" for start, stop in accepted_gap_ranges(profile))


def gap_length_warning(profile: str, length: int, region: str) -> str | None:
    """Describe extrapolation without promoting the profile's evidence boundary."""
    if type(length) is not int or length < 1:
        raise ValueError("gap length must be a positive integer")
    if gap_capability(profile, length).accepted:
        return None
    return (
        f"{region} length {length} is outside the validated range "
        f"{accepted_gap_range_label(profile)} for {profile or 'the selected profile'}; "
        "continuing with experimental reconstruction. Quality is not established "
        "for this length; independent output validation remains mandatory."
    )


def gap_capability(profile: str, length: int) -> GapCapability:
    """Longer inputs require a separately frozen study, never inherited acceptance."""
    if type(length) is not int or length < 1:
        return GapCapability("unsupported", False, False, "invalid-gap-length")
    accepted_ranges = accepted_gap_ranges(profile)
    if any(start <= length <= stop for start, stop in accepted_ranges):
        stratum = "baseline-3-12" if length <= 12 else "research-13-25" if length <= 25 else "research-26-50"
        return GapCapability(stratum, True, True, "accepted-profile-range")
    if 13 <= length <= 25:
        return GapCapability("research-13-25", False, True, "long-gap-stratum-not-promoted")
    if 26 <= length <= 50:
        return GapCapability("research-26-50", False, True, "high-cost-long-gap-stratum-not-promoted")
    return GapCapability("unsupported", False, False, "outside-declared-reconstruction-strata")
