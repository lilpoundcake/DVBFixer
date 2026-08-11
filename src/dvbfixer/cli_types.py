"""Reusable argparse value validators for DVBfixer command parsers."""

from __future__ import annotations

import argparse
import math
from collections.abc import Callable


def _numeric_type(
    converter: Callable[[str], int | float],
    valid: Callable[[int | float], bool],
    label: str,
):
    def parse(value: str) -> int | float:
        try:
            parsed = converter(value)
        except ValueError as exc:
            raise argparse.ArgumentTypeError(f"expected {label}, got {value!r}") from exc
        if not valid(parsed):
            raise argparse.ArgumentTypeError(f"expected {label}, got {value!r}")
        return parsed

    parse.__name__ = label.replace(" ", "_")
    parse._dvbfixer_numeric = True  # type: ignore[attr-defined]
    return parse


positive_int = _numeric_type(int, lambda value: value > 0, "a positive integer")
nonnegative_int = _numeric_type(int, lambda value: value >= 0, "a non-negative integer")
positive_float = _numeric_type(
    float, lambda value: math.isfinite(value) and value > 0, "a positive number"
)
nonnegative_float = _numeric_type(
    float, lambda value: math.isfinite(value) and value >= 0, "a non-negative number"
)


def _require_integer(value: str, *, context: str) -> None:
    try:
        int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"{context} residue number must be an integer, got {value!r}"
        ) from exc


def atom_spec(value: str) -> str:
    """Validate ``CHAIN:RESNUM:ATOM`` while preserving its string representation."""
    parts = value.split(":")
    if len(parts) != 3 or not parts[2]:
        raise argparse.ArgumentTypeError(
            f"expected CHAIN:RESNUM:ATOM, got {value!r}"
        )
    _require_integer(parts[1], context="atom selector")
    return value


def disulfide_spec(value: str) -> str:
    """Validate ``CHAIN1:NUM1:CHAIN2:NUM2`` for ``top --ss``."""
    parts = value.split(":")
    if len(parts) != 4:
        raise argparse.ArgumentTypeError(
            f"expected CHAIN1:NUM1:CHAIN2:NUM2, got {value!r}"
        )
    _require_integer(parts[1], context="first disulfide selector")
    _require_integer(parts[3], context="second disulfide selector")
    return value


_HIS_STATES = frozenset({"HIE", "HID", "HIP", "HSE", "HSD", "HSP"})
_PROTONATION_STATES = frozenset({
    "ASH", "ASPP", "ASPH", "GLH", "GLUP", "GLUH",
    "HIP", "HSP", "HISH", "HIE", "HSE", "HISE", "HID", "HSD", "HISD",
    "CYX", "CYM", "LYN", "LSN", "ASP", "GLU", "HIS", "CYS", "LYS",
})


def histidine_spec(value: str) -> str:
    """Validate ``CHAIN:RESNUM:STATE`` for ``top --his``."""
    parts = value.split(":")
    if len(parts) != 3:
        raise argparse.ArgumentTypeError(
            f"expected CHAIN:RESNUM:STATE, got {value!r}"
        )
    _require_integer(parts[1], context="histidine selector")
    if parts[2].upper() not in _HIS_STATES:
        choices = ", ".join(sorted(_HIS_STATES))
        raise argparse.ArgumentTypeError(
            f"unknown histidine state {parts[2]!r}; choose one of: {choices}"
        )
    return value


def protonation_spec(value: str) -> str:
    """Validate ``all`` or comma-separated ``CHAIN:NUM[:STATE]`` selectors."""
    if value == "all":
        return value
    for selector in value.split(","):
        parts = selector.split(":")
        if len(parts) not in (2, 3):
            raise argparse.ArgumentTypeError(
                f"expected all or CHAIN:NUM[:STATE], got {selector!r}"
            )
        _require_integer(parts[1], context="protonation selector")
        if len(parts) == 3 and parts[2].upper() not in _PROTONATION_STATES:
            raise argparse.ArgumentTypeError(
                f"unknown protonation state {parts[2]!r} in {selector!r}"
            )
    return value
