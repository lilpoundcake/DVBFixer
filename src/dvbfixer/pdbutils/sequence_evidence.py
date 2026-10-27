"""Observed PDB peptide geometry as a tie-breaker, never a sequence substitute."""

from __future__ import annotations

import math
from collections.abc import Iterable, Sequence

ResidueKey = tuple[int | str, str]
Point = tuple[float, float, float]


def _key(number: int | str, insertion_code: str) -> tuple[str, str]:
    try:
        normalized = str(int(number))
    except ValueError:
        normalized = str(number).strip()
    return normalized, insertion_code.strip()


def observed_peptide_continuity(
    lines: Iterable[str], chain: str, residues: Sequence[ResidueKey],
) -> tuple[bool | None, ...]:
    """Classify ordered C–N boundaries using unique observed coordinates in Å.

    Missing/duplicate/altloc/invalid coordinates give no evidence. A plausible
    peptide is 1.0–1.8 Å; a separation >=2 Å or explicit TER is a break. This
    does not change empirical chain splitting or infer missing sequence length.
    Chain IDs and insertion codes are case-sensitive.
    """
    keys = tuple(_key(number, code) for number, code in residues)
    wanted = set(keys)
    points: dict[tuple[tuple[str, str], str], list[Point]] = {}
    segments: dict[tuple[str, str], set[int]] = {}
    segment = 0
    last_chain = ""
    for line in lines:
        if line.startswith("TER"):
            if last_chain == chain:
                segment += 1
            continue
        if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 54:
            continue
        last_chain = line[21]
        if last_chain != chain:
            continue
        key = _key(line[22:26].strip(), line[26])
        name = line[12:16].strip()
        if key not in wanted or name not in {"C", "N"}:
            continue
        segments.setdefault(key, set()).add(segment)
        address = (key, name)
        if line[16:17].strip():
            # A single labelled atom still denotes unresolved alternatives.
            points.setdefault(address, []).extend([(math.nan, math.nan, math.nan)])
            continue
        try:
            point = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            point = (math.nan, math.nan, math.nan)
        points.setdefault(address, []).append(point)
    result: list[bool | None] = []
    for left, right in zip(keys, keys[1:]):
        carbons, nitrogens = points.get((left, "C"), []), points.get((right, "N"), [])
        if len(carbons) != 1 or len(nitrogens) != 1:
            result.append(None)
            continue
        distance = math.dist(carbons[0], nitrogens[0])
        if not math.isfinite(distance):
            result.append(None)
        elif segments.get(left, set()).isdisjoint(segments.get(right, set())):
            result.append(False)
        else:
            result.append(True if 1.0 <= distance <= 1.8 else False if distance >= 2.0 else None)
    return tuple(result)
