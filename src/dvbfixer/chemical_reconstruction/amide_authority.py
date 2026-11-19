"""Pinned ff19SB amide metrics; no wildcard typing or protonation assignment."""
from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass
from pathlib import Path

from dvbfixer.domain.structure_identity import ExactAtomRef

from .contracts import ContractError, Unit, digest_record
from .metric_geometry import GeometryMetric

FFBONDED_SHA256 = "b6fa511c5bd3f11a19484a418d10e32972e2c19788d464466926c7d8fbde1034"
PINNED_SELECTION_DIGEST = "f16d19b955b40d3e2993bc92bc36b4405ea4997f03991b8058d4671d0c42848b"


def _select(content: str, section: str, types: tuple[str, ...]) -> tuple[str, ...]:
    current = ""
    selected = []
    for line in content.splitlines():
        clean = line.split(";", 1)[0].strip()
        if clean.startswith("[") and clean.endswith("]"):
            current = clean[1:-1].strip()
            continue
        tokens = tuple(clean.split())
        if current == section and tokens[:len(types)] in {types, tuple(reversed(types))}:
            selected.append(tokens)
    if len(selected) != 1:
        raise ContractError("amide-parameter-row-not-unique")
    return selected[0]


@dataclass(frozen=True)
class AmideMetricAuthority:
    source: str
    source_sha256: str
    bond_row: tuple[str, ...]
    carbon_angle_row: tuple[str, ...]
    alpha_angle_row: tuple[str, ...]

    def __post_init__(self) -> None:
        if self.source_sha256 != FFBONDED_SHA256 or not isinstance(self.source, str) or not self.source.strip():
            raise ContractError("pinned-amide-parameter-source-required")
        for row, types in ((self.bond_row, ("H", "N")), (self.carbon_angle_row, ("C", "N", "H")),
                           (self.alpha_angle_row, ("XC", "N", "H"))):
            if (not isinstance(row, tuple) or len(row) != len(types)+3 or row[:len(types)] != types
                    or row[len(types)] != "1"):
                raise ContractError("amide-parameter-type-or-function-unqualified")
            values = tuple(float(t) for t in row[-2:])
            if any(not math.isfinite(v) or v <= 0 for v in values):
                raise ContractError("invalid-amide-parameter")
        if float(self.carbon_angle_row[-1]) != float(self.alpha_angle_row[-1]):
            raise ContractError("amide-angle-stiffness-not-equal")
        if any(float(row[-2]) >= 180 for row in (self.carbon_angle_row, self.alpha_angle_row)):
            raise ContractError("invalid-amide-angle")
        if digest_record((self.bond_row, self.carbon_angle_row, self.alpha_angle_row)) != PINNED_SELECTION_DIGEST:
            raise ContractError("amide-selected-rows-pin-mismatch")

    @property
    def digest(self) -> str:
        return digest_record(("ff19sb-amide-metric-authority-v1", self))

    def metrics(self, nitrogen: ExactAtomRef, hydrogen: ExactAtomRef, alpha: ExactAtomRef,
                 carbon: ExactAtomRef, *, atom_types: tuple[str, str, str]) -> tuple[GeometryMetric, ...]:
        if atom_types != ("N", "H", "XC") or (nitrogen.name, hydrogen.name, alpha.name, carbon.name) != ("N", "H", "CA", "C"):
            raise ContractError("explicit-amide-N-H-XC-types-and-ports-required")
        if nitrogen.component != hydrogen.component or nitrogen.component != alpha.component or carbon.component == nitrogen.component:
            raise ContractError("amide-needs-local-N-H-CA-and-external-C")
        rows = ((self.bond_row, (nitrogen, hydrogen), Unit.ANGSTROM, 10.),
                (self.carbon_angle_row, (carbon, nitrogen, hydrogen), Unit.DEGREE, 1.),
                (self.alpha_angle_row, (alpha, nitrogen, hydrogen), Unit.DEGREE, 1.))
        return tuple(GeometryMetric(atoms, float(row[-2])*scale, unit,
                                   f"{self.source}; exact row {' '.join(row)}; authority {self.digest}", self.source_sha256)
                     for row, atoms, unit, scale in rows)


def load_amide_authority(path: Path, *, expected_sha256: str) -> AmideMetricAuthority:
    if expected_sha256 != FFBONDED_SHA256:
        raise ContractError("pinned-amide-parameter-source-required")
    content = path.read_bytes()
    if hashlib.sha256(content).hexdigest() != expected_sha256:
        raise ContractError("amide-parameter-source-pin-mismatch")
    text = content.decode("ascii")
    return AmideMetricAuthority(str(path), expected_sha256, _select(text, "bondtypes", ("H", "N")),
                                 _select(text, "angletypes", ("C", "N", "H")),
                                 _select(text, "angletypes", ("XC", "N", "H")))
