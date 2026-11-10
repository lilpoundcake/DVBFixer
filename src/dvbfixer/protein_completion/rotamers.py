"""Offline, bounded traditional BBDEP single-chi reader and periodic lookup.

No data are bundled, fetched or silently substituted. Only SER and the explicit
free-cysteine category CYH are admitted. Multi-chi/density files are not this
format slice. Integrity checks do not independently establish a data license.
"""
from __future__ import annotations

import hashlib
import json
import math
import re
from collections.abc import Mapping
from dataclasses import asdict, dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Literal

LOOKUP_POLICY = "bbdep-single-chi-periodic-mixture-v1"
MAX_LIBRARY_BYTES = 16 * 1024 * 1024
MAX_ROWS = 37 * 37 * 3
GridKey = tuple[int, int]


class RotamerLibraryError(ValueError):
    """Unavailable, incompatible or invalid local authority; never a fallback."""


@dataclass(frozen=True)
class LibraryAuthority:
    library_id: str
    category: Literal["SER", "CYH"]
    source: str
    kind: Literal["external-local", "synthetic-control"]
    license_label: str
    license_evidence_sha256: str

    def __post_init__(self) -> None:
        if self.category not in {"SER", "CYH"} or self.kind not in {"external-local", "synthetic-control"}:
            raise RotamerLibraryError("unsupported residue category or authority kind")
        for text in (self.library_id, self.source, self.license_label):
            if not isinstance(text, str) or not text.strip() or len(text) > 2048:
                raise RotamerLibraryError("bounded explicit authority metadata required")
        if not isinstance(self.license_evidence_sha256, str) or not re.fullmatch(r"[0-9a-f]{64}", self.license_evidence_sha256):
            raise RotamerLibraryError("license-evidence SHA-256 required; metadata is caller-attested")


@dataclass(frozen=True)
class RotamerRow:
    label: int
    probability: float
    mean_degrees: float
    sigma_degrees: float


@dataclass(frozen=True)
class RotamerEvidence:
    model_digest: str
    category: str
    phi_degrees: float
    psi_degrees: float
    label: int
    probability: float
    mean_degrees: float
    sigma_degrees: float
    lookup_policy: str = LOOKUP_POLICY


def periodic_degrees(value: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise RotamerLibraryError("finite torsion required")
    return (value+180) % 360-180


def _grid(value: str) -> int:
    try:
        angle = int(value)
    except ValueError as exc:
        raise RotamerLibraryError("integer backbone grid required") from exc
    if not -180 <= angle <= 180 or angle % 10:
        raise RotamerLibraryError("expected ten-degree backbone grid")
    return angle


@dataclass(frozen=True)
class SingleChiLibrary:
    authority: LibraryAuthority
    data_sha256: str
    cells: Mapping[GridKey, tuple[RotamerRow, ...]]
    view_sha256: str = field(init=False)

    def __post_init__(self) -> None:
        # Copy before freezing so mutation of a caller's mapping cannot modify
        # this scientific view. File loading validates all rows and grid cells.
        if (not isinstance(self.authority, LibraryAuthority) or not isinstance(self.data_sha256, str)
                or not re.fullmatch(r"[0-9a-f]{64}", self.data_sha256)):
            raise RotamerLibraryError("explicit authority and data SHA-256 required")
        if not isinstance(self.cells, Mapping) or not 0 < len(self.cells) <= 36*36:
            raise RotamerLibraryError("bounded nonempty parsed grid required")
        cells = {}
        for key, values in self.cells.items():
            if (not isinstance(key, tuple) or len(key) != 2
                    or any(isinstance(x, bool) or not isinstance(x, int) or not -180 <= x < 180 or x % 10 for x in key)):
                raise RotamerLibraryError("invalid canonical grid key")
            rows = tuple(values)
            if len(rows) != 3 or any(not isinstance(r, RotamerRow) for r in rows):
                raise RotamerLibraryError("three typed single-chi rows required")
            if any(type(r.label) is not int for r in rows) or {r.label for r in rows} != {1, 2, 3}:
                raise RotamerLibraryError("single-chi labels 1/2/3 required")
            for r in rows:
                if (any(isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x)
                        for x in (r.probability, r.mean_degrees, r.sigma_degrees))
                        or not 0 <= r.probability <= 1 or not -180 <= r.mean_degrees < 180
                        or not 0 < r.sigma_degrees <= 180):
                    raise RotamerLibraryError("invalid parsed rotamer parameters")
            if abs(math.fsum(r.probability for r in rows)-1) > 5e-6:
                raise RotamerLibraryError("non-normalized parsed cell")
            cells[key] = tuple(sorted(rows, key=lambda r: r.label))
        if self.authority.kind == "external-local" and len(cells) != 36*36:
            raise RotamerLibraryError("external library requires complete periodic grid")
        object.__setattr__(self, "cells", MappingProxyType(cells))
        canonical = [(key, [asdict(row) for row in cells[key]]) for key in sorted(cells)]
        object.__setattr__(self, "view_sha256", hashlib.sha256(
            json.dumps(canonical, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest())

    @property
    def model_digest(self) -> str:
        value = {"policy": LOOKUP_POLICY, "authority": asdict(self.authority),
                 "data_sha256": self.data_sha256, "view_sha256": self.view_sha256}
        return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()

    def lookup(self, phi: float, psi: float) -> tuple[RotamerEvidence, ...]:
        """Periodic bilinear probabilities and circular mixture moments.

        Means/variance use weights proportional to grid weight * rotamer
        probability. This is our explicitly versioned interpolation, not the
        authors' KDE/bicubic implementation or analytic probability derivatives.
        """
        phi, psi = periodic_degrees(phi), periodic_degrees(psi)
        left, bottom = math.floor(phi/10)*10, math.floor(psi/10)*10
        u, v = (phi-left)/10, (psi-bottom)/10
        corners = ((left, bottom, (1-u)*(1-v)), (left+10, bottom, u*(1-v)),
                   (left, bottom+10, (1-u)*v), (left+10, bottom+10, u*v))
        support = []
        for x, y, weight in corners:
            if weight == 0:
                continue
            key = (int(periodic_degrees(x)), int(periodic_degrees(y)))
            if key not in self.cells:
                raise RotamerLibraryError(f"missing interpolation cell {key}")
            support.append((self.cells[key], weight))
        choices = []
        for label in (1, 2, 3):
            samples = [(next(r for r in rows if r.label == label), weight) for rows, weight in support]
            probability = math.fsum(row.probability*w for row, w in samples)
            if probability == 0:
                continue  # Impossible label, not an epsilon-probability fabrication.
            weights = [(row, w*row.probability/probability) for row, w in samples]
            sine = math.fsum(w*math.sin(math.radians(row.mean_degrees)) for row, w in weights)
            cosine = math.fsum(w*math.cos(math.radians(row.mean_degrees)) for row, w in weights)
            if math.hypot(sine, cosine) < 1e-8:
                raise RotamerLibraryError("ambiguous circular mean")
            mean = periodic_degrees(math.degrees(math.atan2(sine, cosine)))
            variance = math.fsum(w*(row.sigma_degrees**2+periodic_degrees(row.mean_degrees-mean)**2)
                                 for row, w in weights)
            choices.append(RotamerEvidence(self.model_digest, self.authority.category,
                                            phi, psi, label, probability, mean, math.sqrt(variance)))
        return tuple(sorted(choices, key=lambda row: (-row.probability, row.label)))


def parse_single_chi(data: bytes, *, expected_sha256: str,
                     authority: LibraryAuthority) -> SingleChiLibrary:
    """Strict 17-column traditional single-chi format; no silent skipped rows.

    External files require a complete 36x36 periodic grid. Sparse synthetic
    controls remain explicitly labelled and never acquire external provenance.
    "+180" boundary aliases may duplicate "-180" only with identical numbers.
    """
    if not isinstance(data, bytes) or not 0 < len(data) <= MAX_LIBRARY_BYTES:
        raise RotamerLibraryError("empty or oversized rotamer file")
    sha = hashlib.sha256(data).hexdigest()
    if not isinstance(expected_sha256, str) or not re.fullmatch(r"[0-9a-f]{64}", expected_sha256) or sha != expected_sha256:
        raise RotamerLibraryError("rotamer file checksum mismatch")
    try:
        text = data.decode("ascii")
    except UnicodeDecodeError as exc:
        raise RotamerLibraryError("ASCII rotamer table required") from exc
    cells: dict[GridKey, dict[int, RotamerRow]] = {}
    seen_raw = set()
    count = 0
    for number, line in enumerate(text.splitlines(), 1):
        if len(line) > 4096:
            raise RotamerLibraryError(f"oversized rotamer line {number}")
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        count += 1
        if count > MAX_ROWS:
            raise RotamerLibraryError("too many rotamer rows")
        fields = line.split()
        if len(fields) != 17 or fields[0] != authority.category:
            raise RotamerLibraryError(f"wrong column count or residue category on line {number}")
        phi, psi = _grid(fields[1]), _grid(fields[2])
        try:
            observations = int(fields[3])
            labels = tuple(int(x) for x in fields[4:8])
            values = tuple(float(x) for x in fields[8:])
        except ValueError as exc:
            raise RotamerLibraryError(f"invalid numeric row {number}") from exc
        if observations < 0 or labels[0] not in (1, 2, 3) or labels[1:] != (0, 0, 0):
            raise RotamerLibraryError("only explicit single-chi labels 1/2/3 supported")
        if (any(not math.isfinite(x) for x in values) or not 0 <= values[0] <= 1
                or not -180 <= values[1] <= 180 or not 0 < values[5] <= 180
                or any(values[i] != 0 for i in (2, 3, 4, 6, 7, 8))):
            raise RotamerLibraryError("invalid probability, chi mean or sigma")
        raw_key = phi, psi, labels[0]
        if raw_key in seen_raw:
            raise RotamerLibraryError("duplicate rotamer row")
        seen_raw.add(raw_key)
        row = RotamerRow(labels[0], values[0], periodic_degrees(values[1]), values[5])
        key = int(periodic_degrees(phi)), int(periodic_degrees(psi))
        group = cells.setdefault(key, {})
        if row.label in group and group[row.label] != row:
            raise RotamerLibraryError("inconsistent periodic boundary aliases")
        group[row.label] = row
    if not cells:
        raise RotamerLibraryError("no rotamer rows")
    for rows in cells.values():
        if set(rows) != {1, 2, 3} or abs(math.fsum(r.probability for r in rows.values())-1) > 5e-6:
            raise RotamerLibraryError("incomplete labels or non-normalized cell probabilities")
    if authority.kind == "external-local" and set(cells) != {(p, q) for p in range(-180, 180, 10) for q in range(-180, 180, 10)}:
        raise RotamerLibraryError("external library requires complete periodic grid")
    return SingleChiLibrary(authority, sha, {key: tuple(rows[k] for k in sorted(rows)) for key, rows in cells.items()})


def load_single_chi(path: Path, *, expected_sha256: str,
                    authority: LibraryAuthority) -> SingleChiLibrary:
    """Bounded plain local table only; no URL, gzip expansion or network fallback."""
    if path.is_symlink() or not path.is_file():
        raise RotamerLibraryError("rotamer authority must be an explicit regular file")
    try:
        with path.open("rb") as handle:
            data = handle.read(MAX_LIBRARY_BYTES+1)
    except OSError as exc:
        raise RotamerLibraryError(f"rotamer authority unavailable: {exc}") from exc
    return parse_single_chi(data, expected_sha256=expected_sha256, authority=authority)
