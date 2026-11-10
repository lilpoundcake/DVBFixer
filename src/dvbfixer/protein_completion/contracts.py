"""Immutable observation-bound inputs and generated-only private results."""
from __future__ import annotations

import hashlib
import json
import math
import re
from dataclasses import asdict, dataclass
from typing import Literal

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.protein_completion.rotamers import RotamerEvidence

Point = tuple[float, float, float]
REVISION = "canonical-sidechain-internal-v1"
MAX_CONTEXT_ATOMS = 50_000


@dataclass(frozen=True)
class Observation:
    atom: ExactAtomRef
    element: str
    xyz: Point  # angstrom
    occupancy: float = 1.0

    def __post_init__(self) -> None:
        if not isinstance(self.atom, ExactAtomRef):
            raise ValueError("exact atom identity required")
        if (not isinstance(self.xyz, tuple) or len(self.xyz) != 3
                or any(isinstance(x, bool) or not isinstance(x, (float, int))
                       or not math.isfinite(x) or abs(x) > 1e5 for x in self.xyz)):
            raise ValueError("finite bounded immutable angstrom coordinates required")
        if (isinstance(self.occupancy, bool) or not isinstance(self.occupancy, (float, int))
                or not math.isfinite(self.occupancy) or not 0 < self.occupancy <= 1):
            raise ValueError("occupancy must be in (0, 1]")
        if not isinstance(self.element, str) or not re.fullmatch(r"[A-Z][A-Z]?", self.element):
            raise ValueError("explicit canonical element symbol required")


@dataclass(frozen=True)
class CompletionRequest:
    source_sha256: str  # caller's original bytes; observation binding is separate
    residue: ComponentInstanceRef
    residue_name: str
    observations: tuple[Observation, ...]  # all retained context, including target
    missing_atoms: tuple[str, ...]  # exact requested set, not an instruction to delete
    links: tuple[tuple[ExactAtomRef, ExactAtomRef], ...] = ()  # explicit external links

    def __post_init__(self) -> None:
        if not isinstance(self.source_sha256, str) or not re.fullmatch(r"[0-9a-f]{64}", self.source_sha256):
            raise ValueError("source SHA-256 required")
        if not isinstance(self.residue, ComponentInstanceRef):
            raise ValueError("exact component identity required")
        if not isinstance(self.residue_name, str) or not re.fullmatch(r"[A-Z0-9]{1,8}", self.residue_name):
            raise ValueError("explicit residue name required")
        if (not isinstance(self.observations, tuple) or not self.observations
                or len(self.observations) > MAX_CONTEXT_ATOMS
                or any(not isinstance(o, Observation) for o in self.observations)):
            raise ValueError("bounded immutable observations required")
        atoms = [o.atom for o in self.observations]
        if len(set(atoms)) != len(atoms):
            raise ValueError("duplicate observed atom identity")
        if (not isinstance(self.missing_atoms, tuple) or len(self.missing_atoms) > 32
                or any(not isinstance(n, str) or not re.fullmatch(r"[A-Z][A-Z0-9]{0,3}", n)
                       for n in self.missing_atoms)
                or len(set(self.missing_atoms)) != len(self.missing_atoms)):
            raise ValueError("unique immutable missing-atom names required")
        if not isinstance(self.links, tuple) or len(self.links) > MAX_CONTEXT_ATOMS:
            raise ValueError("bounded immutable links required")
        observed = set(atoms)
        declared = observed | {ExactAtomRef(self.residue, name) for name in self.missing_atoms}
        pairs = set()
        for link in self.links:
            if (not isinstance(link, tuple) or len(link) != 2
                    or any(not isinstance(a, ExactAtomRef) for a in link)
                    or link[0] == link[1] or any(a not in declared for a in link)):
                raise ValueError("link endpoints must be distinct declared atoms")
            pair = frozenset(link)
            if pair in pairs:
                raise ValueError("duplicate external link")
            pairs.add(pair)

    @property
    def digest(self) -> str:
        # Input ordering is retained deliberately: no identity/case normalization.
        payload = {"revision": REVISION, "request": asdict(self)}
        return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":"),
                                         allow_nan=False).encode()).hexdigest()


@dataclass(frozen=True)
class Addition:
    atom: ExactAtomRef
    element: str
    xyz: Point


@dataclass(frozen=True)
class Candidate:
    request_digest: str
    additions: tuple[Addition, ...]
    torsion_degrees: float | None
    rotamer_evidence: RotamerEvidence | None = None


@dataclass(frozen=True)
class Audit:
    accepted: bool
    findings: tuple[str, ...]
    ca_volume_angstrom3: float | None = None
    minimum_nonbonded_clearance_angstrom: float | None = None


@dataclass(frozen=True)
class Attempt:
    candidate: Candidate
    audit: Audit


@dataclass(frozen=True)
class CompletionResult:
    status: Literal["completed", "unchanged", "unsupported", "rejected"]
    request_digest: str
    reason: str
    attempts: tuple[Attempt, ...] = ()  # includes rejected proposals
    selected_index: int | None = None
    revision: str = REVISION
    sampling_policy: str = "three-staggered-torsions-v1; no empirical rotamer probabilities"
    model_digest: str | None = None

    @property
    def selected(self) -> Candidate | None:
        return None if self.selected_index is None else self.attempts[self.selected_index].candidate
