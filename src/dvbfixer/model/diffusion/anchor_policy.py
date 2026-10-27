"""Pure, identity-based research policy for observed-anchor positional restraints.

This policy does not grant publication acceptance or change generated masks.
Chain-normalized B weights and population criteria remain under qualification.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass
from statistics import mean, median
from typing import Literal

from dvbfixer.model.diffusion.contract import AtomIdentity, ResidueIdentity

ANCHOR_POLICY_REVISION = "chain-backbone-b-anchor-restraints-v1"
BACKBONE = frozenset({"N", "CA", "C", "O"})


class NumericalRefinementError(RuntimeError):
    """An explicitly numerical minimization failure, not invalid input/config."""


@dataclass(frozen=True)
class AnchorPositionRestraint:
    atom: AtomIdentity
    reference_angstrom: tuple[float, float, float]
    k_kj_mol_nm2: float
    backbone_b: float | None = None
    chain_reference_b: float | None = None
    fallback_reason: str | None = None

    def __post_init__(self) -> None:
        if len(self.reference_angstrom) != 3 or not all(math.isfinite(x) for x in self.reference_angstrom):
            raise ValueError("anchor reference must contain three finite coordinates")
        if not math.isfinite(self.k_kj_mol_nm2) or self.k_kj_mol_nm2 <= 0:
            raise ValueError("anchor restraint strength must be finite and positive")
        for value in (self.backbone_b, self.chain_reference_b):
            if value is not None and (not math.isfinite(value) or value <= 0):
                raise ValueError("reported B values must be finite and positive")


def atom_residue(atom: AtomIdentity) -> ResidueIdentity:
    return ResidueIdentity(atom.chain, atom.residue_number, atom.insertion_code)


@dataclass(frozen=True)
class AnchorMovementLimits:
    """Explicit original-frame acceptance bounds; intentionally NO default values."""

    backbone_rmsd_angstrom_max: float
    heavy_displacement_angstrom_max: float

    def __post_init__(self) -> None:
        if not all(math.isfinite(x) and x > 0 for x in
                   (self.backbone_rmsd_angstrom_max, self.heavy_displacement_angstrom_max)):
            raise ValueError("anchor movement bounds must be finite and positive")


@dataclass(frozen=True)
class AnchorMovement:
    residue: ResidueIdentity
    backbone_rmsd_angstrom: float
    heavy_displacement_angstrom_max: float
    within_limits: bool


def check_anchor_movements(
    restraints: tuple[AnchorPositionRestraint, ...],
    candidate_coordinates_angstrom: Mapping[AtomIdentity, tuple[float, float, float]],
    limits: AnchorMovementLimits,
) -> tuple[AnchorMovement, ...]:
    """Per-anchor checks, never complex-averaged or Kabsch-aligned.

    Missing/duplicate identities, incomplete backbone or nonfinite coordinates
    fail closed instead of allowing an undefined metric to pass a bound.
    This does not replace independent geometry or outside-mask validation.
    """
    if len({r.atom for r in restraints}) != len(restraints):
        raise ValueError("duplicate anchor restraint identity")
    grouped: dict[ResidueIdentity, list[tuple[str, float]]] = {}
    for restraint in restraints:
        xyz = candidate_coordinates_angstrom.get(restraint.atom)
        if xyz is None or len(xyz) != 3 or not all(math.isfinite(x) for x in xyz):
            raise ValueError("anchor coordinate missing or nonfinite")
        displacement = math.dist(restraint.reference_angstrom, xyz)
        grouped.setdefault(atom_residue(restraint.atom), []).append((restraint.atom.atom_name, displacement))
    result = []
    for r, atoms in sorted(grouped.items()):
        backbone = [d for name, d in atoms if name in BACKBONE]
        if len(backbone) != 4:
            raise ValueError("anchor movement acceptance requires complete backbone")
        rmsd = math.sqrt(mean(d * d for d in backbone))
        maximum = max(d for _name, d in atoms)
        result.append(AnchorMovement(r, rmsd, maximum,
                      rmsd <= limits.backbone_rmsd_angstrom_max and maximum <= limits.heavy_displacement_angstrom_max))
    return tuple(result)


def build_anchor_restraints(
    source_text: str,
    anchors: tuple[ResidueIdentity, ...],
    *,
    observed_atoms: tuple[AtomIdentity, ...],
    base_k_kj_mol_nm2: float = 10000.0,
    bfactor_source: Literal["auto", "experimental", "disabled"] = "auto",
    minimum_chain_residues: int = 8,
) -> tuple[AnchorPositionRestraint, ...]:
    """Restrain existing observed heavy atoms; never synthesize missing atoms.

    Auto trusts only explicit X-ray/neutron headers without confidence markers.
    Unknown/EM/predicted provenance falls back to uniform. The explicit experimental
    override is an operator assertion, not a confidence classifier. Eight complete
    nonconstant backbone samples is a conservative engineering pilot criterion,
    NOT a literature-derived threshold. References include all eligible observed
    residues of a chain, independently of which anchors have been released.
    """
    if bfactor_source not in {"auto", "experimental", "disabled"}:
        raise ValueError("unknown B-factor provenance policy")
    if minimum_chain_residues < 2:
        raise ValueError("chain B reference requires at least two residues")
    if not math.isfinite(base_k_kj_mol_nm2) or base_k_kj_mol_nm2 <= 0:
        raise ValueError("base restraint strength must be finite and positive")
    if len(set(anchors)) != len(anchors) or len(set(observed_atoms)) != len(observed_atoms):
        raise ValueError("duplicate anchor or observed atom identity")
    selected = set(anchors)
    observed = set(observed_atoms)
    records: dict[AtomIdentity, tuple[tuple[float, float, float], float | None]] = {}
    header = "\n".join(line for line in source_text.splitlines()
                       if not line.startswith(("ATOM  ", "HETATM", "ANISOU"))).upper()
    experimental = any(line.startswith("EXPDTA") and any(method in line.upper() for method in
                       ("X-RAY DIFFRACTION", "NEUTRON DIFFRACTION")) for line in source_text.splitlines())
    confidence = any(marker in header for marker in ("PLDDT", "ALPHAFOLD", "PREDICTED ALIGNED ERROR"))
    provenance_reason = ("B weighting disabled" if bfactor_source == "disabled" else
                         "B-factor experimental provenance unavailable" if bfactor_source == "auto" and
                         (not experimental or confidence) else None)
    for line in source_text.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atom = AtomIdentity(line[21], line[22:26].strip(), line[26].strip(), line[12:16].strip())
        if atom not in observed:
            continue
        if atom in records:
            raise ValueError("ambiguous source atom/alternate state")
        element = line[76:78].strip().upper()
        if element in {"H", "D"} or (not element and atom.atom_name.lstrip("0123456789").startswith(("H", "D"))):
            raise ValueError("observed anchor policy requires heavy atoms")
        xyz = tuple(float(line[i:i + 8]) for i in (30, 38, 46))
        if not all(math.isfinite(x) for x in xyz):
            raise ValueError("nonfinite source coordinates")
        try:
            value = float(line[60:66])
            b = value if math.isfinite(value) and value > 0 else None
        except ValueError:
            b = None
        records[atom] = (xyz, b)
    if not observed <= records.keys():
        raise ValueError("source omits requested observed atoms")
    backbone: dict[ResidueIdentity, dict[str, float | None]] = {}
    for atom, (_xyz, b) in records.items():
        if atom.atom_name in BACKBONE:
            backbone.setdefault(atom_residue(atom), {})[atom.atom_name] = b
    means = {r: mean(float(b) for b in values.values() if b is not None)
             for r, values in backbone.items() if set(values) == BACKBONE and all(b is not None for b in values.values())}
    chain_values: dict[str, list[float]] = {}
    for r, b in means.items():
        chain_values.setdefault(r.chain, []).append(b)
    if any(not any(atom_residue(atom) == r for atom in records) for r in anchors):
        raise ValueError("requested anchor is not observed")
    restraints = []
    for atom in observed_atoms:
        r = atom_residue(atom)
        if r not in selected:
            continue
        values = chain_values.get(r.chain, [])
        reason = provenance_reason
        b = means.get(r)
        reference = None
        if reason is None:
            if b is None:
                reason = "anchor backbone B values incomplete or invalid"
            elif len(values) < minimum_chain_residues:
                reason = "insufficient complete observed chain backbone B samples"
            elif max(values) - min(values) <= 1e-6:
                reason = "chain backbone B distribution is constant"
            else:
                reference = median(values)
        weight = max(0.5, min(2.0, reference / b)) if reference is not None and b is not None else 1.0
        restraints.append(AnchorPositionRestraint(atom, records[atom][0], base_k_kj_mol_nm2 * weight,
                                                 b, reference, reason))
    return tuple(restraints)
