"""Private spatial mobility, independent of fatal diagnostic thresholds.

All generated atoms are mobile. Observed residues within a fixed spatial shell,
direct gap anchors and one observed sequence neighbor are bounded, never free.
Masks are frozen from the initial proposal; repeated attempts cannot ratchet
observed reference coordinates or silently expand their allowed displacement.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from dvbfixer.ffutils.geometry import assert_all_l
from dvbfixer.model.diffusion.anchor_acceptance import _structure
from dvbfixer.model.diffusion.anchor_policy import (
    AnchorMovementLimits,
    atom_residue,
    build_anchor_restraints,
    check_anchor_movements,
)
from dvbfixer.model.diffusion.contract import DiffusionContractError, ResidueIdentity
from dvbfixer.model.diffusion.gap_neighbors import source_gap_components as source_gap_components
from dvbfixer.model.diffusion.gap_portfolio import _coordinates, digest
from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

REVISION = "all-generated-plus-bounded-observed-shell-v1"


@dataclass(frozen=True)
class SpatialMobility:
    source_sha256: str
    proposal_sha256: str
    observed_residues: tuple[ResidueIdentity, ...]
    shell_radius_angstrom: float = 5.0
    observed_sequence_padding: int = 1
    revision: str = REVISION


def spatial_mobility(request, source_text, proposal_text):
    """Use coordinates, not WARNING/ERROR severity, to identify first neighbors."""
    if digest(source_text) != request.normalized_pdb.sha256:
        raise DiffusionContractError("spatial relaxation source binding differs")
    assessment = assessment_adapter(request, source_text)(proposal_text)
    if assessment.admission_failures:
        raise DiffusionContractError("spatial relaxation input changed source observations")
    source = _coordinates(source_text)[0]
    proposal = _coordinates(proposal_text)[0]
    if not set(request.generated_atoms) <= proposal.keys():
        raise DiffusionContractError("spatial relaxation proposal omits generated atoms")
    generated = np.asarray([proposal[a] for a in request.generated_atoms])
    selected = {r for g in request.gaps for r in (g.left_anchor, g.right_anchor) if r is not None}
    observed = tuple(request.fixed_atoms)
    # Bound temporary arrays even for large protein complexes.
    for start in range(0, len(observed), 128):
        atoms = observed[start:start + 128]
        points = np.asarray([source[a] for a in atoms])
        distance = np.linalg.norm(points[:, None, :] - generated[None, :, :], axis=2).min(axis=1)
        selected.update(atom_residue(a) for a, d in zip(atoms, distance, strict=True) if d <= 5.0)
    # One layer only, determined before minimization; no recursive whole-chain
    # expansion and no artificial sequence adjacency across a declared gap.
    expanded = set(selected)
    for placement in request.sequence_placements:
        positions = dict(zip(placement.observed_target_indices, placement.observed_residues, strict=True))
        for index, residue in positions.items():
            if residue in selected:
                expanded.update(positions[j] for j in (index - 1, index + 1) if j in positions)
    backbone = {"N", "CA", "C", "O"}
    for residue in expanded:
        if not backbone <= {a.atom_name for a in observed if atom_residue(a) == residue}:
            raise DiffusionContractError("spatial observed shell has incomplete backbone")
    return SpatialMobility(digest(source_text), digest(proposal_text), tuple(sorted(expanded)))


def spatial_restraints(request, source_text, mobility):
    if mobility.revision != REVISION or mobility.source_sha256 != digest(source_text):
        raise DiffusionContractError("spatial mobility policy/source differs")
    return build_anchor_restraints(source_text, mobility.observed_residues,
        observed_atoms=request.fixed_atoms, bfactor_source="auto")


def audit_spatial_relaxation(request, source_text, proposal_text, candidate_text):
    """Recompute the shell, exact record custody, per-residue bounds and all-L."""
    mobility = spatial_mobility(request, source_text, proposal_text)
    restraints = spatial_restraints(request, source_text, mobility)
    allowed_observed = frozenset(r.atom for r in restraints)
    mobile = set(request.generated_atoms) | set(allowed_observed)
    xyz, rows = _coordinates(candidate_text)
    _, original = _coordinates(proposal_text)
    if set(rows) != set(original) or not mobile <= xyz.keys():
        raise DiffusionContractError("spatial candidate changed atom inventory")
    if candidate_text != rewrite_generated_coordinates(proposal_text, {a: xyz[a] for a in mobile}):
        raise DiffusionContractError("spatial candidate changed nonmobile records or metadata")
    limits = AnchorMovementLimits(0.5, 1.25, 0.5)
    movements = check_anchor_movements(restraints, xyz, limits)
    structure = _structure(candidate_text, "spatial-relaxation.pdb")
    assert_all_l(structure.pdb.topology, structure.pdb.positions)
    assessment = assessment_adapter(request, source_text)(candidate_text, allowed_observed)
    passed = assessment.geometry_passed and not assessment.admission_failures and all(m.within_limits for m in movements)
    return mobility, movements, assessment, passed

