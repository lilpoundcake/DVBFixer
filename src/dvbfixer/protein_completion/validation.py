"""Independent measurement-only audit: never reconstructs a candidate to validate it."""
from __future__ import annotations

import math

from dvbfixer.diagnose.geometry import (
    VDW_RADII_ANGSTROM,
    angle_degrees,
    distance_angstrom,
    signed_dihedral_degrees,
)
from dvbfixer.domain.structure_identity import ExactAtomRef
from dvbfixer.model.diffusion.protein_atoms import CANONICAL_HEAVY_ATOMS
from dvbfixer.protein_completion.backbone import backbone_angles
from dvbfixer.protein_completion.branched import BRANCHES, TETRAHEDRAL_ANGLE, named_beta_volume
from dvbfixer.protein_completion.contracts import Audit, Candidate, CompletionRequest
from dvbfixer.protein_completion.geometry import cross, dot, sub
from dvbfixer.protein_completion.rotamers import RotamerLibraryError, SingleChiLibrary

SUPPORTED = frozenset({"ALA", "SER", "CYS", "VAL", "THR"})
# Explicit v1 construction geometry, not an MD force field or calibrated library.
SIDECHAIN = {"SER": ("OG", "O", 1.417, 110.8), "CYS": ("SG", "S", 1.81, 114.4)}
_METALS = frozenset({"MG", "CA", "NA", "K", "ZN", "FE"})


def sidechain_bonds(residue_name: str) -> tuple[tuple[str, str], ...]:
    """The supported kernel's heavy-atom graph, including optional terminal OXT.

    This is connectivity only, not a source of bond orders or force-field terms.
    Unsupported residues refuse instead of receiving a guessed graph.
    """
    if residue_name not in SUPPORTED:
        raise ValueError("unsupported canonical graph")
    edges = (("N", "CA"), ("CA", "C"), ("C", "O"), ("CA", "CB"), ("C", "OXT"))
    if residue_name in BRANCHES:
        branch = BRANCHES[residue_name]
        return edges + (("CB", branch.first), ("CB", branch.second))
    return edges + ((("CB", SIDECHAIN[residue_name][0]),) if residue_name in SIDECHAIN else ())


def admission_reason(request: CompletionRequest) -> str | None:
    """Capability admission separate from a successful geometry result."""
    return _admission_reason(request)


def _admission_reason(request: CompletionRequest, *, anchored_branch: bool = False) -> str | None:
    if anchored_branch:
        branch = BRANCHES.get(request.residue_name)
        if (branch is None or len(request.missing_atoms) != 1
                or request.missing_atoms[0] not in {branch.first, branch.second}):
            return "anchored_branch_requires_one_distal_atom"
    if request.residue_name not in SUPPORTED:
        return "unsupported_residue"
    if request.residue.alternate_location:
        return "alternate_location_requires_selection"
    target = {o.atom.name: o for o in request.observations if o.atom.component == request.residue}
    if any(o.atom.component.model != request.residue.model for o in request.observations):
        return "multiple_models_require_selection"
    if any(o.atom.component.alternate_location or o.occupancy != 1 for o in request.observations):
        return "context_altloc_or_partial_occupancy"
    if not {"N", "CA", "C", "O"} <= target.keys():
        return "missing_backbone_not_supported"
    expected = set(CANONICAL_HEAVY_ATOMS[request.residue_name])
    if set(target)-expected-{"OXT"}:
        return "unexpected_target_atoms_or_observed_hydrogens"
    if any(o.element != name[0] for name, o in target.items()):
        return "target_element_mismatch"
    if set(request.missing_atoms) != expected-set(target):
        return "missing_mask_does_not_match_observations"
    if "CB" in request.missing_atoms and set(target)-{"N", "CA", "C", "O", "OXT"}:
        return "missing_cb_with_observed_distal_atoms_requires_closure"
    if request.residue_name in BRANCHES:
        branch = BRANCHES[request.residue_name]
        pair = {branch.first, branch.second}
        if not anchored_branch and pair & set(request.missing_atoms) and not pair <= set(request.missing_atoms):
            return "partial_beta_branch_requires_observed_anchor_solver"
    for o in request.observations:
        if o.element not in VDW_RADII_ANGSTROM:
            return "unsupported_context_element"
        # Include the reach of an absent Cβ + SG branch, not merely distances
        # to currently observed atoms. Eight angstroms conservatively covers
        # that first-slice reach plus a four-angstrom coordination neighbourhood.
        if o.element in _METALS and math.dist(o.xyz, target["CA"].xyz) < 8:
            return "nearby_metal_requires_coordination_model"
    for a, b in request.links:
        if a.component == request.residue or b.component == request.residue:
            inside, outside = (a, b) if a.component == request.residue else (b, a)
            if (outside.component == request.residue
                    or (inside.name, outside.name) not in {("N", "C"), ("C", "N")}):
                return "nonpeptide_target_link_requires_link_aware_solver"
            observed = {o.atom: o for o in request.observations}
            if (inside not in observed or outside not in observed
                    or observed[inside].element != inside.name
                    or observed[outside].element != outside.name):
                return "invalid_peptide_link"
            if not 1.0 <= math.dist(observed[inside].xyz, observed[outside].xyz) <= 1.8:
                return "invalid_peptide_link_length"
    for atom in (ExactAtomRef(request.residue, "N"), ExactAtomRef(request.residue, "C")):
        if sum(atom in pair for pair in request.links) > 1:
            return "multiple_external_peptide_partners"
    return None


def audit_candidate(request: CompletionRequest, candidate: Candidate, *,
                    library: SingleChiLibrary | None = None, pdb_roundtrip: bool = False) -> Audit:
    """Check exact generated ownership, fixed source geometry, chirality and sterics.

    Original observations are authoritative and never supplied by the candidate.
    Caller-provided context/links must be complete; no claim of file-level audit.
    """
    return _audit_candidate(request, candidate, library=library, pdb_roundtrip=pdb_roundtrip)


def _audit_candidate(request: CompletionRequest, candidate: Candidate, *,
                     library: SingleChiLibrary | None = None, pdb_roundtrip: bool = False,
                     anchored_branch: bool = False, measured_chi: bool = False) -> Audit:
    # Research comparator only: geometry at an arbitrary measured VAL/THR chi,
    # not admission to a generator's staggered policy. Public audits never opt in.
    reason = _admission_reason(request, anchored_branch=anchored_branch)
    if reason:
        return Audit(False, (reason,))
    if candidate.request_digest != request.digest:
        return Audit(False, ("request_digest_mismatch",))
    definition = SIDECHAIN.get(request.residue_name)
    empirical_required = library is not None and definition is not None and definition[0] in request.missing_atoms
    if candidate.rotamer_evidence is not None or empirical_required:
        if library is None or not empirical_required:
            return Audit(False, ("rotamer_authority_required_or_inapplicable",))
        if library.authority.category != {"SER": "SER", "CYS": "CYH"}.get(request.residue_name):
            return Audit(False, ("rotamer_category_mismatch",))
        try:
            # Recompute from authoritative observations and declared links;
            # candidate-supplied angles, probabilities and means are not trusted.
            phi, psi = backbone_angles(request)
            choices = library.lookup(phi, psi)
        except RotamerLibraryError as exc:
            return Audit(False, (f"rotamer_lookup:{exc}",))
        if candidate.rotamer_evidence not in choices:
            return Audit(False, ("rotamer_evidence_mismatch",))
    if not isinstance(candidate.additions, tuple):
        return Audit(False, ("immutable_additions_required",))
    expected = {ExactAtomRef(request.residue, n) for n in request.missing_atoms}
    try:
        actual = [a.atom for a in candidate.additions]
        if set(actual) != expected or len(actual) != len(expected):
            return Audit(False, ("generated_atom_set_mismatch",))
        for addition in candidate.additions:
            if addition.element != addition.atom.name[0]:
                return Audit(False, ("generated_element_mismatch",))
            if (not isinstance(addition.xyz, tuple) or len(addition.xyz) != 3
                    or any(isinstance(x, bool) or not isinstance(x, (int, float))
                           or not math.isfinite(x) or abs(x) > 1e5 for x in addition.xyz)):
                return Audit(False, ("invalid_generated_coordinates",))
    except (AttributeError, TypeError):
        return Audit(False, ("malformed_addition",))
    points = {o.atom: o.xyz for o in request.observations}
    elements = {o.atom: o.element for o in request.observations}
    points.update({a.atom: a.xyz for a in candidate.additions})
    elements.update({a.atom: a.element for a in candidate.additions})
    target = {atom.name: xyz for atom, xyz in points.items() if atom.component == request.residue}
    findings = []

    def bond(a: str, b: str, reference: float, tolerance: float) -> None:
        length = distance_angstrom(target[a], target[b])
        if length is None or abs(length-reference) > tolerance:
            findings.append(f"bond_length:{a}-{b}")

    def angle(a: str, b: str, c: str, reference: float, tolerance: float) -> None:
        value = angle_degrees(target[a], target[b], target[c])
        if value is None or abs(value-reference) > tolerance:
            findings.append(f"bond_angle:{a}-{b}-{c}")

    bond("N", "CA", 1.46, 0.18)
    bond("CA", "C", 1.53, 0.18)
    bond("C", "O", 1.24, 0.16)
    angle("N", "CA", "C", 111, 20)
    angle("CA", "C", "O", 121, 18)
    bond("CA", "CB", 1.53, 0.10)
    angle("N", "CA", "CB", 110.5, 12)
    angle("C", "CA", "CB", 110.1, 12)
    ca = target["CA"]
    volume = dot(cross(sub(target["N"], ca), sub(target["C"], ca)), sub(target["CB"], ca))
    if volume <= 0.5:
        findings.append("non_L_or_nearplanar_ca")
    edges = [(a, b) for a, b in sidechain_bonds(request.residue_name) if a in target and b in target]
    if "OXT" in target:
        bond("C", "OXT", 1.26, 0.18)
        angle("O", "C", "OXT", 125, 20)
    if request.residue_name in SIDECHAIN:
        terminal, _, length, degrees = SIDECHAIN[request.residue_name]
        bond("CB", terminal, length, 0.10)
        angle("CA", "CB", terminal, degrees, 12)
        chi = signed_dihedral_degrees(target["N"], ca, target["CB"], target[terminal])
        if terminal in request.missing_atoms:
            declared = candidate.rotamer_evidence
            allowed = ({declared.mean_degrees} if declared is not None else {-60.0, 60.0, 180.0})
            if (isinstance(candidate.torsion_degrees, bool)
                    or not isinstance(candidate.torsion_degrees, (float, int))
                    or candidate.torsion_degrees not in allowed or chi is None
                    # Rounded PDB readback has 0.001 A precision. Raw proposals
                    # retain the original exact construction/metadata gate.
                    or abs((chi-candidate.torsion_degrees+180) % 360-180) > (0.15 if pdb_roundtrip else 1e-5)):
                findings.append("torsion_metadata_or_grid_mismatch")
    elif request.residue_name in BRANCHES:
        branch = BRANCHES[request.residue_name]
        for branch_atom, length in ((branch.first, branch.first_length), (branch.second, branch.second_length)):
            bond("CB", branch_atom, length, .10)
            angle("CA", "CB", branch_atom, TETRAHEDRAL_ANGLE, 12)
        angle(branch.first, "CB", branch.second, TETRAHEDRAL_ANGLE, 12)
        beta = named_beta_volume(ca, target["CB"], target[branch.first], target[branch.second])
        if branch.expected_named_volume_sign*beta <= .5:
            findings.append("thr_cb_not_R_or_nearplanar" if request.residue_name == "THR"
                            else "val_branch_name_orientation")
        if anchored_branch or branch.first in request.missing_atoms:
            chi = signed_dihedral_degrees(target["N"], ca, target["CB"], target[branch.first])
            if (isinstance(candidate.torsion_degrees, bool)
                    or not isinstance(candidate.torsion_degrees, (float, int))
                    or not math.isfinite(candidate.torsion_degrees)
                    or (not -180 <= candidate.torsion_degrees <= 180 if anchored_branch or measured_chi else
                        candidate.torsion_degrees not in {-60., 60., 180.}) or chi is None
                    or abs((chi-candidate.torsion_degrees+180) % 360-180) > (0.15 if pdb_roundtrip else 1e-5)):
                findings.append("torsion_metadata_or_grid_mismatch")
        elif candidate.torsion_degrees is not None:
            findings.append("unexpected_torsion_metadata")
    elif candidate.torsion_degrees is not None:
        findings.append("unexpected_torsion_metadata")
    graph: dict[ExactAtomRef, set[ExactAtomRef]] = {a: set() for a in points}
    for a, b in [(ExactAtomRef(request.residue, a), ExactAtomRef(request.residue, b))
                 for a, b in edges] + list(request.links):
        graph[a].add(b)
        graph[b].add(a)
    minimum = None
    checked = set()
    for atom in sorted(expected):
        excluded = {atom} | graph[atom]
        # Ignore directly bonded and angle-connected pairs, not entire residues.
        excluded.update(neighbour for a in graph[atom] for neighbour in graph[a])
        for other in sorted(points):
            if other in excluded or frozenset((atom, other)) in checked:
                continue
            checked.add(frozenset((atom, other)))
            clearance = (math.dist(points[atom], points[other])
                         - VDW_RADII_ANGSTROM[elements[atom]] - VDW_RADII_ANGSTROM[elements[other]])
            minimum = clearance if minimum is None else min(minimum, clearance)
            if clearance < -0.6:
                findings.append(f"generated_clash:{atom.name}:{other.component}:{other.name}")
    return Audit(not findings, tuple(findings), volume, minimum)
