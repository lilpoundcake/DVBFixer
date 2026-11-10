"""Private single-branch VAL/THR completion anchored to a fixed partner.

One algebraic solution, not chi sampling or a rotamer distribution. Default
kernel/source/cohort/evidence entry points deliberately retain their refusals.
"""
from __future__ import annotations

from dvbfixer.diagnose.geometry import signed_dihedral_degrees
from dvbfixer.domain.structure_identity import ExactAtomRef
from dvbfixer.protein_completion.branched import BRANCHES, TETRAHEDRAL_ANGLE
from dvbfixer.protein_completion.contracts import (
    Addition,
    Attempt,
    Audit,
    Candidate,
    CompletionRequest,
    CompletionResult,
)
from dvbfixer.protein_completion.geometry import two_angle_tip
from dvbfixer.protein_completion.validation import _admission_reason, _audit_candidate

ANCHORED_REVISION = "canonical-val-thr-observed-branch-anchor-v1"
ANCHORED_POLICY = "observed-partner-two-angle-closure-v1; no empirical rotamer probabilities"


def audit_anchored_branch(request: CompletionRequest, candidate: Candidate, *,
                          pdb_roundtrip: bool = False) -> Audit:
    """Measurement-only geometry, fixed ownership, named stereo, chi and clashes."""
    return _audit_candidate(request, candidate, anchored_branch=True, pdb_roundtrip=pdb_roundtrip)


def complete_anchored_branch(request: CompletionRequest) -> CompletionResult:
    """Generate exactly one missing distal atom without touching retained anchors."""
    reason = _admission_reason(request, anchored_branch=True)
    if reason:
        return CompletionResult("unsupported", request.digest, reason,
                                revision=ANCHORED_REVISION, sampling_policy=ANCHORED_POLICY)
    branch = BRANCHES[request.residue_name]
    missing = request.missing_atoms[0]
    fixed = branch.second if missing == branch.first else branch.first
    points = {o.atom.name: o.xyz for o in request.observations if o.atom.component == request.residue}
    orientation = branch.expected_named_volume_sign*(-1 if missing == branch.first else 1)
    length = branch.first_length if missing == branch.first else branch.second_length
    try:
        tip = two_angle_tip(points["CB"], points["CA"], points[fixed], length=length,
                            first_angle=TETRAHEDRAL_ANGLE, second_angle=TETRAHEDRAL_ANGLE,
                            orientation=orientation)
        proposed = points | {missing: tip}
        chi = signed_dihedral_degrees(*(proposed[n] for n in ("N", "CA", "CB", branch.first)))
        if chi is None:
            raise ValueError("undefined observed-anchor chi")
        candidate = Candidate(request.digest, (Addition(ExactAtomRef(request.residue, missing), missing[0], tip),), chi)
        audit = audit_anchored_branch(request, candidate)
    except ValueError as exc:
        return CompletionResult("rejected", request.digest, f"observed_anchor_geometry:{exc}",
                                revision=ANCHORED_REVISION, sampling_policy=ANCHORED_POLICY)
    return CompletionResult("completed" if audit.accepted else "rejected", request.digest,
                            "observed_anchor_geometry_control" if audit.accepted else "anchored_candidate_failed_geometry",
                            (Attempt(candidate, audit),), 0 if audit.accepted else None,
                            revision=ANCHORED_REVISION, sampling_policy=ANCHORED_POLICY)
