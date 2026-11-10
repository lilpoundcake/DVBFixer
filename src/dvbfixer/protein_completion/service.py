"""Private deterministic canonical-sidechain proposals and independent rejection.

No original atom is returned as generated, moved, reflected or optimized. The
default three-staggered-chi control is distinct from an explicitly provided local
single-chi empirical authority. Neither guarantees the native conformation.
"""
from __future__ import annotations

from typing import TypedDict

from dvbfixer.domain.structure_identity import ExactAtomRef
from dvbfixer.protein_completion.backbone import backbone_angles
from dvbfixer.protein_completion.branched import BRANCHED_REVISION, BRANCHES, TETRAHEDRAL_ANGLE
from dvbfixer.protein_completion.contracts import (
    Addition,
    Attempt,
    Candidate,
    CompletionRequest,
    CompletionResult,
)
from dvbfixer.protein_completion.geometry import extend, l_cb
from dvbfixer.protein_completion.rotamers import (
    LOOKUP_POLICY,
    RotamerEvidence,
    RotamerLibraryError,
    SingleChiLibrary,
)
from dvbfixer.protein_completion.validation import SIDECHAIN, admission_reason, audit_candidate


class _ModelMetadata(TypedDict):
    sampling_policy: str
    model_digest: str | None
    revision: str


def complete_sidechain(request: CompletionRequest, *,
                       library: SingleChiLibrary | None = None) -> CompletionResult:
    """Return all attempts; select only independently passing generated geometry."""
    digest = request.digest
    policy = (LOOKUP_POLICY if library is not None else
              "three-staggered-torsions-v1; no empirical rotamer probabilities")
    metadata: _ModelMetadata = {
        "sampling_policy": policy, "model_digest": library.model_digest if library is not None else None,
        "revision": "canonical-sidechain-rotamer-adapter-v2" if library is not None else "canonical-sidechain-internal-v1",
    }
    if request.residue_name in BRANCHES:
        metadata["revision"] = BRANCHED_REVISION
        if library is not None:
            return CompletionResult("unsupported", digest, "branched_empirical_authority_not_supported", **metadata)
    reason = admission_reason(request)
    if reason:
        return CompletionResult("unsupported", digest, reason, **metadata)
    if not request.missing_atoms:
        # No new atoms does not mean the input structure is geometrically valid.
        audit = audit_candidate(request, Candidate(digest, (), None), library=library)
        return CompletionResult("unchanged" if audit.accepted else "rejected", digest,
                                "no_atoms_requested" if audit.accepted else ";".join(audit.findings), **metadata)
    target = {o.atom.name: o.xyz for o in request.observations if o.atom.component == request.residue}
    try:
        if "CB" not in target:
            target["CB"] = l_cb(target["N"], target["CA"], target["C"])
        terminal = SIDECHAIN.get(request.residue_name)
        branch = BRANCHES.get(request.residue_name)
        empirical: tuple[RotamerEvidence, ...] = ()
        chis: tuple[float | None, ...]
        if library is not None and terminal and terminal[0] in request.missing_atoms:
            if library.authority.category != {"SER": "SER", "CYS": "CYH"}.get(request.residue_name):
                return CompletionResult("unsupported", digest, "rotamer_category_mismatch", **metadata)
            phi, psi = backbone_angles(request)
            empirical = library.lookup(phi, psi)
            chis = tuple(row.mean_degrees for row in empirical)
        else:
            sampling = ((terminal is not None and terminal[0] in request.missing_atoms)
                        or (branch is not None and branch.first in request.missing_atoms))
            chis = (-60.0, 60.0, 180.0) if sampling else (None,)
        attempts = []
        for index, chi in enumerate(chis):
            proposal = dict(target)
            if chi is not None and terminal is not None:
                name, _, length, angle = terminal
                proposal[name] = extend(target["N"], target["CA"], target["CB"], length, angle, chi)
            if chi is not None and branch is not None:
                proposal[branch.first] = extend(target["N"], target["CA"], target["CB"],
                                                branch.first_length, TETRAHEDRAL_ANGLE, chi)
                proposal[branch.second] = extend(target["N"], target["CA"], target["CB"],
                                                 branch.second_length, TETRAHEDRAL_ANGLE,
                                                 chi+branch.second_chi_offset)
            additions = tuple(Addition(ExactAtomRef(request.residue, name), name[0], proposal[name])
                              for name in sorted(request.missing_atoms))
            candidate = Candidate(digest, additions, chi, empirical[index] if empirical else None)
            attempts.append(Attempt(candidate, audit_candidate(request, candidate, library=library)))
    except RotamerLibraryError as exc:
        return CompletionResult("unsupported", digest, f"rotamer_lookup:{exc}", **metadata)
    except ValueError as exc:
        return CompletionResult("rejected", digest, f"construction_geometry:{exc}", **metadata)
    passing = [i for i, attempt in enumerate(attempts) if attempt.audit.accepted]
    # Rank only feasible candidates by worst nonbonded clearance, capped at zero:
    # no invented H-bond/electrostatic energy or rotamer probabilities. Stable
    # enumeration is the tie break; no coordinate-reference scoring occurs here.
    def score(index: int) -> tuple[float, float]:
        clearance = attempts[index].audit.minimum_nonbonded_clearance_angstrom
        evidence = attempts[index].candidate.rotamer_evidence
        # Feasibility was gated first. Within that feasible set prefer the
        # empirical prior; capped clearance breaks ties. This is not a physical
        # posterior or H-bond/solvation/packing energy.
        return (evidence.probability if evidence else 0.0,
                min(0.0, clearance if clearance is not None else 0.0))

    selected = max(passing, key=score) if passing else None
    return CompletionResult("completed" if passing else "rejected", digest,
                            "geometry_candidate_not_native_conformation" if passing else "all_candidates_failed_geometry",
                            tuple(attempts), selected, **metadata)
