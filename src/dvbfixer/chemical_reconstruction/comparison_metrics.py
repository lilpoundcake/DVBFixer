"""Common geometry/native/overlap measurements, never proposal qualification."""
from __future__ import annotations

import math
from dataclasses import dataclass

from dvbfixer.domain.structure_identity import ExactAtomRef

from .authority_geometry import geometry_view
from .candidate_audit import measure_additions
from .contracts import LinkKind, Observation, ReconstructionRequest, digest_record
from .heldout_scoring import measure_heavy_holdout
from .proposals import AuditPolicy, CandidateAudit, GeometryPlan, Strategy


@dataclass(frozen=True)
class CommonHeavyMeasurement:
    request_digest: str
    additions_digest: str
    reference_digest: str
    geometry: CandidateAudit
    atom_count: int
    rmsd_angstrom: float
    maximum_angstrom: float
    severe_overlap_pairs: tuple[tuple[ExactAtomRef, ExactAtomRef, float], ...]
    minimum_eligible_separation_angstrom: float | None
    scope: str = "heavy-exact-name-fixed-frame; 0.7A-overlap-not-vdW; no-constructor-qualification"


def measure_common_heavy(request: ReconstructionRequest, policy: AuditPolicy,
                           additions: tuple[Observation, ...], reference: tuple[Observation, ...]) -> CommonHeavyMeasurement:
    """Measure the same graph/constraints without calling any coordinate builder."""
    # This is a geometry measurement plan, not an external constructor's plan.
    plan = GeometryPlan(request.digest, Strategy.INTERNAL, policy)
    geometry = measure_additions(request, plan, additions)
    count, rmsd, maximum = measure_heavy_holdout(request, additions, reference)
    _, bonds = geometry_view(request)
    adjacent: dict[ExactAtomRef, set[ExactAtomRef]] = {}
    pairs = [tuple(pair) for pair in bonds]
    pairs.extend((c.link.component_atom, c.link.other_atom) for c in request.connections if c.kind == LinkKind.COVALENT)
    for a, b in pairs:
        adjacent.setdefault(a, set()).add(b)
        adjacent.setdefault(b, set()).add(a)
    positions = {o.identity: o.angstrom for o in request.observations}
    positions.update({o.identity: o.angstrom for o in additions})
    distances = {}
    for a in sorted(request.ownership.generated):
        near = adjacent.get(a, set())
        excluded = {a} | near | {n for b in near for n in adjacent.get(b, set())}
        for b in positions:
            if b not in excluded:
                pair = (min(a, b), max(a, b))
                distances[pair] = math.dist(positions[a], positions[b])
    overlaps = tuple((a, b, distance) for (a, b), distance in sorted(distances.items())
                     if distance < policy.minimum_separation_angstrom)
    return CommonHeavyMeasurement(request.digest, digest_record(additions), digest_record(reference), geometry,
        count, rmsd, maximum, overlaps, min(distances.values()) if distances else None)
