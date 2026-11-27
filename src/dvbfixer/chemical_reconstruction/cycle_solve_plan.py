"""Frozen private R1 numerical admission and complete per-start context proofs."""
from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, fields, is_dataclass
from enum import Enum

from dvbfixer.domain.structure_identity import ExactAtomRef

from .contracts import ConstraintKind, _sha, _tuple, digest_record
from .cycle_audit import CycleGeometryAudit, chemical_adjacency, measure_cycle_geometry
from .cycle_initialization import CycleInitializerRecord, audit_cycle_initializer
from .cycle_source import CycleSourceMetadata, audit_cycle_source
from .template_comparison import validate_masked_source


class CycleSolvePlanError(ValueError):
    """Numerical-profile refusal, not an infeasible-chemistry claim."""


@dataclass(frozen=True)
class CycleSolverSettings:
    method: str = "trf"
    jac: str = "2-point"
    tr_solver: str = "exact"
    x_scale: float = 1.0
    diff_step: float = 1e-6
    loss: str = "linear"
    ftol: float = 1e-10
    xtol: float = 1e-10
    gtol: float = 1e-10
    max_nfev: int = 600
    max_optimizer_calls: int = 50000
    movement_cap_angstrom: float = 8.0
    coordinate_bound_angstrom: float = 8.0/math.sqrt(3.)
    maximum_residuals: int = 65536


SETTINGS = CycleSolverSettings()


def record_bytes(value: object) -> bytes:
    """Private exact-field wire payload; NOT a .dat/public loader."""
    def default(item: object) -> object:
        if is_dataclass(item) and not isinstance(item, type):
            return {f.name: getattr(item, f.name) for f in fields(item)}
        if isinstance(item, frozenset):
            return sorted(item)
        if isinstance(item, Enum):
            return item.value
        raise TypeError("unknown private cycle wire value")
    return json.dumps(value, default=default, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False).encode()


@dataclass(frozen=True)
class CycleBranchPlan:
    index: int
    start_digest: str
    pairs: tuple[tuple[int, int], ...]
    omitted_far_count: int
    omitted_far_digest: str
    minimum_omitted_distance_angstrom: float | None
    residual_count: int
    reason: str

    def __post_init__(self) -> None:
        _sha(self.start_digest)
        _sha(self.omitted_far_digest)
        _tuple(self.pairs, tuple, "cycle context pair indices", maximum=SETTINGS.maximum_residuals)
        if (type(self.index) is not int or not 0 <= self.index < 8
                or type(self.omitted_far_count) is not int or self.omitted_far_count < 0
                or type(self.residual_count) is not int or self.residual_count < 0
                or any(len(p) != 2 or any(type(i) is not int for i in p) or not 0 <= p[0] < p[1] for p in self.pairs)):
            raise CycleSolvePlanError("cycle-branch-plan-shape")

    @property
    def ready(self) -> bool:
        return self.reason == "declared-numerical-profile-only"

    @property
    def digest(self) -> str:
        return digest_record(self)


def point_order(initializer: CycleInitializerRecord) -> tuple[ExactAtomRef, ...]:
    problem = initializer.problem
    return tuple(sorted({o.identity for o in problem.source.observations} | problem.source.ownership.generated))


def _branch_plan(initializer: CycleInitializerRecord, index: int) -> CycleBranchPlan:
    problem, start = initializer.problem, initializer.starts[index]
    proof = hashlib.sha256(b"cycle-far-pairs-v1:canonical-row-sha256\n")
    if start.status != "coordinates-initialized":
        return CycleBranchPlan(index, start.digest, (), 0, proof.hexdigest(), None, 0, "cycle-start-failed:"+start.reason)
    points = {o.identity: o.angstrom for o in (*problem.source.observations, *start.additions)}
    order = point_order(initializer)
    indices = {a: i for i, a in enumerate(order)}
    generated = problem.source.ownership.generated
    adjacency = chemical_adjacency(problem.source)
    excluded = {}
    for a in generated:
        near = adjacency.get(a, set())
        excluded[a] = {a} | near | {n for b in near for n in adjacency.get(b, set())}
    pairs: list[tuple[int, int]] = []
    included = omitted = 0
    minimum = None
    # Canonical lexicographic pairs in O(observed*generated), not observed squared.
    for a in order:
        seconds = order[indices[a]+1:] if a in generated else tuple(b for b in sorted(generated) if b > a)
        for b in seconds:
            owner, other = (a, b) if a in generated else (b, a)
            if other in excluded[owner]:
                continue
            distance = math.dist(points[a], points[b])
            if (a not in generated or b not in generated) and distance > SETTINGS.movement_cap_angstrom+problem.policy.minimum_separation_angstrom:
                omitted += 1
                minimum = distance if minimum is None else min(minimum, distance)
                proof.update(digest_record((a, b, distance)).encode("ascii"))
            else:
                included += 1
                if len(pairs) < SETTINGS.maximum_residuals:
                    pairs.append((indices[a], indices[b]))
    chemical_count = len(problem.coverage.metrics)+problem.coverage.constraint_scalar_count
    count = chemical_count+included
    reason = "declared-numerical-profile-only" if count <= SETTINGS.maximum_residuals else "cycle-numerical-residual-budget"
    # A budget refusal retains counts/proof, not a silently truncated usable vector.
    return CycleBranchPlan(index, start.digest, tuple(pairs) if reason == "declared-numerical-profile-only" else (),
                           omitted, proof.hexdigest(), minimum, count, reason)


@dataclass(frozen=True)
class CycleSolvePlan:
    initializer: CycleInitializerRecord
    metadata: CycleSourceMetadata
    initial_audits: tuple[CycleGeometryAudit, ...]
    branches: tuple[CycleBranchPlan, ...]
    settings: CycleSolverSettings = SETTINGS
    revision: str = "cycle-joint-solver-plan-v1"

    def __post_init__(self) -> None:
        if (self.revision != "cycle-joint-solver-plan-v1" or digest_record(self.settings) != digest_record(SETTINGS)
                or audit_cycle_initializer(self.initializer) or self.initializer.status != "initialized-unvalidated"):
            raise CycleSolvePlanError("cycle-solve-plan-initializer-or-settings")
        problem = self.initializer.problem
        if len(problem.source.observations)+len(problem.source.ownership.generated) > 100000:
            raise CycleSolvePlanError("cycle-metadata-topology-total-budget")
        metadata = self.metadata
        if (metadata.problem_digest != problem.digest or metadata.request_digest != problem.source.digest
                or metadata.source_sha256 != problem.source.source_sha256
                or metadata.observations_digest != problem.coverage.observations_digest):
            raise CycleSolvePlanError("cycle-solve-plan-metadata-binding")
        validate_masked_source(metadata.masked_pdb.encode("ascii"), problem.source)
        for q in (*problem.source.constraints, *problem.extra_constraints):
            if q.kind == ConstraintKind.TORSION or q.kind == ConstraintKind.PLANARITY and q.lower not in {None, 0.}:
                raise CycleSolvePlanError("cycle-unqualified-numerical-constraint")
        _tuple(self.initial_audits, CycleGeometryAudit, "cycle initial audits", maximum=8)
        _tuple(self.branches, CycleBranchPlan, "cycle branch plans", maximum=8)
        if len(self.initial_audits) != 8 or len(self.branches) != 8:
            raise CycleSolvePlanError("cycle-solve-plan-eight-branches-required")
        for index, report in enumerate(self.initial_audits):
            start = self.initializer.starts[index]
            if (report.phase != "initial" or report.problem_digest != problem.digest or report.source_metadata_digest != metadata.digest
                    or report.additions_digest != digest_record(start.additions)):
                raise CycleSolvePlanError("cycle-initial-audit-binding")
            if report.global_l_status in {"unavailable", "unmeasurable"}:
                raise CycleSolvePlanError("cycle-initial-audit-"+report.global_l_status)
            if digest_record(self.branches[index]) != digest_record(_branch_plan(self.initializer, index)):
                raise CycleSolvePlanError("cycle-context-pair-proof-binding")
        # Preserve the existing codec's exact byte/node/depth limits. Local import
        # avoids a registry import cycle; parsing does not restore/construct DTOs.
        from .research_decode import parse_research_json
        parse_research_json(record_bytes(self))

    @property
    def digest(self) -> str:
        return digest_record(("cycle-joint-solver-plan-v1", self))

    @property
    def constructor_approved(self) -> bool:
        return False


def compile_cycle_solve_plan(initializer: CycleInitializerRecord, prepared: bytes,
                             metadata: CycleSourceMetadata) -> CycleSolvePlan:
    """External preflight measures starts; COMPLETE prepared bytes are not retained."""
    if initializer.status != "initialized-unvalidated" or audit_cycle_initializer(initializer):
        raise CycleSolvePlanError("cycle-initializer-not-ready-for-numerical-planning")
    problem = initializer.problem
    if len(problem.source.observations)+len(problem.source.ownership.generated) > 100000:
        raise CycleSolvePlanError("cycle-metadata-topology-total-budget")
    audit_cycle_source(problem, prepared, metadata)
    initial = tuple(measure_cycle_geometry(problem, prepared, metadata, start.additions, phase="initial")
                    for start in initializer.starts)
    return CycleSolvePlan(initializer, metadata, initial, tuple(_branch_plan(initializer, i) for i in range(8)))
