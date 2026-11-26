"""Independent R1 coordinate measurements; no optimizer or fabricated legacy plan."""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

from dvbfixer.domain.structure_identity import ExactAtomRef

from .closure_diagnostics import MeasuredTerm, _constraint_value
from .contracts import (
    LinkKind,
    Observation,
    ReconstructionRequest,
    Unit,
    _sha,
    _tuple,
    digest_record,
)
from .cycle_contract import CycleConstructionProblem
from .cycle_initialization import CycleInitializerRecord, audit_cycle_initializer
from .cycle_source import CycleSourceMetadata, audit_cycle_source
from .measurement import angle_degrees, satisfied
from .metadata_chirality import MetadataTopologyError, assert_metadata_all_l

Phase = Literal["initial", "raw", "rounded", "external-witness"]
GlobalLStatus = Literal["passed", "rejected", "unavailable", "unmeasurable", "not-measured"]


def chemical_adjacency(request: ReconstructionRequest) -> dict[ExactAtomRef, set[ExactAtomRef]]:
    """Declared graph/link exclusions only, never an opaque bond guess or crop."""
    adjacent: dict[ExactAtomRef, set[ExactAtomRef]] = {}
    pairs = [(ExactAtomRef(c.identity, b.atom1), ExactAtomRef(c.identity, b.atom2))
             for c in request.components if c.resolved is not None for b in c.resolved.graph.bonds]
    pairs.extend((c.link.component_atom, c.link.other_atom) for c in request.connections if c.kind == LinkKind.COVALENT)
    for a, b in pairs:
        adjacent.setdefault(a, set()).add(b)
        adjacent.setdefault(b, set()).add(a)
    return adjacent


@dataclass(frozen=True)
class CycleGeometryAudit:
    problem_digest: str
    source_metadata_digest: str
    additions_digest: str
    phase: Phase
    findings: tuple[str, ...]
    terms: tuple[MeasuredTerm, ...]
    eligible_pair_count: int
    severe_overlap_pairs: tuple[tuple[ExactAtomRef, ExactAtomRef, float], ...]
    minimum_eligible_separation_angstrom: float | None
    global_l_status: GlobalLStatus
    global_l_detail: str

    def __post_init__(self) -> None:
        for pin in (self.problem_digest, self.source_metadata_digest, self.additions_digest):
            _sha(pin)
        _tuple(self.findings, str, "cycle geometry findings", maximum=1024)
        _tuple(self.terms, MeasuredTerm, "cycle geometry terms", maximum=1024)
        _tuple(self.severe_overlap_pairs, tuple, "cycle severe pairs", maximum=1600000)
        if (self.phase not in {"initial", "raw", "rounded", "external-witness"}
                or self.global_l_status not in {"passed", "rejected", "unavailable", "unmeasurable", "not-measured"}
                or self.findings != tuple(sorted(set(self.findings)))
                or type(self.eligible_pair_count) is not int or self.eligible_pair_count < len(self.severe_overlap_pairs)):
            raise ValueError("cycle-geometry-audit-shape")

    @property
    def passed(self) -> bool:
        return not self.findings and self.global_l_status == "passed" and all(t.passed for t in self.terms) and not self.severe_overlap_pairs

    @property
    def constructor_approved(self) -> Literal[False]:
        return False

    @property
    def md_ready(self) -> Literal[False]:
        return False

    @property
    def source_coverage_verified(self) -> Literal[False]:
        return False


def _global_l(observations: tuple[Observation, ...], metadata: CycleSourceMetadata) -> tuple[GlobalLStatus, str]:
    try:
        from dvbfixer.ffutils.geometry import ChiralityError
    except ModuleNotFoundError as exc:
        if exc.name not in {"openmm", "numpy"} and not (exc.name or "").startswith("openmm."):
            raise
        return "unavailable", str(exc)
    try:
        assert_metadata_all_l(observations, dict(metadata.residue_names))
    except ModuleNotFoundError as exc:
        if exc.name not in {"openmm", "numpy"} and not (exc.name or "").startswith("openmm."):
            raise
        return "unavailable", str(exc)
    except MetadataTopologyError as exc:
        return "unmeasurable", str(exc)
    except ChiralityError as exc:
        return "rejected", str(exc)
    return "passed", "metadata-only-assert_all_l-no-repair"


def measure_cycle_geometry(problem: CycleConstructionProblem, prepared: bytes, metadata: CycleSourceMetadata,
                           additions: tuple[Observation, ...], *, phase: Phase) -> CycleGeometryAudit:
    """Authenticate source metadata first; never extract or score withheld positions."""
    audit_cycle_source(problem, prepared, metadata)
    _tuple(additions, Observation, "cycle measured additions", maximum=16)
    errors = set()
    generated = problem.source.ownership.generated
    if len(additions) != len(generated) or {o.identity for o in additions} != generated:
        errors.add("cycle-generated-mask")
    elements = {ExactAtomRef(c.identity, a.name): a.element for c in problem.source.components
                if c.resolved is not None for a in c.resolved.graph.atoms}
    if any(o.identity not in elements or o.element.upper() != elements[o.identity].upper() for o in additions):
        errors.add("cycle-generated-element")
    if any(o.unit != Unit.ANGSTROM or o.occupancy != 1 for o in additions):
        errors.add("cycle-generated-unit-or-occupancy")
    if errors:
        return CycleGeometryAudit(problem.digest, metadata.digest, digest_record(additions), phase, tuple(sorted(errors)),
                                  (), 0, (), None, "not-measured", "invalid-additions-not-spliced")
    points = {o.identity: o.angstrom for o in problem.source.observations}
    points.update({o.identity: o.angstrom for o in additions})
    terms = []
    for metric in problem.coverage.metrics:
        value = None
        error = None
        try:
            p = tuple(points[a] for a in metric.atoms)
            value = math.dist(p[0], p[1]) if metric.kind == "bond" else angle_degrees(*p)
        except ValueError as exc:
            error = str(exc)
        passed = value is not None and math.isfinite(value) and abs(value-metric.target) <= metric.tolerance
        if not passed:
            errors.add("cycle-"+metric.kind+"-geometry" if error is None else "cycle-degenerate-metric")
        terms.append(MeasuredTerm("cycle-"+metric.kind, "distance" if metric.kind == "bond" else "angle", metric.atoms,
            metric.role, "angstrom" if metric.kind == "bond" else "degree", value, metric.target-metric.tolerance,
            metric.target+metric.tolerance, passed, error=error))
    for index, constraint in enumerate((*problem.source.constraints, *problem.extra_constraints)):
        value = lower = upper = None
        error = None
        try:
            p = tuple(points[a] for a in constraint.atoms)
            value, lower, upper = _constraint_value(constraint, p)
            passed = math.isfinite(value) and satisfied(constraint, p)
        except (KeyError, ValueError) as exc:
            error, passed = str(exc), False
        if not passed:
            errors.add("cycle-declared-constraint" if error is None else "cycle-unmeasurable-constraint")
        terms.append(MeasuredTerm("cycle-constraint", constraint.kind.value, constraint.atoms,
            "generated-involving" if set(constraint.atoms) & generated else "fixed-only", constraint.unit.value,
            value, lower, upper, passed, label=f"{index}:{constraint.label}", error=error))
    adjacency = chemical_adjacency(problem.source)
    pair_count = 0
    overlaps = []
    minimum = None
    for a in sorted(generated):
        near = adjacency.get(a, set())
        excluded = {a} | near | {n for b in near for n in adjacency.get(b, set())}
        for b in sorted(points):
            if b in excluded or b in generated and b < a:
                continue
            pair_count += 1
            distance = math.dist(points[a], points[b])
            minimum = distance if minimum is None else min(minimum, distance)
            if distance < problem.policy.minimum_separation_angstrom:
                overlaps.append((min(a, b), max(a, b), distance))
    if overlaps:
        errors.add("cycle-generated-context-minimum-separation")
    global_status, detail = _global_l((*problem.source.observations, *additions), metadata)
    if global_status != "passed":
        errors.add("cycle-global-L-"+global_status)
    return CycleGeometryAudit(problem.digest, metadata.digest, digest_record(additions), phase, tuple(sorted(errors)),
                              tuple(terms), pair_count, tuple(sorted(overlaps)), minimum, global_status, detail)


def audit_cycle_geometry(problem: CycleConstructionProblem, prepared: bytes, metadata: CycleSourceMetadata,
                         additions: tuple[Observation, ...], report: CycleGeometryAudit) -> None:
    expected = measure_cycle_geometry(problem, prepared, metadata, additions, phase=report.phase)
    if digest_record(report) != digest_record(expected):
        raise ValueError("cycle-saved-geometry-binding")


def round_cycle_additions(additions: tuple[Observation, ...]) -> tuple[Observation, ...]:
    """Declared PDB transport rounding, never a coordinate repair."""
    return tuple(Observation(o.identity, o.element,
        (float(f"{o.position[0]:.3f}"), float(f"{o.position[1]:.3f}"), float(f"{o.position[2]:.3f}")), o.unit, o.occupancy)
        for o in additions)


@dataclass(frozen=True)
class CyclePairAudit:
    initializer_digest: str
    branch_index: int
    raw: CycleGeometryAudit
    rounded: CycleGeometryAudit
    raw_maximum_movement_angstrom: float | None
    rounded_maximum_movement_angstrom: float | None
    rounding_matches: bool
    findings: tuple[str, ...]

    def __post_init__(self) -> None:
        _sha(self.initializer_digest)
        _tuple(self.findings, str, "cycle pair findings", maximum=1024)
        if (type(self.branch_index) is not int or not 0 <= self.branch_index < 8 or type(self.rounding_matches) is not bool
                or self.raw.phase != "raw" or self.rounded.phase != "rounded"
                or self.findings != tuple(sorted(set(self.findings)))):
            raise ValueError("cycle-pair-audit-shape")

    @property
    def passed(self) -> bool:
        return (not self.findings and self.raw.passed and self.rounded.passed and self.rounding_matches
                and self.raw_maximum_movement_angstrom is not None and self.raw_maximum_movement_angstrom <= 8.
                and self.rounded_maximum_movement_angstrom is not None and self.rounded_maximum_movement_angstrom <= 8.)

    @property
    def constructor_approved(self) -> Literal[False]:
        return False

    @property
    def md_ready(self) -> Literal[False]:
        return False

    @property
    def source_coverage_verified(self) -> Literal[False]:
        return False


def measure_cycle_pair(problem: CycleConstructionProblem, prepared: bytes, metadata: CycleSourceMetadata,
                       initializer: CycleInitializerRecord, branch_index: int, raw: tuple[Observation, ...],
                       rounded: tuple[Observation, ...]) -> CyclePairAudit:
    """Independent pair/cap check; does NOT claim a solver converged or select it."""
    if (audit_cycle_initializer(initializer) or initializer.problem.digest != problem.digest
            or initializer.status != "initialized-unvalidated" or type(branch_index) is not int or not 0 <= branch_index < 8):
        raise ValueError("cycle-pair-initializer-binding")
    start = initializer.starts[branch_index]
    if start.status != "coordinates-initialized":
        raise ValueError("cycle-pair-failed-start")
    raw_audit = measure_cycle_geometry(problem, prepared, metadata, raw, phase="raw")
    rounded_audit = measure_cycle_geometry(problem, prepared, metadata, rounded, phase="rounded")
    errors = set()
    if not raw_audit.passed or not rounded_audit.passed:
        errors.add("cycle-raw-or-rounded-geometry")
    matches = digest_record(round_cycle_additions(raw)) == digest_record(rounded)
    if not matches:
        errors.add("cycle-rounding-transport")
    initial = {o.identity: o.angstrom for o in start.additions}

    def movement(values: tuple[Observation, ...]) -> float | None:
        if len(values) != len(initial) or {o.identity for o in values} != initial.keys():
            errors.add("cycle-movement-mask")
            return None
        maximum = max(math.dist(o.angstrom, initial[o.identity]) for o in values)
        if maximum > 8.:
            errors.add("cycle-per-atom-movement-cap")
        return maximum

    raw_movement, rounded_movement = movement(raw), movement(rounded)
    return CyclePairAudit(initializer.digest, branch_index, raw_audit, rounded_audit, raw_movement, rounded_movement,
                          matches, tuple(sorted(errors)))


def audit_cycle_pair(problem: CycleConstructionProblem, prepared: bytes, metadata: CycleSourceMetadata,
                     initializer: CycleInitializerRecord, raw: tuple[Observation, ...], rounded: tuple[Observation, ...],
                     report: CyclePairAudit) -> None:
    expected = measure_cycle_pair(problem, prepared, metadata, initializer, report.branch_index, raw, rounded)
    if digest_record(report) != digest_record(expected):
        raise ValueError("cycle-saved-pair-binding")
