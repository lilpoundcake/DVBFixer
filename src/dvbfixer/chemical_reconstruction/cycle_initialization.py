"""Frozen R1 boundary scaffold, never a geometry-approved parent or native pose.

Only initialization is implemented here. A record is deliberately NOT solver-ready:
initial geometry, numerical plan, independent global-L checks and solver traces are
separate owners. Saved-record checks never call NumPy or a coordinate constructor.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

from dvbfixer.domain.structure_identity import ExactAtomRef

from .contracts import (
    ContractError,
    Observation,
    Unit,
    _number,
    _point,
    _sha,
    _text,
    _tuple,
    digest_record,
)
from .cycle_contract import CycleConstructionProblem
from .measurement import Point, _cross, _dot, _sub, _unit

ALGORITHM = "boundary-harmonic-orbit-v1"
EQUATION_TOLERANCE = 1e-10


class CycleInitializationError(ValueError):
    """An explicit scaffold failure, not permission to replace the start."""


@dataclass(frozen=True)
class CycleObservedFrame:
    origin: Point
    axes: tuple[Point, Point, Point]

    def __post_init__(self) -> None:
        _point(self.origin)
        _tuple(self.axes, tuple, "observed frame axes", maximum=3)
        if len(self.axes) != 3:
            raise CycleInitializationError("cycle-initializer-frame-shape")
        for axis in self.axes:
            _point(axis)

    def local(self, point: Point) -> Point:
        delta = _sub(point, self.origin)
        return _dot(delta, self.axes[0]), _dot(delta, self.axes[1]), _dot(delta, self.axes[2])

    def source(self, point: Point) -> Point:
        return tuple_point(tuple(self.origin[i] + math.fsum(point[j]*self.axes[j][i] for j in range(3)) for i in range(3)))


def tuple_point(values: tuple[float, ...]) -> Point:
    if len(values) != 3 or not all(_number(v) for v in values):
        raise CycleInitializationError("cycle-initializer-nonfinite-point")
    return values[0], values[1], values[2]


@dataclass(frozen=True)
class CycleHarmonicSystem:
    problem_digest: str
    source_digest: str
    observations_digest: str
    graph_digest: str
    constraints_digest: str
    variables: tuple[ExactAtomRef, ...]
    anchors: tuple[ExactAtomRef, ExactAtomRef, ExactAtomRef]
    frame: CycleObservedFrame
    matrix: tuple[tuple[float, ...], ...]
    rhs: tuple[Point, ...]
    length_scale_angstrom: float

    def __post_init__(self) -> None:
        for pin in (self.problem_digest, self.source_digest, self.observations_digest, self.graph_digest, self.constraints_digest):
            _sha(pin)
        _tuple(self.variables, ExactAtomRef, "harmonic variables", maximum=16)
        _tuple(self.anchors, ExactAtomRef, "harmonic anchors", maximum=3)
        _tuple(self.matrix, tuple, "harmonic matrix", maximum=16)
        _tuple(self.rhs, tuple, "harmonic RHS", maximum=16)
        size = len(self.variables)
        if (not size or len(self.matrix) != size or len(self.rhs) != size or len(self.anchors) != 3
                or not isinstance(self.frame, CycleObservedFrame)
                or not _number(self.length_scale_angstrom) or self.length_scale_angstrom <= 0):
            raise CycleInitializationError("cycle-harmonic-system-shape")
        for row in self.matrix:
            if len(row) != size or not all(_number(v) for v in row):
                raise CycleInitializationError("cycle-harmonic-matrix-shape")
        for point in self.rhs:
            tuple_point(point)


def harmonic_system(problem: CycleConstructionProblem) -> CycleHarmonicSystem:
    """Assemble scalar graph weights/boundaries only; no ideal-coordinate access."""
    if not isinstance(problem, CycleConstructionProblem):
        raise CycleInitializationError("cycle-initializer-typed-problem-required")
    points = {o.identity: o.angstrom for o in problem.source.observations}
    origin, second, third = (points[a] for a in problem.anchors)
    try:
        e1 = _unit(_sub(second, origin))
        e3 = _unit(_cross(_sub(second, origin), _sub(third, origin)))
        frame = CycleObservedFrame(origin, (e1, _cross(e3, e1), e3))
    except ValueError as exc:
        raise CycleInitializationError("cycle-initializer-invalid-observed-frame") from exc
    variables = problem.coverage.variables
    indices = {a: i for i, a in enumerate(variables)}
    size = len(variables)
    matrix = [[0.]*size for _ in variables]
    rhs = [[0., 0., 0.] for _ in variables]
    lengths = []
    for metric in problem.coverage.metrics:
        if metric.kind != "bond":
            continue
        a, b = metric.atoms
        if a not in indices and b not in indices:
            continue
        length = metric.target
        try:
            weight = 1. / (length*length)
        except (ZeroDivisionError, OverflowError) as exc:
            raise CycleInitializationError("cycle-initializer-invalid-bond-weight") from exc
        if not math.isfinite(weight) or weight <= 0:
            raise CycleInitializationError("cycle-initializer-invalid-bond-weight")
        lengths.append(length)
        for atom, other in ((a, b), (b, a)):
            if atom not in indices:
                continue
            i = indices[atom]
            matrix[i][i] += weight
            if other in indices:
                matrix[i][indices[other]] -= weight
            else:
                boundary = frame.local(points[other])
                for axis in range(3):
                    rhs[i][axis] += weight*boundary[axis]
    lengths.sort()
    if not lengths:
        raise CycleInitializationError("cycle-initializer-no-generated-bond-lengths")
    middle = len(lengths)//2
    scale = float(lengths[middle] if len(lengths) % 2 else (lengths[middle-1]+lengths[middle])/2)
    return CycleHarmonicSystem(problem.digest, problem.source.digest, problem.coverage.observations_digest,
        problem.coverage.graph_digest, digest_record((*problem.source.constraints, *problem.extra_constraints)),
        variables, problem.anchors, frame, tuple(tuple(row) for row in matrix),
        tuple(tuple_point(tuple(row)) for row in rhs), scale)


def equation_residual(system: CycleHarmonicSystem, harmonic: tuple[Point, ...]) -> float:
    """Measure a SAVED system solution, without a linear solve or initializer."""
    _tuple(harmonic, tuple, "harmonic solution", maximum=16)
    if len(harmonic) != len(system.variables):
        raise CycleInitializationError("cycle-harmonic-solution-shape")
    for point in harmonic:
        tuple_point(point)
    try:
        values = tuple(abs(math.fsum(row[j]*harmonic[j][axis] for j in range(len(harmonic)))-rhs[axis])
                       for row, rhs in zip(system.matrix, system.rhs, strict=True) for axis in range(3))
    except (OverflowError, ValueError) as exc:
        raise CycleInitializationError("cycle-harmonic-nonfinite-equation") from exc
    if not all(math.isfinite(value) for value in values):
        raise CycleInitializationError("cycle-harmonic-nonfinite-equation")
    return max(values)


def offset_position(system: CycleHarmonicSystem, harmonic: Point, atom_index: int, branch_index: int) -> Point:
    """Frozen scalar orbit formula, also used to verify saved start records."""
    sign = 1 if branch_index < 4 else -1
    phi = (branch_index % 4)*math.pi/2
    theta = 2*math.pi*(atom_index+1)/(len(system.variables)+1)
    x, y, z = math.cos(theta), math.sin(theta), (-1)**atom_index/2
    rotated = x, y*math.cos(phi)-z*math.sin(phi), y*math.sin(phi)+z*math.cos(phi)
    point = tuple_point(tuple(harmonic[k]+sign*system.length_scale_angstrom/2*rotated[k] for k in range(3)))
    try:
        return system.frame.source(point)
    except OverflowError as exc:
        raise CycleInitializationError("cycle-initializer-nonfinite-point") from exc


@dataclass(frozen=True)
class CycleStart:
    problem_digest: str
    index: int
    sign: Literal[1, -1]
    quarter_turn: int
    additions: tuple[Observation, ...]
    status: Literal["coordinates-initialized", "failed"] = "coordinates-initialized"
    reason: str = "artificial-start-not-geometry-approved"

    def __post_init__(self) -> None:
        _sha(self.problem_digest)
        _tuple(self.additions, Observation, "cycle initial additions", maximum=16)
        _text(self.reason, "cycle start reason")
        if (type(self.index) is not int or not 0 <= self.index < 8
                or type(self.sign) is not int or self.sign != (1 if self.index < 4 else -1)
                or type(self.quarter_turn) is not int or self.quarter_turn != self.index % 4):
            raise CycleInitializationError("cycle-initializer-branch-label")
        if (self.status == "coordinates-initialized" and (not self.additions or self.reason != "artificial-start-not-geometry-approved")
                or self.status == "failed" and self.additions or self.status not in {"coordinates-initialized", "failed"}):
            raise CycleInitializationError("cycle-initializer-start-status")

    @property
    def digest(self) -> str:
        return digest_record((ALGORITHM, self))


def _element_map(problem: CycleConstructionProblem) -> dict[ExactAtomRef, str]:
    return {ExactAtomRef(c.identity, atom.name): atom.element for c in problem.source.components
            if c.identity in problem.source.targets and c.resolved is not None for atom in c.resolved.graph.atoms}


def _make_start(problem: CycleConstructionProblem, system: CycleHarmonicSystem,
                harmonic: tuple[Point, ...], index: int) -> CycleStart:
    elements = _element_map(problem)
    additions = tuple(Observation(a, elements[a], offset_position(system, harmonic[j], j, index), Unit.ANGSTROM)
                      for j, a in enumerate(system.variables))
    return CycleStart(problem.digest, index, 1 if index < 4 else -1, index % 4, additions)


@dataclass(frozen=True)
class CycleInitializerRecord:
    problem: CycleConstructionProblem
    status: Literal["initialized-unvalidated", "refused", "unavailable"]
    reason: str
    system: CycleHarmonicSystem | None = None
    harmonic: tuple[Point, ...] = ()
    maximum_equation_residual: float | None = None
    starts: tuple[CycleStart, ...] = ()
    algorithm: Literal["boundary-harmonic-orbit-v1"] = "boundary-harmonic-orbit-v1"

    def __post_init__(self) -> None:
        errors = audit_cycle_initializer(self)
        if errors:
            raise CycleInitializationError(",".join(errors))

    @property
    def digest(self) -> str:
        return digest_record(("cycle-initializer-record-v1", self))

    @property
    def coordinates_produced(self) -> bool:
        return any(start.additions for start in self.starts)

    @property
    def geometry_approved(self) -> Literal[False]:
        return False

    @property
    def solver_ready(self) -> Literal[False]:
        return False

    @property
    def md_ready(self) -> Literal[False]:
        return False

    @property
    def source_coverage_verified(self) -> Literal[False]:
        return False


def audit_cycle_initializer(record: CycleInitializerRecord) -> tuple[str, ...]:
    """Constructor-free provenance/equation/formula readback, NOT geometry approval."""
    if not isinstance(record.problem, CycleConstructionProblem):
        return ("cycle-initializer-typed-problem-required",)
    _text(record.reason, "initializer reason")
    _tuple(record.starts, CycleStart, "cycle starts", maximum=8)
    _tuple(record.harmonic, tuple, "saved harmonic solution", maximum=16)
    if record.algorithm != ALGORITHM:
        return ("cycle-initializer-algorithm-binding",)
    try:
        expected = harmonic_system(record.problem)
    except CycleInitializationError as exc:
        return (() if record.status == "refused" and record.reason == str(exc) and record.system is None
                and not record.harmonic and not record.starts and record.maximum_equation_residual is None
                else ("cycle-initializer-assembly-refusal-binding",))
    if not isinstance(record.system, CycleHarmonicSystem) or digest_record(record.system) != digest_record(expected):
        return ("cycle-initializer-system-source-binding",)
    if record.status == "unavailable":
        return (() if record.reason == "cycle-initializer-numpy-unavailable" and not record.harmonic and not record.starts
                and record.maximum_equation_residual is None else ("cycle-initializer-unavailable-binding",))
    if record.status == "refused" and not record.harmonic:
        return (() if record.reason in {"cycle-initializer-singular-system", "cycle-initializer-nonfinite-solution",
                                       "cycle-harmonic-nonfinite-equation"}
                and not record.starts and record.maximum_equation_residual is None else ("cycle-initializer-failure-binding",))
    residual = equation_residual(expected, record.harmonic)
    if digest_record(residual) != digest_record(record.maximum_equation_residual):
        return ("cycle-initializer-equation-metadata",)
    if residual > EQUATION_TOLERANCE:
        return (() if record.status == "refused" and record.reason == "cycle-initializer-equation-residual"
                and not record.starts else ("cycle-initializer-equation-failure-binding",))
    if (record.status != "initialized-unvalidated" or record.reason != "artificial-starts-no-geometry-or-solver-approval"
            or len(record.starts) != 8):
        return ("cycle-initializer-status-or-denominator",)
    elements = _element_map(record.problem)
    for index, start in enumerate(record.starts):
        if start.problem_digest != record.problem.digest or start.index != index:
            return ("cycle-initializer-start-identity-binding",)
        try:
            positions = tuple(offset_position(expected, record.harmonic[j], j, index) for j in range(len(expected.variables)))
            for position in positions:
                _point(position)
        except (CycleInitializationError, ContractError) as exc:
            if start.status != "failed" or start.reason != str(exc) or start.additions:
                return ("cycle-initializer-start-failure-binding",)
            continue
        if (start.status != "coordinates-initialized" or tuple(o.identity for o in start.additions) != expected.variables):
            return ("cycle-initializer-start-identity-binding",)
        for j, addition in enumerate(start.additions):
            if (addition.element != elements[addition.identity] or addition.unit != Unit.ANGSTROM
                    or type(addition.occupancy) is not float or addition.occupancy != 1.
                    or digest_record(addition.position) != digest_record(positions[j])):
                return ("cycle-initializer-start-formula-binding",)
    return ()


def _direct_solve(system: CycleHarmonicSystem) -> tuple[Point, ...]:
    import numpy as np

    try:
        value = np.linalg.solve(np.asarray(system.matrix, dtype=np.float64), np.asarray(system.rhs, dtype=np.float64))
    except np.linalg.LinAlgError as exc:
        raise CycleInitializationError("cycle-initializer-singular-system") from exc
    if value.shape != (len(system.variables), 3) or not np.isfinite(value).all():
        raise CycleInitializationError("cycle-initializer-nonfinite-solution")
    return tuple(tuple_point(tuple(float(v) for v in row)) for row in value)


def initialize_cycle_problem(problem: CycleConstructionProblem) -> CycleInitializerRecord:
    """One declared scaffold; no external starts, fitting, search or optimization."""
    if not isinstance(problem, CycleConstructionProblem):
        raise CycleInitializationError("cycle-initializer-typed-problem-required")
    try:
        system = harmonic_system(problem)
    except CycleInitializationError as exc:
        return CycleInitializerRecord(problem, "refused", str(exc))
    try:
        harmonic = _direct_solve(system)
    except ModuleNotFoundError as exc:
        if exc.name != "numpy":
            raise
        return CycleInitializerRecord(problem, "unavailable", "cycle-initializer-numpy-unavailable", system)
    except CycleInitializationError as exc:
        return CycleInitializerRecord(problem, "refused", str(exc), system)
    try:
        residual = equation_residual(system, harmonic)
    except CycleInitializationError as exc:
        return CycleInitializerRecord(problem, "refused", str(exc), system)
    if residual > EQUATION_TOLERANCE:
        return CycleInitializerRecord(problem, "refused", "cycle-initializer-equation-residual", system, harmonic, residual)
    starts_list = []
    for index in range(8):
        try:
            starts_list.append(_make_start(problem, system, harmonic, index))
        except (CycleInitializationError, ContractError) as exc:
            starts_list.append(CycleStart(problem.digest, index, 1 if index < 4 else -1, index % 4, (), "failed", str(exc)))
    starts = tuple(starts_list)
    return CycleInitializerRecord(problem, "initialized-unvalidated", "artificial-starts-no-geometry-or-solver-approval",
                                  system, harmonic, residual, starts)
