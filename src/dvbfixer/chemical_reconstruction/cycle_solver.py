"""One frozen joint numerical attempt per start; no source/native input or selection."""
from __future__ import annotations

import contextlib
import io
import platform
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from importlib.metadata import PackageNotFoundError, version
from typing import Literal

from .contracts import Observation, _number, _sha, _text, _tuple, digest_record
from .cycle_audit import round_cycle_additions
from .cycle_residuals import (
    CycleNumericalError,
    residual_values,
    squared_residual,
    vector_additions,
)
from .cycle_solve_plan import SETTINGS, CycleSolvePlan, record_bytes


@dataclass(frozen=True)
class CycleSolverStats:
    success: bool
    status: int
    message: str
    nfev: int
    njev: int | None

    def __post_init__(self) -> None:
        _text(self.message, "cycle solver message")
        if (type(self.success) is not bool or type(self.status) is not int or self.status not in {0, 1, 2, 3, 4}
                or self.success != (self.status > 0) or type(self.nfev) is not int or not 1 <= self.nfev <= SETTINGS.max_nfev
                or self.njev is not None and (type(self.njev) is not int or not 0 <= self.njev <= SETTINGS.max_nfev)):
            raise ValueError("cycle-solver-stats-binding")


@dataclass(frozen=True)
class CycleSolution:
    x: tuple[float, ...]
    stats: CycleSolverStats


@dataclass(frozen=True)
class CycleBranchResult:
    index: int
    branch_digest: str
    status: Literal["refused", "unavailable", "numerical-failed", "coordinates-returned", "aborted"]
    reason: str
    stats: CycleSolverStats | None
    returned_x: tuple[float, ...] | None
    raw: tuple[Observation, ...]
    rounded: tuple[Observation, ...]
    initial_residuals: tuple[float, ...] | None
    final_residuals: tuple[float, ...] | None
    optimizer_residual_calls: int
    external_residual_calls: int
    solver_calls: int
    last_finite_x: tuple[float, ...] | None
    last_finite_residuals: tuple[float, ...] | None

    def __post_init__(self) -> None:
        _sha(self.branch_digest)
        _text(self.reason, "cycle branch reason")
        _tuple(self.raw, Observation, "cycle raw coordinates", maximum=16)
        _tuple(self.rounded, Observation, "cycle rounded coordinates", maximum=16)
        if type(self.index) is not int or not 0 <= self.index < 8:
            raise ValueError("cycle-trace-branch-index")
        for values, limit in ((self.returned_x, 48), (self.last_finite_x, 48), (self.initial_residuals, 65536),
                              (self.final_residuals, 65536), (self.last_finite_residuals, 65536)):
            if values is not None:
                _tuple(values, float, "cycle finite numeric trace", maximum=limit)
                if not all(type(v) is float and _number(v) for v in values):
                    raise ValueError("cycle-trace-float-representation")
        if (type(self.optimizer_residual_calls) is not int or not 0 <= self.optimizer_residual_calls <= SETTINGS.max_optimizer_calls+1
                or type(self.external_residual_calls) is not int or not 0 <= self.external_residual_calls <= 2
                or type(self.solver_calls) is not int or self.solver_calls not in {0, 1}
                or (self.last_finite_x is None) != (self.last_finite_residuals is None)
                or self.status not in {"refused", "unavailable", "numerical-failed", "coordinates-returned", "aborted"}):
            raise ValueError("cycle-trace-count-or-status")


@dataclass(frozen=True)
class CycleConstructionResult:
    plan_digest: str
    environment: tuple[tuple[str, str], ...]
    branches: tuple[CycleBranchResult, ...]
    aborted_reason: str | None = None
    revision: str = "cycle-joint-numerical-traces-v1"

    def __post_init__(self) -> None:
        _sha(self.plan_digest)
        _tuple(self.environment, tuple, "cycle numerical environment", maximum=16)
        _tuple(self.branches, CycleBranchResult, "cycle complete branch traces", maximum=8)
        if (self.revision != "cycle-joint-numerical-traces-v1" or tuple(b.index for b in self.branches) != tuple(range(len(self.branches)))
                or self.aborted_reason is None and (len(self.branches) != 8 or any(b.status == "aborted" for b in self.branches))):
            raise ValueError("cycle-result-denominator-or-revision")
        from .research_decode import parse_research_json
        parse_research_json(record_bytes(self))

    @property
    def digest(self) -> str:
        return digest_record(self)

    @property
    def constructor_approved(self) -> Literal[False]:
        return False

    @property
    def md_ready(self) -> Literal[False]:
        return False


class CycleExecutionAborted(RuntimeError):
    def __init__(self, partial: CycleConstructionResult) -> None:
        self.partial = partial
        super().__init__(partial.aborted_reason)


class CycleEvaluationTracker:
    """All callbacks counted, including numerical differences and blocked breach.

    The 50,001st callback is retained as a budget breach but NOT evaluated. Initial
    and final measurement calls have a separate counter and do not consume the
    optimizer ceiling. Last finite evaluation is diagnostic, never a proposal.
    """
    def __init__(self, plan: CycleSolvePlan, index: int) -> None:
        self.plan, self.index = plan, index
        self.optimizer_calls = self.external_calls = 0
        self.last_x: tuple[float, ...] | None = None
        self.last_values: tuple[float, ...] | None = None

    def evaluate(self, x: tuple[float, ...], *, optimizer: bool) -> tuple[float, ...]:
        if optimizer:
            self.optimizer_calls += 1
            if self.optimizer_calls > self.plan.settings.max_optimizer_calls:
                raise CycleNumericalError("cycle-optimizer-residual-call-budget")
        else:
            self.external_calls += 1
        values = residual_values(self.plan, self.index, x)
        squared_residual(values)
        self.last_x, self.last_values = x, values
        return values

    def callback(self, values: Sequence[float]) -> tuple[float, ...]:
        if any(not _number(v) for v in values):
            # Count failed callbacks too; conversion must not hide bool/nonfinite.
            self.optimizer_calls += 1
            if self.optimizer_calls > self.plan.settings.max_optimizer_calls:
                raise CycleNumericalError("cycle-optimizer-residual-call-budget")
            raise CycleNumericalError("cycle-nonfinite-optimizer-vector")
        return self.evaluate(tuple(float(v) for v in values), optimizer=True)


def cycle_environment() -> tuple[tuple[str, str], ...]:
    output = io.StringIO()
    try:
        import numpy as np
    except ModuleNotFoundError as exc:
        if exc.name != "numpy":
            raise
        output.write("numpy-unavailable")
    else:
        with contextlib.redirect_stdout(output):
            np.show_config()
    rows = [("python", platform.python_version()), ("platform", platform.platform()), ("numpy_build", output.getvalue())]
    for name in ("numpy", "scipy", "openmm"):
        try:
            value = version(name)
        except PackageNotFoundError:
            value = "unavailable"
        rows.append((name, value))
    return tuple(rows)


def _solve(function: Callable[[Sequence[float]], tuple[float, ...]], variables: int) -> CycleSolution:
    import numpy as np
    from scipy.optimize import least_squares

    s = SETTINGS
    try:
        value = least_squares(function, [0.]*variables, bounds=(-s.coordinate_bound_angstrom, s.coordinate_bound_angstrom),
            method=s.method, jac=s.jac, tr_solver=s.tr_solver, x_scale=s.x_scale, diff_step=s.diff_step, loss=s.loss,
            ftol=s.ftol, xtol=s.xtol, gtol=s.gtol, max_nfev=s.max_nfev)
    except (np.linalg.LinAlgError, FloatingPointError, OverflowError) as exc:
        raise CycleNumericalError("cycle-solver-numerical-failure:"+str(exc)) from exc
    stats = CycleSolverStats(bool(value.success), int(value.status), str(value.message), int(value.nfev),
                             int(value.njev) if value.njev is not None else None)
    if value.x.shape != (variables,) or not np.isfinite(value.x).all():
        raise CycleNumericalError("cycle-nonfinite-or-incomplete-solver-return")
    return CycleSolution(tuple(float(v) for v in value.x), stats)


def _run_branch(plan: CycleSolvePlan, index: int) -> CycleBranchResult:
    branch = plan.branches[index]
    tracker = CycleEvaluationTracker(plan, index)
    initial = final = None
    stats = None
    returned = None
    raw: tuple[Observation, ...] = ()
    rounded: tuple[Observation, ...] = ()
    solver_calls = 0

    def record(status: Literal["refused", "unavailable", "numerical-failed", "coordinates-returned", "aborted"],
               reason: str) -> CycleBranchResult:
        return CycleBranchResult(index, branch.digest, status, reason, stats, returned, raw, rounded, initial, final,
            tracker.optimizer_calls, tracker.external_calls, solver_calls, tracker.last_x, tracker.last_values)

    if not branch.ready:
        return record("refused", branch.reason)
    try:
        initial = tracker.evaluate(tuple(0. for _ in range(3*len(plan.initializer.problem.coverage.variables))), optimizer=False)
        solver_calls = 1
        solution = _solve(tracker.callback, 3*len(plan.initializer.problem.coverage.variables))
        stats = solution.stats
        if (len(solution.x) != 3*len(plan.initializer.problem.coverage.variables)
                or any(type(v) is not float or not _number(v) for v in solution.x)):
            raise CycleNumericalError("cycle-nonfinite-or-incomplete-solver-return")
        returned = solution.x
        raw = vector_additions(plan, index, returned)
        rounded = round_cycle_additions(raw)
        final = tracker.evaluate(returned, optimizer=False)
        return record("coordinates-returned", "finite-numerical-return-not-geometry-approved")
    except ModuleNotFoundError as exc:
        if exc.name in {"scipy", "numpy"} or (exc.name or "").startswith("scipy."):
            solver_calls = 0
            return record("unavailable", "cycle-solver-dependency-unavailable:"+str(exc))
        raise _BranchAborted(record("aborted", type(exc).__name__+":"+str(exc))) from exc
    except CycleNumericalError as exc:
        return record("numerical-failed", str(exc))
    except Exception as exc:
        # NOT a scientific refusal: propagate with the partial numeric evidence.
        raise _BranchAborted(record("aborted", type(exc).__name__+":"+str(exc))) from exc


class _BranchAborted(RuntimeError):
    def __init__(self, branch: CycleBranchResult) -> None:
        self.branch = branch


def solve_cycle_plan(plan: CycleSolvePlan, *, on_branch: Callable[[CycleBranchResult], None] | None = None) -> CycleConstructionResult:
    """Numerics ONLY: no original/prepared source bytes, witness, scoring or selection.

    A caller can durably retain each branch via on_branch; its return value is never
    used. Unexpected errors abort with a partial result, not eight fabricated trials.
    """
    if not isinstance(plan, CycleSolvePlan):
        raise ValueError("typed-cycle-solve-plan-required")
    environment = cycle_environment()
    branches: list[CycleBranchResult] = []
    try:
        for index in range(8):
            try:
                branch = _run_branch(plan, index)
            except _BranchAborted as exc:
                branches.append(exc.branch)
                raise
            branches.append(branch)
            if on_branch is not None:
                on_branch(branch)
    except Exception as exc:
        reason = branches[-1].reason if isinstance(exc, _BranchAborted) else type(exc).__name__+":"+str(exc)
        raise CycleExecutionAborted(CycleConstructionResult(plan.digest, environment, tuple(branches), reason)) from exc
    return CycleConstructionResult(plan.digest, environment, tuple(branches))
