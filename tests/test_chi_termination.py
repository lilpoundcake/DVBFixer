"""Explicit native diagnostics are not equivalent to a convergence certificate."""
from scripts.research_chi_termination import classify


def test_explicit_max_iteration_warning_is_distinct_from_line_search_abort():
    row=classify('core.optimization.Minimizer: [ WARNING ] LBFGS MAX CYCLES 200 EXCEEDED, BUT FUNC NOT CONVERGED!')
    assert row['explicit_max200_warning'] and not row['explicit_line_search_abort']
    assert not row['termination_certified']


def test_missing_diagnostics_do_not_certify_normal_convergence():
    row=classify('')
    assert not any(row.values())


def test_explicit_line_search_abort_is_not_inferred_from_energy_slopes():
    row=classify('Line search failed even after resetting Hessian; aborting at iter#18')
    assert row['explicit_line_search_abort'] and not row['explicit_max200_warning']
