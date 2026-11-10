"""Verbose native diagnostics cannot silently change another optimizer option."""
import pytest
from scripts.research_chi_checkpoint_verbose import NUMERICAL_GETTERS, check_equal


def test_all_exposed_numerical_options_must_match():
    values=dict.fromkeys(NUMERICAL_GETTERS,0.)
    check_equal(values,values)
    with pytest.raises(ValueError,match='changed numerical'):
        check_equal(values,dict(values,max_iter=400))


def test_missing_numerical_getter_is_not_accepted():
    with pytest.raises(ValueError,match='changed numerical'):
        check_equal({}, {})
