"""Predeclared stopping controls do not introduce native-based profile selection."""
import pytest
from scripts.research_chi_stopping_controls import check_settings, profiles


def test_profiles_change_one_stopping_setting_at_each_sequential_step():
    rows=profiles()
    assert len(rows)==3
    assert {k for k in rows[0] if k!='id' and rows[0][k]!=rows[1][k]}=={'tolerance'}
    assert {k for k in rows[1] if k!='id' and rows[1][k]!=rows[2][k]}=={'gradient_cutoff'}
    assert all(r['type']=='lbfgs_armijo_nonmonotone_atol' for r in rows)


def test_runtime_settings_must_match_plan_not_a_silent_backend_default():
    settings=dict(type='lbfgs_armijo_nonmonotone_atol',tolerance=1e-6,gradient_cutoff=1e-4,
                  max_iter=200,cartesian=False,neighbor_list=False)
    check_settings(settings,settings)
    with pytest.raises(ValueError,match='gradient_cutoff'):
        check_settings(dict(settings,gradient_cutoff=1.),settings)
