"""State controls have fixed enrollment and distinguish platform replay modes."""
import pytest
from scripts.research_legacy_start_controls import FOCAL, SEEDS, platform_properties


def test_declared_control_budget_is_not_an_adaptive_failure_search():
    assert SEEDS==(1,2)
    assert FOCAL==('1CRN-THR2-CB-OG1-CG2','1ENH-THR27-CB-OG1-CG2')
    assert len(SEEDS)*len(FOCAL)*2*2==16


def test_unknown_platform_refuses_rather_than_silently_dropping_properties():
    class Unknown:
        def getName(self):
            return 'unknown'
    with pytest.raises(ValueError,match='undeclared'):
        platform_properties(Unknown(),None)
