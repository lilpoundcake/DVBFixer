"""Prospective complete starting pools preserve import blockers and empty geometry."""
import pytest
from scripts.research_prospective_guard_views import scored_pool_status


def test_empty_pool_and_failed_import_are_not_ready_refinement_starts():
    assert scored_pool_status([],[])=='empty_geometry_pool'
    assert scored_pool_status([dict(id='own-0')],[dict(id='own-0',score_status='import_contract_blocker',import_audit=dict(accepted=False))])=='score_import_blocker'


def test_missing_or_extra_start_inventory_refuses():
    with pytest.raises(ValueError,match='inventory'):
        scored_pool_status([dict(id='own-0')],[])


def test_only_exact_complete_scored_pool_is_ready():
    assert scored_pool_status([dict(id='own-0')],[dict(id='own-0',score_status='scored',import_audit=dict(accepted=True))])=='ready'
