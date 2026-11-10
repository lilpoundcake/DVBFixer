"""Prospective enrollment fixes sources and reuses whole-mask admission rules."""
from scripts.prospective_geometry_cohort import SYSTEMS
from scripts.research_start_geometry import ARMS, SEEDS


def test_prospective_sources_and_controls_are_not_replaced_after_admission():
    assert list(SYSTEMS)==['1FNA','1SHG','1VII']
    assert ARMS==('raw','beta_only','full_generated')
    assert SEEDS==(1,2)
    assert 3*4*len(ARMS)*len(SEEDS)==72
