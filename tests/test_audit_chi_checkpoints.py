"""Independent norm bounds preserve uncertainty instead of certifying a vector."""
import pytest
from scripts.audit_chi_checkpoints import checkpoint_bounds


def test_small_l2_norm_bounds_every_component():
    row=checkpoint_bounds(.9e-4,2,1e-4)
    assert row['certainly_at_or_below_gmax'] and not row['certainly_above_gmax']


def test_intermediate_norm_does_not_determine_max_component():
    row=checkpoint_bounds(1.1e-4,2,1e-4)
    assert not row['certainly_at_or_below_gmax'] and not row['certainly_above_gmax']


def test_large_norm_guarantees_some_component_is_above_cutoff():
    assert checkpoint_bounds(2e-4,2,1e-4)['certainly_above_gmax']
    with pytest.raises(ValueError,match='invalid'):
        checkpoint_bounds(float('nan'),1,1e-4)
