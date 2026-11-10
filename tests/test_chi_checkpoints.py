"""Native L2 checkpoint norms are not silently treated as maximum components."""
import math

import pytest
from scripts.research_chi_checkpoints import norm_bounds


def test_one_dimensional_norm_is_exact_absolute_component():
    assert norm_bounds(.01,1)==dict(maximum_component_lower_bound=.01,maximum_component_upper_bound=.01)


def test_two_dimensional_norm_only_bounds_component_maximum():
    result=norm_bounds(.01,2)
    assert result['maximum_component_lower_bound']==pytest.approx(.01/math.sqrt(2))
    assert result['maximum_component_upper_bound']==.01


@pytest.mark.parametrize('norm,dimensions',[(-1.,1),(float('nan'),1),(1.,0)])
def test_invalid_diagnostic_is_not_convergence_evidence(norm,dimensions):
    with pytest.raises(ValueError,match='norm'):
        norm_bounds(norm,dimensions)
