"""Force-sum auditing rejects dimensional, finite and native-gradient drift."""
import math

import pytest
from scripts.audit_graph_methyl_attribution import check_vector


def test_native_vector_roundoff_is_bounded_but_not_hidden():
    assert check_vector([1., 2., 3.], [1.+1e-10, 2., 3.]) == pytest.approx(1e-10)


@pytest.mark.parametrize('bad', [[1., 2.], [1., 2., math.nan], [1., 2., 4.]])
def test_incomplete_nonfinite_or_wrong_native_closure_is_refused(bad):
    with pytest.raises(ValueError):
        check_vector(bad, [1., 2., 3.])
