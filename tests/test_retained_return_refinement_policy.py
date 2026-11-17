from types import SimpleNamespace

import pytest
from scripts.refine_retained_sampler_returns import release_anchors


def test_anchor_release_requires_remaining_geometry_failure_not_admission_failure():
    assert not release_anchors(SimpleNamespace(admission_failures=(), geometry_passed=True))
    assert release_anchors(SimpleNamespace(admission_failures=(), geometry_passed=False))
    with pytest.raises(ValueError, match="admission failure"):
        release_anchors(SimpleNamespace(admission_failures=("observations-changed",), geometry_passed=False))
