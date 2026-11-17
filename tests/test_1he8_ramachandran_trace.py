import numpy as np
from scripts import analyze_1he8_ramachandran_trace as trace


def test_undefined_dihedral_is_not_reported_as_outlier(monkeypatch):
    def forbidden(*args):
        raise AssertionError("undefined angles must not reach reference lookup")

    monkeypatch.setattr(trace, "_ramachandran_reference_populated", forbidden)
    result = trace.angles(np.zeros((5, 3)))
    assert not result["defined"]
    assert not result["gross_outlier"]


def test_angles_reuse_canonical_reference_and_are_scale_invariant(monkeypatch):
    calls = []
    monkeypatch.setattr(trace, "_ramachandran_reference_populated", lambda phi, psi: calls.append((phi, psi)) or True)
    points = np.array([[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [1., 1., 1.], [2., 1., 1.]])
    first = trace.angles(points)
    second = trace.angles(points*10)
    assert first == second
    assert first["defined"] and not first["gross_outlier"]
    assert len(calls) == 2
