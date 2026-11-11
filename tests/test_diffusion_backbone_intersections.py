"""Synthetic diagnostic controls, not a calibrated threading benchmark."""
from dataclasses import replace

import numpy as np
import pytest

from dvbfixer.model.diffusion.backbone_intersections import (
    BackboneSegment,
    screen_segments,
    segment_distance,
)
from dvbfixer.model.diffusion.contract import AtomIdentity


def segment(chain, start, stop, insertion=""):
    return BackboneSegment(AtomIdentity(chain, "1", insertion, "N"), AtomIdentity(chain, "1", insertion, "CA"), start, stop)


def test_true_3d_intersection_and_source_ownership():
    a = segment("A", (-1, 0, 0), (1, 0, 0))
    b = segment("a", (0, -1, 0), (0, 1, 0), "B")
    result = screen_segments((a, b), frozenset((a.first,)))
    finding, = result.alarms
    assert finding.category == "3d-interior-intersection"
    assert finding.ownership == "generated-observed"
    assert finding.first_fraction == finding.second_fraction == .5
    assert not finding.topological_defect_established
    assert not result.threading_tested and not result.publication_gate
    assert screen_segments((a, b), frozenset()).alarms[0].ownership == "source-only"


def test_crossing_in_projection_is_not_a_3d_intersection():
    a = segment("A", (-1, 0, 0), (1, 0, 0))
    b = segment("B", (0, -1, 2), (0, 1, 2))
    assert not screen_segments((a, b), frozenset()).alarms
    close = replace(b, start=(0, -1, .3), stop=(0, 1, .3))
    alarm, = screen_segments((a, close), frozenset()).alarms
    assert alarm.category == "3d-near-contact"
    assert alarm.distance_angstrom == pytest.approx(.3)


def test_parallel_and_endpoint_minima_are_supported():
    a = segment("A", (0, 0, 0), (1, 0, 0))
    b = segment("B", (2, .2, 0), (3, .2, 0))
    assert segment_distance(a, b)[0] == pytest.approx(np.sqrt(1.04))
    c = replace(b, start=(.5, .2, 0), stop=(1.5, .2, 0))
    assert segment_distance(a, c)[0] == pytest.approx(.2)


def test_rigid_transform_invariance():
    a = segment("A", (-1, 0, 0), (1, 0, 0))
    b = segment("B", (0, -1, .2), (0, 1, .2))
    rotation, _ = np.linalg.qr(np.random.default_rng(7).normal(size=(3, 3)))
    def transform(s):
        return replace(s, start=tuple(np.asarray(s.start) @ rotation + 19), stop=tuple(np.asarray(s.stop) @ rotation + 19))
    original, = screen_segments((a, b), frozenset()).alarms
    moved, = screen_segments((transform(a), transform(b)), frozenset()).alarms
    assert original.distance_angstrom == pytest.approx(moved.distance_angstrom, abs=1e-12)
    assert original.first_fraction == pytest.approx(moved.first_fraction)


def test_only_bonded_neighbors_are_excluded_not_entire_chains():
    a = segment("A", (-1, 0, 0), (1, 0, 0))
    b = segment("A", (0, -1, 0), (0, 1, 0), "B")
    assert len(screen_segments((a, b), frozenset()).alarms) == 1
    explicit_link = frozenset((frozenset((a.second, b.first)),))
    assert not screen_segments((a, b), frozenset(), bonded_atoms=explicit_link).alarms


def test_degenerate_and_inconsistent_segments_fail_closed():
    with pytest.raises(ValueError, match="zero-length"):
        segment("A", (0, 0, 0), (0, 0, 0))
    a = segment("A", (0, 0, 0), (1, 0, 0))
    with pytest.raises(ValueError, match="inconsistent"):
        screen_segments((a, replace(a, start=(0, 1, 0))), frozenset())


def test_budget_failure_is_not_a_clean_result():
    segments = tuple(segment(str(i), (0, i*.01, 0), (1, i*.01, 0)) for i in range(3))
    with pytest.raises(ValueError, match="candidate budget"):
        screen_segments(segments, frozenset(), max_candidate_pairs=1)
