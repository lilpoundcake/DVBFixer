import json
from dataclasses import asdict, replace

import numpy as np
import pytest

from dvbfixer.model.diffusion.backbone_intersections import BackboneSegment
from dvbfixer.model.diffusion.contract import AtomIdentity
from dvbfixer.model.diffusion.loop_surface_probe import probe_loop_surface

SQUARE = np.array([[-1., -1, 0], [1, -1, 0], [1, 1, 0], [-1, 1, 0]])


def path(start, stop):
    return BackboneSegment(AtomIdentity("a", "3", "B", "N"), AtomIdentity("a", "3", "B", "CA"), start, stop)


def test_surface_piercing_is_not_a_proven_threading_defect():
    result = probe_loop_surface(SQUARE, (path((.2, .1, -1), (.2, .1, 1)),))
    hit, = result.piercings
    assert hit.fraction == pytest.approx(.5)
    assert not hit.topological_defect_established
    assert result.closure_ambiguous and result.open_path_end_sensitive
    assert not result.publication_gate and not result.disulfide_topology_tested
    json.dumps(asdict(result))


def test_projected_crossing_without_surface_piercing_and_open_end_sensitivity():
    assert not probe_loop_surface(SQUARE, (path((-.5, 0, 1), (.5, 0, 1)),)).piercings
    assert not probe_loop_surface(SQUARE, (path((.2, .1, .2), (.2, .1, 1)),)).piercings


def test_fan_edge_hits_are_deduplicated_and_marked_ambiguous():
    result = probe_loop_surface(SQUARE, (path((0, 0, -1), (0, 0, 1)),))
    assert len(result.piercings) == 1
    assert result.piercings[0].surface_edge_ambiguous


def test_surface_probe_rigid_invariance():
    segment = path((.2, .1, -1), (.2, .1, 1))
    matrix, _ = np.linalg.qr(np.random.default_rng(7).normal(size=(3, 3)))
    moved = replace(segment, start=tuple(np.asarray(segment.start) @ matrix + 5), stop=tuple(np.asarray(segment.stop) @ matrix + 5))
    original, = probe_loop_surface(SQUARE, (segment,)).piercings
    transformed, = probe_loop_surface(SQUARE @ matrix + 5, (moved,)).piercings
    assert original.fraction == pytest.approx(transformed.fraction)
    assert np.asarray(original.point) @ matrix + 5 == pytest.approx(transformed.point)


def test_surface_probe_refuses_degeneracy_and_unbounded_work():
    with pytest.raises(ValueError, match="degenerate"):
        probe_loop_surface([[0, 0, 0], [0, 0, 0], [1, 0, 0]], ())
    with pytest.raises(ValueError, match="budget"):
        probe_loop_surface(SQUARE, (path((0, 0, -1), (0, 0, 1)),), max_tests=1)
