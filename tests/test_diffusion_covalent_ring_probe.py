import json
from dataclasses import asdict

import numpy as np
import pytest

from dvbfixer.model.diffusion.contract import AtomIdentity
from dvbfixer.model.diffusion.covalent_ring_probe import probe_covalent_ring


def fixture():
    ring = tuple(AtomIdentity('A', '972', 'B', f'C{i}') for i in range(6))
    points = {a: (float(np.cos(i*np.pi/3)), float(np.sin(i*np.pi/3)), 0.) for i,a in enumerate(ring)}
    pair = (AtomIdentity('a', '1044', '', 'CA'), AtomIdentity('a', '1044', '', 'CB'))
    points.update({pair[0]: (.2, .1, -1.), pair[1]: (.2, .1, 1.)})
    edges = frozenset([frozenset((a,b)) for a,b in zip(ring, ring[1:]+ring[:1])]+[frozenset(pair)])
    return ring, pair, points, edges


def test_real_closed_ring_and_covalent_bond_piercing_is_nonfatal():
    ring, pair, points, edges = fixture()
    result = probe_covalent_ring(ring, points, edges, frozenset(ring+pair))
    finding, = result.piercings
    assert finding.bond == pair
    assert finding.ownership == 'generated-generated'
    assert finding.fraction == pytest.approx(.5)
    assert not result.publication_gate and not result.topology_or_trajectory_qualified
    json.dumps(asdict(result))


def test_topology_is_required_and_nonplanar_result_is_not_a_clean_pass():
    ring, pair, points, edges = fixture()
    with pytest.raises(ValueError, match='closed cycle'):
        probe_covalent_ring(ring, points, edges-{frozenset((ring[0], ring[1]))}, frozenset())
    points[ring[0]] = (1., 0., 2.)
    result = probe_covalent_ring(ring, points, edges, frozenset())
    assert result.undefined_reason is not None
    assert not result.piercings


def test_near_contact_without_piercing_and_local_bonded_neighbors_are_excluded():
    ring, pair, points, edges = fixture()
    points[pair[1]] = (.2, .1, -.2)
    assert not probe_covalent_ring(ring, points, edges, frozenset()).piercings
    points[pair[1]] = (.2, .1, 1.)
    edges |= {frozenset((ring[0], pair[0]))}
    assert not probe_covalent_ring(ring, points, edges, frozenset()).piercings


def test_ring_probe_is_rigid_invariant_and_preserves_case_and_insertion_identity():
    ring, pair, points, edges = fixture()
    matrix, _ = np.linalg.qr(np.random.default_rng(7).normal(size=(3,3)))
    moved = {a: tuple(np.asarray(xyz) @ matrix+5) for a, xyz in points.items()}
    first = probe_covalent_ring(ring, points, edges, frozenset(ring))
    second = probe_covalent_ring(ring, moved, edges, frozenset(ring))
    assert len(first.piercings) == len(second.piercings) == 1
    assert first.piercings[0].fraction == pytest.approx(second.piercings[0].fraction)
    assert second.piercings[0].bond[0].chain == 'a'
    assert all(a.insertion_code == 'B' for a in second.ring)


def test_invalid_coordinates_and_work_budgets_fail_closed():
    ring, pair, points, edges = fixture()
    points[pair[0]] = (0., 0., float('nan'))
    with pytest.raises(ValueError, match='finite'):
        probe_covalent_ring(ring, points, edges, frozenset())
    _, _, points, _ = fixture()
    with pytest.raises(ValueError, match='budget'):
        probe_covalent_ring(ring, points, edges, frozenset(), max_tests=1)
