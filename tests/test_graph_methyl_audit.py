"""Independent grid, HA and H-contact descriptors are not acceptance gates."""
import math

import pytest
from scripts.audit_graph_methyl_control import (
    h_close_pairs,
    ha_geometry,
    independently_place,
    select_grid,
)


def test_independent_tetrahedral_frame_preserves_radii_and_signed_phase():
    points = [[0., 0., 0.], [0., 0., 1.], [1., 0., 2.], [.01, .01, -.1], [.02, .01, -.11], [.01, .02, -.09]]
    group = dict(parent=0, axis=1, frame=2, hydrogens=[3, 4, 5])
    for phase in range(0, 360, 15):
        hs = independently_place(points, group, phase)
        for old, actual in zip(points[3:], hs, strict=True):
            assert math.hypot(*old) == pytest.approx(math.hypot(*actual), abs=1e-12)
            assert actual[2]/math.hypot(*actual) == pytest.approx(-1/3, abs=1e-12)
        signed = math.degrees(math.atan2(hs[0][1], hs[0][0]))
        assert abs((signed-phase+180)%360-180) < 1e-10


def test_independent_grid_ties_use_shortest_signed_phase():
    grid = [dict(phase=a, energy_kj_mol=100.) for a in range(0, 360, 15)]
    grid[1]['energy_kj_mol'] = grid[-1]['energy_kj_mol'] = 1.
    assert select_grid(grid) == 345
    grid[0]['energy_kj_mol'] = 1.+1e-9
    assert select_grid(grid) == 0
    with pytest.raises(ValueError):
        select_grid(grid[1:])
    grid[0]['energy_kj_mol'] = math.nan
    with pytest.raises(ValueError):
        select_grid(grid)


def test_HA_descriptor_uses_actual_bonds_not_residue_name_guesses():
    rows = [dict(key=['a', '82', 'A', 'GLN', n], element=e) for n, e in
        [('CA', 'C'), ('N', 'N'), ('C', 'C'), ('CB', 'C'), ('HA', 'H')]]
    points = [[0., 0., 0.], [.1, 0., 0.], [0., .1, 0.], [0., 0., .1], [.03, .03, .03]]
    bonds = [[0, i] for i in range(1, 5)]
    result = ha_geometry(rows, bonds, points)
    assert result['HA_count'] == 1 and len(result['all_acute']) == 1
    assert result['all_acute'][0]['atom'] == ['a', '82', 'A', 'GLN', 'HA']
    points[4] = [-.03, -.03, -.03]
    assert not ha_geometry(rows, bonds, points)['all_acute']


def test_close_H_descriptor_excludes_bond_and_angle_pairs_not_other_contacts():
    pytest.importorskip('scipy')
    rows = [dict(key=['a', '1', '', 'ALA', str(i)], element=e) for i, e in enumerate(['C', 'H', 'H', 'H'])]
    points = [[0., 0., 0.], [.01, 0., 0.], [-.01, 0., 0.], [.04, 0., 0.]]
    result = h_close_pairs(rows, [[0, 1], [0, 2]], points)
    assert result['below_1_angstrom'] == 3
    assert all('3' in [k[-1] for k in p['atoms']] for p in result['closest_pairs'])
