"""PDB handoff and analytic constraint controls; no scientific sampling."""
import pytest
from scripts.audit_post_h_1vcu import penalty
from scripts.research_post_h_1vcu import patch_prepared


def source():
    names = ['CA', 'N', 'C', 'CB', 'HA']
    points = [[0., 0., 0.], [0.14, 0., 0.01], [-0.07, 0.121, 0.01], [-0.07, -0.121, 0.01], [0., 0., 0.11]]
    atoms = [dict(identity=dict(component=dict(chain_id='A', sequence_number='0', insertion_code=''), name=n),
                  residue_name='SER', element='H' if n == 'HA' else 'N' if n == 'N' else 'C') for n in names]
    lines = ['HEADER    preserved\n']
    for i, (name, point) in enumerate(zip(names, points, strict=True), 1):
        xyz = ''.join(f'{10*x:8.3f}' for x in point)
        lines.append(f'ATOM  {i:5d} {name:>4s} SER A   0    {xyz}  1.00  0.00           {atoms[i-1]["element"]:>2s}\n')
    lines.append('END\n')
    return ''.join(lines).encode(), dict(atoms=atoms, initial_nm=points,
        proposed_nm=points[:4]+[[0., 0., -0.11]], moved_indices=[4], bonds=[[0, i] for i in range(1, 5)])


def test_only_target_coordinate_columns_change_and_baseline_is_exact():
    raw, saved = source()
    assert patch_prepared(raw, saved, False) == raw
    altered = patch_prepared(raw, saved, True)
    before, after = raw.splitlines(keepends=True), altered.splitlines(keepends=True)
    assert before[:5] == after[:5] and before[6:] == after[6:]
    assert before[5][:30]+before[5][54:] == after[5][:30]+after[5][54:]
    assert float(after[5][46:54]) == -1.1


def test_handoff_refuses_geometry_or_identity_mismatch():
    raw, saved = source()
    saved['initial_nm'][1][0] += 0.01
    with pytest.raises(ValueError, match='handoff'):
        patch_prepared(raw, saved, True)


def test_penalty_cartesian_derivative_matches_finite_difference():
    import numpy as np

    points = np.array([[0.1, 0.2, 0.3], [0.21, 0.18, 0.33], [-0.2, 0.4, 0.5]])
    constraints = [(0, 1, 0.10), (1, 2, 0.45)]
    energy, gradient, error, allowance = penalty(points, constraints, 2000.)
    assert energy > 0 and error > 0 and np.all(allowance >= 0)
    for i in range(3):
        for j in range(3):
            left, right = points.copy(), points.copy()
            left[i, j] -= 1e-6
            right[i, j] += 1e-6
            numeric = (penalty(right, constraints, 2000.)[0]-penalty(left, constraints, 2000.)[0])/2e-6
            assert gradient[i, j] == pytest.approx(numeric, abs=1e-6)
    assert np.sum(gradient, axis=0) == pytest.approx([0., 0., 0.], abs=1e-12)


def test_zero_distance_constraint_refuses():
    with pytest.raises(ValueError, match='degenerate'):
        penalty([[0., 0., 0.], [0., 0., 0.]], [(0, 1, 0.1)], 1.)
