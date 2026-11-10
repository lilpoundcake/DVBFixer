"""Independent outward/tangential frame and static force/gradient controls."""
import math

import pytest
from scripts.audit_ha_minimizer_path import (
    force_readback,
    outward_frame,
    projection,
    sphere_tangent,
)


def template():
    radius = math.sqrt(8/9)
    points = [[0., 0., 0.], [radius, 0., 1/3], [-radius/2, radius*math.sqrt(3)/2, 1/3],
        [-radius/2, -radius*math.sqrt(3)/2, 1/3], [.1, 0., -.1]]
    return [dict(key=['A', '16', '', 'GLN', name], element=name[0], xyz=point)
            for name, point in zip(('CA', 'N', 'C', 'CB', 'HA'), points, strict=True)]


def test_outward_and_sphere_tangent_are_geometric_and_translation_invariant():
    rows = template()
    centre, outward = outward_frame(rows)
    assert outward == pytest.approx([0., 0., -1.])
    point = rows[-1]['xyz']
    cosine, tangent = sphere_tangent(centre, point, outward)
    assert cosine == pytest.approx(1/math.sqrt(2))
    assert math.hypot(*tangent) == pytest.approx(1.)
    assert projection(tangent, point) == pytest.approx(0, abs=1e-12)
    moved = [dict(r, xyz=[x+3 for x in r['xyz']]) for r in rows]
    new_centre, new_outward = outward_frame(moved)
    new_cosine, new_tangent = sphere_tangent(new_centre, moved[-1]['xyz'], new_outward)
    assert new_cosine == pytest.approx(cosine)
    assert new_tangent == pytest.approx(tangent)


def test_degenerate_sphere_tangent_is_not_assigned_an_arbitrary_direction():
    assert sphere_tangent([0, 0, 0], [0, 0, -1], [0, 0, -1]) == (1., None)
    assert projection([1, 2, 3], None) is None
    with pytest.raises(ValueError):
        sphere_tangent([0, 0, 0], [0, 0, 0], [0, 0, -1])


@pytest.mark.parametrize('fault', [None, 'gradient', 'constraint', 'energy'])
def test_native_gradient_sign_force_sum_and_energy_checks(fault):
    openmm = pytest.importorskip('openmm')
    rows = template()
    system = openmm.System()
    for _ in range(4):
        system.addParticle(0)
    system.addParticle(1)
    for expression in ('2*x', '3*z'):
        force = openmm.CustomExternalForce(expression)
        force.addParticle(4, [])
        system.addForce(force)
    gradient = [0.]*15
    gradient[12:15] = [2., 0., 3.]
    stats = {'system energy': -.01, 'restraint energy': 0., 'max constraint error': 0.}
    if fault == 'gradient':
        gradient[-1] = -3.
    elif fault == 'constraint':
        stats['restraint energy'] = 1.
    elif fault == 'energy':
        stats['system energy'] = 10.
    node = dict(label='callback', ordinal=0, iteration=0,
        points_nm=[[x/10 for x in row['xyz']] for row in rows],
        callback=dict(gradient_kj_mol_nm=gradient, statistics=stats))
    xml = openmm.XmlSerializer.serialize(system)
    if fault:
        with pytest.raises(ValueError):
            force_readback(xml, rows, [node])
    else:
        result = force_readback(xml, rows, [node])[0]
        assert result['total_HA_force_kj_mol_nm'] == pytest.approx([-2., 0., -3.])
        assert result['outward_projection_kj_mol_nm'] == pytest.approx(3.)
        assert result['native_gradient_maximum_error_kj_mol_nm'] == pytest.approx(0.)
        assert len(result['force_components']) == 2
