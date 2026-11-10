"""Static force readback math is separate from optimizer paths and chirality repair."""
import math

import pytest
from scripts.audit_stereo_force_projection import (
    focal_centres,
    grouped_system,
    guard_parameters,
    projected_force,
    resolve,
    volume_gradient,
)


def test_signed_volume_gradient_matches_all_cartesian_finite_differences():
    points = [[0.2, 0.3, 0.4], [1.0, 0.2, 0.1], [-0.3, 0.8, 0.2], [0.1, -0.4, 0.9]]
    _, gradient, norm = volume_gradient(points)
    assert norm > 0
    for i in range(4):
        for j in range(3):
            positive = [p[:] for p in points]
            negative = [p[:] for p in points]
            positive[i][j] += 1e-6
            negative[i][j] -= 1e-6
            finite = (volume_gradient(positive)[0]-volume_gradient(negative)[0]) / 2e-6
            assert finite == pytest.approx(gradient[i][j], abs=1e-9)
    assert all(sum(g[j] for g in gradient) == pytest.approx(0) for j in range(3))
    assert projected_force(gradient, [0, 1, 2, 3], gradient, norm) == pytest.approx(norm)
    assert projected_force([[-x for x in g] for g in gradient], [0, 1, 2, 3], gradient, norm) == pytest.approx(-norm)


@pytest.mark.parametrize('points', [
    [[0, 0, 0]]*4,
    [[math.nan, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
    [[0, 0]]*4,
])
def test_undefined_gradient_is_refused(points):
    with pytest.raises(ValueError):
        volume_gradient(points)


def test_val_has_no_beta_centre_and_identity_resolution_is_exact():
    case = dict(target=dict(chain_id='d', sequence_number='82', insertion_code='A'), residue='VAL')
    centres = focal_centres(case)
    assert [c['label'] for c in centres] == ['target_ca', 'context_ca']
    rows = [dict(key=['d', '82', 'A', 'VAL', name]) for name in centres[0]['names']]
    assert resolve(rows, centres[0]) == [0, 1, 2, 3]
    with pytest.raises(ValueError, match='duplicate'):
        resolve(rows+rows[:1], centres[0])
    with pytest.raises(ValueError, match='missing'):
        resolve([dict(key=['D', *r['key'][1:]]) for r in rows], centres[0])


def test_force_group_assignment_changes_no_coefficients_or_original_system():
    openmm = pytest.importorskip('openmm')
    system = openmm.System()
    system.addParticle(12)
    force = openmm.CustomExternalForce('0.5*k*x*x')
    force.addGlobalParameter('k', 2)
    force.addParticle(0, [])
    force.setForceGroup(7)
    system.addForce(force)
    original = openmm.XmlSerializer.serialize(system)
    grouped, groups = grouped_system(system)
    assert groups[0]['original_group'] == 7
    assert grouped.getForce(0).getForceGroup() == 0
    assert openmm.XmlSerializer.serialize(system) == original


def test_canonical_ca_guard_projection_matches_analytic_force():
    openmm = pytest.importorskip('openmm')
    from openmm.app import Element, Topology
    from openmm.unit import Quantity, kilojoules_per_mole, nanometer

    from dvbfixer.ffutils.geometry import build_ca_chirality_force

    top = Topology()
    residue = top.addResidue('ALA', top.addChain('d'), id='82', insertionCode='A')
    for name in ('N', 'CA', 'C', 'CB'):
        top.addAtom(name, Element.getBySymbol(name[0]), residue)
    initial = [(0.1, 0, 0), (0, 0, 0), (0, 0.1, 0), (0, 0, 0.1)]
    force, count = build_ca_chirality_force(top, Quantity(initial, nanometer))
    assert count == 1
    system = openmm.System()
    for _ in initial:
        system.addParticle(12)
    system.addForce(force)
    grouped, _ = grouped_system(system)
    indices = [1, 0, 2, 3]
    parameters = guard_parameters(grouped, indices)
    assert parameters is not None
    points = [*initial[:3], (0, 0, 0.01)]
    volume, gradient, norm = volume_gradient([points[i] for i in indices])
    integrator = openmm.VerletIntegrator(0.001)
    context = openmm.Context(grouped, integrator, openmm.Platform.getPlatformByName('Reference'))
    try:
        context.setPositions(points*nanometer)
        forces = context.getState(forces=True).getForces().value_in_unit(kilojoules_per_mole/nanometer)
        observed = projected_force(forces, indices, gradient, norm)
        expected = parameters['k'] * (parameters['qmin']-volume/parameters['vref_nm3']) / parameters['vref_nm3'] * norm
        assert observed > 0
        assert observed == pytest.approx(expected, rel=1e-10)
    finally:
        del context
        del integrator
