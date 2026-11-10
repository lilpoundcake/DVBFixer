"""Read-only native iteration reporting, without changing the optimizer result."""
import math

import pytest
from scripts.research_ha_minimizer_path import reporter_type


def test_recorder_copies_values_and_never_requests_early_stop():
    openmm = pytest.importorskip('openmm')
    recorder = reporter_type(openmm)(1)
    points, gradient = [0., 1., 2.], [3., 4., 5.]
    stats = {'system energy': 7.}
    assert recorder.report(3, points, gradient, stats) is False
    points[0] = 10.
    gradient[0] = 11.
    stats['system energy'] = 20.
    row = recorder.records[0]
    assert row['positions_nm'] == [0., 1., 2.]
    assert row['gradient_kj_mol_nm'] == [3., 4., 5.]
    assert row['statistics'] == {'system energy': 7.}
    assert row['ordinal'] == 0 and row['iteration'] == 3


@pytest.mark.parametrize('fault', ['length', 'nonfinite', 'bound'])
def test_bad_or_unbounded_reports_fail_instead_of_stopping_with_success(fault):
    openmm = pytest.importorskip('openmm')
    recorder = reporter_type(openmm)(1, maximum=0 if fault == 'bound' else 1)
    points = [0., 1.] if fault == 'length' else [0., 1., math.nan] if fault == 'nonfinite' else [0., 1., 2.]
    with pytest.raises(ValueError):
        recorder.report(0, points, [0., 0., 0.], {'system energy': 1.})


def test_real_reporter_does_not_change_native_reference_minimum():
    openmm = pytest.importorskip('openmm')
    from openmm.unit import nanometer

    system = openmm.System()
    system.addParticle(0)
    system.addParticle(1)
    force = openmm.HarmonicBondForce()
    force.addBond(0, 1, .1, 1000.)
    system.addForce(force)
    finals = []
    recorder = reporter_type(openmm)(2)
    for observed in (False, True):
        integrator = openmm.VerletIntegrator(0.)
        context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('Reference'))
        try:
            context.setPositions([[0., 0., 0.], [.2, .01, 0.]]*nanometer)
            if observed:
                openmm.LocalEnergyMinimizer.minimize(context, 1., 50, recorder)
            else:
                openmm.LocalEnergyMinimizer.minimize(context, 1., 50)
            finals.append([list(p) for p in context.getState(positions=True).getPositions().value_in_unit(nanometer)])
        finally:
            del context
            del integrator
    assert finals[0] == finals[1]
    assert recorder.records
    assert all(r['positions_nm'][:3] == [0., 0., 0.] for r in recorder.records)
