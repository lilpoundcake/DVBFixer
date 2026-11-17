import json

import numpy as np
from openmm import (
    Context,
    CustomExternalForce,
    LocalEnergyMinimizer,
    Platform,
    System,
    VerletIntegrator,
)
from openmm.unit import nanometer
from scripts.trace_private_minimization import PassiveReporter


def test_passive_callback_does_not_change_reference_minimum(tmp_path):
    system = System()
    system.addParticle(12)
    force = CustomExternalForce("k*(x*x+y*y+z*z)")
    force.addGlobalParameter("k", 1000)
    force.addParticle(0)
    system.addForce(force)
    results = []
    reporter = PassiveReporter(tmp_path / "trace")
    for callback in (None, reporter):
        integrator = VerletIntegrator(.001)
        context = Context(system, integrator, Platform.getPlatformByName("Reference"))
        context.setPositions(np.array([[1., 2., 3.]])*nanometer)
        LocalEnergyMinimizer.minimize(context, 1., 100, callback)
        results.append(context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(nanometer))
        del context, integrator
    np.testing.assert_array_equal(results[0], results[1])
    assert reporter.calls > 0
    rows = [json.loads(line) for line in (tmp_path / "trace/callbacks.jsonl").read_text().splitlines()]
    assert all(not row["requested_termination"] for row in rows)
    assert all(row["coordinates_finite"] and row["gradient_finite"] for row in rows)
