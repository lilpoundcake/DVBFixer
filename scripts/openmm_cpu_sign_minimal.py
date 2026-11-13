"""Portable synthetic public-API reproducer, no dvbfixer/private inputs required.

Default: describe only (no Context/minimization). Explicit execution:
    python openmm_cpu_sign_minimal.py --run-native --platform OpenCL

The compact wrapper was packaged AFTER the frozen factorial study, not counted
as an additional executed native trial. It neither patches nor installs OpenMM.
"""
from __future__ import annotations

import argparse
import json
import math
import struct


def make_system():
    import openmm

    system = openmm.System()
    for _ in range(3):
        system.addParticle(12.)
    bond = openmm.HarmonicBondForce()
    bond.addBond(0, 1, 0.12, 10000.)
    system.addForce(bond)
    force = openmm.CustomExternalForce('0.5*ks*(x*x+3*y*y+z*z)')
    force.addGlobalParameter('ks', 5e7)
    force.addParticle(2, [])
    system.addForce(force)
    system.addConstraint(0, 1, 0.1)
    return system


def measure(x, native, statistics, use_opencl_allowance=False):
    """Scalar-only analysis of synthetic points; correct sign is not redefined."""
    p = [x[i:i+3] for i in (0, 3, 6)]
    delta = [b-a for a, b in zip(p[0], p[1], strict=True)]
    r = math.sqrt(sum(v*v for v in delta))
    if r <= 1e-15:
        raise ValueError('degenerate pair')
    strength = statistics['restraint strength']
    physical = [-10000.*(r-0.12)*v/r for v in delta]+[10000.*(r-0.12)*v/r for v in delta]
    physical += [5e7*p[2][i]*scale for i, scale in enumerate((1., 3., 1.))]
    gp = [-strength*(r-0.1)*v/r for v in delta]+[strength*(r-0.1)*v/r for v in delta]+[0., 0., 0.]
    correct = [a+b for a, b in zip(physical, gp, strict=True)]
    source = [a-b for a, b in zip(physical, gp, strict=True)]
    source32 = [struct.unpack('f', struct.pack('f', v))[0] for v in source]
    allowance = [32*2**-23*abs(strength)*(r+0.1)]*6+[0., 0., 0.] if use_opencl_allowance else [0.]*9

    def comparison(expected):
        errors = [abs(a-b) for a, b in zip(expected, native, strict=True)]
        bounds = [1e-3+2e-5*max(1., abs(a), abs(b)) for a, b in zip(expected, native, strict=True)]
        return dict(strict_passed=all(e <= b for e, b in zip(errors, bounds, strict=True)),
            original_constraint_allowance_passed=all(e <= b+c for e, b, c in zip(errors, bounds, allowance, strict=True)),
            maximum_error=max(errors))

    return dict(physical=comparison(correct), source_sign=comparison(source),
        source_sign_float32=comparison(source32),
        particle0_x_physical=correct[0], particle0_x_source_sign=source[0],
        particle0_x_reported=native[0], maximum_reported_gradient=max(abs(v) for v in native),
        warning='Matching an incorrect source sign is NOT a physically correct gradient.')


def run(backend):
    import openmm
    import openmm.version

    integrator = openmm.VerletIntegrator(0.001)
    integrator.setConstraintTolerance(1e-5)
    selected = openmm.Platform.getPlatformByName(backend)
    properties = dict(DeviceIndex='0', OpenCLPlatformIndex='0', Precision='single',
        UseCpuPme='false', DisablePmeStream='false') if backend == 'OpenCL' else {}
    system = make_system()
    context = openmm.Context(system, integrator, selected, properties)
    context.setPositions([openmm.Vec3(-0.0500025, 0., 0.), openmm.Vec3(0.0500025, 0., 0.), openmm.Vec3(10., 10., 0.)])

    class Reporter(openmm.MinimizationReporter):
        def __init__(self):
            super().__init__()
            self.records = []

        def report(self, iteration, x, grad, args):
            if len(self.records) >= 256:
                raise RuntimeError('diagnostic callback budget exceeded')
            self.records.append(dict(iteration=iteration, **measure(list(x), list(grad), dict(args), backend == 'OpenCL')))
            return False

    reporter, returned = Reporter(), False
    try:
        openmm.LocalEnergyMinimizer.minimize(context, 10., 1, reporter)
        returned = True
    finally:
        print(json.dumps(dict(runtime_version=openmm.version.version, revision=openmm.version.git_revision,
            platform=backend, returned=returned, callbacks=reporter.records,
            limits='Synthetic diagnostic, not convergence, source/binary attestation or a dependency fix.'), indent=2, allow_nan=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-native', action='store_true')
    parser.add_argument('--platform', choices=('OpenCL', 'CPU', 'Reference'), default='OpenCL')
    args = parser.parse_args()
    if args.run_native:
        run(args.platform)
    else:
        print('Description only; no Context or minimizer call. Use --run-native to execute one synthetic diagnostic.')
