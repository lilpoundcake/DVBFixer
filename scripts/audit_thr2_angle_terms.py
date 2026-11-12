"""All three retained whole starts: analytic per-angle Cbeta force attribution.

Saved-only analysis; no optimizer, force changes, H placement or sampling.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from scripts.audit_stereo_force_projection import grouped_system, projected_force, volume_gradient
from scripts.audit_thr2_native_capture import read, verify_files
from scripts.research_dunbrack_prior import sha, write
from scripts.research_thr2_native_capture import check


def angle_forces(points, theta0, strength):
    """Cartesian forces for E=0.5*k*(theta-theta0)^2; nm/radians/kJ/mol."""
    u = [x-y for x, y in zip(points[0], points[1], strict=True)]
    v = [x-y for x, y in zip(points[2], points[1], strict=True)]
    a, b = math.sqrt(sum(x*x for x in u)), math.sqrt(sum(x*x for x in v))
    if min(a, b) <= 1e-12:
        raise ValueError('undefined zero-length angle arm')
    cosine = max(-1., min(1., sum(x*y for x, y in zip(u, v, strict=True))/(a*b)))
    sine = math.sqrt(max(0., 1-cosine*cosine))
    if sine <= 1e-10:
        raise ValueError('singular analytic angle derivative; retain failure')
    theta = math.acos(cosine)
    scale = strength*(theta-theta0)/sine
    first = [scale*(v[i]/(a*b)-cosine*u[i]/a**2) for i in range(3)]
    third = [scale*(u[i]/(a*b)-cosine*v[i]/b**2) for i in range(3)]
    middle = [-x-y for x, y in zip(first, third, strict=True)]
    return theta, [first, middle, third]


def audit(plan_path, root, output):
    import openmm
    from openmm.unit import kilojoules_per_mole, nanometer, radian

    parent = read(plan_path)
    check(parent)
    execution = read(root/'evidence.json')
    if execution['plan_sha256'] != sha(plan_path.read_bytes()) or [r['label'] for r in execution['native']] != parent['protocol']['labels']:
        raise ValueError('complete three-start enrollment differs')
    reports, inputs = [], {}
    for result in execution['native']:
        source = root/result['label']
        verify_files(source, result)
        inputs.update({str(source/name):digest for name, digest in result['artifacts'].items()})
        call = next(c for c in read(source/'capture.json')['calls'] if c['stage'] == 'whole')
        prefix = call['before']['prefix']
        atoms = read(source/f'{prefix}-atoms.json.gz')['atoms']
        keys = {tuple(row['key']):i for i, row in enumerate(atoms)}
        ids = [keys[('A', '2', '', 'THR', name)] for name in ('CB', 'CA', 'OG1', 'CG2')]
        state = openmm.XmlSerializer.deserialize((source/f'{prefix}-state.xml').read_text())
        xyz = state.getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
        volume, gradient, length = volume_gradient([xyz[i] for i in ids])
        system = openmm.XmlSerializer.deserialize((source/f'{prefix}-system.xml').read_text())
        grouped, groups = grouped_system(system)
        native = [f for f in system.getForces() if isinstance(f, openmm.HarmonicAngleForce)]
        if len(native) != 1:
            raise ValueError('exactly one native harmonic angle force required')
        force = native[0]
        terms, total_forces = [], [[0., 0., 0.] for _ in atoms]
        for index in range(force.getNumAngles()):
            a, b, c, theta0, k = force.getAngleParameters(index)
            indices = [a, b, c]
            theta, forces = angle_forces([xyz[i] for i in indices], theta0.value_in_unit(radian),
                k.value_in_unit(kilojoules_per_mole/radian**2))
            for i, vector in zip(indices, forces, strict=True):
                for j in range(3):
                    total_forces[i][j] += vector[j]
            projection = sum(sum(x*y for x, y in zip(vector, gradient[ids.index(i)], strict=True))
                for i, vector in zip(indices, forces, strict=True) if i in ids)/length
            if set(ids).intersection(indices):
                terms.append(dict(index=index, atoms=[atoms[i]['key'] for i in indices],
                    hydrogens=[atoms[i]['key'] for i in indices if atoms[i]['element'] in ('H', 'D')],
                    theta_degrees=math.degrees(theta), equilibrium_degrees=math.degrees(theta0.value_in_unit(radian)),
                    cb_projection_kj_mol_nm=projection))
        integrator = openmm.VerletIntegrator(0.)
        context = openmm.Context(grouped, integrator, openmm.Platform.getPlatformByName('Reference'))
        try:
            context.setState(state)
            group = next(g['index'] for g in groups if g['type'] == 'HarmonicAngleForce')
            forces = context.getState(forces=True, groups=1 << group).getForces(asNumpy=True).value_in_unit(kilojoules_per_mole/nanometer)
            error = max(abs(forces[i][j]-total_forces[i][j]) for i in range(len(atoms)) for j in range(3))
            if error > 1e-7:
                raise ValueError(f'analytic force fails independent native readback: {error}')
            native_projection = projected_force(forces, ids, gradient, length)
        finally:
            del context
        by_h = {}
        for term in terms:
            # Multi-H methyl angles are one joint term, never counted twice.
            key = '+'.join('/'.join(k) for k in sorted(term['hydrogens'])) if term['hydrogens'] else 'heavy-only'
            by_h[key] = by_h.get(key, 0.)+term['cb_projection_kj_mol_nm']
        reports.append(dict(label=result['label'], cb_volume_angstrom3=1000*volume,
            native_angle_projection_kj_mol_nm=native_projection,
            analytic_angle_projection_kj_mol_nm=sum(t['cb_projection_kj_mol_nm'] for t in terms),
            maximum_cartesian_force_error_kj_mol_nm=float(error), by_hydrogen=by_h, terms=terms))
    write(output, dict(schema='dvbfixer.thr2-angle-terms.v1', plan_sha256=sha(plan_path.read_bytes()),
        execution_sha256=sha((root/'evidence.json').read_bytes()), inputs=inputs,
        audit_script_sha256=sha(Path(__file__).read_bytes()), starts=reports,
        extra_minimization=False, production_changed=False,
        limits='initial local harmonic-angle force decomposition, not an entire-trajectory causal intervention'))
    print(json.dumps([{k:v for k, v in r.items() if k != 'terms'} for r in reports], indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('plan', 'root', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan.resolve(), args.root.resolve(), args.output.resolve())
