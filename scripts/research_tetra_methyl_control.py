"""Private generated-only tetrahedral methyl-H initial-condition control.

ALL24 original H starts x unchanged/tetrahedral_zero/tetrahedral_energy_selected.
Same H radii and native physics. No heavy rewrite, new seed or public repair.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path

from scripts.audit_legacy_start_controls import coordinate_delta, load_atoms
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import reporter_type, verify
from scripts.research_methyl_start_control import ANGLES, methyl_indices, selected_angle

ARMS = ('unchanged', 'tetrahedral_zero', 'tetrahedral_energy_selected')


def methyl_hch_angles(rows, points):
    indices = methyl_indices(rows)
    centre = points[indices[0]]
    vectors = [[x-y for x, y in zip(points[i], centre)] for i in indices[2:]]
    if min(math.hypot(*v) for v in vectors) <= 1e-12:
        raise ValueError('undefined methyl H-C-H geometry')
    return [math.degrees(math.acos(max(-1., min(1., sum(x*y for x, y in zip(vectors[i], vectors[j]))/(math.hypot(*vectors[i])*math.hypot(*vectors[j])))))) for i, j in ((0, 1), (0, 2), (1, 2))]


def tetrahedral_methyl(rows, points, degrees):
    if degrees not in ANGLES or len(points) != len(rows):
        raise ValueError('undeclared tetrahedral phase or particle inventory')
    indices = methyl_indices(rows)
    ca = [i for i, row in enumerate(rows) if row['key'] == ['A', '24', '', 'THR', 'CA']]
    if len(ca) != 1 or any(len(p) != 3 or not all(math.isfinite(x) for x in p) for p in points):
        raise ValueError('canonical CA frame or finite coordinates missing')
    centre = points[indices[0]]
    axis = [x-y for x, y in zip(points[indices[1]], centre)]
    length = math.hypot(*axis)
    if length <= 1e-12:
        raise ValueError('undefined CG2-to-CB axis')
    axis = [x/length for x in axis]
    reference = [x-y for x, y in zip(points[ca[0]], centre)]
    axial = sum(x*y for x, y in zip(reference, axis))
    first = [x-axial*y for x, y in zip(reference, axis)]
    norm = math.hypot(*first)
    if norm <= 1e-12:
        raise ValueError('CA projection cannot define a methyl azimuth frame')
    first = [x/norm for x in first]
    second = [axis[1]*first[2]-axis[2]*first[1], axis[2]*first[0]-axis[0]*first[2], axis[0]*first[1]-axis[1]*first[0]]
    result = [list(p) for p in points]
    for ordinal, index in enumerate(indices[2:]):
        radius = math.dist(points[index], centre)
        if radius <= 1e-12:
            raise ValueError('undefined original methyl C-H radius')
        azimuth = math.radians(degrees+120*ordinal)
        radial = math.sqrt(8/9)
        direction = [-a/3+radial*(math.cos(azimuth)*u+math.sin(azimuth)*v) for a, u, v in zip(axis, first, second)]
        result[index] = [x+radius*d for x, d in zip(centre, direction)]
    return result


def plan(rigid_plan_path, rigid_run, output):
    import openmm

    parent_raw = rigid_plan_path.read_bytes()
    parent = json.loads(parent_raw)
    verify(parent)
    run_raw = (rigid_run/'evidence.json').read_bytes()
    run = json.loads(run_raw)
    if run['plan_sha256'] != sha(parent_raw) or len(run['cases']) != 48 or any(r['status'] != 'completed' for r in run['cases']):
        raise ValueError('complete preceding rigid control not bound')
    inputs = {**parent['inputs'], str(rigid_plan_path): sha(parent_raw), str(rigid_run/'evidence.json'): sha(run_raw)}
    for row in run['cases']:
        folder = rigid_run/row['id']/f"seed-{row['seed']}"/row['source_arm']/row['arm']
        for name, digest in row['artifacts'].items():
            p = folder/name
            if sha(p.read_bytes()) != digest:
                raise ValueError('preceding rigid-control artifact changed')
            inputs[str(p)] = digest
    cases = []
    for case in parent['cases']:
        rows = load_atoms(Path(case['source'])/'h-internal-before.json.gz')
        points = json.loads(gzip.decompress(Path(case['path']).read_bytes()))['initial_nm']
        for angle in ANGLES:
            tetrahedral_methyl(rows, points, angle)
        cases.append(dict(case, original_methyl_HCH_degrees=methyl_hch_angles(rows, points),
            canonical_methyl_HCH_degrees=methyl_hch_angles(rows, tetrahedral_methyl(rows, points, 0))))
    sources = {**parent['source_hashes']}
    sources[str(Path(__file__).resolve())] = sha(Path(__file__).read_bytes())
    write(output, dict(schema='dvbfixer.tetra-methyl-start-plan.v1', cases=cases, inputs=inputs,
        source_hashes=sources, python=sys.executable, openmm_version=openmm.__version__,
        policy=dict(arms=ARMS, angles_degrees=ANGLES, energy_tie_kj_mol=1e-8,
            per_cell_seconds=900, total_seconds=5400, tolerance_kj_mol_nm=1., max_iterations=50,
            platform='Reference', maximum_callbacks=512,
            timing='separate prospective control after negative rigid-rotation result, before tetrahedral H sampling',
            enrollment='ALL24 original H starts, three matched arms, no new H-placement seeds',
            intervention='only generated THR24 HG21/HG22/HG23; preserve each original C-H radius, canonical sp3 angular geometry around fixed CG2-to-CB axis, CA-projection azimuth frame',
            selection='canonical phase zero is a separate arm; selected arm uses exactly 24 predeclared phases and original initial System energy only, never final geometry/native accuracy',
            invariants='all nonmethyl starting coordinates, original H System parameters/masses and Reference/1.0/50 call unchanged; exact unchanged endpoint required',
            retention='all 72 cells and candidate grids/callbacks/failures, no retry or favorable-case filtering',
            limits='H-stage only; no whole-system minimization, force ablation, H regeneration, seed/iteration tuning, production or repair qualification'),
        production_changed=False, native_coordinates_read=False))


def worker(plan_path, index, arm, output):
    import openmm
    from openmm.unit import angstrom, dalton, kilojoules_per_mole, nanometer

    raw = plan_path.read_bytes()
    record = json.loads(raw)
    verify(record)
    if openmm.__version__ != record['openmm_version']:
        raise ValueError('OpenMM version changed')
    case = record['cases'][index]
    source = Path(case['source'])
    rows = load_atoms(source/'h-internal-before.json.gz')
    points = json.loads(gzip.decompress(Path(case['path']).read_bytes()))['initial_nm']
    system = openmm.XmlSerializer.deserialize((source/'h-internal-system.xml').read_text())
    serialized = openmm.XmlSerializer.serialize(system)
    indices = methyl_indices(rows)
    if any(system.getParticleMass(i).value_in_unit(dalton) <= 0 for i in indices[2:]):
        raise ValueError('generated methyl hydrogens are not mobile')
    integrator = openmm.VerletIntegrator(0.)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('Reference'))
    recorder = reporter_type(openmm)(len(rows), 512)
    try:
        context.setPositions(points*nanometer)
        original_energy = context.getState(energy=True).getPotentialEnergy().value_in_unit(kilojoules_per_mole)
        candidates, energies = [], {}
        for angle in ANGLES:
            candidate = tetrahedral_methyl(rows, points, angle)
            context.setPositions(candidate*nanometer)
            energy = context.getState(energy=True).getPotentialEnergy().value_in_unit(kilojoules_per_mole)
            candidates.append(dict(angle_degrees=angle, energy_kj_mol=energy, H_points_nm=[candidate[i] for i in indices[2:]]))
            energies[angle] = energy
        selected = selected_angle(energies)
        applied = None if arm == 'unchanged' else 0 if arm == 'tetrahedral_zero' else selected
        initial = [list(p) for p in points] if applied is None else tetrahedral_methyl(rows, points, applied)
        context.setPositions(initial*nanometer)
        initial_rows = [dict(r, xyz=list(p)) for r, p in zip(rows, context.getState(positions=True).getPositions().value_in_unit(angstrom), strict=True)]
        openmm.LocalEnergyMinimizer.minimize(context, 1., 50, recorder)
        state = context.getState(positions=True, energy=True)
        final_points = [list(p) for p in state.getPositions().value_in_unit(nanometer)]
        final_rows = [dict(r, xyz=list(p)) for r, p in zip(rows, state.getPositions().value_in_unit(angstrom), strict=True)]
        if openmm.XmlSerializer.serialize(system) != serialized:
            raise ValueError('tetrahedral initialization changed original System')
        payload = json.dumps(dict(initial_nm=initial, final_nm=final_points, initial_atoms=initial_rows,
            final_atoms=final_rows, callbacks=recorder.records), sort_keys=True, allow_nan=False).encode()
        with (output/'path.json.gz').open('xb') as stream:
            stream.write(gzip.compress(payload, mtime=0))
        difference = coordinate_delta(load_atoms(source/'h-internal-after.json.gz'), final_rows)
        write(output/'trace.json', dict(case=case, arm=arm, plan_sha256=sha(raw), candidates=candidates,
            selected_degrees=selected, applied_degrees=applied, original_energy_kj_mol=original_energy,
            callbacks=len(recorder.records), baseline_endpoint_difference=difference,
            serialized_system_sha256=sha(serialized.encode()),
            final_energy_kj_mol=state.getPotentialEnergy().value_in_unit(kilojoules_per_mole),
            whole_system_minimization=False, production_changed=False, native_coordinates_read=False))
        if arm == 'unchanged' and any(not difference[k]['exact'] for k in ('heavy', 'hydrogen')):
            raise ValueError('tetrahedral unchanged control does not reproduce original endpoint')
    finally:
        del context
        del integrator


def run(plan_path, output):
    raw = plan_path.read_bytes()
    record = json.loads(raw)
    verify(record)
    output.mkdir(parents=True, exist_ok=False)
    (output/'executed-plan.json').write_bytes(raw)
    environment = dict(os.environ, PYTHONPATH=os.pathsep.join((str(Path('src').resolve()), str(Path('.').resolve()))))
    started, results = time.monotonic(), []
    for index, case in enumerate(record['cases']):
        for arm in ARMS:
            folder = (output/case['id']/f"seed-{case['seed']}"/case['source_arm']/arm).resolve()
            folder.mkdir(parents=True, exist_ok=False)
            row = dict(id=case['id'], seed=case['seed'], source_arm=case['source_arm'], arm=arm)
            remaining = record['policy']['total_seconds']-(time.monotonic()-started)
            if remaining <= 0:
                row.update(status='budget_not_executed', returncode=None)
            else:
                command = [record['python'], str(Path(__file__).resolve()), 'worker', '--plan', str(plan_path.resolve()),
                    '--index', str(index), '--arm', arm, '--output', str(folder)]
                with (folder/'run.log').open('xb') as log:
                    try:
                        child = subprocess.run(command, env=environment, stdout=log, stderr=subprocess.STDOUT,
                            timeout=min(remaining, record['policy']['per_cell_seconds']), check=False)
                        row.update(status='completed' if child.returncode == 0 else 'execution_blocker', returncode=child.returncode)
                    except subprocess.TimeoutExpired:
                        row.update(status='timeout', returncode=None)
            row['artifacts'] = {p.name: sha(p.read_bytes()) for p in sorted(folder.iterdir()) if p.is_file()}
            results.append(row)
            write(folder/'result.json', row)
            print(json.dumps({k: v for k, v in row.items() if k != 'artifacts'}), flush=True)
    verify(record)
    write(output/'evidence.json', dict(schema='dvbfixer.tetra-methyl-start-run.v1', cases=results,
        plan_sha256=sha(raw), production_changed=False, native_coordinates_read=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('plan')
    p.add_argument('--rigid-plan', type=Path, required=True)
    p.add_argument('--rigid-run', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    for name in ('run', 'worker'):
        p = sub.add_parser(name)
        p.add_argument('--plan', type=Path, required=True)
        p.add_argument('--output', type=Path, required=True)
        if name == 'worker':
            p.add_argument('--index', type=int, required=True)
            p.add_argument('--arm', choices=ARMS, required=True)
    args = parser.parse_args()
    if args.command == 'plan':
        plan(args.rigid_plan, args.rigid_run, args.output)
    elif args.command == 'run':
        run(args.plan, args.output)
    else:
        worker(args.plan, args.index, args.arm, args.output)
