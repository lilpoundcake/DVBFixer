"""Private paired initial methyl-H rotation control on ALL24 native H starts.

Fixed 24-angle grid, initial native energy only. Same forces/Reference/50 updates.
No observed coordinate rewrite, H regeneration, force ablation or public repair.
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

ANGLES = tuple(range(0, 360, 15))
ARMS = ('unchanged', 'energy_selected_methyl')
TIE_KJ_MOL = 1e-8


def methyl_indices(rows):
    keys = {tuple(r['key']): i for i, r in enumerate(rows)}
    if len(keys) != len(rows):
        raise ValueError('duplicate atom identity')
    names = ('CG2', 'CB', 'HG21', 'HG22', 'HG23')
    expected = [('A', '24', '', 'THR', n) for n in names]
    if any(k not in keys for k in expected):
        raise ValueError('exact THR24 methyl identity missing')
    indices = [keys[k] for k in expected]
    if any(rows[i]['element'] != ('C' if n in ('CG2', 'CB') else 'H') for n, i in zip(names, indices, strict=True)):
        raise ValueError('unexpected methyl element')
    return indices


def rotate_methyl(rows, points, degrees):
    if degrees not in ANGLES or len(points) != len(rows):
        raise ValueError('undeclared angle or particle inventory')
    indices = methyl_indices(rows)
    if any(len(p) != 3 or not all(math.isfinite(x) for x in p) for p in points):
        raise ValueError('invalid methyl-control coordinates')
    result = [list(p) for p in points]
    if degrees == 0:
        return result
    origin, end = (points[i] for i in indices[:2])
    axis = [x-y for x, y in zip(end, origin)]
    norm = math.hypot(*axis)
    if norm <= 1e-12:
        raise ValueError('undefined methyl rotation axis')
    axis = [x/norm for x in axis]
    theta = math.radians(degrees)
    cosine, sine = math.cos(theta), math.sin(theta)
    for i in indices[2:]:
        vector = [x-y for x, y in zip(points[i], origin)]
        dot = sum(x*y for x, y in zip(axis, vector))
        cross = [axis[1]*vector[2]-axis[2]*vector[1], axis[2]*vector[0]-axis[0]*vector[2], axis[0]*vector[1]-axis[1]*vector[0]]
        result[i] = [o+v*cosine+c*sine+a*dot*(1-cosine) for o, v, c, a in zip(origin, vector, cross, axis)]
    return result


def selected_angle(energies):
    if set(energies) != set(ANGLES) or any(not math.isfinite(e) for e in energies.values()):
        raise ValueError('complete finite declared energy grid required')
    minimum = min(energies.values())
    eligible = [a for a in ANGLES if energies[a] <= minimum+TIE_KJ_MOL]
    return min(eligible, key=lambda a: (min(a, 360-a), a if a <= 180 else a-360))


def plan(path_plan_path, path_root, output):
    import openmm

    parent_raw = path_plan_path.read_bytes()
    parent = json.loads(parent_raw)
    verify(parent)
    run_raw = (path_root/'evidence.json').read_bytes()
    run = json.loads(run_raw)
    if run['plan_sha256'] != sha(parent_raw) or len(run['cases']) != 24:
        raise ValueError('complete native path cohort not bound')
    if [(c['id'], c['seed'], c['source_arm']) for c in parent['cases']] != [(c['id'], c['seed'], c['source_arm']) for c in run['cases']]:
        raise ValueError('native path cohort identity differs')
    inputs = {**parent['inputs'], str(path_plan_path): sha(parent_raw), str(path_root/'evidence.json'): sha(run_raw)}
    cases = []
    for case, result in zip(parent['cases'], run['cases'], strict=True):
        if result['status'] != 'completed':
            raise ValueError('native path source has execution blockers')
        folder = path_root/case['id']/f"seed-{case['seed']}"/case['source_arm']
        for name, digest in result['artifacts'].items():
            p = folder/name
            if sha(p.read_bytes()) != digest:
                raise ValueError('native path source artifact changed')
            inputs[str(p)] = digest
        source = Path(case['source'])
        rows = load_atoms(source/'h-internal-before.json.gz')
        indices = methyl_indices(rows)
        observed = {tuple(r['key']) for r in load_atoms(source/'h-input.json.gz')}
        if any(tuple(rows[i]['key']) in observed for i in indices[2:]):
            raise ValueError('methyl H is an original observation')
        system = openmm.XmlSerializer.deserialize((source/'h-internal-system.xml').read_text())
        if system.getNumConstraints() or system.usesPeriodicBoundaryConditions() or any(system.isVirtualSite(i) for i in range(system.getNumParticles())):
            raise ValueError('constraints/PBC/virtual sites outside this control')
        cases.append(dict(id=case['id'], seed=case['seed'], source_arm=case['source_arm'], source=str(source), path=str(folder/'path.json.gz')))
    sources = {**parent['source_hashes']}
    for name in ('scripts/research_methyl_start_control.py', 'scripts/research_ha_minimizer_path.py',
                 'scripts/audit_legacy_start_controls.py', 'scripts/research_dunbrack_prior.py'):
        sources[str(Path(name).resolve())] = sha(Path(name).read_bytes())
    write(output, dict(schema='dvbfixer.methyl-start-control-plan.v1', cases=cases, inputs=inputs,
        source_hashes=sources, python=sys.executable, openmm_version=openmm.__version__,
        policy=dict(arms=ARMS, angles_degrees=ANGLES, energy_tie_kj_mol=TIE_KJ_MOL,
            tolerance_kj_mol_nm=1., max_iterations=50, platform='Reference',
            per_cell_seconds=900, total_seconds=3600, maximum_callbacks=512,
            enrollment='ALL24 exact native H Systems/starts, two matched arms, no new H seeds',
            intervention='only generated THR24 HG21/HG22/HG23 rigidly rotate together around fixed CG2-to-CB axis; preserve all other initial coordinates, C-H lengths and internal geometry',
            selection='enumerate fixed 24-angle grid, choose minimum original whole-System initial potential energy; near-ties choose shortest rotation then signed angle; never score final/native accuracy',
            invariants='all original forces/parameters/masses unchanged; no H regeneration, force ablation, integration or whole-system minimization; unchanged control must reproduce exact endpoint',
            audits='all grid energies/selection/geometry invariants, pair/HA initial forces and native callback path, all final H geometry and exact heavy immobility, all failures retained',
            limits='post-attribution mechanism control on one protein, not production selection, general rotamer search or repair qualification'),
        production_changed=False, native_coordinates_read=False))


def worker(plan_path, index, arm, output):
    import openmm
    from openmm.unit import angstrom, dalton, kilojoules_per_mole, nanometer

    raw = plan_path.read_bytes()
    record = json.loads(raw)
    verify(record)
    if openmm.__version__ != record['openmm_version']:
        raise ValueError('OpenMM version differs')
    case = record['cases'][index]
    source = Path(case['source'])
    rows = load_atoms(source/'h-internal-before.json.gz')
    path = json.loads(gzip.decompress(Path(case['path']).read_bytes()))
    points = path['initial_nm']
    indices = methyl_indices(rows)
    system = openmm.XmlSerializer.deserialize((source/'h-internal-system.xml').read_text())
    serialized = openmm.XmlSerializer.serialize(system)
    if any(system.getParticleMass(i).value_in_unit(dalton) <= 0 for i in indices[2:]):
        raise ValueError('generated methyl H is not mobile')
    integrator = openmm.VerletIntegrator(0.)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('Reference'))
    recorder = reporter_type(openmm)(len(rows), record['policy']['maximum_callbacks'])
    try:
        energies = {}
        candidates = []
        for angle in ANGLES:
            candidate = rotate_methyl(rows, points, angle)
            context.setPositions(candidate*nanometer)
            energy = context.getState(energy=True).getPotentialEnergy().value_in_unit(kilojoules_per_mole)
            energies[angle] = energy
            candidates.append(dict(angle_degrees=angle, energy_kj_mol=energy, H_points_nm=[candidate[i] for i in indices[2:]]))
        choice = selected_angle(energies)
        applied = 0 if arm == 'unchanged' else choice
        initial_points = rotate_methyl(rows, points, applied)
        context.setPositions(initial_points*nanometer)
        before = context.getState(positions=True).getPositions().value_in_unit(angstrom)
        before_rows = [dict(r, xyz=list(p)) for r, p in zip(rows, before, strict=True)]
        openmm.LocalEnergyMinimizer.minimize(context, 1., 50, recorder)
        state = context.getState(positions=True, energy=True)
        final_points = [list(p) for p in state.getPositions().value_in_unit(nanometer)]
        final_rows = [dict(r, xyz=list(p)) for r, p in zip(rows, state.getPositions().value_in_unit(angstrom), strict=True)]
        if openmm.XmlSerializer.serialize(system) != serialized:
            raise ValueError('initial methyl rotation changed System')
        payload = json.dumps(dict(initial_nm=initial_points, final_nm=final_points, initial_atoms=before_rows,
            final_atoms=final_rows, callbacks=recorder.records), sort_keys=True, allow_nan=False).encode()
        with (output/'path.json.gz').open('xb') as stream:
            stream.write(gzip.compress(payload, mtime=0))
        baseline = coordinate_delta(load_atoms(source/'h-internal-after.json.gz'), final_rows)
        write(output/'trace.json', dict(case=case, arm=arm, plan_sha256=sha(raw),
            candidates=candidates, selected_degrees=choice, applied_degrees=applied,
            callbacks=len(recorder.records), baseline_endpoint_difference=baseline,
            serialized_system_sha256=sha(serialized.encode()), final_energy_kj_mol=state.getPotentialEnergy().value_in_unit(kilojoules_per_mole),
            whole_system_minimization=False, production_changed=False, native_coordinates_read=False))
        if arm == 'unchanged' and any(not baseline[k]['exact'] for k in ('heavy', 'hydrogen')):
            raise ValueError('unchanged H replay differs from original endpoint')
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
    write(output/'evidence.json', dict(schema='dvbfixer.methyl-start-control-run.v1', cases=results,
        plan_sha256=sha(raw), production_changed=False, native_coordinates_read=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('plan')
    p.add_argument('--path-plan', type=Path, required=True)
    p.add_argument('--path-root', type=Path, required=True)
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
        plan(args.path_plan, args.path_root, args.output)
    elif args.command == 'run':
        run(args.plan, args.output)
    else:
        worker(args.plan, args.index, args.arm, args.output)
