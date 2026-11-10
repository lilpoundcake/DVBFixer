"""Private paired saved-System replay with one generated HA start intervention.

ALL24 starts, unchanged/outward_HA arms. No H regeneration, new heavy atoms,
parameter/seed search, emergency reflection or public candidate publication.
The serialized original --no-solvent physics and native two-phase budget stay
fixed. Raw failures/timeouts and observation movement remain auditable.
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

from scripts.audit_legacy_start_controls import load_atoms
from scripts.audit_stereo_force_projection import historical_verify, verify
from scripts.research_beta_wall import restraint
from scripts.research_dunbrack_prior import sha, write

ARMS = ('unchanged', 'outward_HA')


def outward_ha(rows, points_nm):
    """Place only GLN16 HA opposite the sum of its three unit heavy bonds.

    Preserve its original CA-H distance. No empirical table, native/reference
    coordinate, energy-directed orientation choice or heavy-coordinate rewrite.
    """
    mapped = {tuple(row['key']): i for i, row in enumerate(rows)}
    if len(mapped) != len(rows) or len(points_nm) != len(rows):
        raise ValueError('ambiguous/incomplete saved atom inventory')
    names = ('CA', 'N', 'C', 'CB', 'HA')
    keys = [('A', '16', '', 'GLN', name) for name in names]
    if any(key not in mapped for key in keys):
        raise ValueError('exact GLN16 context atoms missing')
    indices = [mapped[key] for key in keys]
    for i, name in zip(indices, names, strict=True):
        expected = {'H', 'D'} if name == 'HA' else {name[0]}
        if rows[i]['element'] not in expected:
            raise ValueError('GLN16 atom element differs')
    selected = [points_nm[i] for i in indices]
    if any(len(p) != 3 or not all(math.isfinite(x) for x in p) for p in selected):
        raise ValueError('nonfinite or malformed GLN16 coordinates')
    centre = selected[0]
    bonds = [tuple(x-y for x, y in zip(p, centre)) for p in selected[1:4]]
    lengths = [math.hypot(*v) for v in bonds]
    radius = math.dist(selected[4], centre)
    if min(*lengths, radius) <= 1e-12:
        raise ValueError('undefined heavy/HA bond vector')
    units = [tuple(x/r for x in v) for v, r in zip(bonds, lengths, strict=True)]
    direction = tuple(-sum(v[j] for v in units) for j in range(3))
    norm = math.hypot(*direction)
    if norm <= 1e-12:
        raise ValueError('undefined outward heavy-vector sum')
    direction = tuple(x/norm for x in direction)
    if any(sum(x*y for x, y in zip(direction, u)) >= 0 for u in units):
        raise ValueError('construction is not outward to every heavy neighbour')
    point = tuple(x+radius*d for x, d in zip(centre, direction))
    result = [list(p) for p in points_nm]
    result[indices[4]] = list(point)
    if any(result[i] != list(p) for i, p in enumerate(points_nm) if i != indices[4]):
        raise ValueError('HA construction moved another coordinate')
    if not math.isclose(math.dist(point, centre), radius, rel_tol=1e-12, abs_tol=1e-14):
        raise ValueError('HA bond length changed')
    return result, indices[4]


def initial_points(rows, state):
    from openmm.unit import angstrom, nanometer

    saved = state.getPositions().value_in_unit(angstrom)
    if len(saved) != len(rows) or any(list(point) != row['xyz'] for point, row in zip(saved, rows, strict=True)):
        raise ValueError('saved State differs from captured starting coordinates')
    return [[float(x) for x in p] for p in state.getPositions().value_in_unit(nanometer)]


def require_generated_ha(source):
    for line in source.read_text().splitlines():
        if line.startswith(('ATOM  ', 'HETATM')) and (line[21], line[22:26].strip(), line[26].strip(), line[12:16].strip()) == ('A', '16', '', 'HA'):
            raise ValueError('GLN16 HA is an original observation, not generated')


def plan(force_plan_path, angle_plan_path, output):
    import openmm

    force_raw = force_plan_path.read_bytes()
    force = json.loads(force_raw)
    verify(force)
    historical_verify(json.loads(Path(force['parent_plan']).read_bytes()), Path(force['historical_root']))
    angle_raw = angle_plan_path.read_bytes()
    angle = json.loads(angle_raw)
    verify(angle)
    if angle['parent_plan'] != str(force_plan_path):
        raise ValueError('angle/source force plans differ')
    cohort = json.loads(Path(force['parent_plan']).read_bytes())
    cases_by_id = {case['id']: case for case in cohort['cases']}
    cases = []
    inputs = {**force['inputs'], str(force_plan_path): sha(force_raw), str(angle_plan_path): sha(angle_raw)}
    for execution in force['executions']:
        if execution['status'] != 'completed':
            raise ValueError('saved cohort has inherited execution blockers')
        case = cases_by_id[execution['id']]
        folder = Path(force['run_root'])/execution['id']/f"seed-{execution['seed']}"/execution['arm']
        required = ['phase-1-before-system.xml', 'initial-integrator.xml', 'initial-state.xml',
            'phase-1-before.json.gz', 'phase-1-after.json.gz', 'phase-2-after.json.gz',
            'phase-2-before-system.xml', 'topology.json', 'trace.json']
        owned = {str(folder/name): sha((folder/name).read_bytes()) for name in required}
        source = Path(case['original_masked_source'])
        require_generated_ha(source)
        owned[str(source)] = sha(source.read_bytes())
        inputs.update(owned)
        trace = json.loads((folder/'trace.json').read_bytes())
        if len(trace['h_calls']) != 1 or trace['h_calls'][0]['platform'] != 'Reference':
            raise ValueError('generated-H source provenance differs')
        if [p['phase'] for p in trace['phases']] != [1, 2] or any(
                p[side]['platform']['name'] != 'Reference' for p in trace['phases'] for side in ('before', 'after')):
            raise ValueError('original two-phase Reference capture differs')
        rows = load_atoms(folder/'phase-1-before.json.gz')
        state = openmm.XmlSerializer.deserialize((folder/'initial-state.xml').read_text())
        outward_ha(rows, initial_points(rows, state))
        cases.append(dict(id=execution['id'], seed=execution['seed'], source_arm=execution['arm'],
            source=str(folder), target=case['target'], residue=case['residue'],
            original_masked_source=str(source), inputs=owned))
    sources = {str(p): sha(p.read_bytes()) for p in Path('src/dvbfixer').rglob('*.py')}
    for name in ('scripts/research_ha_start_control.py', 'scripts/research_beta_wall.py',
                 'scripts/audit_legacy_start_controls.py', 'scripts/research_dunbrack_prior.py'):
        sources[name] = sha(Path(name).read_bytes())
    write(output, dict(schema='dvbfixer.generated-ha-start-control-plan.v1', cases=cases,
        inputs=inputs, source_hashes=sources, python=sys.executable, openmm_version=openmm.__version__,
        source_git_head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        policy=dict(arms=ARMS, source_cases=24, maximum_replay_cells=48,
            phases=2, platform='Reference', tolerance_kj_mol_nm=10., max_iterations_per_phase=1000,
            timeout_seconds_per_cell=1800, total_wall_budget_seconds=7200,
            source='ALL24 frozen starts; observed heavy source immutable before replay, generated GLN16 HA only may change',
            intervention='ALL starts, not only inward/failing cases; deterministic negative sum of CA-to-N/C/CB unit heavy bonds, original CA-HA length',
            invariants='same serialized System/integrator/State, all non-HA starting coordinates exact; HA must be unrestrained; all parameters/anchors/vrefs unchanged',
            physics='original no-solvent nonperiodic two native phases; only original positional weights divided10 in phase2',
            retention='raw states and all failures; timeout/budget nonexecutions separate; no resume/adoption or additional retries',
            audit='unchanged replay vs bound source endpoints, complete target/global raw/rounded heavy gate and original-observation motion; HA angles reported separately',
            limits='retrospective H-start intervention, not original THR2 replay; no native reference/RMSD, parameter tuning, H regeneration, reflection or public repair'),
        production_changed=False, native_coordinates_read=False))


def worker(plan_path, case_index, arm, output):
    import openmm
    from openmm.unit import angstrom, kilojoules_per_mole, nanometer

    raw = plan_path.read_bytes()
    record = json.loads(raw)
    case = record['cases'][case_index]
    verify(dict(inputs=case['inputs'], source_hashes=record['source_hashes']))
    if openmm.__version__ != record['openmm_version']:
        raise ValueError('OpenMM version changed since enrollment')
    folder = Path(case['source'])
    rows = load_atoms(folder/'phase-1-before.json.gz')
    system = openmm.XmlSerializer.deserialize((folder/'phase-1-before-system.xml').read_text())
    integrator = openmm.XmlSerializer.deserialize((folder/'initial-integrator.xml').read_text())
    initial = openmm.XmlSerializer.deserialize((folder/'initial-state.xml').read_text())
    if system.usesPeriodicBoundaryConditions() or any(system.isVirtualSite(i) for i in range(system.getNumParticles())):
        raise ValueError('periodic/virtual-site system outside the declared no-solvent replay')
    if len(rows) != system.getNumParticles():
        raise ValueError('particle/source atom inventory differs')
    original_system = openmm.XmlSerializer.serialize(system)
    require_generated_ha(Path(case['original_masked_source']))
    original_points = initial_points(rows, initial)
    corrected, index = outward_ha(rows, original_points)
    positional = restraint(system)
    if any(i == index and p[0] != 0 for i, p in (positional.getParticleParameters(j) for j in range(positional.getNumParticles()))):
        raise ValueError('HA has an original positional restraint')
    points = original_points if arm == 'unchanged' else corrected
    allowed = [] if arm == 'unchanged' else [rows[index]['key']]
    changed = [rows[i]['key'] for i, (a, b) in enumerate(zip(original_points, points, strict=True)) if a != b]
    if any(key not in allowed for key in changed):
        raise ValueError('undeclared starting coordinate change')
    if openmm.XmlSerializer.serialize(system) != original_system:
        raise ValueError('H-start construction changed System parameters')
    write(output/'input-change.json', dict(arm=arm, allowed_atoms=allowed, changed_atoms=changed,
        HA_index=index, HA_movement_angstrom=10*math.dist(original_points[index], points[index]),
        non_HA_coordinates_exact=True, original_system_exact=True,
        original_state_sha256=sha((folder/'initial-state.xml').read_bytes())))
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('Reference'))
    try:
        context.setState(initial)
        if arm == 'outward_HA':
            context.setPositions(points*nanometer)
        frames = []
        for stage in ('initial', 'phase1', 'phase2'):
            if stage == 'phase2':
                for i in range(positional.getNumParticles()):
                    particle, parameters = positional.getParticleParameters(i)
                    positional.setParticleParameters(i, particle, [parameters[0]/10, *parameters[1:]])
                positional.updateParametersInContext(context)
                expected = openmm.XmlSerializer.deserialize((folder/'phase-2-before-system.xml').read_text())
                if openmm.XmlSerializer.serialize(system) != openmm.XmlSerializer.serialize(expected):
                    raise ValueError('phase2 parameters differ from bound original System')
            if stage != 'initial':
                openmm.LocalEnergyMinimizer.minimize(context, 10*kilojoules_per_mole/nanometer, 1000)
            system_path = output/f'{stage}-system.xml'
            with system_path.open('x') as stream:
                stream.write(openmm.XmlSerializer.serialize(system))
            state = context.getState(positions=True, energy=True)
            xyz = state.getPositions(asNumpy=True).value_in_unit(angstrom).tolist()
            material = [dict(row, xyz=point) for row, point in zip(rows, xyz, strict=True)]
            energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
            payload = json.dumps(dict(atoms=material, energy_kj_mol=energy), sort_keys=True, allow_nan=False).encode()
            path = output/f'{stage}.json.gz'
            with path.open('xb') as stream:
                stream.write(gzip.compress(payload, mtime=0))
            frames.append(dict(stage=stage, sha256=sha(path.read_bytes()),
                system_sha256=sha(system_path.read_bytes()), energy_kj_mol=energy))
        write(output/'trace.json', dict(case_index=case_index, arm=arm, frames=frames,
            platform=context.getPlatform().getName(), plan_sha256=sha(raw),
            integrator_unchanged=True, no_integration=True, no_H_regeneration=True,
            emergency_reflection=False, production_changed=False, native_coordinates_read=False))
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
    started = time.monotonic()
    results = []
    for i, case in enumerate(record['cases']):
        for arm in record['policy']['arms']:
            folder = output/case['id']/f"source-seed-{case['seed']}"/case['source_arm']/arm
            folder.mkdir(parents=True, exist_ok=False)
            row = dict(id=case['id'], seed=case['seed'], source_arm=case['source_arm'], arm=arm)
            remaining = record['policy']['total_wall_budget_seconds']-(time.monotonic()-started)
            if remaining <= 0:
                row.update(status='budget_not_executed', returncode=None)
            else:
                argv = [record['python'], str(Path(__file__).resolve()), 'worker', '--plan', str(plan_path.resolve()),
                    '--case-index', str(i), '--arm', arm, '--output', str(folder.resolve())]
                with (folder/'run.log').open('xb') as log:
                    try:
                        process = subprocess.run(argv, env=environment, stdout=log, stderr=subprocess.STDOUT,
                            timeout=min(remaining, record['policy']['timeout_seconds_per_cell']), check=False)
                        row.update(status='completed' if process.returncode == 0 else 'execution_blocker', returncode=process.returncode)
                    except subprocess.TimeoutExpired:
                        row.update(status='timeout', returncode=None)
            row['artifacts'] = {p.name: sha(p.read_bytes()) for p in sorted(folder.iterdir()) if p.is_file()}
            write(folder/'result.json', row)
            results.append(row)
            print(json.dumps({k: row[k] for k in ('id', 'seed', 'source_arm', 'arm', 'status')}), flush=True)
    verify(record)
    write(output/'evidence.json', dict(schema='dvbfixer.generated-ha-start-control-run.v1', cases=results,
        plan_sha256=sha(raw), production_changed=False, native_coordinates_read=False,
        source_observations_restored=False, scientific_qualification=False))
    print(json.dumps(dict(enrolled=len(results), completed=sum(r['status']=='completed' for r in results)), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('plan')
    p.add_argument('--force-plan', type=Path, required=True)
    p.add_argument('--angle-plan', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    for name in ('run', 'worker'):
        p = sub.add_parser(name)
        p.add_argument('--plan', type=Path, required=True)
        p.add_argument('--output', type=Path, required=True)
        if name == 'worker':
            p.add_argument('--case-index', type=int, required=True)
            p.add_argument('--arm', choices=ARMS, required=True)
    args = parser.parse_args()
    if args.command == 'plan':
        plan(args.force_plan, args.angle_plan, args.output)
    elif args.command == 'worker':
        worker(args.plan, args.case_index, args.arm, args.output)
    else:
        run(args.plan, args.output)
