"""Private whole-system continuation of all three frozen tetra-methyl H arms.

Existing post-H helper, original saved Systems/State/integrator and two native
Reference/1000-iteration phases. No solvent, new H, reflection or public output.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import subprocess
import sys
import time
from pathlib import Path

from scripts.audit_legacy_start_controls import coordinate_delta, load_atoms
from scripts.audit_stereo_force_projection import historical_verify
from scripts.research_beta_wall import restraint
from scripts.research_dunbrack_prior import sha, write
from scripts.research_h_stage_transport import post_h_repair
from scripts.research_ha_minimizer_path import verify


def plan(h_plan_path, h_run, force_path, preflight_path, output):
    raw = h_plan_path.read_bytes()
    parent = json.loads(raw)
    verify(parent)
    force = json.loads(force_path.read_bytes())
    verify(force)
    historical_verify(json.loads(Path(force['parent_plan']).read_bytes()), Path(force['historical_root']))
    checked = json.loads(preflight_path.read_bytes())
    if checked['plan_sha256'] != sha(raw) or checked['exact_original_starts'] != 24 or checked['enrolled'] != 24 or any(not r['nm_start_exact'] for r in checked['cases']):
        raise ValueError('exact post-H nm transport preflight not complete')
    if checked['script_sha256'] != sha(Path('scripts/research_h_stage_transport.py').read_bytes()) or checked['production_geometry_sha256'] != sha(Path('src/dvbfixer/ffutils/geometry.py').read_bytes()):
        raise ValueError('post-H transport implementation changed after preflight')
    run_raw = (h_run/'evidence.json').read_bytes()
    run = json.loads(run_raw)
    expected = [(c['id'], c['seed'], c['source_arm'], a) for c in parent['cases'] for a in parent['policy']['arms']]
    if run['plan_sha256'] != sha(raw) or [(r['id'], r['seed'], r['source_arm'], r['arm']) for r in run['cases']] != expected or len(expected) != 72 or any(r['status'] != 'completed' for r in run['cases']):
        raise ValueError('complete H-arm cohort differs')
    inputs = {**parent['inputs'], **force['inputs'], str(h_plan_path): sha(raw), str(h_run/'evidence.json'): sha(run_raw),
        str(force_path): sha(force_path.read_bytes()), str(preflight_path): sha(preflight_path.read_bytes())}
    sources = {**parent['source_hashes'], **force['source_hashes']}
    for p in Path('src/dvbfixer').rglob('*.py'):
        sources[str(p)] = sha(p.read_bytes())
    for name in ('scripts/research_methyl_whole_control.py', 'scripts/research_h_stage_transport.py',
                 'scripts/research_beta_wall.py', 'scripts/audit_legacy_start_controls.py', 'scripts/research_dunbrack_prior.py'):
        sources[str(Path(name).resolve())] = sha(Path(name).read_bytes())
    science = {c['id']: c for c in force['cases']}
    original = {(c['id'], c['seed'], c['arm']): c for c in force['executions']}
    cases = []
    for c in parent['cases']:
        base = original[(c['id'], c['seed'], c['source_arm'])]
        if base['status'] != 'completed':
            raise ValueError('original whole-system case has blocker')
        folder = Path(force['run_root'])/c['id']/f"seed-{c['seed']}"/c['source_arm']
        for name, digest in base['artifacts'].items():
            p = folder/name
            if sha(p.read_bytes()) != digest:
                raise ValueError('original whole-system artifact changed')
            inputs[str(p)] = digest
        source_pdb = Path(science[c['id']]['original_masked_source'])
        inputs[str(source_pdb)] = sha(source_pdb.read_bytes())
        cases.append(dict(id=c['id'], seed=c['seed'], source_arm=c['source_arm'], source=c['source'],
            original_whole=str(folder), target=science[c['id']]['target'], residue=science[c['id']]['residue'],
            observed_source=str(source_pdb)))
    for row in run['cases']:
        folder = h_run/row['id']/f"seed-{row['seed']}"/row['source_arm']/row['arm']
        for name, digest in row['artifacts'].items():
            p = folder/name
            if sha(p.read_bytes()) != digest:
                raise ValueError('frozen H-arm artifact changed')
            inputs[str(p)] = digest
    write(output, dict(schema='dvbfixer.methyl-whole-control-plan.v1', cases=cases, inputs=inputs,
        source_hashes=sources, h_run=str(h_run), arms=parent['policy']['arms'], python=sys.executable,
        openmm_version=parent['openmm_version'],
        policy=dict(cells=72, phases=2, platform='Reference', tolerance_kj_mol_nm=10., max_iterations_per_phase=1000,
            per_cell_seconds=1800, total_seconds=7200,
            enrollment='ALL24 source contexts x ALL3 predeclared H arms; no new seed/angle/iteration search',
            bridge='existing post-H repair exactly once on reconstructed bound topology; exact original nm starts preflighted',
            physics='original no-solvent saved System/State/integrator; only original positional k divided10 in phase2, original anchors/vrefs unchanged',
            invariants='heavy coordinates remain exact through H transport, original whole-system parameters unchanged; unchanged arm must reproduce both saved endpoints exactly',
            audit='actual-frame raw/rounded target and whole-context heavy geometry, all H geometry, original-observation motion, baseline endpoint differences and every failure',
            limits='no emergency reflection, coordinate restoration, native accuracy, public candidate/PDB publication or production qualification'),
        production_changed=False, native_coordinates_read=False))


def worker(plan_path, index, arm, output):
    import openmm
    from openmm.unit import angstrom, kilojoules_per_mole, nanometer

    raw = plan_path.read_bytes()
    record = json.loads(raw)
    verify(record)
    if openmm.__version__ != record['openmm_version']:
        raise ValueError('whole-system OpenMM version differs')
    case = record['cases'][index]
    source = Path(case['source'])
    original = Path(case['original_whole'])
    h_folder = Path(record['h_run'])/case['id']/f"seed-{case['seed']}"/case['source_arm']/arm
    h_data = json.loads(gzip.decompress((h_folder/'path.json.gz').read_bytes()))
    rows = load_atoms(source/'h-internal-before.json.gz')
    bonds = json.loads((source/'topology.json').read_bytes())['bonds']
    points, repaired = post_h_repair(rows, bonds, h_data['final_nm'])
    system = openmm.XmlSerializer.deserialize((original/'phase-1-before-system.xml').read_text())
    integrator = openmm.XmlSerializer.deserialize((original/'initial-integrator.xml').read_text())
    initial_state = openmm.XmlSerializer.deserialize((original/'initial-state.xml').read_text())
    if system.usesPeriodicBoundaryConditions() or system.getNumParticles() != len(rows):
        raise ValueError('original no-solvent particle contract differs')
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('Reference'))
    frames = []
    differences = []
    try:
        context.setState(initial_state)
        context.setPositions(points*nanometer)
        positional = restraint(system)
        for stage in ('initial', 'phase1', 'phase2'):
            if stage == 'phase2':
                for i in range(positional.getNumParticles()):
                    particle, parameters = positional.getParticleParameters(i)
                    positional.setParticleParameters(i, particle, [parameters[0]/10, *parameters[1:]])
                positional.updateParametersInContext(context)
            expected_path = original/('phase-2-before-system.xml' if stage == 'phase2' else 'phase-1-before-system.xml')
            canonical_expected = openmm.XmlSerializer.serialize(openmm.XmlSerializer.deserialize(expected_path.read_text()))
            actual_system = openmm.XmlSerializer.serialize(system)
            if actual_system != canonical_expected:
                raise ValueError('whole-system parameters/anchors/reference changed')
            if stage != 'initial':
                openmm.LocalEnergyMinimizer.minimize(context, 10*kilojoules_per_mole/nanometer, 1000)
            state = context.getState(positions=True, energy=True)
            material = [dict(r, xyz=list(p)) for r, p in zip(rows, state.getPositions().value_in_unit(angstrom), strict=True)]
            energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
            payload = json.dumps(dict(atoms=material, energy_kj_mol=energy), sort_keys=True, allow_nan=False).encode()
            path = output/f'{stage}.json.gz'
            with path.open('xb') as stream:
                stream.write(gzip.compress(payload, mtime=0))
            with (output/f'{stage}-system.xml').open('x') as stream:
                stream.write(actual_system)
            reference = load_atoms(original/('phase-1-before.json.gz' if stage == 'initial' else f"phase-{1 if stage == 'phase1' else 2}-after.json.gz"))
            delta = coordinate_delta(reference, material)
            differences.append(delta)
            frames.append(dict(stage=stage, energy_kj_mol=energy, atom_sha256=sha(path.read_bytes()),
                system_sha256=sha(actual_system.encode()), baseline_difference=delta))
        write(output/'trace.json', dict(case=case, arm=arm, plan_sha256=sha(raw), frames=frames,
            post_H_repairs=repaired, original_parameters_exact=True, emergency_reflection=False,
            restored_observations=False, production_changed=False, native_coordinates_read=False))
        if arm == 'unchanged' and any(not d[k]['exact'] for d in differences for k in ('heavy', 'hydrogen')):
            raise ValueError('whole-system unchanged control not exact; retain artifacts and blocker')
    finally:
        del context
        del integrator


def run(plan_path, output):
    raw = plan_path.read_bytes()
    record = json.loads(raw)
    verify(record)
    output.mkdir(parents=True, exist_ok=False)
    (output/'executed-plan.json').write_bytes(raw)
    env = dict(os.environ, PYTHONPATH=os.pathsep.join((str(Path('src').resolve()), str(Path('.').resolve()))))
    started, results = time.monotonic(), []
    for index, case in enumerate(record['cases']):
        for arm in record['arms']:
            folder = (output/case['id']/f"seed-{case['seed']}"/case['source_arm']/arm).resolve()
            folder.mkdir(parents=True, exist_ok=False)
            row = dict(id=case['id'], seed=case['seed'], source_arm=case['source_arm'], arm=arm)
            remaining = record['policy']['total_seconds']-(time.monotonic()-started)
            if remaining <= 0:
                row.update(status='budget_not_executed', returncode=None)
            else:
                argv = [record['python'], str(Path(__file__).resolve()), 'worker', '--plan', str(plan_path.resolve()),
                    '--index', str(index), '--arm', arm, '--output', str(folder)]
                with (folder/'run.log').open('xb') as log:
                    try:
                        process = subprocess.run(argv, env=env, stdout=log, stderr=subprocess.STDOUT,
                            timeout=min(remaining, record['policy']['per_cell_seconds']), check=False)
                        row.update(status='completed' if process.returncode == 0 else 'execution_blocker', returncode=process.returncode)
                    except subprocess.TimeoutExpired:
                        row.update(status='timeout', returncode=None)
            row['artifacts'] = {p.name: sha(p.read_bytes()) for p in sorted(folder.iterdir()) if p.is_file()}
            results.append(row)
            write(folder/'result.json', row)
            print(json.dumps({k: v for k, v in row.items() if k != 'artifacts'}), flush=True)
    verify(record)
    write(output/'evidence.json', dict(schema='dvbfixer.methyl-whole-control-run.v1', cases=results,
        plan_sha256=sha(raw), production_changed=False, native_coordinates_read=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('plan')
    p.add_argument('--h-plan', type=Path, required=True)
    p.add_argument('--h-run', type=Path, required=True)
    p.add_argument('--force-plan', type=Path, required=True)
    p.add_argument('--preflight', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    for name in ('run', 'worker'):
        p = sub.add_parser(name)
        p.add_argument('--plan', type=Path, required=True)
        p.add_argument('--output', type=Path, required=True)
        if name == 'worker':
            p.add_argument('--index', type=int, required=True)
            p.add_argument('--arm', choices=('unchanged', 'tetrahedral_zero', 'tetrahedral_energy_selected'), required=True)
    args = parser.parse_args()
    if args.command == 'plan':
        plan(args.h_plan, args.h_run, args.force_plan, args.preflight, args.output)
    elif args.command == 'run':
        run(args.plan, args.output)
    else:
        worker(args.plan, args.index, args.arm, args.output)
