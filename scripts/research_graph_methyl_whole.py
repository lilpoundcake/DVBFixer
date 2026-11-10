"""Private downstream native-phase control after the graph H study.

Independent saved-State native baselines precede ALL three-arm continuations.
No solvent, new seeds, force edits, safety-net changes or candidate PDB output.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import subprocess
import time
from pathlib import Path

from scripts.preflight_graph_methyl_transport import mapped_start
from scripts.research_beta_wall import restraint
from scripts.research_dunbrack_prior import sha, write
from scripts.research_graph_methyl_control import ARMS
from scripts.research_h_stage_transport import topology_from_rows
from scripts.research_ha_minimizer_path import verify


def plan(h_plan, h_root, h_audit, transport_preflight, output):
    raw = h_plan.read_bytes()
    parent = json.loads(raw)
    verify(parent)
    run_raw = (h_root/'evidence.json').read_bytes()
    run = json.loads(run_raw)
    audit_raw = h_audit.read_bytes()
    audit = json.loads(audit_raw)
    preflight_raw = transport_preflight.read_bytes()
    preflight = json.loads(preflight_raw)
    if audit['plan_sha256'] != sha(raw) or audit['run_sha256'] != sha(run_raw) or run['plan_sha256'] != sha(raw):
        raise ValueError('upstream H result/audit binding differs')
    if audit['summary']['statuses'] != {'completed': 30, 'refused_source': 18} or audit['summary']['arms']['unchanged']['exact_original_endpoints'] != 10:
        raise ValueError('complete exact upstream H study not accepted')
    if preflight['captured'] != 10 or preflight['exact_starts'] != 10:
        raise ValueError('original exact H/full transport not accepted')
    if preflight['script_sha256'] != sha(Path('scripts/preflight_graph_methyl_transport.py').read_bytes()):
        raise ValueError('transport implementation changed')
    inputs = {**parent['inputs'], str(h_plan): sha(raw), str(h_root/'evidence.json'): sha(run_raw),
        str(h_audit): sha(audit_raw), str(transport_preflight): sha(preflight_raw)}
    for row in run['cases']:
        folder = h_root/row['id']/f"seed-{row['seed']}"/row['arm']
        for name, digest in row['artifacts'].items():
            if sha((folder/name).read_bytes()) != digest:
                raise ValueError('H intervention artifacts changed')
            inputs[str(folder/name)] = digest
    sources = dict(parent['source_hashes'])
    for name in ('scripts/research_graph_methyl_whole.py', 'scripts/preflight_graph_methyl_transport.py',
                 'scripts/research_h_stage_transport.py', 'scripts/research_beta_wall.py'):
        sources[name] = sha(Path(name).read_bytes())
    write(output, dict(schema='dvbfixer.graph-methyl-whole-plan.v1', cases=parent['cases'], inputs=inputs,
        source_hashes=sources, h_root=str(h_root), python=parent['python'], openmm_version=parent['openmm_version'],
        policy=dict(arms=ARMS, per_cell_seconds=1800, total_seconds=14400, baseline_cells=10, continuation_cells=30,
            platform='Reference', no_solvent=True, tolerance_kj_mol_nm=10., maximum_iterations_per_phase=1000,
            baseline='original saved whole-State, native Simulation.minimizeEnergy twice; full untouched forces/anchors and original tenfold phase2 restraint-weight reduction',
            continuations='all30 H endpoints, original H-label production repair once, object-index map to full labels, same saved System/State/integrator; unchanged must equal all baseline nm states exactly',
            audits='raw and PDB-rounded global CA/THR-beta signed volumes/0.5A3 shape descriptors, H angles/bonds/close pairs, original heavy motion, all failures; no target-mask accuracy claim',
            limits='two native phases only BEFORE emergency safety net, no reflection, H regeneration, retuned force/phase/seed/budget, public candidate/PDB or scientific qualification'),
        production_changed=False))


def worker(plan_path, index, arm, baseline_root, output):
    import openmm
    from openmm.app import Simulation
    from openmm.unit import kilojoules_per_mole, nanometer

    raw = plan_path.read_bytes()
    record = json.loads(raw)
    verify(record)
    case = record['cases'][index]
    if case['status'] != 'captured' or arm not in ('native_baseline', *ARMS) or openmm.__version__ != record['openmm_version']:
        raise ValueError('undeclared native continuation case/arm/runtime')
    source = Path(case['source'])
    whole = json.loads((source/'whole-topology.json').read_bytes())
    system = openmm.XmlSerializer.deserialize((source/'whole-system.xml').read_text())
    initial_xml = openmm.XmlSerializer.serialize(system)
    integrator = openmm.XmlSerializer.deserialize((source/'whole-integrator.xml').read_text())
    original_state = openmm.XmlSerializer.deserialize((source/'whole-state.xml').read_text())
    original_nm = [list(p) for p in original_state.getPositions().value_in_unit(nanometer)]
    repairs = None
    initial = original_nm
    if arm != 'native_baseline':
        h_source = Path(record['h_root'])/case['id']/f"seed-{case['seed']}"/arm
        path = json.loads(gzip.decompress((h_source/'path.json.gz').read_bytes()))
        h_topology = json.loads((source/'H-topology.json').read_bytes())
        initial, repairs = mapped_start(h_topology['atoms'], h_topology['bonds'], path['final_nm'], whole['H_to_whole_indices'])
        if arm == 'unchanged' and initial != original_nm:
            raise ValueError('original unchanged whole nm start differs')
    topology = topology_from_rows(whole['atoms'], whole['bonds'])
    simulation = Simulation(topology, system, integrator, openmm.Platform.getPlatformByName('Reference'))
    frames = []
    try:
        simulation.context.setState(original_state)
        simulation.context.setPositions(initial*nanometer)
        positional = restraint(system)
        for stage in ('initial', 'phase1', 'phase2'):
            if stage == 'phase2':
                for i in range(positional.getNumParticles()):
                    particle, parameters = positional.getParticleParameters(i)
                    positional.setParticleParameters(i, particle, [parameters[0]/10, *parameters[1:]])
                positional.updateParametersInContext(simulation.context)
            if stage != 'initial':
                if arm == 'native_baseline':
                    simulation.minimizeEnergy(maxIterations=1000)
                else:
                    openmm.LocalEnergyMinimizer.minimize(simulation.context, 10*kilojoules_per_mole/nanometer, 1000)
            state = simulation.context.getState(positions=True, energy=True)
            points = [list(p) for p in state.getPositions().value_in_unit(nanometer)]
            xml = openmm.XmlSerializer.serialize(system)
            if stage == 'initial' and xml != initial_xml:
                raise ValueError('continuation changed initial System')
            payload = dict(positions_nm=points, energy_kj_mol=state.getPotentialEnergy().value_in_unit(kilojoules_per_mole))
            raw_frame = gzip.compress(json.dumps(payload, sort_keys=True, allow_nan=False).encode(), mtime=0)
            with (output/f'{stage}.json.gz').open('xb') as stream:
                stream.write(raw_frame)
            with (output/f'{stage}-system.xml').open('x') as stream:
                stream.write(xml)
            frame = dict(stage=stage, coordinate_sha256=sha(raw_frame), system_sha256=sha(xml.encode()))
            if arm != 'native_baseline':
                baseline = baseline_root/case['id']/f"seed-{case['seed']}"
                expected_system = (baseline/f'{stage}-system.xml').read_text()
                expected_points = json.loads(gzip.decompress((baseline/f'{stage}.json.gz').read_bytes()))['positions_nm']
                if xml != expected_system:
                    raise ValueError('continuation force/anchor/phase reduction differs from baseline')
                frame['exact_baseline_nm_state'] = points == expected_points
                if arm == 'unchanged' and points != expected_points:
                    raise ValueError('unchanged native baseline endpoint mismatch; retain partial evidence')
            frames.append(frame)
        write(output/'trace.json', dict(case=case, arm=arm, plan_sha256=sha(raw), frames=frames,
            post_H_repairs=repairs, emergency_reflection=False, observations_restored=False, production_changed=False))
    finally:
        del simulation
        del integrator


def run(plan_path, output, baseline_root=None):
    raw = plan_path.read_bytes()
    record = json.loads(raw)
    verify(record)
    baseline_mode = baseline_root is None
    baseline_binding = None
    if not baseline_mode:
        evidence = (baseline_root/'evidence.json').read_bytes()
        baseline = json.loads(evidence)
        if baseline['plan_sha256'] != sha(raw) or baseline['mode'] != 'native_baseline':
            raise ValueError('independent baseline plan/mode differs')
        if [(c['id'], c['seed']) for c in baseline['cases']] != [(c['id'], c['seed']) for c in record['cases']]:
            raise ValueError('independent baseline enrollment differs')
        for cell in baseline['cases']:
            if cell['status'] not in {'completed', 'refused_source'}:
                raise ValueError('independent baseline has blocker; do not start continuations')
            folder = baseline_root/cell['id']/f"seed-{cell['seed']}"
            for name, digest in cell['artifacts'].items():
                if sha((folder/name).read_bytes()) != digest:
                    raise ValueError('independent baseline artifact changed')
        baseline_binding = sha(evidence)
    output.mkdir(parents=True, exist_ok=False)
    (output/'executed-plan.json').write_bytes(raw)
    env = dict(os.environ, PYTHONPATH=os.pathsep.join((str(Path('src').resolve()), str(Path('.').resolve()))))
    results, started = [], time.monotonic()
    for index, case in enumerate(record['cases']):
        for arm in ('native_baseline',) if baseline_mode else ARMS:
            folder = output/case['id']/f"seed-{case['seed']}"
            if not baseline_mode:
                folder = folder/arm
            folder = folder.resolve()
            folder.mkdir(parents=True, exist_ok=False)
            row = dict(id=case['id'], seed=case['seed'], arm=arm, status=case['status'], returncode=None)
            if case['status'] == 'captured':
                remaining = record['policy']['total_seconds']-(time.monotonic()-started)
                if remaining <= 0:
                    row['status'] = 'budget_not_executed'
                else:
                    command = [record['python'], str(Path(__file__).resolve()), 'worker', '--plan', str(plan_path.resolve()),
                        '--index', str(index), '--arm', arm, '--output', str(folder)]
                    if not baseline_mode:
                        command += ['--baseline-root', str(baseline_root.resolve())]
                    with (folder/'run.log').open('xb') as log:
                        try:
                            process = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT,
                                timeout=min(remaining, record['policy']['per_cell_seconds']), check=False)
                            row.update(status='completed' if process.returncode == 0 else 'execution_blocker', returncode=process.returncode)
                        except subprocess.TimeoutExpired:
                            row['status'] = 'timeout'
            row['artifacts'] = {p.name: sha(p.read_bytes()) for p in sorted(folder.iterdir()) if p.is_file()}
            results.append(row)
            write(folder/'result.json', row)
            print(json.dumps({k: v for k, v in row.items() if k != 'artifacts'}), flush=True)
    verify(record)
    write(output/'evidence.json', dict(schema='dvbfixer.graph-methyl-whole-run.v1', cases=results,
        mode='native_baseline' if baseline_mode else 'continuations', plan_sha256=sha(raw),
        baseline_run_sha256=baseline_binding, production_changed=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('plan')
    p.add_argument('--h-plan', type=Path, required=True)
    p.add_argument('--h-root', type=Path, required=True)
    p.add_argument('--h-audit', type=Path, required=True)
    p.add_argument('--transport-preflight', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    for name in ('run', 'worker'):
        p = sub.add_parser(name)
        p.add_argument('--plan', type=Path, required=True)
        p.add_argument('--output', type=Path, required=True)
        p.add_argument('--baseline-root', type=Path)
        if name == 'worker':
            p.add_argument('--index', type=int, required=True)
            p.add_argument('--arm', choices=['native_baseline', *ARMS], required=True)
    args = parser.parse_args()
    if args.command == 'plan':
        plan(args.h_plan, args.h_root, args.h_audit, args.transport_preflight, args.output)
    elif args.command == 'run':
        run(args.plan, args.output, args.baseline_root)
    else:
        worker(args.plan, args.index, args.arm, args.baseline_root, args.output)
