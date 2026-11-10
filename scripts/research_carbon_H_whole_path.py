"""Observational replay of ALL24 THR24 native whole paths, no new science arms.

Same saved initial nm coordinates, original Systems/State/integrator, two native
10/1000 phases. Reporter never changes arrays, requests exit or retunes physics.
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

from scripts.carbon_H_capture_adapter import read_json
from scripts.research_beta_wall import restraint
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import reporter_type, verify


def plan(attribution_plan, attribution_audit, whole_plan, output):
    import openmm
    from openmm.unit import dalton

    raw = attribution_plan.read_bytes()
    parent = json.loads(raw)
    verify(parent)
    audit_raw = attribution_audit.read_bytes()
    audit = json.loads(audit_raw)
    whole_raw = whole_plan.read_bytes()
    whole = json.loads(whole_raw)
    verify(whole)
    if (audit['plan_sha256'] != sha(raw) or audit['summary']['paths'] != 24
            or audit['summary']['states'] != 72 or len(parent['cases']) != 24
            or whole['openmm_version'] != openmm.__version__):
        raise ValueError('complete original matched attribution/whole source required')
    inputs = {**parent['inputs'], **whole['inputs'], str(attribution_plan):sha(raw),
        str(attribution_audit):sha(audit_raw), str(whole_plan):sha(whole_raw)}
    cases = []
    for case in parent['cases']:
        original = next(c for c in whole['cases'] if (c['id'], c['seed'], c['conditioning'])
            == (case['id'], case['seed'], case['conditioning']))
        system = openmm.XmlSerializer.deserialize((Path(case['source'])/'initial-system.xml').read_text())
        if system.getNumConstraints() or any(system.getParticleMass(i).value_in_unit(dalton) <= 0 for i in range(system.getNumParticles())):
            raise ValueError('unconstrained all-mobile native whole path required; no subset exclusion')
        cases.append(dict(case, integrator=original['whole']['integrator']))
    sources = {name:sha(Path(name).read_bytes()) for name in (
        'scripts/research_carbon_H_whole_path.py', 'scripts/research_ha_minimizer_path.py',
        'scripts/research_beta_wall.py', 'scripts/carbon_H_capture_adapter.py')}
    write(output, dict(schema='dvbfixer.carbon-H-whole-path-plan.v1', cases=cases, inputs=inputs,
        source_hashes=sources, python=sys.executable, openmm_version=openmm.__version__,
        policy=dict(platform='Reference', no_solvent=True, tolerance_kj_mol_nm=10., maximum_iterations_per_phase=1000,
            maximum_callbacks_per_phase=2048, per_cell_seconds=1800, total_seconds=14400,
            enrollment='ALL24 original THR24 factorial paths: two old seeds/three conditioning contexts/four original arms, no new independent samples',
            start='original whole snapshot nm start already contains original H-label repair; never run repair again or regenerate H',
            physics='same native System/State/integrator, original anchors, original phase2 positional weights divided10; no force/phase/strength changes',
            observation='original native minimizer reporter, no array mutation/additional objective or early-exit request; keep coordinates/gradients/objective at every callback',
            controls='every initial and BOTH nm endpoints must equal the frozen original whole path exactly; same per-phase serialized Systems',
            limits='observational computational replays, no retry/iteration/seed/force/phase search, reflection, candidate/PDB or production qualification'),
        production_changed=False))


def worker(plan_path, index, output):
    import openmm
    from openmm.unit import kilojoules_per_mole, nanometer

    raw = plan_path.read_bytes()
    record = json.loads(raw)
    verify(record)
    if openmm.__version__ != record['openmm_version']:
        raise ValueError('native path runtime differs')
    case = record['cases'][index]
    source = Path(case['source'])
    system = openmm.XmlSerializer.deserialize((source/'initial-system.xml').read_text())
    state = openmm.XmlSerializer.deserialize(Path(case['state']).read_text())
    integrator = openmm.XmlSerializer.deserialize(Path(case['integrator']).read_text())
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('Reference'))
    phases = []
    try:
        context.setState(state)
        initial = read_json(source/'initial.json.gz')['positions_nm']
        context.setPositions(initial*nanometer)
        positional = restraint(system)
        for phase in ('phase1', 'phase2'):
            if phase == 'phase2':
                for i in range(positional.getNumParticles()):
                    particle, parameters = positional.getParticleParameters(i)
                    positional.setParticleParameters(i, particle, [parameters[0]/10, *parameters[1:]])
                positional.updateParametersInContext(context)
            xml = openmm.XmlSerializer.serialize(system)
            if xml != (source/f'{phase}-system.xml').read_text():
                raise ValueError('native path force/anchor parameters differ')
            before = [list(p) for p in context.getState(positions=True).getPositions().value_in_unit(nanometer)]
            reporter = reporter_type(openmm)(len(before), record['policy']['maximum_callbacks_per_phase'])
            path = output/f'{phase}-path.json.gz'
            returned = False
            try:
                openmm.LocalEnergyMinimizer.minimize(context, 10*kilojoules_per_mole/nanometer, 1000, reporter)
                returned = True
            finally:
                final = [list(p) for p in context.getState(positions=True).getPositions().value_in_unit(nanometer)]
                payload = dict(initial_nm=before, final_nm=final, callbacks=reporter.records, native_call_returned=returned)
                with path.open('xb') as stream:
                    stream.write(gzip.compress(json.dumps(payload, sort_keys=True, allow_nan=False).encode(), mtime=0))
            exact = final == read_json(source/f'{phase}.json.gz')['positions_nm']
            phases.append(dict(phase=phase, path_sha256=sha(path.read_bytes()), callbacks=len(reporter.records),
                exact_original_nm_endpoint=exact, system_sha256=sha(xml.encode())))
            if openmm.XmlSerializer.serialize(system) != xml or not exact:
                raise ValueError('native reporter replay changed System/endpoint; retain blocker, no retry')
        write(output/'trace.json', dict(case_id=[case['id'], case['seed'], case['conditioning'], case['arm']],
            plan_sha256=sha(raw), phases=phases, production_changed=False,
            reporter_requests_early_exit=False, new_coordinates_or_physics=False))
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
    results, started = [], time.monotonic()
    for index, case in enumerate(record['cases']):
        folder = (output/case['id']/f"seed-{case['seed']}"/case['conditioning']/case['arm']).resolve()
        folder.mkdir(parents=True, exist_ok=False)
        row = dict(id=case['id'], seed=case['seed'], conditioning=case['conditioning'], arm=case['arm'],
            status='budget_not_executed', returncode=None)
        remaining = record['policy']['total_seconds']-(time.monotonic()-started)
        if remaining > 0:
            command = [record['python'], str(Path(__file__).resolve()), 'worker', '--plan', str(plan_path.resolve()),
                '--index', str(index), '--output', str(folder)]
            with (folder/'run.log').open('xb') as log:
                try:
                    process = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT,
                        timeout=min(remaining, record['policy']['per_cell_seconds']), check=False)
                    row.update(status='completed' if process.returncode == 0 else 'execution_blocker', returncode=process.returncode)
                except subprocess.TimeoutExpired:
                    row['status'] = 'timeout'
        row['artifacts'] = {p.name:sha(p.read_bytes()) for p in sorted(folder.iterdir()) if p.is_file()}
        results.append(row)
        write(folder/'result.json', row)
        print(json.dumps({k:v for k, v in row.items() if k != 'artifacts'}), flush=True)
    verify(record)
    write(output/'evidence.json', dict(schema='dvbfixer.carbon-H-whole-path-run.v1', cases=results,
        plan_sha256=sha(raw), production_changed=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('plan')
    for name in ('attribution-plan', 'attribution-audit', 'whole-plan', 'output'):
        p.add_argument('--'+name, type=Path, required=True)
    for name in ('run', 'worker'):
        p = sub.add_parser(name)
        p.add_argument('--plan', type=Path, required=True)
        p.add_argument('--output', type=Path, required=True)
        if name == 'worker':
            p.add_argument('--index', type=int, required=True)
    args = parser.parse_args()
    if args.command == 'plan':
        plan(args.attribution_plan, args.attribution_audit, args.whole_plan, args.output)
    elif args.command == 'run':
        run(args.plan, args.output)
    else:
        worker(args.plan, args.index, args.output)
