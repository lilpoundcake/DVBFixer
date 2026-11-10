"""Private Reference H-context collection; stop before first whole minimization.

Actual ForceField topology objects bind internal H labels to post-H/full labels.
This captures original native H placement only, not the graph intervention arms.
"""
from __future__ import annotations

import argparse
import functools
import inspect
import json
import os
import random
import subprocess
import sys
import time
from pathlib import Path
from unittest.mock import patch

from scripts.graph_methyl_geometry import discover
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify
from scripts.research_legacy_start_controls import compressed, snapshot_atoms


class CapturedBoundary(Exception):
    """Intentional non-error exit before executing any whole minimization."""


def labels(topology):
    return [dict(key=[a.residue.chain.id, str(a.residue.id), a.residue.insertionCode.strip(), a.residue.name, a.name],
        element=a.element.symbol if a.element else None) for a in topology.atoms()]


def object_lineage(original_atoms, final_topology):
    final = list(final_topology.atoms())
    index = {id(a): i for i, a in enumerate(final)}
    if len(index) != len(final) or len(final) != len(original_atoms) or any(id(a) not in index for a in original_atoms):
        raise ValueError('native topology object lineage lost; never guess label aliases')
    return [index[id(a)] for a in original_atoms]


def plan(source_path, prepared, output):
    import openmm
    from openmm.app import Modeller

    raw = source_path.read_bytes()
    source = json.loads(raw)
    verify(source)
    evidence_raw = (prepared/'evidence.json').read_bytes()
    evidence = json.loads(evidence_raw)
    if evidence['source_plan_sha256'] != sha(raw) or [c['id'] for c in evidence['cases']] != [c['id'] for c in source['cases']]:
        raise ValueError('complete preparation source cohort differs')
    inputs = {**source['inputs'], str(source_path): sha(raw), str(prepared/'evidence.json'): sha(evidence_raw)}
    cases = []
    for case, result in zip(source['cases'], evidence['cases'], strict=True):
        folder = prepared/case['id']
        for name, digest in result['artifacts'].items():
            path = folder/name
            if sha(path.read_bytes()) != digest:
                raise ValueError('prepared source artifact changed')
            inputs[str(path)] = digest
        cases.append(dict(id=case['id'], source_status=case['status'], preparation_status=result['status'],
            prepared=str(folder/'prepared.pdb'), dat=str(folder/'prepared.dat')))
    sources = dict(source['source_hashes'])
    for name in ('scripts/capture_methyl_diversity.py', 'scripts/graph_methyl_geometry.py',
                 'scripts/research_legacy_start_controls.py'):
        sources[name] = sha(Path(name).read_bytes())
    inputs[inspect.getfile(Modeller)] = sha(Path(inspect.getfile(Modeller)).read_bytes())
    write(output, dict(schema='dvbfixer.methyl-diversity-capture-plan.v1', cases=cases, inputs=inputs,
        source_hashes=sources, runtime=source['runtime'], openmm_version=openmm.__version__, python=sys.executable,
        policy=dict(seeds=[1, 2], per_cell_seconds=1800, total_seconds=7200,
            capture='original H System/Integrator/nm start/endpoint, native topology object lineage and first full-System start; no whole minimization',
            settings='minimize --rebuild-h --no-solvent --platform Reference --atom-naming standard; H seed1/2 and Reference explicitly forced',
            graph='bind generic implementation before capture; discover only, no intervention arms/phase grids evaluated'),
        production_changed=False, new_H_control_executed=False))


def worker(plan_path, label, seed, folder):
    import openmm
    from openmm.app import ForceField, Modeller, Simulation
    from openmm.unit import dalton, nanometer
    from pdbfixer import PDBFixer

    from dvbfixer.cli import main

    record = json.loads(plan_path.read_bytes())
    verify(record)
    if seed not in record['policy']['seeds'] or openmm.__version__ != record['openmm_version']:
        raise ValueError('capture seed/runtime differs')
    case = next(c for c in record['cases'] if c['id'] == label)
    if case['preparation_status'] != 'prepared':
        raise ValueError('source/preparation refusal must not execute')
    for name, path in (('prepared.pdb', case['prepared']), ('prepared.dat', case['dat'])):
        with (folder/name).open('xb') as stream:
            stream.write(Path(path).read_bytes())
    native_h, native_f = Modeller.addHydrogens, ForceField.createSystem
    native_min, native_sim, native_fix = openmm.LocalEnergyMinimizer.minimize, Simulation.minimizeEnergy, PDBFixer.__init__
    reference = openmm.Platform.getPlatformByName('Reference')
    active, systems, h_events = False, [], []
    internal_atoms = []

    def fixer_init(fixer, *args, **kwargs):
        native_fix(fixer, *args, **kwargs)
        fixer.platform = reference

    def create_system(forcefield, topology, *args, **kwargs):
        system = native_f(forcefield, topology, *args, **kwargs)
        if active:
            systems.append(dict(pointer=int(system.this), topology=topology, atoms=list(topology.atoms()),
                rows=labels(topology), bonds=[[a.index, b.index] for a, b in topology.bonds()]))
        return system

    def observe(context, *args, **kwargs):
        nonlocal internal_atoms
        if not active or args != (1.0, 50) or kwargs or h_events:
            raise ValueError('unexpected native H minimizer call/settings')
        system = context.getSystem()
        matches = [s for s in systems if s['pointer'] == int(system.this)]
        if len(matches) != 1 or context.getPlatform().getName() != 'Reference':
            raise ValueError('H Context lacks unique actual topology/System lineage')
        observed = matches[0]
        internal_atoms = observed['atoms']
        rows, bonds = observed['rows'], observed['bonds']
        if system.getNumParticles() != len(rows):
            raise ValueError('H particle/topology inventory differs')
        xml = openmm.XmlSerializer.serialize(system)
        (folder/'H-system.xml').write_text(xml)
        (folder/'H-integrator.xml').write_text(openmm.XmlSerializer.serialize(context.getIntegrator()))
        before = [list(p) for p in context.getState(positions=True).getPositions().value_in_unit(nanometer)]
        write(folder/'H-topology.json', dict(atoms=rows, bonds=bonds))
        result = native_min(context, *args, **kwargs)
        after = [list(p) for p in context.getState(positions=True).getPositions().value_in_unit(nanometer)]
        if openmm.XmlSerializer.serialize(system) != xml:
            raise ValueError('H observer changed System')
        compressed(folder/'H-path.json.gz', dict(initial_nm=before, final_nm=after))
        mobile = [i for i in range(len(rows)) if system.getParticleMass(i).value_in_unit(dalton) > 0]
        h_events.append(dict(rows=rows, bonds=bonds, mobile=mobile))
        return result

    @functools.wraps(native_h)
    def hydrogen(model, *args, **kwargs):
        nonlocal active
        if active or h_events:
            raise ValueError('second/reentrant native H call')
        bound = inspect.signature(native_h).bind(model, *args, **kwargs)
        bound.apply_defaults()
        if bound.arguments['platform'] is not None:
            raise ValueError('unexpected H platform argument')
        original = snapshot_atoms(model.topology, model.positions)
        compressed(folder/'H-input.json.gz', dict(atoms=original))
        saved = random.getstate()
        random.seed(seed)
        active = True
        try:
            with patch.object(openmm.LocalEnergyMinimizer, 'minimize', observe):
                result = native_h(model, *args, **dict(kwargs, platform=reference))
        finally:
            active = False
            random.setstate(saved)
        if len(h_events) != 1 or random.getstate() != saved:
            raise ValueError('H observer/random-state restoration differs')
        event = h_events[0]
        original_keys = {tuple(r['key']) for r in original}
        if len(original_keys) != len(original):
            raise ValueError('ambiguous pre-H identities')
        new_h = [i for i, row in enumerate(event['rows']) if row['element'] in {'H', 'D'} and tuple(row['key']) not in original_keys]
        event['new_h'] = new_h
        event['post_H_lineage'] = object_lineage(internal_atoms, model.topology)
        event['post_H_labels'] = labels(model.topology)
        event['graph'] = discover(event['rows'], event['bonds'], set(new_h), set(event['mobile']))
        write(folder/'H-ownership.json', event)
        compressed(folder/'H-return.json.gz', dict(atoms=snapshot_atoms(model.topology, model.positions)))
        return result

    def boundary(simulation, *args, **kwargs):
        if len(h_events) != 1:
            raise ValueError('whole boundary without exactly one H context')
        lineage = object_lineage(internal_atoms, simulation.topology)
        state = simulation.context.getState(positions=True, velocities=True, parameters=True, integratorParameters=True)
        for name, obj in (('system', simulation.system), ('integrator', simulation.integrator), ('state', state)):
            (folder/f'whole-{name}.xml').write_text(openmm.XmlSerializer.serialize(obj))
        write(folder/'whole-topology.json', dict(atoms=labels(simulation.topology),
            bonds=[[a.index, b.index] for a, b in simulation.topology.bonds()], H_to_whole_indices=lineage))
        compressed(folder/'whole-start.json.gz', dict(atoms=snapshot_atoms(simulation.topology, state.getPositions()),
            positions_nm=[list(p) for p in state.getPositions().value_in_unit(nanometer)]))
        raise CapturedBoundary()

    stopped = False
    saved_argv = sys.argv[:]
    try:
        sys.argv = ['dvbfixer', 'minimize', str(folder/'prepared.pdb'), '--dat', str(folder/'prepared.dat'),
            '--rebuild-h', '--no-solvent', '--platform', 'Reference', '--atom-naming', 'standard', '-o', str(folder/'not-published.pdb')]
        with (patch.object(Modeller, 'addHydrogens', hydrogen), patch.object(ForceField, 'createSystem', create_system),
                patch.object(Simulation, 'minimizeEnergy', boundary), patch.object(PDBFixer, '__init__', fixer_init)):
            try:
                main()
            except CapturedBoundary:
                stopped = True
    finally:
        sys.argv = saved_argv
    if (not stopped or Modeller.addHydrogens is not native_h or ForceField.createSystem is not native_f
            or openmm.LocalEnergyMinimizer.minimize is not native_min or Simulation.minimizeEnergy is not native_sim
            or PDBFixer.__init__ is not native_fix or (folder/'not-published.pdb').exists()):
        raise ValueError('capture boundary/hooks/output isolation failed')
    write(folder/'capture.json', dict(id=label, seed=seed, plan_sha256=sha(plan_path.read_bytes()),
        hooks_restored=True, whole_minimization_executed=False, H_calls=1,
        groups=len(h_events[0]['graph']['groups']), production_changed=False, new_H_control_executed=False))


def run(plan_path, root):
    raw = plan_path.read_bytes()
    record = json.loads(raw)
    verify(record)
    root.mkdir(parents=True, exist_ok=False)
    (root/'executed-plan.json').write_bytes(raw)
    env = dict(os.environ, PATH=str(Path(record['python']).parent)+os.pathsep+os.environ.get('PATH', ''),
        PYTHONPATH=os.pathsep.join((str(Path('src').resolve()), str(Path('.').resolve()))))
    started, results = time.monotonic(), []
    for case in record['cases']:
        for seed in record['policy']['seeds']:
            folder = (root/case['id']/f'seed-{seed}').resolve()
            folder.mkdir(parents=True, exist_ok=False)
            row = dict(id=case['id'], seed=seed, status=case['preparation_status'], returncode=None)
            if case['preparation_status'] == 'prepared':
                remaining = record['policy']['total_seconds']-(time.monotonic()-started)
                if remaining <= 0:
                    row['status'] = 'budget_not_executed'
                else:
                    command = [record['python'], str(Path(__file__).resolve()), 'worker', '--plan', str(plan_path.resolve()),
                        '--id', case['id'], '--seed', str(seed), '--output', str(folder)]
                    with (folder/'run.log').open('xb') as log:
                        try:
                            child = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT,
                                timeout=min(remaining, record['policy']['per_cell_seconds']), check=False)
                            row.update(status='captured' if child.returncode == 0 else 'execution_blocker', returncode=child.returncode)
                        except subprocess.TimeoutExpired:
                            row['status'] = 'timeout'
            row['artifacts'] = {p.name: sha(p.read_bytes()) for p in sorted(folder.iterdir()) if p.is_file()}
            results.append(row)
            write(folder/'result.json', row)
            print(json.dumps({k: v for k, v in row.items() if k != 'artifacts'}), flush=True)
    verify(record)
    write(root/'evidence.json', dict(schema='dvbfixer.methyl-diversity-capture-run.v1', cases=results,
        plan_sha256=sha(raw), production_changed=False, new_H_control_executed=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('plan')
    p.add_argument('--source-plan', type=Path, required=True)
    p.add_argument('--prepared', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    for name in ('run', 'worker'):
        p = sub.add_parser(name)
        p.add_argument('--plan', type=Path, required=True)
        p.add_argument('--output', type=Path, required=True)
        if name == 'worker':
            p.add_argument('--id', required=True)
            p.add_argument('--seed', type=int, choices=[1, 2], required=True)
    args = parser.parse_args()
    if args.command == 'plan':
        plan(args.source_plan, args.prepared, args.output)
    elif args.command == 'run':
        run(args.plan, args.output)
    else:
        worker(args.plan, args.id, args.seed, args.output)
