"""Versioned capture fix: preserve ForceField signature for OpenMM ArgTracker.

The original collector and failed run remain immutable. No sampling settings,
seeds, source chemistry or force parameters are changed by this amendment.
"""
from __future__ import annotations

import argparse
import functools
import inspect
import json
import random
import sys
from pathlib import Path
from unittest.mock import patch

from scripts import capture_methyl_diversity as v1
from scripts.branched_template_baseline import preflight
from scripts.graph_methyl_geometry import discover
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify
from scripts.research_legacy_start_controls import compressed, snapshot_atoms


def forcefield_observer(native, active, systems):
    @functools.wraps(native)
    def observe(forcefield, topology, *args, **kwargs):
        system = native(forcefield, topology, *args, **kwargs)
        if active():
            systems.append(dict(pointer=int(system.this), topology=topology, atoms=list(topology.atoms()),
                rows=v1.labels(topology), bonds=[[a.index, b.index] for a, b in topology.bonds()]))
        return system
    return observe


def plan(original_plan, failed_root, output):
    raw = original_plan.read_bytes()
    record = json.loads(raw)
    verify(record)
    failed_raw = (failed_root/'evidence.json').read_bytes()
    failed = json.loads(failed_raw)
    expected = [(c['id'], s) for c in record['cases'] for s in record['policy']['seeds']]
    if failed['plan_sha256'] != sha(raw) or [(c['id'], c['seed']) for c in failed['cases']] != expected:
        raise ValueError('complete original failed capture cohort not bound')
    inputs = dict(record['inputs'])
    inputs[str(original_plan)] = sha(raw)
    inputs[str(failed_root/'evidence.json')] = sha(failed_raw)
    for case in failed['cases']:
        folder = failed_root/case['id']/f"seed-{case['seed']}"
        for name, digest in case['artifacts'].items():
            path = folder/name
            if sha(path.read_bytes()) != digest:
                raise ValueError('failed capture artifacts changed')
            inputs[str(path)] = digest
        if case['status'] == 'execution_blocker':
            if any(name in case['artifacts'] for name in ('H-system.xml', 'H-path.json.gz', 'whole-state.xml')):
                raise ValueError('amendment must not silently resample an already executed H stage')
            if 'flexibleConstraints' not in (folder/'run.log').read_text():
                raise ValueError('undeclared failed capture reason')
        elif case['status'] != 'refused_source':
            raise ValueError('undeclared original capture status')
    record.update(schema='dvbfixer.methyl-diversity-capture-plan.v2', inputs=inputs)
    record['source_hashes'] = dict(record['source_hashes'])
    record['source_hashes']['scripts/capture_methyl_diversity_v2.py'] = sha(Path(__file__).read_bytes())
    record['amendment'] = dict(original_plan_sha256=sha(raw), failed_run_sha256=sha(failed_raw),
        reason='functools.wraps preserves native createSystem signature/defaults used by OpenMM ArgTracker; original ten blockers precede H minimization',
        source_or_scientific_protocol_change=False)
    write(output, record)


def worker(plan_path, label, seed, folder):
    import openmm
    from openmm.app import ForceField, Modeller, Simulation
    from openmm.unit import dalton, nanometer
    from pdbfixer import PDBFixer

    from dvbfixer.cli import main

    record = json.loads(plan_path.read_bytes())
    verify(record)
    if (seed not in record['policy']['seeds'] or openmm.__version__ != record['openmm_version']
            or preflight() != record['runtime']):
        raise ValueError('capture seed/runtime/resources differ')
    case = next(c for c in record['cases'] if c['id'] == label)
    if case['preparation_status'] != 'prepared':
        raise ValueError('source/preparation refusal must not execute')
    for name, path in (('prepared.pdb', case['prepared']), ('prepared.dat', case['dat'])):
        with (folder/name).open('xb') as stream:
            stream.write(Path(path).read_bytes())
    native_h, native_f = Modeller.addHydrogens, ForceField.createSystem
    native_min, native_sim, native_fix = openmm.LocalEnergyMinimizer.minimize, Simulation.minimizeEnergy, PDBFixer.__init__
    reference = openmm.Platform.getPlatformByName('Reference')
    active, systems, events, internal_atoms = False, [], [], []

    def fixer_init(fixer, *args, **kwargs):
        native_fix(fixer, *args, **kwargs)
        fixer.platform = reference

    create_system = forcefield_observer(native_f, lambda: active, systems)
    if inspect.signature(create_system) != inspect.signature(native_f):
        raise ValueError('observer lost native ForceField signature')

    def observe(context, *args, **kwargs):
        nonlocal internal_atoms
        if not active or args != (1.0, 50) or kwargs or events:
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
        events.append(dict(rows=rows, bonds=bonds, mobile=mobile))
        return result

    @functools.wraps(native_h)
    def hydrogen(model, *args, **kwargs):
        nonlocal active
        if active or events:
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
        if len(events) != 1 or random.getstate() != saved:
            raise ValueError('H observer/random-state restoration differs')
        event = events[0]
        original_keys = {tuple(r['key']) for r in original}
        if len(original_keys) != len(original):
            raise ValueError('ambiguous pre-H identities')
        new_h = [i for i, row in enumerate(event['rows']) if row['element'] in {'H', 'D'} and tuple(row['key']) not in original_keys]
        event.update(new_h=new_h, post_H_lineage=v1.object_lineage(internal_atoms, model.topology),
            post_H_labels=v1.labels(model.topology))
        event['graph'] = discover(event['rows'], event['bonds'], set(new_h), set(event['mobile']))
        write(folder/'H-ownership.json', event)
        compressed(folder/'H-return.json.gz', dict(atoms=snapshot_atoms(model.topology, model.positions)))
        return result

    def boundary(simulation, *args, **kwargs):
        if len(events) != 1:
            raise ValueError('whole boundary without exactly one H context')
        lineage = v1.object_lineage(internal_atoms, simulation.topology)
        state = simulation.context.getState(positions=True, velocities=True, parameters=True, integratorParameters=True)
        for name, obj in (('system', simulation.system), ('integrator', simulation.integrator), ('state', state)):
            (folder/f'whole-{name}.xml').write_text(openmm.XmlSerializer.serialize(obj))
        write(folder/'whole-topology.json', dict(atoms=v1.labels(simulation.topology),
            bonds=[[a.index, b.index] for a, b in simulation.topology.bonds()], H_to_whole_indices=lineage))
        compressed(folder/'whole-start.json.gz', dict(atoms=snapshot_atoms(simulation.topology, state.getPositions()),
            positions_nm=[list(p) for p in state.getPositions().value_in_unit(nanometer)]))
        raise v1.CapturedBoundary()

    stopped, saved_argv = False, sys.argv[:]
    try:
        sys.argv = ['dvbfixer', 'minimize', str(folder/'prepared.pdb'), '--dat', str(folder/'prepared.dat'),
            '--rebuild-h', '--no-solvent', '--platform', 'Reference', '--atom-naming', 'standard', '-o', str(folder/'not-published.pdb')]
        with (patch.object(Modeller, 'addHydrogens', hydrogen), patch.object(ForceField, 'createSystem', create_system),
                patch.object(Simulation, 'minimizeEnergy', boundary), patch.object(PDBFixer, '__init__', fixer_init)):
            try:
                main()
            except v1.CapturedBoundary:
                stopped = True
    finally:
        sys.argv = saved_argv
    if (not stopped or Modeller.addHydrogens is not native_h or ForceField.createSystem is not native_f
            or openmm.LocalEnergyMinimizer.minimize is not native_min or Simulation.minimizeEnergy is not native_sim
            or PDBFixer.__init__ is not native_fix or (folder/'not-published.pdb').exists()):
        raise ValueError('capture boundary/hooks/output isolation failed')
    write(folder/'capture.json', dict(id=label, seed=seed, plan_sha256=sha(plan_path.read_bytes()),
        hooks_restored=True, whole_minimization_executed=False, H_calls=1,
        groups=len(events[0]['graph']['groups']), production_changed=False, new_H_control_executed=False))


def run(plan_path, output):
    # Borrow only the unchanged bounded subprocess orchestrator. Point its child
    # entry at this version; the original file and failed plan stay unchanged.
    with patch.object(v1, '__file__', str(Path(__file__).resolve())):
        v1.run(plan_path, output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('plan')
    p.add_argument('--original-plan', type=Path, required=True)
    p.add_argument('--failed-root', type=Path, required=True)
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
        plan(args.original_plan, args.failed_root, args.output)
    elif args.command == 'run':
        run(args.plan, args.output)
    else:
        worker(args.plan, args.id, args.seed, args.output)
