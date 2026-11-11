"""Saved-state native-force decomposition, not an intervention or trajectory replay.

Enrollment: all initial/callback/returned nodes in all 18 bounded replays. Only
measurement force-group metadata changes, with exact restoration verification.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from scripts.audit_stereo_force_projection import (
    grouped_system,
    projected_force,
    volume_gradient,
)
from scripts.audit_thr2_native_capture import read, verify_files
from scripts.research_dunbrack_prior import sha, write
from scripts.research_thr2_native_capture import check


def freeze(capture_plan, root, output):
    plan = read(capture_plan)
    check(plan)
    evidence = read(root/'evidence.json')
    inputs = {str(capture_plan):sha(capture_plan.read_bytes()),
        str(root/'evidence.json'):sha((root/'evidence.json').read_bytes())}
    for row in evidence['native']:
        folder = root/row['label']
        verify_files(folder, row)
        inputs.update({str(folder/name):digest for name, digest in row['artifacts'].items()})
    for row in evidence['replays']:
        folder = root/'replays'/row['label']/f"call-{row['number']}"/row['mode']
        verify_files(folder, row)
        inputs.update({str(folder/name):digest for name, digest in row['artifacts'].items()})
    write(output, dict(schema='dvbfixer.thr2-force-terms-plan.v1', capture_plan=str(capture_plan), root=str(root),
        inputs=inputs, source_hashes={str(Path(__file__).resolve()):sha(Path(__file__).read_bytes()),
            str(Path('scripts/audit_thr2_native_capture.py').resolve()):sha(Path('scripts/audit_thr2_native_capture.py').read_bytes())},
        enrollment='ALL18 original saved replays, all initial/callback/returned states; no outcome selection',
        measurement='same platform/properties and State parameters; force-group metadata only; CA and CB signed-volume unit-gradient force projections',
        publication='private every-node terms; public initial/first sign-crossing pair/last-callback/returned and complete per-path extrema',
        limits='local projections are not optimizer steps, integrated work, energy-barrier certificates or original uncaptured causality',
        extra_minimization=False, production_changed=False))


def audit(plan_path, private_output, public_output):
    import openmm
    from openmm.unit import kilojoules_per_mole, nanometer

    plan = read(plan_path)
    for name, digest in {**plan['inputs'], **plan['source_hashes']}.items():
        if sha(Path(name).read_bytes()) != digest:
            raise ValueError(f'saved-force input/source changed: {name}')
    parent = read(Path(plan['capture_plan']))
    check(parent)
    root = Path(plan['root'])
    results, public, errors = [], [], []
    for row in read(root/'evidence.json')['replays']:
        report = {key:row[key] for key in ('label', 'number', 'mode', 'status')}
        source = root/row['label']
        folder = root/'replays'/row['label']/f"call-{row['number']}"/row['mode']
        if not (folder/'path.json.gz').exists():
            results.append(report)
            public.append(report)
            continue
        call = next(c for c in read(source/'capture.json')['calls'] if c['number'] == row['number'])
        prefix = call['before']['prefix']
        atoms = read(source/f'{prefix}-atoms.json.gz')['atoms']
        keys = {tuple(a['key']):i for i, a in enumerate(atoms)}
        target = ('A', '2', '', 'THR')
        indices = {name:[keys[(*target, atom)] for atom in names] for name, names in (
            ('ca', ('CA', 'N', 'C', 'CB')), ('cb', ('CB', 'CA', 'OG1', 'CG2')))}
        original = openmm.XmlSerializer.deserialize((source/f'{prefix}-system.xml').read_text())
        system, groups = grouped_system(original)
        integrator = openmm.XmlSerializer.deserialize((source/f'{prefix}-integrator.xml').read_text())
        platform = call['before']['platform']
        name = 'Reference' if row['mode'] == 'Reference' else platform['name']
        properties = {} if row['mode'] == 'Reference' else platform['constructor_properties']
        context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName(name), properties)
        context.setState(openmm.XmlSerializer.deserialize((source/f'{prefix}-state.xml').read_text()))
        path = read(folder/'path.json.gz')
        entries = [('initial', path['initial_nm'])]
        count = len(atoms)
        entries += [(f"callback-{cb['ordinal']}", [cb['positions_nm'][3*i:3*i+3] for i in range(count)]) for cb in path['callbacks']]
        entries.append(('returned', path['final_nm']))
        nodes = []
        try:
            for label, points in entries:
                context.setPositions(points*nanometer)
                centres = {}
                for centre, ids in indices.items():
                    volume, gradient, length = volume_gradient([points[i] for i in ids])
                    centres[centre] = dict(volume_angstrom3=1000*volume, gradient=gradient, length=length)
                terms = []
                for group in groups:
                    state = context.getState(energy=True, forces=True, groups=1 << group['index'])
                    forces = state.getForces(asNumpy=True).value_in_unit(kilojoules_per_mole/nanometer)
                    terms.append(dict(index=group['index'], type=group['type'], name=group['name'],
                        energy_kj_mol=state.getPotentialEnergy().value_in_unit(kilojoules_per_mole),
                        **{centre+'_projection_kj_mol_nm':projected_force(forces, indices[centre], c['gradient'], c['length'])
                            for centre, c in centres.items()}))
                all_state = context.getState(energy=True, forces=True)
                forces = all_state.getForces(asNumpy=True).value_in_unit(kilojoules_per_mole/nanometer)
                total = {centre+'_projection_kj_mol_nm':projected_force(forces, indices[centre], c['gradient'], c['length'])
                    for centre, c in centres.items()}
                total['energy_kj_mol'] = all_state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
                err = {key:abs(value-sum(t[key] for t in terms)) for key, value in total.items()}
                if not all(math.isfinite(v) for v in [*total.values(), *err.values()]):
                    raise ValueError('nonfinite force decomposition')
                errors.append(err)
                nodes.append(dict(stage=label, volumes_angstrom3={k:v['volume_angstrom3'] for k, v in centres.items()},
                    total=total, terms=terms, sum_error=err))
        finally:
            del context
        selected = {0, len(nodes)-2, len(nodes)-1}
        for i in range(1, len(nodes)):
            if nodes[i-1]['volumes_angstrom3']['cb']*nodes[i]['volumes_angstrom3']['cb'] < 0:
                selected.update((i-1, i))
                break
        bounds = {key:dict(minimum=min(n['total'][key] for n in nodes), maximum=max(n['total'][key] for n in nodes))
            for key in ('ca_projection_kj_mol_nm', 'cb_projection_kj_mol_nm')}
        report.update(stage=call['stage'], node_count=len(nodes), total_projection_bounds=bounds, groups=groups)
        results.append(dict(report, nodes=nodes))
        public.append(dict(report, selected_nodes=[nodes[i] for i in sorted(selected)]))
    summary = dict(paths=len(results), nodes=sum(r.get('node_count', 0) for r in results),
        maximum_sum_errors={key:max(e[key] for e in errors) for key in errors[0]} if errors else {})
    common = dict(schema='dvbfixer.thr2-force-terms-audit.v1', plan_sha256=sha(plan_path.read_bytes()),
        summary=summary, extra_minimization=False, production_changed=False, limits=plan['limits'])
    write(private_output, dict(common, paths=results))
    write(public_output, dict(common, paths=public, private_audit_sha256=sha(private_output.read_bytes())))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('plan')
    for name in ('capture-plan', 'root', 'output'):
        p.add_argument('--'+name, type=Path, required=True)
    p = sub.add_parser('audit')
    for name in ('plan', 'private-output', 'public-output'):
        p.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'plan':
        freeze(args.capture_plan.resolve(), args.root.resolve(), args.output.resolve())
    else:
        audit(args.plan.resolve(), args.private_output.resolve(), args.public_output.resolve())
