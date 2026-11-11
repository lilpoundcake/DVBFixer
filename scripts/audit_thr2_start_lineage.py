"""Saved-only matched-start audit: distinguish H geometry from heavy/FF changes."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from scripts.audit_thr2_native_capture import delta, read, verify_files
from scripts.research_dunbrack_prior import sha, write
from scripts.research_thr2_native_capture import check


def audit(plan_path, root, output):
    import openmm
    from openmm.unit import nanometer

    plan = read(plan_path)
    check(plan)
    execution = read(root/'evidence.json')
    if execution['plan_sha256'] != sha(plan_path.read_bytes()):
        raise ValueError('native execution binding differs')
    if [r['label'] for r in execution['native']] != plan['protocol']['labels']:
        raise ValueError('native label enrollment differs')
    reference, comparisons, bindings = None, [], {}
    for result in execution['native']:
        folder = root/result['label']
        verify_files(folder, result)
        capture = read(folder/'capture.json')
        whole = next(c for c in capture['calls'] if c['stage'] == 'whole')
        prefix = whole['before']['prefix']
        files = {name:folder/f'{prefix}-{name}.xml' for name in ('system', 'integrator', 'state')}
        atom_file = folder/f'{prefix}-atoms.json.gz'
        topology_file = folder/f"call-{whole['number']}-topology.json"
        for path in [*files.values(), atom_file, topology_file, folder/'capture.json']:
            bindings[str(path)] = sha(path.read_bytes())
        state = openmm.XmlSerializer.deserialize(files['state'].read_text())
        rows = read(atom_file)['atoms']
        keys = [row['key'] for row in rows]
        if len({tuple(k) for k in keys}) != len(keys):
            raise ValueError('duplicate native atom identity')
        xyz = state.getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
        if len(xyz) != len(keys):
            raise ValueError('State/topology inventory differs')
        data = dict(keys=keys, elements=[row['element'] for row in rows], xyz=xyz,
            system=files['system'].read_bytes(), integrator=files['integrator'].read_bytes(),
            topology=topology_file.read_bytes(), parameters=dict(state.getParameters()))
        if reference is None:
            reference = data
        if data['keys'] != reference['keys'] or data['elements'] != reference['elements']:
            raise ValueError('matched-start complete atom identity differs')
        comparison = dict(label=result['label'], reference_label=plan['protocol']['labels'][0],
            system_xml_sha256=sha(data['system']), integrator_xml_sha256=sha(data['integrator']),
            exact_system_xml=data['system'] == reference['system'],
            exact_integrator_xml=data['integrator'] == reference['integrator'],
            exact_bond_inventory=data['topology'] == reference['topology'],
            exact_context_parameters=data['parameters'] == reference['parameters'])
        for name, hydrogen in (('heavy', False), ('hydrogen', True)):
            indices = [i for i, element in enumerate(data['elements']) if (element in ('H', 'D')) == hydrogen]
            comparison[name] = dict(count=len(indices), **delta([reference['xyz'][i] for i in indices], [xyz[i] for i in indices]))
        comparisons.append(comparison)
    write(output, dict(schema='dvbfixer.thr2-start-lineage.v1', plan_sha256=sha(plan_path.read_bytes()),
        execution_sha256=sha((root/'evidence.json').read_bytes()),
        audit_script_sha256=sha(Path(__file__).read_bytes()), inputs=bindings, comparisons=comparisons,
        extra_minimization=False, production_changed=False,
        limits='matched saved-start evidence, not isolation of Python RNG from native H minimization numerical variability'))
    print(json.dumps(comparisons, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('plan', 'root', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan.resolve(), args.root.resolve(), args.output.resolve())
