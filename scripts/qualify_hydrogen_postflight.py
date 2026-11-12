"""Matched saved/local-graph baseline versus candidate; no H sampling or optimizer."""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
from pathlib import Path

from scripts.audit_thr2_native_capture import read
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify
from scripts.research_legacy_start_controls import compressed, snapshot_atoms

ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT/'tmp/hydrogen-postflight-contract-01'


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ValueError('cannot load explicitly bound qualification module')
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def freeze(output):
    import openmm.version

    cases = load('hydrogen_postflight_cases', ROOT/'tests/test_hydrogen_postflight.py')
    baseline = EVIDENCE/'baseline_geometry.py'
    protocol = EVIDENCE/'protocol-before-implementation.md'
    with protocol.open('xb') as stream:
        stream.write((ROOT/'docs/plans/hydrogen-postflight-contract.md').read_bytes())
    inputs = {str(p):sha(p.read_bytes()) for p in (baseline, EVIDENCE/'baseline.json', ROOT/'tests/test_hydrogen_postflight.py', protocol)}
    native = ROOT/'tmp/pdbfixer-research/ser-hg-capture-run-01'
    raw = (native/'evidence.json').read_bytes()
    execution = json.loads(raw)
    inputs[str(native/'evidence.json')] = sha(raw)
    for row in execution['cases']:
        for name, digest in row['artifacts'].items():
            p = native/row['label']/name
            if sha(p.read_bytes()) != digest:
                raise ValueError('saved native SER artifact changed')
            inputs[str(p)] = digest
    write(output, dict(schema='dvbfixer.hydrogen-postflight-qualification-plan.v1',
        baseline=str(baseline), cases_source=str(ROOT/'tests/test_hydrogen_postflight.py'),
        native_root=str(native), inputs=inputs,
        source_hashes={str(Path(__file__).resolve()):sha(Path(__file__).read_bytes())},
        runtime_revision=openmm.version.git_revision,
        cases=[dict(resname=r, distance=d) for r in cases.CHEMISTRIES for d in cases.DISTANCES],
        native_labels=[r['label'] for r in execution['cases']],
        policy=dict(upper_oh_angstrom=1.2, equality_is_failure=True, local_graphs=48, saved_native_states=6,
            heavy_and_valid_sibling_motion=False, native_sampling=False, force_field_change=False,
            retention='every old/new outcome, errors and private coordinates; no candidate retry or threshold search',
            limits='48 local graph controls are not independent physical structures; saved SER raw states are one synthetic two-residue input')))


def reconstruct(rows, bonds):
    from openmm import Vec3, unit
    from openmm.app import Topology, element

    topology = Topology()
    chains, residues, atoms = {}, {}, []
    for row in rows:
        chain, resid, icode, resname, atom = row['key']
        if chain not in chains:
            chains[chain] = topology.addChain(chain)
        key = chain, resid, icode, resname
        if key not in residues:
            residues[key] = topology.addResidue(resname, chains[chain], resid, icode)
        atoms.append(topology.addAtom(atom, element.Element.getBySymbol(row['element']), residues[key]))
    for a, b in bonds:
        topology.addBond(atoms[a], atoms[b])
    positions = [Vec3(*(x/10 for x in row['xyz'])) for row in rows]*unit.nanometer
    return topology, positions


def metrics(module, topology, positions):
    from openmm import unit

    points = positions.value_in_unit(unit.angstrom)
    distances = []
    for a, b in topology.bonds():
        if a.element.symbol == 'O' and b.element.symbol == 'H':
            parent, hydrogen = a, b
        elif b.element.symbol == 'O' and a.element.symbol == 'H':
            parent, hydrogen = b, a
        else:
            continue
        distances.append(math.dist(points[parent.index], points[hydrogen.index]))
    return dict(oh_distances_angstrom=distances, all_oh_below_upper_bound=all(d < 1.2 for d in distances),
        findings=[dict(hydrogen=f.h_atom.name, parent=f.parent.name, reason=f.reason)
            for f in module.detect_misplaced_hydrogens(topology, positions)])


def run(plan_path, output, public_output):
    import openmm.version

    from dvbfixer.ffutils import geometry as candidate

    plan = read(plan_path)
    verify(plan)
    if openmm.version.git_revision != plan['runtime_revision']:
        raise ValueError('qualification runtime differs')
    candidate_path = Path(candidate.__file__).resolve()
    candidate_raw = candidate_path.read_bytes()
    baseline = load('hydrogen_postflight_baseline', Path(plan['baseline']))
    cases = load('hydrogen_postflight_cases', Path(plan['cases_source']))
    output.mkdir(parents=True, exist_ok=False)
    (output/'executed-plan.json').write_bytes(plan_path.read_bytes())
    enrollment = []
    for case in plan['cases']:
        top, pos = cases.oh_graph(case['resname'], case['distance'])
        enrollment.append(dict(kind='local-graph', **case, rows=snapshot_atoms(top, pos),
            bonds=[[a.index, b.index] for a, b in top.bonds()]))
    native = Path(plan['native_root'])
    for label in plan['native_labels']:
        source = native/label
        enrollment.append(dict(kind='saved-native-SER', label=label,
            rows=read(source/'repair-1-before-atoms.json.gz')['atoms'], bonds=read(source/'repair-1-before-topology.json')['bonds']))
    reports = []
    for index, case in enumerate(enrollment):
        folder = output/f'case-{index:02d}'
        folder.mkdir()
        compressed(folder/'input.json.gz', dict(atoms=case['rows'], bonds=case['bonds']))
        report = {k:v for k, v in case.items() if k not in ('rows', 'bonds')}
        report['index'] = index
        report['arms'] = {}
        for label, module in (('baseline', baseline), ('candidate', candidate)):
            topology, positions = reconstruct(case['rows'], case['bonds'])
            before = snapshot_atoms(topology, positions)
            outcome = dict(before=metrics(module, topology, positions), status='completed')
            try:
                outcome['repairs'] = module.repair_misplaced_hydrogens(topology, positions)
            except Exception as error:
                outcome.update(status='error', error_type=type(error).__name__, error=str(error))
            after = snapshot_atoms(topology, positions)
            compressed(folder/f'{label}.json.gz', dict(atoms=after))
            outcome['after'] = metrics(module, topology, positions)
            outcome['identity_unchanged'] = [a['key'] for a in before] == [a['key'] for a in after]
            outcome['heavy_coordinates_unchanged'] = all(a['xyz'] == b['xyz'] for a, b in zip(before, after, strict=True) if a['element'] != 'H')
            outcome['changed_atom_keys'] = [a['key'] for a, b in zip(before, after, strict=True) if a['xyz'] != b['xyz']]
            if label == 'candidate' and case['kind'] == 'local-graph':
                outcome['valid_input_unchanged'] = before == after if case['distance'] < 1.2 else None
                outcome['unselected_atoms_unchanged'] = all(a['xyz'] == b['xyz'] for i, (a, b) in enumerate(zip(before, after, strict=True)) if i != 1)
            report['arms'][label] = outcome
        reports.append(report)
    candidates = [r['arms']['candidate'] for r in reports]
    passed = all(a['status'] == 'completed' and a['after']['all_oh_below_upper_bound'] and a['heavy_coordinates_unchanged']
        and a['identity_unchanged'] and a.get('unselected_atoms_unchanged', True) and a.get('valid_input_unchanged') is not False for a in candidates)
    summary = dict(local_graphs=48, saved_native_states=6, candidate_errors=sum(a['status'] != 'completed' for a in candidates),
        candidate_oh_bound_pass=sum(a['after']['all_oh_below_upper_bound'] for a in candidates),
        baseline_oh_bound_pass=sum(r['arms']['baseline']['after']['all_oh_below_upper_bound'] for r in reports),
        candidate_heavy_and_identity_preserved=all(a['heavy_coordinates_unchanged'] and a['identity_unchanged'] for a in candidates),
        qualification_passed=passed)
    write(public_output, dict(schema='dvbfixer.hydrogen-postflight-qualification.v1', plan_sha256=sha(plan_path.read_bytes()),
        candidate_source_sha256=sha(candidate_raw), baseline_source_sha256=sha(Path(plan['baseline']).read_bytes()),
        qualification_script_sha256=sha(Path(__file__).read_bytes()), summary=summary, cases=reports,
        new_native_sampling=False, force_field_changed=False, limits=plan['policy']['limits']))
    verify(plan)
    if candidate_path.read_bytes() != candidate_raw:
        raise ValueError('candidate source changed during qualification')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('plan')
    p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('run')
    for name in ('plan', 'output', 'public-output'):
        p.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'plan':
        freeze(args.output.resolve())
    else:
        run(args.plan.resolve(), args.output.resolve(), args.public_output.resolve())
