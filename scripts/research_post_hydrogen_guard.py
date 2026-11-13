"""Bounded saved-return proposals: after H-only, no new native H placement."""
from __future__ import annotations

import argparse
import time
from dataclasses import asdict, replace
from pathlib import Path

from scripts.audit_thr2_native_capture import read, verify_files
from scripts.hydrogen_initialization_guard import GuardRefusal, propose, validate
from scripts.qualify_hydrogen_initialization_guard import admit, atoms_from_rows
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify
from scripts.research_legacy_start_controls import compressed
from scripts.research_thr2_hb_only import check_plan

ROOT = Path(__file__).resolve().parents[1]


def bind_return(initial, positions, returned_atoms=None, returned_bonds=None):
    """Retain the original new-H ownership, never infer it from relaxed geometry."""
    if returned_atoms is not None and atoms_from_rows(returned_atoms) != initial.atoms:
        raise GuardRefusal('return atom identity/order differs')
    if returned_bonds is not None:
        before = sorted(tuple(sorted(b)) for b in initial.bonds)
        after = sorted(tuple(sorted(b)) for b in returned_bonds)
        if before != after:
            raise GuardRefusal('return bond graph differs')
    if len(positions) != len(initial.atoms):
        raise GuardRefusal('return coordinate inventory differs')
    returned = replace(initial, positions_nm=tuple(tuple(p) for p in positions))
    validate(returned)
    if any(p != q for a, p, q in zip(initial.atoms, initial.positions_nm, returned.positions_nm, strict=True)
           if a.identity not in initial.new_hydrogens):
        raise GuardRefusal('native H-only return moved an observation')
    return returned


def load_return(case):
    import openmm
    from openmm.unit import nanometer

    source = Path(case['source'])
    prefix = case['initial_prefix']
    original = read(source/f'{prefix}-atoms.json.gz')
    before = read(source/f"h-{case['h_call']}-before.json.gz")
    bonds = read(source/case['bonds'])['bonds']
    system = openmm.XmlSerializer.deserialize((source/f'{prefix}-system.xml').read_text())
    positions = openmm.XmlSerializer.deserialize((source/f'{prefix}-state.xml').read_text()).getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
    initial = admit(original['atoms'], before, positions, bonds, system)
    returned = openmm.XmlSerializer.deserialize(Path(case['returned_state']).read_text()).getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
    if case.get('returned_path'):
        path = read(Path(case['returned_path']))
        if not path['returned'] or path['final_nm'] != returned or path['initial_nm'] != positions:
            raise GuardRefusal('unchanged control return/State/source mismatch')
    snapshot = read(Path(case['returned_atoms'])) if case.get('returned_atoms') else {}
    return bind_return(initial, returned, snapshot.get('atoms'), snapshot.get('bonds'))


def check(plan):
    verify(plan)
    check_plan(read(Path(plan['control_plan'])))


def freeze(output):
    initial_path = ROOT/'docs/research/hydrogen-initialization-guard-static-plan.json'
    control_path = ROOT/'docs/research/thr2-hb-only-plan.json'
    control_evidence_path = ROOT/'tmp/pdbfixer-research/thr2-hb-only-run-01/evidence.json'
    initial, control = read(initial_path), read(control_path)
    verify(initial)
    check_plan(control)
    inputs = {**initial['inputs'], **initial['source_hashes'], **control['inputs'], **control['source_hashes']}
    for p in (initial_path, control_path, control_evidence_path, ROOT/'docs/plans/post-hydrogen-guard.md'):
        inputs[str(p)] = sha(p.read_bytes())
    old_controls = [c for c in read(control_evidence_path)['cases'] if c['arm'] == 'unchanged']
    cases = []
    for source_index, source in enumerate(initial['cases']):
        if len(source['slots']) != 1:
            raise ValueError('original complete one-return-per-case inventory required')
        slot = source['slots'][0]
        base = dict(source_index=source_index, kind=source['kind'], label=source['label'],
                    source=source['folder'], h_call=slot['h_call'], initial_prefix=slot['prefix'],
                    bonds=slot['bonds'], h_calls=source['h_calls'], requested_ff=source['requested_ff'])
        if source['kind'] == 'original-thr2':
            for mode in ('captured', 'Reference'):
                row, = [c for c in old_controls if c['label'] == source['label'] and c['mode'] == mode]
                folder = control_evidence_path.parent/f"path-{row['index']:02d}"
                verify_files(folder, row)
                inputs.update({str(folder/name):digest for name, digest in row['artifacts'].items()})
                cases.append(dict(base, mode=mode, control_folder=str(folder),
                    returned_path=str(folder/'phase-1/path.json.gz'),
                    returned_state=str(folder/'phase-1/returned-state.xml')))
        else:
            capture = read(Path(source['folder'])/'capture.json')
            native, = capture['native']
            prefix = native['after']['prefix']
            cases.append(dict(base, mode='captured',
                returned_state=str(Path(source['folder'])/f'{prefix}-state.xml'),
                returned_atoms=str(Path(source['folder'])/f'{prefix}-atoms.json.gz')))
    if len(cases) != 18:
        raise ValueError('all six matched THR2 and twelve real returns required')
    files = [Path(__file__).resolve(), ROOT/'scripts/audit_post_hydrogen_guard.py',
             ROOT/'tests/test_post_hydrogen_guard.py']
    write(output, dict(schema='dvbfixer.private-post-h-static-plan.v1', cases=cases,
        inputs=inputs, source_hashes={str(p):sha(p.read_bytes()) for p in files},
        control_plan=str(control_path),
        policy=dict(maximum_seconds=120, evaluations=18, contexts=0, new_h_sampling=False,
                    extra_h_relaxation=False, rule='unchanged private-selective-carbon-h-v1'),
        production_changed=False))


def run(plan_path, root, output):
    plan = read(plan_path)
    check(plan)
    root.mkdir(parents=True, exist_ok=False)
    (root/'executed-plan.json').write_bytes(plan_path.read_bytes())
    started, cases = time.monotonic(), []
    for index, case in enumerate(plan['cases']):
        report = dict(index=index, status='not-executed-budget')
        if time.monotonic()-started <= plan['policy']['maximum_seconds']:
            try:
                value = load_return(case)
                report.update(admitted=True, atoms=len(value.atoms), new_h=len(value.new_hydrogens))
                result = propose(value)
                changed = [i for i, (p, q) in enumerate(zip(value.positions_nm, result.positions_nm, strict=True)) if p != q]
                if {value.atoms[i].identity for i in changed} != set(result.moved) or not set(result.moved) <= value.new_hydrogens:
                    raise AssertionError('moved/observation mask differs')
                name = f'case-{index:02d}.json.gz'
                compressed(root/name, dict(atoms=[asdict(a) for a in value.atoms], initial_nm=value.positions_nm,
                    bonds=value.bonds, new_indices=[i for i, a in enumerate(value.atoms) if a.identity in value.new_hydrogens],
                    proposed_nm=result.positions_nm, moved_indices=changed))
                report.update(status='proposed', proposal_file=name, moved=[asdict(i) for i in result.moved],
                              supported=len(result.supported), unsupported=len(result.unsupported), observed_exact=True)
            except GuardRefusal as error:
                report.update(status='refused', error=str(error))
            except Exception as error:
                report.update(status='execution-error', error_type=type(error).__name__, error=str(error))
        cases.append(report)
    check(plan)
    write(output, dict(schema='dvbfixer.private-post-h-static-evidence.v1',
        plan_sha256=sha(plan_path.read_bytes()), cases=cases,
        artifacts={p.name:sha(p.read_bytes()) for p in sorted(root.glob('*.json.gz'))},
        contexts=0, minimizations=0, new_h_sampling=False, production_changed=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('plan', 'run'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--plan', type=Path)
    parser.add_argument('--root', type=Path)
    args = parser.parse_args()
    if args.command == 'plan':
        freeze(args.output.resolve())
    else:
        run(args.plan.resolve(), args.root.resolve(), args.output.resolve())
