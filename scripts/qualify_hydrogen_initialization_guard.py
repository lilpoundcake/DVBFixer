"""Frozen saved-initializer static cohort; no Context, optimizer or new H draws."""
from __future__ import annotations

import argparse
import math
import time
from dataclasses import asdict
from pathlib import Path

from scripts.audit_thr2_native_capture import read, verify_files
from scripts.hydrogen_initialization_guard import Atom, GuardRefusal, Input, propose, validate
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify
from scripts.research_legacy_start_controls import compressed
from scripts.research_thr2_hb_only import transport

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.ffutils.variants import VARIANT_TO_PARENT

ROOT = Path(__file__).resolve().parents[1]


def atoms_from_rows(rows):
    """Use explicit topology occurrence indices, or require unique legacy keys."""
    metadata = all('chain_index' in r and 'residue_index' in r for r in rows)
    if not metadata and any('chain_index' in r or 'residue_index' in r for r in rows):
        raise GuardRefusal('partial occurrence metadata')
    if not metadata and len({tuple(r['key']) for r in rows}) != len(rows):
        raise GuardRefusal('ambiguous legacy atom keys')
    occurrences, counts, atoms = {}, {}, []
    for row in rows:
        chain, resid, icode, resname, name = row['key']
        base = chain, resid, icode
        group = (row['chain_index'], row['residue_index']) if metadata else base
        if group not in occurrences:
            counts[base] = counts.get(base, 0)+1
            occurrences[group] = base, counts[base]
        if occurrences[group][0] != base:
            raise GuardRefusal('inconsistent occurrence metadata')
        component = ComponentInstanceRef(1, chain, resid, icode, occurrence=occurrences[group][1])
        atoms.append(Atom(ExactAtomRef(component, name), resname, row['element']))
    return tuple(atoms)


def admit(rows, before, positions, bonds, system):
    """Cross-check positive mass with newly absent H; observations must be pinned."""
    from openmm.unit import dalton

    atoms, observed = atoms_from_rows(rows), atoms_from_rows(before['atoms'])
    old = {a.identity:a for a in observed}
    if len(old) != len(observed) or len(set(a.identity for a in atoms)) != len(atoms):
        raise GuardRefusal('duplicate observation/initializer identity')
    if not set(old) <= {a.identity for a in atoms}:
        raise GuardRefusal('observed identity disappeared')
    if (system.getNumParticles() != len(atoms) or len(positions) != len(atoms)
            or system.getNumConstraints() or any(system.isVirtualSite(i) for i in range(len(atoms)))):
        raise GuardRefusal('unsupported native System inventory/constraints/virtual sites')
    new = set()
    for i, atom in enumerate(atoms):
        mass = system.getParticleMass(i).value_in_unit(dalton)
        if not math.isfinite(mass) or mass < 0:
            raise GuardRefusal('invalid mass')
        absent = atom.identity not in old
        if (mass > 0) != absent or (absent and atom.element != 'H'):
            raise GuardRefusal('native mobility differs from new-H ownership')
        if absent:
            new.add(atom.identity)
        else:
            previous = old[atom.identity]
            if (previous.element != atom.element or
                    VARIANT_TO_PARENT.get(previous.residue_name, previous.residue_name)
                    != VARIANT_TO_PARENT.get(atom.residue_name, atom.residue_name)):
                raise GuardRefusal('observed atom chemistry changed')
    # XML State positions are authoritative; table xyz is only redundant readback.
    if any(math.dist(p, [v/10 for v in row['xyz']]) > 1e-12 for p, row in zip(positions, rows, strict=True)):
        raise GuardRefusal('saved row/State coordinates differ')
    by_identity = {a.identity:i for i, a in enumerate(atoms)}
    original = before.get('positions_nm', [[v/10 for v in row['xyz']] for row in before['atoms']])
    transport(original, [positions[by_identity[a.identity]] for a in observed], 'captured')
    value = Input(atoms, tuple(tuple(p) for p in positions), tuple(tuple(b) for b in bonds), frozenset(new))
    validate(value)
    return value


def freeze(output):
    inputs, cases = {}, []
    old_root = ROOT/'tmp/pdbfixer-research/thr2-native-capture-run-01'
    real_root = ROOT/'tmp/hydrogen-postflight-contract-01/glycam-fallback-01/native-run-01'
    for root, kind in ((old_root, 'original-thr2'), (real_root, 'fixed-real-cohort')):
        evidence = root/'evidence.json'
        inputs[str(evidence)] = sha(evidence.read_bytes())
        execution = read(evidence)
        for row in execution['native' if kind == 'original-thr2' else 'cases']:
            folder = root/row['label'] if kind == 'original-thr2' else root/f"case-{row['index']:02d}"
            verify_files(folder, row)
            inputs.update({str(folder/name):digest for name, digest in row['artifacts'].items()})
            capture = read(folder/'capture.json')
            calls = [capture['calls'][0]] if kind == 'original-thr2' else capture['native']
            h_calls = []
            for index, h in enumerate(capture['h_calls'], 1):
                h_calls.append(dict(number=index, returned=h['returned'],
                    forcefield_provided=h.get('forcefield_provided', 'not-recorded'),
                    forcefield_xmls=h.get('forcefield_xmls', 'not-recorded'),
                    error_type=h.get('error_type'), error=h.get('error')))
            slots = [dict(prefix=c['before']['prefix'], h_call=c.get('h_call', 1),
                          bonds='call-1-topology.json' if kind == 'original-thr2' else f"native-{c['number']}-bonds.json") for c in calls]
            cases.append(dict(kind=kind, label=row['label'], folder=str(folder),
                requested_ff=row.get('requested_ff', 'not-recorded'), h_calls=h_calls, slots=slots))
    if len(cases) != 15:
        raise ValueError('Complete three-original/twelve-real enrollment required')
    paths = [Path(__file__).resolve(), ROOT/'scripts/hydrogen_initialization_guard.py',
             ROOT/'scripts/audit_hydrogen_initialization_guard.py',
             ROOT/'tests/test_hydrogen_initialization_guard.py',
             ROOT/'src/dvbfixer/domain/structure_identity.py', ROOT/'src/dvbfixer/ffutils/variants.py',
             ROOT/'src/dvbfixer/model/diffusion/protein_atoms.py', ROOT/'scripts/research_thr2_hb_only.py',
             ROOT/'scripts/audit_thr2_native_capture.py', ROOT/'scripts/research_dunbrack_prior.py',
             ROOT/'scripts/research_ha_minimizer_path.py', ROOT/'scripts/research_legacy_start_controls.py']
    contract = ROOT/'docs/plans/hydrogen-initialization-guard.md'
    inputs[str(contract)] = sha(contract.read_bytes())
    write(output, dict(schema='dvbfixer.private-initializer-static-plan.v1', cases=cases, inputs=inputs,
        source_hashes={str(p):sha(p.read_bytes()) for p in paths},
        policy=dict(maximum_seconds=120, contexts=0, minimizations=0, new_h_sampling=False,
            retention='all admitted/refused/no-native slots and unsuccessful H calls, no retries or source replacement'),
        production_changed=False))


def run(plan_path, output, public_output):
    import openmm
    from openmm.unit import nanometer

    plan = read(plan_path)
    verify(plan)
    output.mkdir(parents=True, exist_ok=False)
    (output/'executed-plan.json').write_bytes(plan_path.read_bytes())
    started, reports = time.monotonic(), []
    for index, case in enumerate(plan['cases']):
        folder = Path(case['folder'])
        report = dict(index=index, **{k:v for k, v in case.items() if k != 'folder'}, results=[])
        for slot_index, slot in enumerate(case['slots']):
            result = dict(slot=slot_index, prefix=slot['prefix'], status='not-executed-budget')
            if time.monotonic()-started <= plan['policy']['maximum_seconds']:
                try:
                    prefix = slot['prefix']
                    rows = read(folder/f'{prefix}-atoms.json.gz')['atoms']
                    before = read(folder/f"h-{slot['h_call']}-before.json.gz")
                    bonds = read(folder/slot['bonds'])['bonds']
                    system = openmm.XmlSerializer.deserialize((folder/f'{prefix}-system.xml').read_text())
                    state = openmm.XmlSerializer.deserialize((folder/f'{prefix}-state.xml').read_text())
                    positions = state.getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
                    value = admit(rows, before, positions, bonds, system)
                    result.update(admitted=True, atoms=len(value.atoms), new_h=len(value.new_hydrogens))
                    proposed = propose(value)
                    changed = [i for i, (a, b) in enumerate(zip(value.positions_nm, proposed.positions_nm, strict=True)) if a != b]
                    if {value.atoms[i].identity for i in changed} != set(proposed.moved):
                        raise AssertionError('moved-mask mismatch')
                    if not set(proposed.moved) <= value.new_hydrogens:
                        raise AssertionError('observed atom moved')
                    compressed(output/f'case-{index:02d}-slot-{slot_index}.json.gz',
                        dict(atoms=[asdict(a) for a in value.atoms], initial_nm=value.positions_nm,
                             bonds=value.bonds, new_indices=[i for i, a in enumerate(value.atoms) if a.identity in value.new_hydrogens],
                             proposed_nm=proposed.positions_nm, moved_indices=changed))
                    result.update(status='proposed', moved=[asdict(i) for i in proposed.moved],
                        supported=len(proposed.supported), unsupported=len(proposed.unsupported), observed_exact=True,
                        proposal_file=f'case-{index:02d}-slot-{slot_index}.json.gz')
                except GuardRefusal as error:
                    result.update(status='refused', error=str(error))
                except Exception as error:
                    result.update(status='execution-error', error_type=type(error).__name__, error=str(error))
            report['results'].append(result)
        reports.append(report)
    verify(plan)
    write(public_output, dict(schema='dvbfixer.private-initializer-static-evidence.v1',
        plan_sha256=sha(plan_path.read_bytes()), cases=reports,
        artifacts={p.name:sha(p.read_bytes()) for p in sorted(output.glob('*.json.gz'))},
        contexts=0, minimizations=0, new_h_sampling=False, production_changed=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('plan', 'run'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--plan', type=Path)
    parser.add_argument('--public-output', type=Path)
    args = parser.parse_args()
    if args.command == 'plan':
        freeze(args.output.resolve())
    else:
        run(args.plan.resolve(), args.output.resolve(), args.public_output.resolve())
