"""Two frozen raw returns: ordinary repair then strict guard, offline only."""
from __future__ import annotations

import argparse
import math
import random
import time
from dataclasses import asdict, replace
from pathlib import Path

from scripts.audit_post_h_1vcu import target_geometry
from scripts.audit_thr2_native_capture import read
from scripts.hydrogen_initialization_guard import GuardRefusal, propose
from scripts.post_hydrogen_public_api import ownership
from scripts.qualify_glycam_fallback import restore
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify
from scripts.research_legacy_start_controls import compressed
from scripts.research_post_h_1vcu import check as check_parent
from scripts.research_thr2_hb_only import transport

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.ffutils.geometry import repair_misplaced_hydrogens
from dvbfixer.ffutils.variants import VARIANT_TO_PARENT
from dvbfixer.model.diffusion.protein_atoms import CANONICAL_HEAVY_ATOMS

ROOT = Path(__file__).resolve().parents[1]


def input_from_table(table):
    from scripts.hydrogen_initialization_guard import Atom, Input, validate

    atoms = tuple(Atom(ExactAtomRef(ComponentInstanceRef(**a['identity']['component']), a['identity']['name']),
                       a['residue_name'], a['element']) for a in table['atoms'])
    value = Input(atoms, tuple(tuple(p) for p in table['positions_nm']),
                  tuple(tuple(b) for b in table['bonds']), frozenset())
    validate(value)
    return value


def reconstruct(table, metadata):
    rows, residues, chains = [], {}, {}
    for atom, point in zip(table['atoms'], table['positions_nm'], strict=True):
        c = atom['identity']['component']
        if c['model'] != 1 or c['occurrence'] != 1 or c['alternate_location']:
            raise ValueError('offline 1VCU transport requires unambiguous one-model components')
        chain = c['chain_id']
        key = chain, c['sequence_number'], c['insertion_code']
        if chain not in chains:
            chains[chain] = len(chains)
        if key not in residues:
            residues[key] = len(residues)
        rows.append(dict(key=[*key, atom['residue_name'], atom['identity']['name']], element=atom['element'],
                         chain_index=chains[chain], residue_index=residues[key], xyz=[10*x for x in point]))
    kinds = {tuple(sorted((a, b))):(kind, order) for a, b, kind, order in metadata}
    saved = dict(atoms=rows, positions_nm=table['positions_nm'], bonds=table['bonds'],
        bond_metadata=[dict(type=kinds[tuple(sorted((a, b)))][0], order=kinds[tuple(sorted((a, b)))][1])
                       for a, b in table['bonds']])
    return restore(saved)


def distance_findings(value):
    adjacent = [[] for _ in value.atoms]
    for a, b in value.bonds:
        adjacent[a].append(b)
        adjacent[b].append(a)
    findings = []
    for h, atom in enumerate(value.atoms):
        if atom.identity not in value.new_hydrogens:
            continue
        parent, = adjacent[h]
        center = value.atoms[parent]
        heavy = [i for i in adjacent[parent] if value.atoms[i].element not in ('H', 'D')]
        hs = [i for i in adjacent[parent] if value.atoms[i].element in ('H', 'D')]
        inventory = CANONICAL_HEAVY_ATOMS.get(VARIANT_TO_PARENT.get(center.residue_name, center.residue_name), ())
        if (center.element != 'C' or len(heavy) != 3 or hs != [h] or center.identity.name not in inventory
                or any(value.atoms[i].identity.component != center.identity.component for i in [h]+heavy)
                or any(value.atoms[i].identity.name not in inventory for i in heavy)):
            continue
        radius = math.dist(value.positions_nm[h], value.positions_nm[parent])
        distances = [math.dist(value.positions_nm[i], value.positions_nm[parent]) for i in heavy]
        if not 0.05 <= radius <= 0.15 or any(not 0.10 <= r <= 0.20 for r in distances):
            findings.append(dict(atom=asdict(atom.identity), residue=atom.residue_name,
                                 radius_nm=radius, heavy_bond_lengths_nm=distances))
    return findings


def freeze(output):
    parent_path = ROOT/'docs/research/post-hydrogen-1vcu-plan.json'
    parent = read(parent_path)
    check_parent(parent)
    root = ROOT/'tmp/pdbfixer-research/post-hydrogen-1vcu-public-api-01'
    inputs = {**parent['inputs'], **parent['source_hashes']}
    cases = []
    for i, label in enumerate(('baseline', 'guarded')):
        journal = root/f'case-{i:02d}/h-1-journal.json.gz'
        data = read(journal)
        if data['stage'] != 'guard-refusal' or data['raw'] is None:
            raise ValueError('both actual retained guard-refusal raw returns are required')
        cases.append(dict(arm=label, journal=str(journal),
            baseline_whole=str(root/'case-00/whole-1-before-atoms.json.gz') if i == 0 else None))
        inputs[str(journal)] = sha(journal.read_bytes())
    for p in (parent_path, root/'evidence.json', root/'case-00/whole-1-before-atoms.json.gz',
              ROOT/'docs/plans/post-hydrogen-after-repair.md', ROOT/'src/dvbfixer/ffutils/geometry.py'):
        inputs[str(p)] = sha(p.read_bytes())
    write(output, dict(schema='dvbfixer.private-post-h-shadow-repair-plan.v1', parent_plan=str(parent_path),
        inputs=inputs, source_hashes={str(p):sha(p.read_bytes()) for p in
            (Path(__file__).resolve(), ROOT/'tests/test_post_h_after_repair.py')}, cases=cases,
        policy=dict(states=2, maximum_seconds=120, contexts=0, optimizations=0, native_h_sampling=0,
                    core_rule_unchanged=True, native_gate_requires_both_admitted=True)))


def run(plan_path, root, output):
    import numpy as np
    from openmm.unit import nanometer

    plan = read(plan_path)
    verify(plan)
    check_parent(read(Path(plan['parent_plan'])))
    root.mkdir(parents=True, exist_ok=False)
    (root/'executed-plan.json').write_bytes(plan_path.read_bytes())
    started, reports = time.monotonic(), []
    for index, case in enumerate(plan['cases']):
        result = dict(index=index, arm=case['arm'], status='not-executed-budget', underlying_api_completed=True,
                      caller_adapter_committed=case['arm'] == 'baseline')
        if time.monotonic()-started <= plan['policy']['maximum_seconds']:
            payload = read(Path(case['journal']))
            before, raw = input_from_table(payload['before']), input_from_table(payload['raw'])
            value, _ = ownership(before, raw, frozenset(), payload['reference_requested'])
            result['raw_distance_findings'] = distance_findings(value)
            result['raw_target'] = target_geometry(payload['raw']['atoms'], raw.positions_nm)
            topology, positions = reconstruct(payload['raw'], payload['raw_bond_metadata'])
            py_rng, np_rng = random.getstate(), np.random.get_state()
            repaired = None
            try:
                count = repair_misplaced_hydrogens(topology, positions)
                repaired = replace(value, positions_nm=tuple(tuple(float(x) for x in p) for p in positions.value_in_unit(nanometer)))
                changes = [i for i, (a, b) in enumerate(zip(value.positions_nm, repaired.positions_nm, strict=True)) if a != b]
                result.update(ordinary_repair_returned=True, ordinary_repairs=count,
                    ordinary_changed_atoms=[asdict(value.atoms[i].identity) for i in changes],
                    heavy_preserved=all(value.atoms[i].element in ('H', 'D') for i in changes),
                    post_repair_distance_findings=distance_findings(repaired),
                    post_repair_target=target_geometry(payload['raw']['atoms'], repaired.positions_nm))
                if case['baseline_whole']:
                    whole = read(Path(case['baseline_whole']))
                    from scripts.qualify_hydrogen_initialization_guard import atoms_from_rows
                    ids = [a.identity for a in atoms_from_rows(whole['atoms'])]
                    if ids != [a.identity for a in repaired.atoms] or whole['bonds'] != [list(b) for b in repaired.bonds]:
                        raise ValueError('shadow/native whole identity or graph differs')
                    result['baseline_handoff_maximum_rounding_nm'] = transport(
                        [list(p) for p in repaired.positions_nm], whole['positions_nm'], 'captured')
                guard = propose(repaired)
                result.update(status='proposed', moved=[asdict(a) for a in guard.moved],
                    proposed_target=target_geometry(payload['raw']['atoms'], guard.positions_nm))
                compressed(root/f'case-{index:02d}-proposed.json.gz', dict(atoms=payload['raw']['atoms'],
                    bonds=value.bonds, raw_nm=value.positions_nm, repaired_nm=repaired.positions_nm,
                    proposed_nm=guard.positions_nm, new_hydrogens=[asdict(a) for a in sorted(value.new_hydrogens)]))
            except Exception as error:
                result.update(status='refused' if isinstance(error, GuardRefusal) else 'repair-or-handoff-error',
                              error_type=type(error).__name__, error=str(error))
            finally:
                current_np = np.random.get_state()
                result['python_rng_unchanged'] = random.getstate() == py_rng
                result['numpy_rng_unchanged'] = (np_rng[0] == current_np[0] and np.array_equal(np_rng[1], current_np[1])
                                                and np_rng[2:] == current_np[2:])
                compressed(root/f'case-{index:02d}-repair.json.gz', dict(raw=payload['raw'],
                    repaired_nm=repaired.positions_nm if repaired else None, result=result))
        reports.append(result)
    verify(plan)
    gate = all(r['status'] == 'proposed' and r.get('heavy_preserved') and r.get('python_rng_unchanged')
               and r.get('numpy_rng_unchanged') for r in reports)
    write(output, dict(schema='dvbfixer.private-post-h-shadow-repair-evidence.v1',
        plan_sha256=sha(plan_path.read_bytes()), cases=reports, native_gate_both_admitted=gate,
        artifacts={p.name:sha(p.read_bytes()) for p in sorted(root.glob('*.json.gz'))},
        new_contexts=0, new_minimizations=0, native_h_sampling=0,
        limits='Offline new timing control of both refused raw returns, not a retry or a pass for the failed public-API pair.'))


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
