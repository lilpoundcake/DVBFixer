"""Independent retained-snapshot replay and post-selection held-out reporting."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import statistics
from collections import Counter
from pathlib import Path

from scripts.research_rosetta_repack import identity, observed_atoms, write_json
from scripts.rosetta_fixed_contract import choose_scored_own, contacts, import_audit


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def json_normal(value):
    return json.loads(json.dumps(value, allow_nan=False))


def rmsd(case: dict, source: bytes, reference: bytes) -> float:
    actual = {identity(a): a['xyz'] for a in observed_atoms(source)}
    truth = {identity(a): a['xyz'] for a in observed_atoms(reference)}
    names = [tuple(k) for k in case['generated']]
    return math.sqrt(sum(math.dist(actual[k], truth[k])**2 for k in names)/len(names))


def source_records_preserved(source: bytes, candidate: bytes, generated: list) -> bool:
    generated = {tuple(k) for k in generated}
    retained = []
    for line in candidate.splitlines(keepends=True):
        if line.startswith(b'ATOM  ') and identity(observed_atoms(line)[0]) in generated:
            continue
        retained.append(line)
    return b''.join(retained) == source


def summarize(plan_root: Path, run_root: Path, output: Path) -> None:
    plan_path, raw_path = plan_root/'plan.json', run_root/'evidence.json'
    plan, raw = json.loads(plan_path.read_bytes()), json.loads(raw_path.read_bytes())
    if raw['plan_sha256'] != sha(plan_path) or (run_root/'executed-plan.json').read_bytes() != plan_path.read_bytes():
        raise ValueError('executed scoring plan changed')
    for name, digest in plan['artifacts'].items():
        if sha(plan_root/name) != digest:
            raise ValueError('scoring plan artifact changed')
    for name, digest in plan['source_hashes'].items():
        if sha(Path(name)) != digest:
            raise ValueError('scoring implementation changed after plan')
    cases, matched, changed = {}, [], []
    validated_snapshots = 0
    statuses, contact_counts = Counter(), Counter()
    topology_limited = []
    for case in plan['cases']:
        name = case['id']
        row = raw['results'][name]
        if row['status'] == 'shared_admission_refused':
            cases[name] = row
            continue
        folder = run_root/name
        own = [c for c in row['cells'] if c['arm'] == 'own']
        selected = choose_scored_own(own)
        if selected != row['selection'] or selected != json.loads((folder/'selection-before-native.json').read_bytes()):
            raise ValueError('selection differs from independent feasible-pool replay')
        reference = (plan_root/case['reference_reporting_only']).read_bytes()
        cells, extras, virtual_coordinates, context_hydrogens = [], [], [], []
        for cell in row['cells']:
            if not source_records_preserved((plan_root/name/'source.pdb').read_bytes(),
                                            (plan_root/cell['pdb']).read_bytes(), case['generated']):
                raise ValueError('original source records changed')
            cell_folder = folder/cell['id']
            if json.loads((cell_folder/'score.json').read_bytes()) != cell:
                raise ValueError('cell score record changed')
            snapshot_path = cell_folder/'pose-snapshot.json.gz'
            if snapshot_path.exists():
                if sha(snapshot_path) != cell['snapshot_sha256']:
                    raise ValueError('scorer snapshot changed')
                snapshot = json.loads(gzip.decompress(snapshot_path.read_bytes()))
                validated_snapshots += 1
                physical = snapshot['heavy_atoms']
                expected = observed_atoms((plan_root/cell['pdb']).read_bytes())
                virtual = snapshot.get('virtual_atoms', [])
                if any(not a['virtual_by_rosetta_type'] for a in virtual):
                    raise ValueError('virtual atom lacks type declaration')
                if {identity(a) for a in virtual} & {identity(a) for a in physical}:
                    raise ValueError('physical and virtual inventories overlap')
                # Re-derive the bounded extra-heavy allowance, not a name-only bypass.
                declared_outgoing = {tuple(a if a[-1] == 'C' else b)[:4] for a, b in case['expected_external']
                                     if {a[-1], b[-1]} == {'C', 'N'}}
                terminal = {(r['chain'], r['resid'], r['icode']) for r in snapshot['residues']
                            if 'CtermProteinFull' in r['name']}
                allowed = [identity(a) for a in physical if a['atom'] == 'OXT'
                           and identity(a)[:3] in terminal and identity(a)[:4] not in declared_outgoing
                           and identity(a) not in {identity(e) for e in expected}]
                audit = import_audit(expected, physical, case['expected_external'], snapshot['external'],
                                     case['target_internal'], snapshot['target_internal'], allowed)
                if json_normal(audit) != cell['import_audit']:
                    raise ValueError('independent physical/topology import audit differs')
                if cell['score_status'] == 'scored' and not audit['accepted']:
                    raise ValueError('scored cell failed independent import audit')
                extras.append(sorted((identity(a), a['xyz']) for a in physical if identity(a) in {tuple(k) for k in allowed}))
                virtual_coordinates.append(sorted((identity(a), a['xyz']) for a in virtual))
                context_hydrogens.append(sorted(((h['chain'], h['resid'], h['icode'], h['atom']), h['xyz'])
                                                for h in snapshot['hydrogens']
                                                if (h['chain'], h['resid'], h['icode']) != tuple(case['target_key'][:3])))
                if cell['score_status'] == 'scored':
                    if not math.isclose(sum(t['raw']*t['weight'] for t in cell['score_terms'].values()), cell['total_score'],
                                        rel_tol=1e-10, abs_tol=1e-8):
                        raise ValueError('score-term sum differs from total')
            candidate_rows = observed_atoms((plan_root/cell['pdb']).read_bytes())
            for atom in candidate_rows:
                atom['element'] = atom['atom'][0]  # this protocol admits canonical protein heavy names only
            findings = contacts(candidate_rows, case['generated'], case['target_internal']+case['expected_external'], plan['radii'])
            old_clashes = [f for f in cell['common_audit']['rounded']['findings'] if f.startswith('generated_clash:')]
            if len(findings) != len(old_clashes):
                raise ValueError('topology description changed old rounded clash denominator')
            statuses[cell['score_status']] += 1
            contact_counts.update(f['topology'] for f in findings)
            cells.append(dict(id=cell['id'], arm=cell['arm'], index=cell['index'], score_status=cell['score_status'],
                              old_source_accepted=cell['old_source_accepted'], common_accepted=cell['common_accepted'],
                              total_score=cell['total_score'],
                              reporting_only_rounded_exact_name_rmsd_angstrom=rmsd(case, (plan_root/cell['pdb']).read_bytes(), reference),
                              rounded_old_clashes_with_topology=findings,
                              score_record_sha256=sha(cell_folder/'score.json'),
                              snapshot_sha256=sha(snapshot_path) if snapshot_path.exists() else None,
                              import_audit=cell.get('import_audit'),
                              hydrogen_count=cell.get('hydrogen_count'), virtual_atom_count=cell.get('virtual_atom_count')))
        same_extras = not extras or all(e == extras[0] for e in extras)
        same_virtual = not virtual_coordinates or all(v == virtual_coordinates[0] for v in virtual_coordinates)
        same_context_h = not context_hydrogens or all(h == context_hydrogens[0] for h in context_hydrogens)
        if not same_extras or not same_virtual or not same_context_h:
            topology_limited.append(name)
        previous = case['old_selected_index']
        new = selected['selected_index']
        reporting = {}
        if previous is not None and new is not None:
            old_cell = next(c for c in cells if c['arm'] == 'own' and c['index'] == previous)
            new_cell = next(c for c in cells if c['arm'] == 'own' and c['index'] == new)
            reporting = dict(old_rounded_rmsd=old_cell['reporting_only_rounded_exact_name_rmsd_angstrom'],
                             scored_rounded_rmsd=new_cell['reporting_only_rounded_exact_name_rmsd_angstrom'])
            if same_extras and same_virtual and same_context_h:
                matched.append(dict(id=name, mode=case['mode'], **reporting))
                if new != previous:
                    changed.append(dict(id=name, previous=previous, selected=new, **reporting))
        cases[name] = dict(status='evaluated', old_selected_index=previous, selection=selected,
                           old_gate_feasible_own_count=sum(c['old_source_accepted'] and c['common_accepted'] for c in own),
                           extra_heavy_coordinates_identical=same_extras, virtual_coordinates_identical=same_virtual,
                           scorer_owned_context_hydrogens_identical=same_context_h,
                           reporting=reporting, cells=cells)
    aggregate = dict(planned_cases=len(cases), shared_refusals=sum(c['status'] == 'shared_admission_refused' for c in cases.values()),
                     independently_validated_snapshots=validated_snapshots,
                     cell_statuses=dict(statuses), retained_clash_topologies=dict(contact_counts),
                     own_selection_statuses=dict(Counter(c['selection']['status'] for c in cases.values() if c['status'] == 'evaluated')),
                     selections_changed=changed, context_coordinate_comparison_blockers=topology_limited)
    for mode in ('paired', 'anchored'):
        rows = [r for r in matched if r['mode'] == mode]
        aggregate[mode] = dict(matched_count=len(rows),
                               old_median_rounded_rmsd=statistics.median(r['old_rounded_rmsd'] for r in rows) if rows else None,
                               scored_median_rounded_rmsd=statistics.median(r['scored_rounded_rmsd'] for r in rows) if rows else None)
    write_json(output, dict(schema='dvbfixer.rosetta-fixed-scoring-summary.v2',
                            plan_path=str(plan_path), plan_sha256=sha(plan_path), raw_path=str(raw_path), raw_sha256=sha(raw_path),
                            summarizer_sha256=sha(Path(__file__)), cases=cases, aggregate=aggregate,
                            environment=json.loads((run_root/'environment.json').read_bytes()),
                            frozen_source_hashes=plan['source_hashes'], frozen_artifact_hashes=plan['artifacts'],
                            hydrogen_protocol=plan['hydrogen_protocol'], virtual_atom_policy=plan['virtual_atom_policy'],
                            protocol_history=[
                                dict(path='tmp/pdbfixer-research/rosetta-fixed-plan-01/execution-error.json',
                                     status='planning_serialization_blocker', scoring_executed=False),
                                dict(path='tmp/pdbfixer-research/rosetta-fixed-run-02.log',
                                     status='missing_pythonpath_invocation_blocker', scoring_executed=False),
                                dict(path='tmp/pdbfixer-research/rosetta-fixed-run-03/evidence.json',
                                     sha256=(sha(Path('tmp/pdbfixer-research/rosetta-fixed-run-03/evidence.json'))
                                             if Path('tmp/pdbfixer-research/rosetta-fixed-run-03/evidence.json').is_file() else None),
                                     status='all_128_import_inventory_blocked_before_scoring',
                                     correction='version-2 physical/typed-virtual separation; no physical atom exception by name')],
                            scientific_qualification=False, prospective_accuracy_benchmark=False,
                            score_scope='full ref2015 after fixed-heavy import and fresh unoptimized H; no calibration or new steric gate',
                            topology_scope=plan['topology_scope']))
    print(json.dumps(aggregate, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    summarize(args.plan, args.run, args.output)


if __name__ == '__main__':
    main()
