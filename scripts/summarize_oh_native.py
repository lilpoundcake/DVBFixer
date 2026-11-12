"""Scalar coverage and retained upstream-blocker inventory; never sample H."""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from scripts.audit_oh_native import read
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify


def duplicates(saved):
    atoms = saved['atoms']
    counts = Counter(tuple(row['key']) for row in atoms)
    neighbors = {i:[] for i in range(len(atoms))}
    for a, b in saved['bonds']:
        neighbors[a].append(b)
        neighbors[b].append(a)
    return [dict(key=list(key), count=count,
        occurrences=[dict(index=i, bonded_neighbors=[atoms[j]['key'] for j in neighbors[i]])
            for i, row in enumerate(atoms) if tuple(row['key']) == key])
        for key, count in sorted(counts.items()) if count > 1]


def summarize(plan_path, run, audit_path, output):
    plan = json.loads(plan_path.read_bytes())
    verify(plan)
    audit = json.loads(audit_path.read_bytes())
    execution = json.loads((run/'evidence.json').read_bytes())
    if audit['plan_sha256'] != sha(plan_path.read_bytes()) or audit['evidence_sha256'] != sha((run/'evidence.json').read_bytes()):
        raise ValueError('audit/execution binding differs')
    coverage, case_rows, findings = Counter(), [], []
    native_ff, fallback = [], []
    for row in audit['cases']:
        folder = run/f"case-{row['index']:02d}"
        captured = json.loads((folder/'capture.json').read_bytes())
        for name, digest in execution['cases'][row['index']]['artifacts'].items():
            if sha((folder/name).read_bytes()) != digest:
                raise ValueError('retained source artifact changed')
        routes = [h for h in captured['h_calls'] if h['returned']]
        if len(routes) != 1:
            raise ValueError('unexpected successful H-call count')
        h = routes[0]
        route = 'native-forcefield' if h['forcefield_provided'] else 'forcefield-free'
        (native_ff if h['forcefield_provided'] else fallback).append(row['index'])
        xmls = h['forcefield_xmls'] or []
        family = ('forcefield-free' if not h['forcefield_provided'] else
                  'amber-glycam' if any('GLYCAM' in name for name in xmls) else
                  'charmm' if any('charmm' in name for name in xmls) else 'amber')
        for repair in row['repairs']:
            for key, count in repair['candidate']['measurement']['oh_counts'].items():
                coverage[family+'/'+key] += count
        failure_details = [dict(number=h['number'], error_type=h.get('error_type'),
                                error=h.get('error')) for h in captured['h_calls'] if not h['returned']]
        saved_duplicates = {}
        for path in sorted(folder.glob('h-*-*.json.gz')):
            duplicate = duplicates(read(path))
            if duplicate:
                saved_duplicates[path.name] = duplicate
        candidate_changes = [r['candidate']['changed_atom_keys'] for r in row['repairs']]
        baseline_remaining = sum(len(r['baseline']['measurement']['violations']) for r in row['repairs'])
        before_violation_count = sum(len(r['before']['violations']) for r in row['repairs'])
        case_rows.append(dict(index=row['index'], label=row['label'], requested_ff=row['requested_ff'],
            status=row['status'], successful_h_route=route, actual_family=family, actual_xmls=xmls,
            failed_h_attempts=failure_details, duplicate_h_snapshot_inventory=saved_duplicates,
            candidate_changed_atoms=sum(len(keys) for keys in candidate_changes),
            baseline_remaining_oh_violations=baseline_remaining,
            candidate_remaining_oh_violations=sum(len(r['candidate']['measurement']['violations']) for r in row['repairs']),
            input_oh_violations=before_violation_count,
            partial_output_not_success=row['output_present'] and row['status'] != 'completed',
            dat_valid=row['dat']['valid']))
        for repair in row['repairs']:
            for finding in repair['before']['violations']:
                findings.append(dict(case=row['index'], actual_family=family,
                    corrected=repair['candidate']['status'] == 'completed' and repair['candidate']['measurement']['passed'],
                    **finding))
    rows = [r for case in audit['cases'] for r in case['repairs']]
    summary = dict(audit['summary'], native_forcefield_h_return_cases=native_ff,
        forcefield_free_h_return_cases=fallback,
        added_contract_violations_repaired=sum(r['corrected'] for r in findings),
        native_forcefield_added_contract_violations_repaired=sum(r['corrected'] and r['actual_family'] != 'forcefield-free' for r in findings),
        all_heavy_valid_h_bonds_identity_preserved=all(all(r['candidate'][key] for key in
            ('heavy_preserved', 'valid_h_preserved', 'bonds_preserved', 'identity_preserved')) for r in rows),
        new_refusals=sum(r['new_refusal'] for r in rows),
        short_oh_warnings=sum(len(r['candidate']['measurement']['short_oh_warnings']) for r in rows),
        cross_residue_warnings=sum(len(r['candidate']['measurement']['cross_residue_warnings']) for r in rows),
        duplicate_atom_identity_cases=[r['index'] for r in case_rows if r['duplicate_h_snapshot_inventory']],
        independent_ligand_oh_postflight_coverage=False, native_charmm_glycan_qualification=False)
    write(output, dict(schema='dvbfixer.oh-native-summary.v1', summary=summary,
        plan_sha256=sha(plan_path.read_bytes()), execution_sha256=sha((run/'evidence.json').read_bytes()),
        audit_sha256=sha(audit_path.read_bytes()), script_sha256=sha(Path(__file__).read_bytes()),
        cases=case_rows, oh_coverage_counts=dict(sorted(coverage.items())), corrected_findings=findings,
        limits='Counts describe saved H states, not independent structures. Duplicate-identity partial outputs can look OH-clean after PDB parser atom merging; failed case remains invalid. No native reruns or main merge.'))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('plan', 'run', 'audit', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    summarize(args.plan.resolve(), args.run.resolve(), args.audit.resolve(), args.output.resolve())
