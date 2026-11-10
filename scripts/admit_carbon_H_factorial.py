"""Freeze capture transport, source refusals and CH2/CH3 geometry admission.

No Context or arm minimization. Coordinates are read only for admission and the
fixed constructors; no energies, phase searches or outcome selection.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from scripts.audit_graph_methyl_whole import global_shape
from scripts.carbon_H_capture_adapter import adapt
from scripts.graph_methyl_geometry import canonical as methyl_geometry
from scripts.graph_methyl_geometry import discover as methyl_groups
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify
from scripts.research_methylene_geometry import canonical as methylene_geometry
from scripts.research_methylene_geometry import discover as methylene_groups


def admit(enrollment_path, output):
    import openmm

    raw = enrollment_path.read_bytes()
    enrollment = json.loads(raw)
    verify(dict(inputs=enrollment['inputs'], source_hashes={}))
    if (enrollment['summary']['contexts'] != 40 or len(enrollment['cases']) != 40
            or enrollment['arms_evaluated'] or enrollment['coordinates_read']):
        raise ValueError('unchanged metadata-only complete enrollment required')
    inputs = {**enrollment['inputs'], str(enrollment_path): sha(raw)}
    parents = []
    for name in enrollment['inputs']:
        if not name.endswith('.json'):
            continue
        parent = json.loads(Path(name).read_bytes())
        if parent.get('schema') not in {'dvbfixer.tetra-methyl-start-plan.v1', 'dvbfixer.methyl-diversity-capture-plan.v2'}:
            continue
        parents.append(parent['schema'])
        verify(parent)
        if parent['openmm_version'] != openmm.__version__:
            raise ValueError('capture runtime differs')
        inputs.update(parent['inputs'])
    if sorted(parents) != ['dvbfixer.methyl-diversity-capture-plan.v2', 'dvbfixer.tetra-methyl-start-plan.v1']:
        raise ValueError('both original capture provenance plans required')
    # Bind every retained capture artifact, including refusals, to its original
    # manifest before loading. Legacy source artifacts are bound by its parent.
    capture_root = Path(enrollment['cases'][24]['system']).parents[2]
    capture_run = json.loads((capture_root/'evidence.json').read_bytes())
    for case in capture_run['cases']:
        folder = capture_root/case['id']/f"seed-{case['seed']}"
        for name, digest in case['artifacts'].items():
            if sha((folder/name).read_bytes()) != digest:
                raise ValueError('native capture artifact changed')
            inputs[str(folder/name)] = digest
    cases = []
    for case in enrollment['cases']:
        row = dict(case)
        if case['status'] == 'refused_source':
            declaration = next(c for c in capture_run['cases'] if (c['id'], c['seed']) == (case['id'], case['seed']))
            row.update(status='refused_source', inherited_capture_record=declaration)
            cases.append(row)
            continue
        # Resource/transport errors are blockers, never geometry exclusions.
        data = adapt(case)
        for key in ('atoms', 'pre_H', 'topology', 'baseline_path', 'system'):
            inputs[case[key]] = sha(Path(case[key]).read_bytes())
        source = Path(case['system']).parent
        for name in (('phase-1-before.json.gz', 'h-internal-after.json.gz')
                if case['context_kind'] == 'geometry-reconstruction-stress' else ('H-ownership.json',)):
            inputs[str(source/name)] = sha((source/name).read_bytes())
        row.update(native_graph_exact=True, native_observations_exact=True,
            particle_count=len(data['rows']), newly_mobile_H=len(data['new_h']),
            initial_heavy_shape=global_shape([dict(r, xyz=[10*x for x in p])
                for r, p in zip(data['rows'], data['initial_nm'], strict=True)]),
            replay_labels=['unchanged is an exact computational replay, not a new independent sample'])
        if case['context_kind'] == 'prepared-protein-view':
            row['replay_labels'].append('all_methyl_zero repeats the prior source-diverse graph-zero H arm')
            if case['id'] == '1TIM':
                row['replay_labels'].append('all_methylene/combined_zero intersect prior all-CH2 controls from unchanged/graph-zero sources')
        try:
            methyl = methyl_groups(data['rows'], data['bonds'], set(data['new_h']), set(data['mobile']))
            methylene = methylene_groups(data['rows'], data['bonds'], set(data['new_h']), set(data['mobile']))
            row.update(methyl=methyl, methylene=methylene)
            for group in methyl['groups']:
                methyl_geometry(data['initial_nm'], group, 0)
            for group in methylene['groups']:
                methylene_geometry(data['initial_nm'], group)
            methyl_H = {i for g in methyl['groups'] for i in g['hydrogens']}
            methylene_H = {i for g in methylene['groups'] for i in g['hydrogens']}
            if methyl_H & methylene_H:
                raise ValueError('overlapping factorial H ownership')
            row['status'] = 'admitted'
        except ValueError as error:
            row.update(status='refused_geometry', refusal=str(error))
        cases.append(row)
    sources = {name: sha(Path(name).read_bytes()) for name in (
        'scripts/admit_carbon_H_factorial.py', 'scripts/carbon_H_capture_adapter.py',
        'scripts/graph_methyl_geometry.py', 'scripts/research_methylene_geometry.py',
        'scripts/audit_graph_methyl_whole.py')}
    verify(dict(inputs=inputs, source_hashes=sources))
    summary = dict(contexts=len(cases), statuses=dict(Counter(c['status'] for c in cases)),
        potential_arm_cells=4*sum(c['status'] == 'admitted' for c in cases))
    write(output, dict(schema='dvbfixer.carbon-H-factorial-admission.v1', cases=cases,
        inputs=inputs, source_hashes=sources, summary=summary, enrollment_sha256=sha(raw),
        openmm_version=openmm.__version__, policy=enrollment['policy'],
        arms_evaluated=False, production_changed=False, scientific_qualification=False))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--enrollment', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    admit(args.enrollment, args.output)
