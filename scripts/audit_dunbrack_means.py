"""Independent saved mean-transport/geometry/import audit before native reporting.

Tracked publications exclude per-case library means/weights/sigmas/table samples.
No generation, minimization, packing, library query or score evaluation here.
"""
from __future__ import annotations

import argparse
import json
import statistics
from dataclasses import asdict
from pathlib import Path

from scripts.audit_ncac_expansion import accuracy
from scripts.audit_ncac_rotamers import check_rotation
from scripts.audit_prospective_ncac import snapshot_check
from scripts.compare_branched_completion import measure
from scripts.research_dunbrack_prior import sha, write
from scripts.research_rosetta_joint import verify
from scripts.rosetta_fixed_contract import choose_scored_own
from scripts.summarize_rosetta_fixed_scoring import source_records_preserved

from dvbfixer.diagnose.geometry import signed_dihedral_degrees
from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.pdb_component_reconstruction import _component_ref, _parse_atom_lines
from dvbfixer.protein_completion.branched_source import plan_branched_source
from dvbfixer.protein_completion.contracts import Addition


def audit(plan_path,generation_root,score_run_root,output):
    raw=plan_path.read_bytes()
    plan=json.loads(raw)
    for p,d in {**plan['inputs'],**plan['source_hashes']}.items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('frozen mean protocol/input changed')
    if (generation_root/'executed-plan.json').read_bytes()!=raw:
        raise ValueError('executed mean plan differs')
    geometry_raw=(generation_root/'geometry-before-scoring.json').read_bytes()
    geometry=json.loads(geometry_raw)
    if geometry['plan_sha256']!=sha(raw) or geometry['native_coordinates_read']:
        raise ValueError('generation/native-order binding differs')
    score_root=generation_root/'score-plan'
    score_plan_raw=(score_root/'plan.json').read_bytes()
    score_plan=json.loads(score_plan_raw)
    verify(score_plan,score_root)
    scored_raw=(score_run_root/'evidence.json').read_bytes()
    scored=json.loads(scored_raw)
    if scored['plan_sha256']!=sha(score_plan_raw) or (score_run_root/'executed-plan.json').read_bytes()!=score_plan_raw:
        raise ValueError('executed score plan differs')
    declarations={c['id']:c for c in score_plan['cases']}
    if set(geometry['results'])!={c['id'] for c in plan['cases']} or set(scored['results'])!=set(declarations):
        raise ValueError('case/view inventory differs')
    ledger,replay={},{}
    snapshots,feasible_total=0,0
    for case in plan['cases']:
        name=case['id']
        result=geometry['results'][name]
        if case['status']!='generation_ready':
            if result['status']!='parent_not_executed' or name+'-means' in declarations:
                raise ValueError('inherited nonexecution generated/scored')
            ledger[name]=dict(status='parent_not_executed',group=case['group'])
            continue
        source=Path(case['source_path']).read_bytes()
        request=plan_branched_source(source,ComponentInstanceRef(**case['target'])).request
        declaration=declarations[name+'-means']
        if len(result['candidates'])!=3:
            raise ValueError('complete three-well proposal pool required')
        feasible,materialized,audits=[],[],[]
        for i,(candidate,sample) in enumerate(zip(result['candidates'],case['samples'],strict=True)):
            if candidate['index']!=i or candidate['well']!=sample['well'] or candidate['raw_weight']!=sample['probability']:
                raise ValueError('candidate/well/raw-weight binding differs')
            check_rotation(case['parent_points'],candidate['points'],candidate['rotation_degrees'])
            delta=(sample['mean_degrees']-case['parent_chi_degrees']+180.)%360.-180.
            if delta!=candidate['rotation_degrees']:
                raise ValueError('declared mean rotation differs')
            branch='CG1' if request.residue_name=='VAL' else 'OG1'
            chi=signed_dihedral_degrees(*(candidate['points'][n] for n in ('N','CA','CB',branch)))
            if chi is None or abs((chi-sample['mean_degrees']+180.)%360.-180.)>1e-8:
                raise ValueError('actual chi does not match private API mean')
            added=tuple(Addition(ExactAtomRef(request.residue,n),n[0],tuple(candidate['points'][n])) for n in request.missing_atoms)
            common,pdb=measure(source,request,added,anchored=False)
            if json.loads(json.dumps(asdict(common)))!=candidate['common_audit']:
                raise ValueError('independent common geometry differs')
            filename=f'{name}-means/own-{i}.pdb'
            if pdb!=(score_root/filename).read_bytes() or not source_records_preserved(source,pdb,declaration['generated']):
                raise ValueError('mean materialization changed source records')
            if common.accepted:
                feasible.append(i)
            audits.append(dict(index=i,common_accepted=common.accepted,observed_preserved=common.observed_preserved,
                               raw_findings=common.raw['findings']))
            materialized.append((added,pdb))
        eligible=[c for c in result['candidates'] if c['index'] in feasible]
        prior=sorted(eligible,key=lambda c:(-c['raw_weight'],c['index']))[0]['index'] if eligible else None
        if prior!=result['prior_selected_index']:
            raise ValueError('independent prior choice differs')
        scored_row=scored['results'][name+'-means']
        if [c['index'] for c in declaration['cells']]!=feasible or [c['index'] for c in scored_row['cells']]!=feasible:
            raise ValueError('complete feasible score pool differs')
        contexts=[]
        for c,d in zip(scored_row['cells'],declaration['cells'],strict=True):
            if any(c[k]!=v for k,v in d.items()):
                raise ValueError('scored mean declaration changed')
            contexts.append(snapshot_check(score_root,score_run_root,declaration,c))
            snapshots+=1
        selection=choose_scored_own(scored_row['cells'])
        if selection!=scored_row['selection'] or selection!=json.loads((score_run_root/(name+'-means')/'selection-before-native.json').read_bytes()):
            raise ValueError('independent fixed-heavy choice differs')
        feasible_total+=len(feasible)
        ledger[name]=dict(status='evaluated',group=case['group'],candidates=audits,prior_selected_index=prior,
            scored_selection=selection,context_coordinates_identical=not contexts or all(c==contexts[0] for c in contexts),
            cb_and_source_fixed=True,library_mean_transport_verified=True)
        replay[name]=(request,source,materialized)
    selection_path=output.with_name(output.stem+'-selection.json')
    if output.exists() or selection_path.exists():
        raise FileExistsError('fresh mean selection/audit destinations required')
    write(selection_path,dict(schema='dvbfixer.dunbrack-means-selection.v1',cases=ledger,native_coordinates_read=False,
        plan_sha256=sha(raw),geometry_sha256=sha(geometry_raw),score_run_sha256=sha(scored_raw),
        validated_snapshots=snapshots,common_feasible_candidates=feasible_total,table_data_included=False,
        audit_script_sha256=sha(Path(__file__).read_bytes())))
    # No native/reference data opened before this all-case publication boundary.
    old_plan=json.loads(Path('tmp/pdbfixer-research/ncac-expansion-plan-01.json').read_bytes())
    new_plan=json.loads(Path('tmp/pdbfixer-research/ncac-new-system-masks-01/plan.json').read_bytes())
    old_parent=json.loads(Path('docs/research/ncac-relative-rotamers-evidence.json').read_bytes())
    new_parent=json.loads(Path('docs/research/ncac-new-system-evidence.json').read_bytes())
    originals={c['id']:c for c in old_plan['cases']+new_plan['cases']}
    reports={}
    for name,(request,source,pool) in replay.items():
        group=ledger[name]['group']
        original=originals[name]
        mask_root=Path(old_plan['mask_root']) if group=='old_eleven' else Path('tmp/pdbfixer-research/ncac-new-system-masks-01')
        mask_plan=json.loads((mask_root/'plan.json').read_bytes())
        reference=(mask_root/original['reference']).read_bytes()
        if sha(reference)!=mask_plan['artifacts'][original['reference']]:
            raise ValueError('native reference changed')
        truth={a.name:a.position for a in _parse_atom_lines(reference.decode().splitlines(keepends=True)) if _component_ref(a)==request.residue}
        prior=ledger[name]['prior_selected_index']
        score=ledger[name]['scored_selection']['selected_index']
        parent=old_parent['cases'][name] if group=='old_eleven' else new_parent['cases'][name]['methods']['ncac']
        parent_scored=parent['selected_rmsd'] if group=='old_eleven' else parent['selected_accuracy']['rounded']['exact_name_rmsd_angstrom']
        reports[name]=dict(group=group,prior_index=prior,scored_index=score,
            prior_accuracy=None if prior is None else accuracy(request,*pool[prior],truth),
            scored_accuracy=None if score is None else accuracy(request,*pool[score],truth),
            parent_relative_scored_rmsd=parent_scored,context_coordinates_identical=ledger[name]['context_coordinates_identical'])
    aggregate={}
    for group in ('old_eleven','new_inputs_eleven'):
        subset=[r for r in reports.values() if r['group']==group and r['prior_index'] is not None and r['scored_index'] is not None and r['context_coordinates_identical']]
        aggregate[group]=dict(planned=sum(c['group']==group for c in plan['cases']),generated=sum(c['group']==group and c['status']=='generation_ready' for c in plan['cases']),
            mutual_selected=len(subset),prior_vs_score_disagreements=sum(r['prior_index']!=r['scored_index'] for r in subset),
            prior_rounded_exact_median=statistics.median(r['prior_accuracy']['rounded']['exact_name_rmsd_angstrom'] for r in subset) if subset else None,
            scored_rounded_exact_median=statistics.median(r['scored_accuracy']['rounded']['exact_name_rmsd_angstrom'] for r in subset) if subset else None,
            parent_relative_scored_median=statistics.median(r['parent_relative_scored_rmsd'] for r in subset) if subset else None)
    write(output,dict(schema='dvbfixer.dunbrack-means-audit.v1',cases=reports,aggregate=aggregate,
        selection_ledger_name=selection_path.name,selection_ledger_sha256=sha(selection_path.read_bytes()),
        table_data_included=False,production_changed=False,scientific_qualification=False,
        audit_script_sha256=sha(Path(__file__).read_bytes())))
    print(json.dumps(aggregate,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--generation',type=Path,required=True)
    parser.add_argument('--score-run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    audit(args.plan,args.generation,args.score_run,args.output)
