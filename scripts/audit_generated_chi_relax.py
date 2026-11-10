"""Independent mobility/raw-rounded geometry audit, then native reporting.

No PyRosetta import, coordinate generation, minimization or score evaluation.
Final score belongs to the raw target-chi-refined pose, not a fresh PDB import.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import statistics
from dataclasses import asdict
from pathlib import Path

from scripts.audit_ncac_expansion import accuracy
from scripts.compare_branched_completion import measure
from scripts.research_dunbrack_prior import sha, write
from scripts.research_rosetta_repack import identity

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.pdb_component_reconstruction import _component_ref, _parse_atom_lines
from dvbfixer.protein_completion.branched_source import plan_branched_source
from dvbfixer.protein_completion.contracts import Addition

BACKBONE_H={'H','H1','H2','H3','HA','HA2','HA3','HXT'}


def mobility_check(before,after,case,bound):
    a={tuple(r['key']):r for r in before}
    b={tuple(r['key']):r for r in after}
    if len(a)!=len(before) or len(b)!=len(after) or a.keys()!=b.keys():
        raise ValueError('refined atom inventory differs')
    distal={tuple(k) for k in case['generated'] if k[-1]!='CB'}
    protected,mobile=[],[]
    for key,row in a.items():
        other=b[key]
        expected=(not row['is_virtual'] and (key in distal if not row['is_hydrogen'] else
            key[:4]==tuple(case['target_key']) and not row['is_backbone'] and key[-1] not in BACKBONE_H))
        if row['mobile']!=expected or any(other[k]!=row[k] for k in ('mobile','is_hydrogen','is_virtual','is_backbone')):
            raise ValueError('unexpected mobility or atom-type change')
        if not all(math.isfinite(x) for x in other['xyz']):
            raise ValueError('nonfinite refined coordinate')
        delta=math.dist(row['xyz'],other['xyz'])
        (mobile if expected else protected).append(delta)
    maximum=max(protected,default=0.)
    if maximum>bound:
        raise ValueError('protected atom moved; no restore/fallback allowed')
    return dict(protected_max_motion_angstrom=maximum,mobile_max_motion_angstrom=max(mobile,default=0.))


def snapshot_atoms(snapshot):
    physical=snapshot['heavy_atoms']
    residue_names={identity(a)[:3]:identity(a)[3] for a in physical}
    result={identity(a):a['xyz'] for a in physical}
    for a in snapshot['hydrogens']:
        key=(a['chain'],a['resid'],a['icode'],residue_names[(a['chain'],a['resid'],a['icode'])],a['atom'])
        if key in result:
            raise ValueError('duplicate H/physical inventory')
        result[key]=a['xyz']
    for a in snapshot['virtual_atoms']:
        key=identity(a)
        if key in result or not a['virtual_by_rosetta_type']:
            raise ValueError('ambiguous virtual inventory')
        result[key]=a['xyz']
    return result


def audit(plan_path,run_root,output):
    raw=plan_path.read_bytes()
    plan=json.loads(raw)
    for p,d in {**plan['inputs'],**plan['source_hashes']}.items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('frozen relaxation input changed')
    if (run_root/'executed-plan.json').read_bytes()!=raw:
        raise ValueError('executed relaxation protocol differs')
    run_raw=(run_root/'evidence.json').read_bytes()
    run=json.loads(run_raw)
    if run['plan_sha256']!=sha(raw) or run['native_coordinates_read'] or run['weights_modified']:
        raise ValueError('run binding/native/weight contract differs')
    if set(run['results'])!={c['id'] for c in plan['cases']}:
        raise ValueError('relaxation pool omitted')
    ledger,replay={},{}
    for case in plan['cases']:
        name=case['id']
        source=Path(case['source']).read_bytes()
        key=case['case']['target_key']
        request=plan_branched_source(source,ComponentInstanceRef(1,key[0],str(key[1]),key[2].strip())).request
        results=run['results'][name]['cells']
        if [c['index'] for c in results]!=[c['declaration']['index'] for c in case['cells']]:
            raise ValueError('feasible-start inventory differs')
        reports,candidates=[],{}
        for original,result in zip(case['cells'],results,strict=True):
            folder=run_root/name/original['declaration']['id']
            if json.loads((folder/'result.json').read_bytes())!=result:
                raise ValueError('retained refinement result differs')
            if result['coordinate_restore_performed']:
                raise ValueError('coordinate restoration forbidden')
            if result['status']!='refined':
                reports.append(dict(index=result['index'],status=result['status'],error_type=result['error_type'],error=result['error']))
                continue
            if not result['original_pose_exact'] or not result['graph_unchanged'] or not result['only_target_chi_enabled'] or result['convergence_claim']:
                raise ValueError('refinement scope/replay claims differ')
            old=Path(original['snapshot']).read_bytes()
            if old!=(folder/'pose-snapshot.json.gz').read_bytes():
                raise ValueError('initial complete pose not byte-identical')
            before_raw=(folder/'before-atoms.json').read_bytes()
            after_raw=(folder/'after-atoms.json').read_bytes()
            if sha(before_raw)!=result['before_sha256'] or sha(after_raw)!=result['after_sha256']:
                raise ValueError('atom trace digest differs')
            before=json.loads(before_raw)['atoms']
            after=json.loads(after_raw)['atoms']
            expected=snapshot_atoms(json.loads(gzip.decompress(old)))
            if {tuple(a['key']):a['xyz'] for a in before}!=expected:
                raise ValueError('mobility trace not bound to original all-atom pose')
            moves=mobility_check(before,after,case['case'],plan['protected_numerical_bound_angstrom'])
            for k,v in moves.items():
                if not math.isclose(v,result['movement'][k],rel_tol=0.,abs_tol=1e-12):
                    raise ValueError('reported mobility differs')
            if result['initial_score']!=original['saved']['total_score'] or not math.isfinite(result['final_raw_score']):
                raise ValueError('score binding/nonfinite result differs')
            if set(result['score_terms'])!=set(original['saved']['score_terms']) or any(t['weight']!=original['saved']['score_terms'][n]['weight'] for n,t in result['score_terms'].items()):
                raise ValueError('refinement score weights changed')
            if not math.isclose(sum(t['raw']*t['weight'] for t in result['score_terms'].values()),result['final_raw_score'],rel_tol=1e-10,abs_tol=1e-8):
                raise ValueError('refined score terms do not sum to total')
            points={a['key'][-1]:a['xyz'] for a in after if a['key'][:4]==key and not a['is_hydrogen'] and not a['is_virtual']}
            added=tuple(Addition(ExactAtomRef(request.residue,n),n[0],tuple(points[n])) for n in request.missing_atoms)
            common,pdb=measure(source,request,added,anchored=False)
            start_pdb=Path(plan['score_plan_root'])/original['declaration']['pdb']
            initial_by_key={tuple(a['key']):a for a in before}
            mobile_heavy=max((math.dist(a['xyz'],initial_by_key[tuple(a['key'])]['xyz']) for a in after
                              if a['mobile'] and not a['is_hydrogen']),default=0.)
            mobile_h=max((math.dist(a['xyz'],initial_by_key[tuple(a['key'])]['xyz']) for a in after
                          if a['mobile'] and a['is_hydrogen']),default=0.)
            reports.append(dict(index=result['index'],status='audited',common_audit=asdict(common),movement=moves,
                                mobile_heavy_max_motion_angstrom=mobile_heavy,mobile_h_max_motion_angstrom=mobile_h,
                                rounded_heavy_output_equals_start=pdb==start_pdb.read_bytes(),
                                initial_score=result['initial_score'],final_raw_score=result['final_raw_score']))
            candidates[result['index']]=(added,pdb)
        blockers=[r for r in reports if r['status']!='audited']
        feasible=[r for r in reports if r['status']=='audited' and r['common_audit']['accepted']]
        selected=min(feasible,key=lambda r:(r['final_raw_score'],r['index']))['index'] if feasible and not blockers else None
        ledger[name]=dict(group=case['group'],candidates=reports,selected_index=selected,
            selection_status='execution_blocker_in_pool' if blockers else 'selected' if selected is not None else 'empty_geometry_feasible_pool')
        replay[name]=(request,source,candidates)
    selection_path=output.with_name(output.stem+'-selection.json')
    if output.exists() or selection_path.exists():
        raise FileExistsError('fresh relaxation selection/report destinations required')
    write(selection_path,dict(schema='dvbfixer.generated-chi-relax-selection.v1',cases=ledger,
        inherited_nonexecutions=plan['inherited_nonexecutions'],native_coordinates_read=False,
        plan_sha256=sha(raw),run_sha256=sha(run_raw),audit_script_sha256=sha(Path(__file__).read_bytes()),
        score_semantics=plan['score_semantics'],table_data_included=False))
    # All strict complete-pool choices committed before native reporting.
    mean_protocol=json.loads(Path('tmp/pdbfixer-research/dunbrack-means-plan-01.json').read_bytes())
    old_plan=json.loads(Path('tmp/pdbfixer-research/ncac-expansion-plan-01.json').read_bytes())
    new_plan=json.loads(Path('tmp/pdbfixer-research/ncac-new-system-masks-01/plan.json').read_bytes())
    originals={c['id']:c for c in old_plan['cases']+new_plan['cases']}
    parent=json.loads(Path('docs/research/dunbrack-means-evidence.json').read_bytes())
    reports={}
    for name,(request,source,candidates) in replay.items():
        group=ledger[name]['group']
        mask_root=Path(old_plan['mask_root']) if group=='old_eleven' else Path('tmp/pdbfixer-research/ncac-new-system-masks-01')
        masks=json.loads((mask_root/'plan.json').read_bytes())
        original=originals[name]
        reference=(mask_root/original['reference']).read_bytes()
        if sha(reference)!=masks['artifacts'][original['reference']]:
            raise ValueError('native reporting reference differs')
        truth={a.name:a.position for a in _parse_atom_lines(reference.decode().splitlines(keepends=True)) if _component_ref(a)==request.residue}
        candidate_reports={str(i):accuracy(request,*data,truth) for i,data in candidates.items()}
        selected=ledger[name]['selected_index']
        reports[name]=dict(group=group,selected_index=selected,
            selected_accuracy=None if selected is None else candidate_reports[str(selected)],candidate_accuracy=candidate_reports,
            parent_mean_scored_rmsd=parent['cases'][name]['scored_accuracy']['rounded']['exact_name_rmsd_angstrom'],
            parent_selected_index=parent['cases'][name]['scored_index'],
            selected_rounded_output_equals_start=None if selected is None else next(c['rounded_heavy_output_equals_start'] for c in ledger[name]['candidates'] if c['index']==selected))
    aggregate={}
    for group in ('old_eleven','new_inputs_eleven'):
        subset=[r for r in reports.values() if r['group']==group and r['selected_index'] is not None]
        aggregate[group]=dict(planned=sum(c['group']==group for c in mean_protocol['cases']),executed_pools=sum(c['group']==group for c in ledger.values()),
            selected=len(subset),refined_rounded_exact_median=statistics.median(r['selected_accuracy']['rounded']['exact_name_rmsd_angstrom'] for r in subset) if subset else None,
            changed_selected_indices=sum(r['selected_index']!=r['parent_selected_index'] for r in subset),
            selected_rounded_outputs_unchanged=sum(r['selected_rounded_output_equals_start'] for r in subset),
            refined_generated_val_symmetry_median=statistics.median(r['selected_accuracy']['rounded']['generated_only_val_symmetry_rmsd_angstrom'] for r in subset) if subset else None,
            matched_parent_mean_scored_median=statistics.median(r['parent_mean_scored_rmsd'] for r in subset) if subset else None)
    write(output,dict(schema='dvbfixer.generated-chi-relax-audit.v1',cases=reports,aggregate=aggregate,
        selection_ledger_name=selection_path.name,selection_ledger_sha256=sha(selection_path.read_bytes()),
        score_semantics=plan['score_semantics'],production_changed=False,table_data_included=False,scientific_qualification=False,
        audit_script_sha256=sha(Path(__file__).read_bytes())))
    print(json.dumps(aggregate,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    audit(args.plan,args.run,args.output)
