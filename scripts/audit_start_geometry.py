"""Saved-only geometry-preparation factorial audit, no native-reference scoring."""
from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

from scripts.audit_legacy_comparator import atoms, motion
from scripts.audit_legacy_start_controls import coordinate_delta, load_atoms
from scripts.audit_legacy_stereo_trace import determinant, raw_rows, state_report
from scripts.research_dunbrack_prior import sha, write
from scripts.research_legacy_comparator import verify

from dvbfixer.diagnose.geometry import signed_dihedral_degrees
from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.ffutils.dat import DatRecord
from dvbfixer.protein_completion.branched import BRANCHES
from dvbfixer.protein_completion.contracts import (
    Addition,
    Candidate,
    CompletionRequest,
    Observation,
)
from dvbfixer.protein_completion.validation import _audit_candidate


def geometry(rows,bonds,target,residue,rounded=False):
    material=[dict(a,xyz=[float(f'{x:.3f}') for x in a['xyz']]) if rounded else a for a in rows]
    identities=[ExactAtomRef(ComponentInstanceRef(1,a['key'][0],a['key'][1],a['key'][2]),a['key'][-1]) for a in material]
    selected=ComponentInstanceRef(**target)
    branch=BRANCHES[residue]
    missing=('CB',branch.first,branch.second)
    observations=[]
    additions=[]
    for a,ref in zip(material,identities,strict=True):
        if a['element'] in {'H','D',None}:
            continue
        if ref.component==selected and ref.name in missing:
            additions.append(Addition(ref,a['element'],tuple(a['xyz'])))
        else:
            observations.append(Observation(ref,a['element'],tuple(a['xyz'])))
    heavy={o.atom for o in observations}|{a.atom for a in additions}
    links=tuple((identities[a],identities[b]) for a,b in bonds if identities[a] in heavy and identities[b] in heavy
        and identities[a].component!=identities[b].component)
    request=CompletionRequest(sha(json.dumps(material,sort_keys=True,allow_nan=False).encode()),selected,residue,
        tuple(observations),missing,links)
    points={o.atom.name:o.xyz for o in observations if o.atom.component==selected}
    points.update({a.atom.name:a.xyz for a in additions})
    chi=signed_dihedral_degrees(*(points[n] for n in ('N','CA','CB',branch.first)))
    audit=_audit_candidate(request,Candidate(request.digest,tuple(additions),chi),measured_chi=True,pdb_roundtrip=rounded)
    stereo=state_report(material,target)
    expected=dict(CA='S',CB='R') if residue=='THR' else dict(CA='S')
    grouped={}
    for a in material:
        if a['element'] not in {'H','D',None}:
            grouped.setdefault(tuple(a['key'][:4]),{})[a['key'][-1]]=a['xyz']
    global_failures=[]
    for key,xyz in grouped.items():
        if {'N','CA','C','CB'}<=xyz.keys() and determinant(xyz,'CA','N','C','CB')<=.5:
            global_failures.append(dict(residue=key,centre='CA',reason='non_L_or_nearplanar'))
        if key[-1]=='THR' and {'CA','CB','OG1','CG2'}<=xyz.keys() and determinant(xyz,'CB','CA','OG1','CG2')<=.5:
            global_failures.append(dict(residue=key,centre='CB',reason='non_R_or_nearplanar'))
    return dict(accepted=audit.accepted and not global_failures and stereo['target']['target_cip']==expected,
        target_geometry=asdict(audit),stereo=stereo,global_stereo_failures=global_failures,
        scope='actual-frame heavy geometry, no original observation-preservation or H-geometry claim')


def heavy_map(rows):
    mapped={tuple(a['key'][:3]+[a['key'][-1]]):a for a in rows if a['element'] not in {'H','D',None}}
    if len(mapped)!=sum(a['element'] not in {'H','D',None} for a in rows):
        raise ValueError('duplicate heavy identity')
    return mapped


def input_audit(case,folder,arm):
    original=Path(case['prepared']).read_bytes()
    modified=(folder/'prepared.pdb').read_bytes()
    target=ComponentInstanceRef(**case['target'])
    change=json.loads((folder/'input-change.json').read_bytes())
    if change['source_sha256']!=sha(original) or change['input_sha256']!=sha(modified):
        raise ValueError('source geometry rewrite digest differs')
    allowed=set(change['allowed_atoms'])
    expected=set() if arm=='raw' or case['residue']=='VAL' else {'OG1','CG2'} if arm=='beta_only' else set(case['mask'])
    if allowed!=expected:
        raise ValueError('geometry arm changed undeclared atoms')
    ownership=DatRecord.load(case['dat']).added_keys()
    for name in allowed:
        if (target.chain_id,target.sequence_number,target.insertion_code,name) not in ownership:
            raise ValueError('source geometry changed an observed atom')
    changed=[]
    for before,after in zip(original.decode().splitlines(keepends=True),modified.decode().splitlines(keepends=True),strict=True):
        if before==after:
            continue
        if not before.startswith(('ATOM  ','HETATM')):
            raise ValueError('source rewrite changed a non-atom record')
        key=(before[21],before[22:26].strip(),before[26].strip(),before[12:16].strip())
        if key[:3]!=(target.chain_id,target.sequence_number,target.insertion_code) or key[-1] not in allowed or before[:30]+before[54:]!=after[:30]+after[54:]:
            raise ValueError('source rewrite changed more than dat-owned target coordinates')
        changed.append(key[-1])
    if sorted(changed)!=sorted(change['changed_atom_names']):
        raise ValueError('source rewrite change inventory differs')
    branch=BRANCHES[case['residue']]
    before={n:a.position for (r,n),a in atoms(original).items() if r==target}
    after={n:a.position for (r,n),a in atoms(modified).items() if r==target}
    for name,xyz in change['unrounded_points'].items():
        if name not in allowed and tuple(xyz)!=before[name]:
            raise ValueError('unrounded construction moved an immutable target anchor')
        if any(abs(x-y)>.000500001 for x,y in zip(xyz,after[name],strict=True)):
            raise ValueError('rewritten source differs from declared construction beyond PDB rounding')
    original_chi=signed_dihedral_degrees(*(before[n] for n in ('N','CA','CB',branch.first)))
    unrounded_chi=signed_dihedral_degrees(*(change['unrounded_points'][n] for n in ('N','CA','CB',branch.first)))
    rounded_chi=signed_dihedral_degrees(*(after[n] for n in ('N','CA','CB',branch.first)))
    if original_chi is None or unrounded_chi is None or rounded_chi is None:
        raise ValueError('source geometry chi became undefined')
    error=abs((unrounded_chi-original_chi+180)%360-180)
    if error>1e-8:
        raise ValueError('source unrounded chi was not preserved')
    if (folder/'prepared.dat').read_bytes()!=Path(case['dat']).read_bytes():
        raise ValueError('source geometry altered dat restraint ownership')
    return dict(changed_atom_names=changed,unrounded_chi_error_degrees=error,
        rounded_chi_error_degrees=abs((rounded_chi-original_chi+180)%360-180),
        observed_heavy_records_preserved=True,val_identity_control=case['residue']=='VAL')


def audit(plan_path,root,output):
    raw=plan_path.read_bytes()
    plan=json.loads(raw)
    verify(plan)
    run_raw=(root/'evidence.json').read_bytes()
    run=json.loads(run_raw)
    if run['plan_sha256']!=sha(raw) or (root/'executed-plan.json').read_bytes()!=raw or run['coordinates_restored']:
        raise ValueError('source geometry execution contract differs')
    expected=[(c['id'],s,a) for c in plan['cases'] for s in plan['protocol']['seeds'] for a in plan['protocol']['arms']]
    if [(r['id'],r['seed'],r['arm']) for r in run['cases']]!=expected:
        raise ValueError('source geometry complete enrollment differs')
    original_plan_path=Path('tmp/pdbfixer-research/legacy-comparator-plan-01.json')
    if str(original_plan_path) not in plan['inputs']:
        raise ValueError('original masked-source enrollment not bound')
    original_plan=json.loads(original_plan_path.read_bytes())
    verify(original_plan)
    observed_sources={c['id']:Path(c['source']) for c in original_plan['cases']}
    cases={c['id']:c for c in plan['cases']}
    results=[]
    endpoints={}
    for row in run['cases']:
        folder=root/row['id']/f"seed-{row['seed']}"/row['arm']
        for name,digest in row['artifacts'].items():
            if sha((folder/name).read_bytes())!=digest:
                raise ValueError('source geometry retained artifact changed')
        case=cases[row['id']]
        report=dict(id=row['id'],seed=row['seed'],arm=row['arm'],status=row['status'])
        if row['status']!='completed':
            report.update(error=row.get('error'),returncode=row.get('returncode'))
            results.append(report)
            continue
        rewritten=input_audit(case,folder,row['arm'])
        trace=json.loads((folder/'trace.json').read_bytes())
        if not trace['hooks_restored'] or trace['plan_sha256']!=sha(raw) or len(trace['h_calls'])!=1:
            raise ValueError('source geometry controlled H reconstruction differs')
        if (trace['id'],trace['seed'],trace['arm'])!=(row['id'],row['seed'],row['arm']):
            raise ValueError('source geometry trace identity differs')
        for h in trace['h_calls']:
            if h['platform']!='Reference' or h['seed']!=row['seed'] or not h['python_random_state_restored']:
                raise ValueError('source geometry H seed/platform restoration differs')
        bonds=json.loads((folder/'topology.json').read_bytes())['bonds']
        phases=[]
        initial=load_atoms(folder/'phase-1-before.json.gz')
        for phase in trace['phases']:
            if phase['before']['platform']['name']!='Reference' or phase['after']['platform']['name']!='Reference':
                raise ValueError('source geometry native platform differs')
            rows=load_atoms(folder/f"phase-{phase['phase']}-after.json.gz")
            phases.append(dict(phase=phase['phase'],raw=geometry(rows,bonds,case['target'],case['residue']),
                rounded=geometry(rows,bonds,case['target'],case['residue'],rounded=True),movement=coordinate_delta(initial,rows)))
        final_raw=(folder/'minimized.pdb').read_bytes()
        final_map=heavy_map(raw_rows(final_raw))
        if final_map.keys()!=heavy_map(initial).keys():
            raise ValueError('final source geometry heavy identities differ')
        # Rebind actual final PDB heavy coordinates; H remains outside diagnostic scope.
        actual=[dict(a,xyz=final_map[tuple(a['key'][:3]+[a['key'][-1]])]['xyz']) if a['element'] not in {'H','D',None} else a for a in rows]
        final=geometry(actual,bonds,case['target'],case['residue'],rounded=True)
        endpoints[(row['id'],row['seed'],row['arm'])]=actual
        report.update(source_rewrite=rewritten,initial_geometry=geometry(initial,bonds,case['target'],case['residue']),
            phases=phases,final_geometry=final,
            original_observation_motion=motion(observed_sources[row['id']].read_bytes(),final_raw),
            all_prepared_heavy_motion=motion(Path(case['prepared']).read_bytes(),final_raw),
            final_geometry_accepted=final['accepted'],
            emergency_chirality_message=('Emergency chirality' in (folder/'minimize.log').read_text()))
        results.append(report)
    controls=[]
    for case in plan['cases']:
        if case['residue']!='VAL':
            continue
        for seed in plan['protocol']['seeds']:
            reference=(case['id'],seed,'raw')
            for arm in ('beta_only','full_generated'):
                alternative=(case['id'],seed,arm)
                if reference in endpoints and alternative in endpoints:
                    controls.append(dict(id=case['id'],seed=seed,arm=arm,delta=coordinate_delta(endpoints[reference],endpoints[alternative])))
    summaries={}
    for arm in plan['protocol']['arms']:
        rows=[r for r in results if r['arm']==arm]
        completed=[r for r in rows if r['status']=='completed']
        thr=[r for r in completed if cases[r['id']]['residue']=='THR']
        summaries[arm]=dict(enrolled=len(rows),completed=len(completed),execution_blockers=len(rows)-len(completed),
            final_geometry_passes=sum(r['final_geometry_accepted'] for r in completed),
            THR_completed=len(thr),THR_final_geometry_passes=sum(r['final_geometry_accepted'] for r in thr),
            THR_final_CB_R=sum(r['final_geometry']['stereo']['target']['target_cip'].get('CB')=='R' for r in thr))
    write(output,dict(schema='dvbfixer.private-start-geometry-audit.v1',cases=results,summary=summaries,VAL_identity_controls=controls,
        plan_sha256=sha(raw),run_sha256=sha(run_raw),audit_script_sha256=sha(Path(__file__).read_bytes()),
        native_coordinates_read=False,native_accuracy_reported=False,profile_selected=False,production_changed=False,
        original_failure_causally_explained=False,
        interpretation='source geometry intervention with downstream H/anchor/CA-reference mediation, not fixed-potential isolation or production qualification'))
    print(json.dumps(summaries,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    audit(args.plan,args.run,args.output)
