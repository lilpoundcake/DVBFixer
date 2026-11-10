"""Saved-only actual native-phase stereo readback, no native reference metrics."""
from __future__ import annotations

import argparse
import gzip
import json
import math
from pathlib import Path

import numpy as np
from scripts.audit_legacy_comparator import atoms, motion
from scripts.branched_stereo_oracle import rdkit_labels
from scripts.research_dunbrack_prior import sha, write
from scripts.research_legacy_comparator import verify


def determinant(points,center,first,second,third):
    xyz={n:np.asarray(v,dtype=float) for n,v in points.items()}
    if any(not np.isfinite(v).all() for v in xyz.values()):
        raise ValueError('nonfinite stereo coordinates')
    value=float(np.linalg.det(np.array([xyz[first]-xyz[center],xyz[second]-xyz[center],xyz[third]-xyz[center]])))
    if not math.isfinite(value):
        raise ValueError('nonfinite stereo determinant')
    return value


def state_report(rows,target):
    grouped={}
    identity=set()
    for atom in rows:
        key=tuple(atom['key'])
        if key in identity:
            raise ValueError('duplicate native-context atom identity')
        identity.add(key)
        if atom['element'] in {'H','D',None}:
            continue
        grouped.setdefault(key[:4],{})[key[-1]]=atom['xyz']
    ca,thr=[],[]
    selected=None
    for key,points in grouped.items():
        if {'N','CA','C','CB'}<=points.keys():
            ca.append(dict(residue=key,volume_angstrom3=determinant(points,'CA','N','C','CB')))
        if key[-1]=='THR' and {'CA','CB','OG1','CG2'}<=points.keys():
            thr.append(dict(residue=key,volume_angstrom3=determinant(points,'CB','CA','OG1','CG2')))
        if key[:3]==(target['chain_id'],target['sequence_number'],target['insertion_code']):
            if selected is not None:
                raise ValueError('ambiguous target in context snapshot')
            selected=dict(residue=key,ca_volume_angstrom3=determinant(points,'CA','N','C','CB'),
                target_cip=rdkit_labels(points,key[-1]))
            if key[-1]=='THR':
                selected['beta_volume_angstrom3']=determinant(points,'CB','CA','OG1','CG2')
    if selected is None:
        raise ValueError('target absent from context snapshot')
    return dict(target=selected,negative_ca=[r for r in ca if r['volume_angstrom3']<0],
        negative_thr_beta=[r for r in thr if r['volume_angstrom3']<0],
        protein_ca_centres=len(ca),thr_beta_centres=len(thr),atom_count=len(rows))


def raw_rows(raw):
    return [dict(key=[component.chain_id,component.sequence_number,component.insertion_code,atom.residue_name,name],
                 element=atom.element,xyz=atom.position) for (component,name),atom in atoms(raw).items()]


def audit(plan_path,run_root,output):
    plan_raw=plan_path.read_bytes()
    plan=json.loads(plan_raw)
    verify(plan)
    if (run_root/'executed-plan.json').read_bytes()!=plan_raw:
        raise ValueError('native-phase executed plan differs')
    run_raw=(run_root/'evidence.json').read_bytes()
    run=json.loads(run_raw)
    if run['plan_sha256']!=sha(plan_raw) or run['native_coordinates_read']:
        raise ValueError('native-phase run binding differs')
    if [c['id'] for c in plan['cases']]!=[c['id'] for c in run['cases']]:
        raise ValueError('native-phase complete inventory differs')
    results=[]
    for case,executed in zip(plan['cases'],run['cases'],strict=True):
        folder=run_root/case['id']
        for name,d in executed['artifacts'].items():
            if sha((folder/name).read_bytes())!=d:
                raise ValueError('native-phase artifact changed')
        row=dict(id=case['id'],status=executed['status'])
        if executed['status']!='completed':
            results.append(row)
            continue
        if (folder/'prepared.pdb').read_bytes()!=Path(case['prepared']).read_bytes() or (folder/'prepared.dat').read_bytes()!=Path(case['dat']).read_bytes():
            raise ValueError('native-phase source copy differs')
        trace=json.loads((folder/'trace.json').read_bytes())
        if not trace['hook_restored'] or not trace['args_forwarded_unchanged'] or trace['coordinates_restored'] or trace['plan_sha256']!=sha(plan_raw):
            raise ValueError('phase observer contract differs')
        if [(e['call'],e['stage']) for e in trace['events']]!=[(i,s) for i in range(1,trace['calls']+1) for s in ('before','after')]:
            raise ValueError('native-phase snapshot inventory differs')
        states=[]
        first=None
        for event in trace['events']:
            raw=(folder/event['snapshot']).read_bytes()
            if sha(raw)!=event['sha256']:
                raise ValueError('native-phase snapshot digest differs')
            snapshot=json.loads(gzip.decompress(raw))['atoms']
            if first is None:
                first=snapshot
            states.append(dict(call=event['call'],stage=event['stage'],**state_report(snapshot,case['target'])))
        if first is None:
            raise ValueError('no native minimizeEnergy call recorded')
        first_heavy={tuple(a['key'][:3]+[a['key'][-1]]):a['xyz'] for a in first if a['element'] not in {'H','D',None}}
        initial=raw_rows(Path(case['prepared']).read_bytes())
        input_heavy={tuple(a['key'][:3]+[a['key'][-1]]):a['xyz'] for a in initial}
        if first_heavy.keys()!=input_heavy.keys():
            raise ValueError('prephase heavy inventory differs from prepared input')
        prephase=max(math.dist(v,first_heavy[k]) for k,v in input_heavy.items())
        actual=(folder/'minimized.pdb').read_bytes()
        baseline=Path(case['baseline']).read_bytes()
        row.update(calls=trace['calls'],states=states,prepared_state=state_report(initial,case['target']),
            final_state=state_report(raw_rows(actual),case['target']),
            heavy_motion_before_first_native_call_angstrom=prephase,
            exact_baseline_pdb_repeat=actual==baseline,baseline_endpoint_difference=motion(baseline,actual))
        results.append(row)
    write(output,dict(schema='dvbfixer.legacy-stereo-trace-audit.v1',cases=results,
        summary=dict(enrolled=len(results),completed=sum(r['status']=='completed' for r in results),
            native_calls=sum(r.get('calls',0) for r in results),
            exact_baseline_pdb_repeats=sum(r.get('exact_baseline_pdb_repeat',False) for r in results)),
        plan_sha256=sha(plan_raw),run_sha256=sha(run_raw),audit_script_sha256=sha(Path(__file__).read_bytes()),
        native_coordinates_read=False,native_accuracy_reported=False,profile_selected=False,
        production_changed=False,causal_mechanism_certified=False))
    print(json.dumps(dict(cases=len(results),calls=sum(r.get('calls',0) for r in results),
        exact_baseline_repeats=sum(r.get('exact_baseline_pdb_repeat',False) for r in results)),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    audit(args.plan,args.run,args.output)
