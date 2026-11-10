"""Independent saved-state H/phase/replay audit, no optimization/native metrics."""
from __future__ import annotations

import argparse
import gzip
import json
import math
from pathlib import Path

from scripts.audit_legacy_stereo_trace import state_report
from scripts.research_dunbrack_prior import sha, write
from scripts.research_legacy_comparator import verify

AMENDMENT=Path('tmp/pdbfixer-research/legacy-start-audit-amendment-01.json')


def load_atoms(path):
    return json.loads(gzip.decompress(path.read_bytes()))['atoms']


def coordinate_delta(first,second):
    """Complete immutable identities; never compare just their intersection."""
    a={tuple(r['key']):r for r in first}
    b={tuple(r['key']):r for r in second}
    if len(a)!=len(first) or len(b)!=len(second) or a.keys()!=b.keys():
        raise ValueError('complete coordinate comparison identity differs')
    result={}
    for label,is_h in (('heavy',False),('hydrogen',True)):
        values=[]
        for key,row in a.items():
            if row['element']!=b[key]['element']:
                raise ValueError('coordinate comparison element differs')
            if row['element'] is None or (row['element'] in {'H','D'})!=is_h:
                continue
            values.append(math.dist(row['xyz'],b[key]['xyz']))
        result[label]=dict(count=len(values),maximum_angstrom=max(values,default=None),
            rms_angstrom=math.sqrt(sum(v*v for v in values)/len(values)) if values else None,
            exact=all(v==0 for v in values))
    return result


def system_inventory(system_xml,topology,rows,target):
    import openmm

    system=openmm.XmlSerializer.deserialize(system_xml)
    indices=[i for i,r in enumerate(rows) if r['key'][:3]==[target['chain_id'],target['sequence_number'],target['insertion_code']] and r['key'][-1]=='CB']
    if len(indices)!=1 or system.getNumParticles()!=len(rows):
        raise ValueError('captured system/atom-map target inventory differs')
    cb=indices[0]
    adjacent={i:set() for i in range(len(rows))}
    for a,b in topology['bonds']:
        adjacent[a].add(b)
        adjacent[b].add(a)
    forces=[]
    cb_impropers=[]
    compound=[]
    for force in system.getForces():
        forces.append(dict(type=type(force).__name__,name=force.getName()))
        if isinstance(force,openmm.PeriodicTorsionForce):
            for i in range(force.getNumTorsions()):
                a,b,c,d,*_=force.getTorsionParameters(i)
                particles=(a,b,c,d)
                if cb in particles and set(particles)-{cb}<=adjacent[cb]:
                    cb_impropers.append([rows[j]['key'] for j in particles])
        if isinstance(force,openmm.CustomCompoundBondForce):
            for i in range(force.getNumBonds()):
                particles,_=force.getBondParameters(i)
                if cb in particles:
                    compound.append(dict(expression=force.getEnergyFunction(),particles=[rows[j]['key'] for j in particles]))
    return dict(forces=forces,target_cb_centered_periodic_impropers=cb_impropers,
        compound_terms_containing_target_cb=compound,
        interpretation='force inventory only, not a force-isolation or barrier certificate')


def capture_audit(plan_path,root,output):
    raw=plan_path.read_bytes()
    plan=json.loads(raw)
    verify(plan)
    amendment_raw=AMENDMENT.read_bytes()
    amendment=json.loads(amendment_raw)
    for name,digest in amendment['inputs'].items():
        if sha(Path(name).read_bytes())!=digest:
            raise ValueError('explicit heavy-coordinate reporting amendment changed')
    executed_raw=(root/'evidence.json').read_bytes()
    executed=json.loads(executed_raw)
    if executed['plan_sha256']!=sha(raw) or (root/'executed-plan.json').read_bytes()!=raw:
        raise ValueError('start capture execution binding differs')
    expected=[(c['id'],s) for c in plan['cases'] for s in plan['protocol']['seeds']]
    if [(r['id'],r['seed']) for r in executed['cases']]!=expected:
        raise ValueError('complete H-start enrollment differs')
    by_case={c['id']:c for c in plan['cases']}
    results=[]
    h_after={}
    starts={}
    for record in executed['cases']:
        case=by_case[record['id']]
        folder=root/record['id']/f"seed-{record['seed']}"
        for name,digest in record['artifacts'].items():
            if sha((folder/name).read_bytes())!=digest:
                raise ValueError('captured full state artifact changed')
        row=dict(id=record['id'],seed=record['seed'],status=record['status'])
        if record['status']!='completed':
            results.append(row)
            continue
        capture=json.loads((folder/'capture.json').read_bytes())
        if not capture['hooks_restored'] or capture['plan_sha256']!=sha(raw) or capture['seed']!=record['seed'] or capture['id']!=record['id']:
            raise ValueError('native start-control hook contract differs')
        h=[]
        for call in capture['h_calls']:
            if not call['python_random_state_restored'] or call['seed']!=record['seed']:
                raise ValueError('H-control Python random state not restored')
            before=load_atoms(folder/f"h-{call['call']}-before.json.gz")
            after=load_atoms(folder/f"h-{call['call']}-after.json.gz")
            delta=coordinate_delta([r for r in before if r['element'] not in {'H','D'}],
                                   [r for r in after if r['element'] not in {'H','D'}])
            h.append(dict(call=call['call'],heavy_delta=delta['heavy'],
                          hydrogen_count=sum(r['element'] in {'H','D'} for r in after)))
            h_after[(record['id'],record['seed'],call['call'])]=after
        phases=[]
        for phase in capture['phases']:
            sides={}
            for side in ('before','after'):
                prefix=phase[side]['prefix']
                rows=load_atoms(folder/f'{prefix}-atoms.json.gz')
                sides[side]=dict(stereo=state_report(rows,case['target']),platform=phase[side]['platform'])
                if phase['phase']==1 and side=='before':
                    starts[(record['id'],record['seed'])]=rows
                    inventory=system_inventory((folder/f'{prefix}-system.xml').read_text(),
                        json.loads((folder/f'{prefix}-topology.json').read_bytes()),rows,case['target'])
            if phase['before']['platform']!=phase['after']['platform']:
                raise ValueError('platform changed inside native phase')
            phases.append(dict(phase=phase['phase'],settings=phase['settings'],**sides))
        row.update(hydrogen_calls=h,phases=phases,force_inventory=inventory)
        results.append(row)
    comparisons=[]
    first,second=plan['protocol']['seeds']
    for case in plan['cases']:
        keys=[k for k in h_after if k[:2]==(case['id'],first)]
        for key in keys:
            other=(case['id'],second,key[2])
            if other in h_after:
                comparisons.append(dict(id=case['id'],hydrogen_call=key[2],
                    post_h_delta=coordinate_delta(h_after[key],h_after[other]),
                    native_phase1_start_delta=coordinate_delta(starts[(case['id'],first)],starts[(case['id'],second)])))
        if not keys and (case['id'],first) in starts and (case['id'],second) in starts:
            comparisons.append(dict(id=case['id'],hydrogen_call=None,
                note='keep-H branch: seed labels did not intervene in H placement; numerical repeat control only',
                native_phase1_start_delta=coordinate_delta(starts[(case['id'],first)],starts[(case['id'],second)])))
    write(output,dict(schema='dvbfixer.legacy-start-capture-audit.v2',cases=results,seed_comparisons=comparisons,
        reporting_amendment_sha256=sha(amendment_raw),
        heavy_coordinate_policy='report exact/nonexact numerical motion verbatim; no observation-preservation or scientific geometry acceptance claim',
        summary=dict(enrolled=len(expected),completed=sum(r['status']=='completed' for r in results),
            native_phases=sum(len(r.get('phases',[])) for r in results)),
        plan_sha256=sha(raw),run_sha256=sha(executed_raw),audit_script_sha256=sha(Path(__file__).read_bytes()),
        native_coordinates_read=False,native_accuracy_reported=False,production_changed=False,
        profile_selected=False,full_causal_attribution=False))


def replay_audit(plan_path,root,output):
    raw=plan_path.read_bytes()
    plan=json.loads(raw)
    verify(plan)
    run_raw=(root/'evidence.json').read_bytes()
    run=json.loads(run_raw)
    if run['plan_sha256']!=sha(raw) or (root/'executed-plan.json').read_bytes()!=raw:
        raise ValueError('saved native replay plan differs')
    parent=json.loads(Path(plan['parent_plan']).read_bytes())
    targets={c['id']:c['target'] for c in parent['cases']}
    p=plan['protocol']
    expected=[(c,s,m,r) for c in p['focal_cases'] for s in p['seeds'] for m in p['platforms'] for r in range(1,p['repeats']+1)]
    if [(r['id'],r['seed'],r['mode'],r['repeat']) for r in run['cases']]!=expected:
        raise ValueError('native replay complete enrollment differs')
    results=[]
    endpoints={}
    for row in run['cases']:
        folder=root/row['id']/f"seed-{row['seed']}"/row['mode']/f"repeat-{row['repeat']}"
        for name,digest in row['artifacts'].items():
            if sha((folder/name).read_bytes())!=digest:
                raise ValueError('native replay artifact changed')
        retained=dict(row)
        if row['status']!='completed':
            results.append(retained)
            continue
        source=Path(plan['capture_root'])/row['id']/f"seed-{row['seed']}"
        capture=json.loads((source/'capture.json').read_bytes())
        phase=next(x for x in capture['phases'] if x['phase']==p['phase'])
        if row['settings']!=phase['settings']:
            raise ValueError('saved-state native replay optimizer settings changed')
        template=load_atoms(source/f"{phase['before']['prefix']}-atoms.json.gz")
        original_after=load_atoms(source/f"{phase['after']['prefix']}-atoms.json.gz")
        sides={}
        values={}
        for stage in ('before','after'):
            saved=json.loads(gzip.decompress((folder/f'{stage}.json.gz').read_bytes()))
            if len(saved['positions_angstrom'])!=len(template) or not math.isfinite(saved['energy_kj_mol']):
                raise ValueError('native replay finite/particle inventory differs')
            rows=[dict(atom,xyz=xyz) for atom,xyz in zip(template,saved['positions_angstrom'],strict=True)]
            if any(not all(math.isfinite(x) for x in a['xyz']) for a in rows):
                raise ValueError('nonfinite replay coordinate')
            sides[stage]=dict(stereo=state_report(rows,targets[row['id']]),energy_kj_mol=saved['energy_kj_mol'])
            values[stage]=rows
        retained.update(states=sides,initial_coordinate_delta=coordinate_delta(template,values['before']),
            endpoint_delta_vs_captured_phase=coordinate_delta(original_after,values['after']))
        endpoints[(row['id'],row['seed'],row['mode'],row['repeat'])]=values['after']
        results.append(retained)
    repeats=[]
    platform_pairs=[]
    for case_id in p['focal_cases']:
        for seed in p['seeds']:
            for mode in p['platforms']:
                a,b=((case_id,seed,mode,i) for i in (1,2))
                if a in endpoints and b in endpoints:
                    repeats.append(dict(id=case_id,seed=seed,mode=mode,delta=coordinate_delta(endpoints[a],endpoints[b])))
            a,b=((case_id,seed,m,1) for m in p['platforms'])
            if a in endpoints and b in endpoints:
                platform_pairs.append(dict(id=case_id,seed=seed,delta=coordinate_delta(endpoints[a],endpoints[b])))
    write(output,dict(schema='dvbfixer.legacy-start-replay-audit.v1',cases=results,repeated_starts=repeats,
        cross_platform_pairs=platform_pairs,summary=dict(enrolled=len(expected),completed=sum(r['status']=='completed' for r in results)),
        plan_sha256=sha(raw),run_sha256=sha(run_raw),audit_script_sha256=sha(Path(__file__).read_bytes()),
        native_coordinates_read=False,native_accuracy_reported=False,profile_selected=False,
        production_changed=False,full_causal_attribution=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    for name in ('capture','replay'):
        p=sub.add_parser(name)
        p.add_argument('--plan',type=Path,required=True)
        p.add_argument('--run',type=Path,required=True)
        p.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.command=='capture':
        capture_audit(args.plan,args.run,args.output)
    else:
        replay_audit(args.plan,args.run,args.output)
