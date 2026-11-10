"""Private library-mean branch transport, not a production completion backend.

Use previously queried per-case means; never export the library, move CB/known
atoms, change geometry gates or read native coordinates during generation.
"""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict
from pathlib import Path

from scripts.research_dunbrack_prior import sha, write


def rotate_distal(points,degrees):
    if not math.isfinite(degrees):
        raise ValueError('finite rotation required')
    distal=set(points)-{'N','CA','C','O','CB'}
    if not {'N','CA','C','O','CB'}<=points.keys() or distal not in ({'CG1','CG2'},{'OG1','CG2'}):
        raise ValueError('whole VAL/THR with fixed backbone and CB required')
    axis=[points['CB'][i]-points['CA'][i] for i in range(3)]
    norm=math.sqrt(sum(x*x for x in axis))
    if norm<1e-8:
        raise ValueError('degenerate CA-CB axis')
    u=[x/norm for x in axis]
    c,s=math.cos(math.radians(degrees)),math.sin(math.radians(degrees))
    output={n:list(v) for n,v in points.items()}
    for n in distal:
        v=[points[n][i]-points['CB'][i] for i in range(3)]
        projection=sum(a*b for a,b in zip(u,v,strict=True))
        cross=[u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]]
        output[n]=[points['CB'][i]+c*v[i]+s*cross[i]+(1-c)*projection*u[i] for i in range(3)]
    return output


def plan(output):
    parent_path=Path('tmp/pdbfixer-research/dunbrack-prior-plan-02.json')
    parent_raw=parent_path.read_bytes()
    parent=json.loads(parent_raw)
    query_root=Path('tmp/pdbfixer-research/dunbrack-prior-run-02')
    queries_raw=(query_root/'selection-before-native.json').read_bytes()
    queries=json.loads(queries_raw)
    inputs=dict(parent['inputs'])
    inputs[str(parent_path)]=sha(parent_raw)
    for name in ('selection-before-native.json','environment.json','weight-semantics.json'):
        inputs[str(query_root/name)]=sha((query_root/name).read_bytes())
    score_settings=Path('tmp/pdbfixer-research/ncac-rotamers-score-plan-01/plan.json')
    inputs[str(score_settings)]=sha(score_settings.read_bytes())
    if queries['plan_sha256']!=sha(parent_raw) or queries['native_coordinates_read']:
        raise ValueError('prior query binding differs')
    original_geo=Path('tmp/pdbfixer-research/ncac-rotamers-score-plan-01/geometry.json')
    originals=json.loads(original_geo.read_bytes())
    cases=[]
    for c in parent['cases']:
        row=queries['results'][c['id']]
        if c['status']!='query_ready':
            cases.append(dict(c,status='parent_not_executed'))
            continue
        if c['group']=='old_eleven':
            candidate=originals['cases'][c['id']]['candidates'][0]
        else:
            p=Path(c['source_path']).parent/'ncac-proposals.json'
            inputs[str(p)]=sha(p.read_bytes())
            candidate=json.loads(p.read_bytes())['candidates'][0]
        cases.append(dict(c,status='generation_ready',parent_points=candidate['points'],samples=row['samples'],
                          parent_chi_degrees=c['candidates'][0]['chi_degrees']))
    hashes=dict(parent['source_hashes'])
    for p in ('scripts/research_dunbrack_means.py','scripts/prospective_ncac_cohort.py'):
        hashes[p]=sha(Path(p).read_bytes())
    for p,d in {**inputs,**hashes}.items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('prior source/query binding changed')
    write(output,dict(schema='dvbfixer.dunbrack-means-plan.v1',cases=cases,inputs=inputs,source_hashes=hashes,
        protocol='three API mean-chi proposals, one per well in well order; rotate generated distal pair around saved CA-CB; CB/backbone/source observations fixed; no optimization or sigma perturbation',
        geometry='unchanged comparator raw/rounded/target/global/steric checks, retain all refusals',
        selection='separate max-raw-family-weight and unchanged ref2015 minimum over complete geometry-feasible mean pool; stable original well index tie; missing feasible score blocks',
        native_policy='native accuracy reporting only after all independent choices commit; not part of generation/scoring',
        scope='retrospective on two already evaluated cohorts; inherited seven nonexecutions retained',
        production_enabled=False,table_redistribution=False))


def generate(plan_path,output):
    from scripts.compare_branched_completion import measure
    from scripts.prospective_ncac_cohort import score_case

    from dvbfixer.diagnose.geometry import signed_dihedral_degrees
    from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
    from dvbfixer.protein_completion.branched_source import plan_branched_source
    from dvbfixer.protein_completion.contracts import Addition

    raw=plan_path.read_bytes()
    record=json.loads(raw)
    for p,d in {**record['inputs'],**record['source_hashes']}.items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('frozen mean-proposal input changed')
    output.mkdir(parents=True,exist_ok=False)
    (output/'executed-plan.json').write_bytes(raw)
    score_root=output/'score-plan'
    score_root.mkdir()
    results,cases,artifacts={},{},{}
    for case in record['cases']:
        name=case['id']
        if case['status']!='generation_ready':
            results[name]=dict(status='parent_not_executed',group=case['group'])
            continue
        source=Path(case['source_path']).read_bytes()
        request=plan_branched_source(source,ComponentInstanceRef(**case['target'])).request
        rows,candidates=[],[]
        for i,sample in enumerate(case['samples']):
            if sample['well']!=i+1:
                raise ValueError('complete well-order sample set required')
            delta=(sample['mean_degrees']-case['parent_chi_degrees']+180.)%360.-180.
            proposed=rotate_distal(case['parent_points'],delta)
            tip='CG1' if request.residue_name=='VAL' else 'OG1'
            actual=signed_dihedral_degrees(*(proposed[n] for n in ('N','CA','CB',tip)))
            if actual is None or abs((actual-sample['mean_degrees']+180.)%360.-180.)>1e-8:
                raise ValueError('transport does not reproduce declared library mean chi')
            added=tuple(Addition(ExactAtomRef(request.residue,n),n[0],tuple(proposed[n])) for n in request.missing_atoms)
            common,pdb=measure(source,request,added,anchored=False)
            rows.append((i,common,pdb))
            candidates.append(dict(index=i,well=sample['well'],raw_weight=sample['probability'],rotation_degrees=delta,
                points=proposed,common_audit=asdict(common),mean_chi_transport_verified=True))
        view=score_case(name,'means',source,request,rows,score_root,artifacts)
        cases[name]=view
        feasible=[c for c in candidates if c['common_audit']['accepted']]
        selected=max(feasible,key=lambda c:(c['raw_weight'],-c['index']))['index'] if feasible else None
        results[name]=dict(status='generated',group=case['group'],candidates=candidates,prior_selected_index=selected)
    write(output/'geometry-before-scoring.json',dict(schema='dvbfixer.dunbrack-means-geometry.v1',results=results,
        plan_sha256=sha(raw),native_coordinates_read=False))
    parent=json.loads(Path('tmp/pdbfixer-research/ncac-rotamers-score-plan-01/plan.json').read_bytes())
    settings={k:parent[k] for k in ('expected_rosetta_revision','score_function','seed','radii','init_flags',
        'hydrogen_protocol','extra_heavy_policy','virtual_atom_policy','topology_scope','selection_policy','native_policy')}
    write(score_root/'plan.json',dict(schema='dvbfixer.rosetta-fixed-scoring-plan.v2',cases=list(cases.values()),
        artifacts=artifacts,source_hashes=record['source_hashes'],**settings,retrospective=True,public_capability=False))
    print(json.dumps(dict(planned=len(results),generated=sum(r['status']=='generated' for r in results.values()),
        proposals=sum(len(r.get('candidates',[])) for r in results.values()),
        feasible=sum(len(c['cells']) for c in cases.values())),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('plan')
    p.add_argument('--output',type=Path,required=True)
    p=sub.add_parser('generate')
    p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.command=='plan':
        plan(args.output)
    else:
        generate(args.plan,args.output)
