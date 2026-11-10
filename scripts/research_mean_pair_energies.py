"""Frozen residue-pair accounting on ALL saved mean prior/score disagreements.

Repeat each original fixed-heavy state exactly before reading its EnergyGraph.
No new candidates/weights, packing, H optimization, coordinate restore or native
input. Pair arithmetic explains model ranking, not native/physical causality.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from unittest.mock import patch

from scripts.research_dunbrack_prior import sha, write

PLAN_ROOT=Path('tmp/pdbfixer-research/dunbrack-means-generation-01/score-plan')
RUN_ROOT=Path('tmp/pdbfixer-research/dunbrack-means-score-run-01')
GEOMETRY=Path('tmp/pdbfixer-research/dunbrack-means-generation-01/geometry-before-scoring.json')
PAIR_TERMS=('fa_atr','fa_rep','fa_sol')


def disagreements(geometry,scores):
    result=[]
    for name,row in geometry['results'].items():
        if row['status']!='generated':
            continue
        prior=row['prior_selected_index']
        selected=scores['results'][name+'-means']['selection']['selected_index']
        if prior is not None and selected is not None and prior!=selected:
            result.append(dict(id=name,prior_index=prior,scored_index=selected))
    return result


def closure(values,expected):
    if not math.isfinite(expected) or any(not math.isfinite(v) for v in values):
        raise ValueError('nonfinite pair accounting')
    total=math.fsum(values)
    bound=1e-5+1e-7*math.fsum(abs(v) for v in values)
    return dict(pair_sum=total,recorded_total=expected,residual=total-expected,
                numerical_bound=bound,closed=abs(total-expected)<=bound)


def plan(output):
    geometry=json.loads(GEOMETRY.read_bytes())
    scores=json.loads((RUN_ROOT/'evidence.json').read_bytes())
    original=json.loads((PLAN_ROOT/'plan.json').read_bytes())
    cases=[]
    inputs={str(p):sha(p.read_bytes()) for p in (GEOMETRY,RUN_ROOT/'evidence.json',PLAN_ROOT/'plan.json',
        RUN_ROOT/'executed-plan.json',RUN_ROOT/'environment.json')}
    for row in disagreements(geometry,scores):
        identifier=row['id']+'-means'
        case=next(c for c in original['cases'] if c['id']==identifier)
        scored=scores['results'][identifier]
        cells=[]
        for role,index in [('prior',row['prior_index']),('selected',row['scored_index'])]:
            cell=next(c for c in case['cells'] if c['index']==index)
            saved=next(c for c in scored['cells'] if c['index']==index)
            if saved['score_status']!='scored' or not saved['import_audit']['accepted']:
                raise ValueError('saved disagreement not scored/import-valid')
            for p in (PLAN_ROOT/cell['pdb'],RUN_ROOT/identifier/cell['id']/'score.json',
                      RUN_ROOT/identifier/cell['id']/'pose-snapshot.json.gz'):
                inputs[str(p)]=sha(p.read_bytes())
            cells.append(dict(role=role,declaration=cell,saved=saved,
                              original_snapshot=str(RUN_ROOT/identifier/cell['id']/'pose-snapshot.json.gz')))
        cases.append(dict(id=row['id'],case=case,cells=cells))
    hashes=dict(original['source_hashes'])
    for name in ('scripts/research_mean_pair_energies.py','scripts/research_rosetta_fixed_scoring.py'):
        hashes[name]=sha(Path(name).read_bytes())
    for name,d in hashes.items():
        if sha(Path(name).read_bytes())!=d:
            raise ValueError('original scoring implementation changed')
    write(output,dict(schema='dvbfixer.mean-pair-accounting-plan.v1',cases=cases,inputs=inputs,
        source_hashes=hashes,plan_root=str(PLAN_ROOT),init_flags=original['init_flags'],
        expected_revision=original['expected_rosetta_revision'],score_function=original['score_function'],
        pair_terms=PAIR_TERMS,onebody_term='fa_dun',
        enrollment='ALL nonempty prior-versus-scored mean-pool disagreements, no native or score-margin filtering',
        replay='score_cell outputs and complete compressed pose snapshot must reproduce original bytes/values exactly before EnergyGraph interpretation',
        accounting='enumerate every residue pair edge; sum target-containing and other edges separately; fa_dun onebody per residue',
        numerical_policy='pair-cache accounting closed if residual <= 1e-5 + 1e-7*sum(abs(terms)); otherwise retain incomplete accounting, never adjust score or geometry',
        native_coordinates_read=False,weights_modified=False,production_enabled=False))


def all_atom_state(pose):
    info=pose.pdb_info()
    rows=[]
    for i in range(1,pose.total_residue()+1):
        residue=pose.residue(i)
        for a in range(1,residue.natoms()+1):
            xyz=residue.xyz(a)
            rows.append([info.chain(i),info.number(i),info.icode(i),residue.name(),residue.atom_name(a),xyz.x,xyz.y,xyz.z])
    return rows


def run(plan_path,output):
    from scripts.research_rosetta_fixed_scoring import score_cell

    raw=plan_path.read_bytes()
    record=json.loads(raw)
    for name,d in {**record['inputs'],**record['source_hashes']}.items():
        if sha(Path(name).read_bytes())!=d:
            raise ValueError('frozen pair-accounting input changed')
    import pyrosetta

    if record['expected_revision'] not in pyrosetta.version():
        raise ValueError('unexpected PyRosetta version')
    pyrosetta.init(record['init_flags'])
    sf=pyrosetta.create_score_function(record['score_function'])
    from pyrosetta.rosetta.core import scoring

    terms={name:getattr(scoring,name) for name in record['pair_terms']}
    output.mkdir(parents=True,exist_ok=False)
    (output/'executed-plan.json').write_bytes(raw)
    results={}
    for row in record['cases']:
        cells={}
        for cell in row['cells']:
            folder=output/row['id']/cell['role']
            folder.mkdir(parents=True)
            original_import=pyrosetta.pose_from_pdb
            captured=[]

            def capture(path,original_import=original_import,captured=captured):
                pose=original_import(path)
                captured.append(pose)
                return pose

            try:
                with patch.object(pyrosetta,'pose_from_pdb',capture):
                    repeated=score_cell(pyrosetta,sf,Path(record['plan_root']),folder,row['case'],cell['declaration'])
                if pyrosetta.pose_from_pdb is not original_import or len(captured)!=1:
                    raise ValueError('pose capture not restored or inventory differs')
                if any(cell['saved'].get(k)!=v for k,v in repeated.items()):
                    raise ValueError('exact repeated score metadata/values differ')
                if (folder/'pose-snapshot.json.gz').read_bytes()!=Path(cell['original_snapshot']).read_bytes():
                    raise ValueError('complete pose snapshot does not reproduce exactly')
                pose=captured[0]
                before=all_atom_state(pose)
                info=pose.pdb_info()
                residue_ids={i:[info.chain(i),info.number(i),info.icode(i),pose.residue(i).name3()] for i in range(1,pose.total_residue()+1)}
                matches=[i for i,key in residue_ids.items() if key==row['case']['target_key']]
                if len(matches)!=1:
                    raise ValueError('target pose identity not unique')
                target=matches[0]
                graph=pose.energies().energy_graph()
                pairs=[]
                for i in range(1,pose.total_residue()+1):
                    for j in range(i+1,pose.total_residue()+1):
                        edge=graph.find_energy_edge(i,j)
                        if edge is None:
                            continue
                        energy=edge.fill_energy_map()
                        pairs.append(dict(first=residue_ids[i],second=residue_ids[j],contains_target=target in (i,j),
                            weighted_terms={name:float(energy[t])*sf.get_weight(t) for name,t in terms.items()}))
                onebody=[]
                for i in range(1,pose.total_residue()+1):
                    emap=pose.energies().onebody_energies(i)
                    onebody.append(dict(residue=residue_ids[i],is_target=i==target,
                                        weighted_fa_dun=float(emap[scoring.fa_dun])*sf.get_weight(scoring.fa_dun)))
                accounting={name:closure([p['weighted_terms'][name] for p in pairs],
                    cell['saved']['score_terms']['ScoreType.'+name]['raw']*cell['saved']['score_terms']['ScoreType.'+name]['weight']) for name in terms}
                accounting['fa_dun']=closure([r['weighted_fa_dun'] for r in onebody],
                    cell['saved']['score_terms']['ScoreType.fa_dun']['raw']*cell['saved']['score_terms']['ScoreType.fa_dun']['weight'])
                if all_atom_state(pose)!=before:
                    raise ValueError('pair accounting moved an atom')
                evidence=dict(role=cell['role'],original_snapshot_sha256=sha(Path(cell['original_snapshot']).read_bytes()),
                    repeated_snapshot_sha256=sha((folder/'pose-snapshot.json.gz').read_bytes()),exact_score_repeat=True,
                    atom_coordinates_unchanged=True,residue_pairs=pairs,fa_dun_onebody=onebody,accounting=accounting)
                write(folder/'pair-accounting.json',evidence)
                cells[cell['role']]=dict(path=str(folder/'pair-accounting.json'),sha256=sha((folder/'pair-accounting.json').read_bytes()),
                    accounting_closed=all(r['closed'] for r in accounting.values()))
            except Exception as error:
                write(output/'execution-error.json',dict(case=row['id'],role=cell['role'],error_type=type(error).__name__,error=str(error)))
                raise
        results[row['id']]=cells
    write(output/'evidence.json',dict(schema='dvbfixer.mean-pair-accounting-run.v1',results=results,
        plan_sha256=sha(raw),native_coordinates_read=False,weights_modified=False))
    print(json.dumps(dict(cases=len(results),exact_repeated_states=2*len(results),
        accounting_closed=all(c['accounting_closed'] for r in results.values() for c in r.values())),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('plan')
    p.add_argument('--output',type=Path,required=True)
    p=sub.add_parser('run')
    p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.command=='plan':
        plan(args.output)
    else:
        run(args.plan,args.output)
