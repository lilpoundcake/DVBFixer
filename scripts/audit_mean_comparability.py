"""Saved-only cross-arm context and generated-VAL-symmetry readback.

No pose import, new scores, coordinates, library query or selector. Enumerate all
mean views and corresponding parent-relative views, not just score disagreements.
"""
from __future__ import annotations

import argparse
import gzip
import json
import statistics
from pathlib import Path

from scripts.research_dunbrack_prior import sha, write
from scripts.research_rosetta_repack import identity
from scripts.summarize_rosetta_fixed_scoring import source_records_preserved

MEANS=Path('tmp/pdbfixer-research/dunbrack-means-generation-01/score-plan')
MEAN_RUN=Path('tmp/pdbfixer-research/dunbrack-means-score-run-01')
PARENTS={
    'old_eleven':(Path('tmp/pdbfixer-research/ncac-rotamers-score-plan-01'),Path('tmp/pdbfixer-research/ncac-rotamers-score-run-01')),
    'new_inputs_eleven':(Path('tmp/pdbfixer-research/ncac-new-system-run-01/score-plan'),Path('tmp/pdbfixer-research/ncac-new-system-score-run-01')),
}


def context(snapshot,generated,target):
    """Exact scorer-owned/observed context inventories, excluding generated target."""
    generated={tuple(k) for k in generated}
    physical=snapshot['heavy_atoms']
    if len({identity(a) for a in physical})!=len(physical):
        raise ValueError('duplicate physical atom identity')
    virtual=snapshot['virtual_atoms']
    if any(not a['virtual_by_rosetta_type'] for a in virtual):
        raise ValueError('untyped virtual atom')
    return dict(heavy=sorted((identity(a),a['xyz']) for a in physical if identity(a) not in generated),
        virtual=sorted((identity(a),a['xyz']) for a in virtual),
        hydrogens=sorted(((a['chain'],a['resid'],a['icode'],a['atom']),a['xyz']) for a in snapshot['hydrogens']
                         if (a['chain'],a['resid'],a['icode'])!=tuple(target[:3])),
        residues=snapshot['residues'],external=snapshot['external'],target_internal=snapshot['target_internal'])


def collect(plan_root,run_root,case_id,inputs):
    for p in (plan_root/'plan.json',run_root/'executed-plan.json',run_root/'evidence.json'):
        inputs[str(p)]=sha(p.read_bytes())
    plan=json.loads((plan_root/'plan.json').read_bytes())
    run=json.loads((run_root/'evidence.json').read_bytes())
    if run['plan_sha256']!=sha((plan_root/'plan.json').read_bytes()) or (run_root/'executed-plan.json').read_bytes()!=(plan_root/'plan.json').read_bytes():
        raise ValueError('score plan execution binding differs')
    case=next(c for c in plan['cases'] if c['id']==case_id)
    row=run['results'][case_id]
    if [c['id'] for c in case['cells']]!=[c['id'] for c in row['cells']]:
        raise ValueError('complete feasible pool differs')
    source_path=plan_root/case_id/'source.pdb'
    inputs[str(source_path)]=sha(source_path.read_bytes())
    cells=[]
    for cell in row['cells']:
        if cell['score_status']!='scored' or not cell['import_audit']['accepted']:
            raise ValueError('parent-certified scored pool contains blocker')
        pdb_path=plan_root/cell['pdb']
        if sha(pdb_path.read_bytes())!=plan['artifacts'][cell['pdb']]:
            raise ValueError('candidate PDB changed')
        snapshot=run_root/case_id/cell['id']/'pose-snapshot.json.gz'
        score=run_root/case_id/cell['id']/'score.json'
        if sha(snapshot.read_bytes())!=cell['snapshot_sha256'] or json.loads(score.read_bytes())!=cell:
            raise ValueError('retained scored state changed')
        for p in (pdb_path,snapshot,score):
            inputs[str(p)]=sha(p.read_bytes())
        cells.append(dict(id=cell['id'],pdb=str(pdb_path),snapshot=str(snapshot),score=str(score)))
    return dict(case_id=case_id,source=str(source_path),generated=case['generated'],target_key=case['target_key'],cells=cells)


def plan(output):
    mean_protocol=Path('tmp/pdbfixer-research/dunbrack-means-plan-01.json')
    original=json.loads(mean_protocol.read_bytes())
    inputs={str(mean_protocol):sha(mean_protocol.read_bytes())}
    cases=[]
    for c in original['cases']:
        if c['status']!='generation_ready':
            cases.append(dict(id=c['id'],group=c['group'],status='inherited_nonexecution'))
            continue
        p_root,r_root=PARENTS[c['group']]
        parent_id=c['id'] if c['group']=='old_eleven' else c['id']+'-ncac'
        cases.append(dict(id=c['id'],group=c['group'],status='readback_ready',
            mean=collect(MEANS,MEAN_RUN,c['id']+'-means',inputs),
            parent=collect(p_root,r_root,parent_id,inputs)))
    for p in ('docs/research/dunbrack-means-evidence.json','docs/research/ncac-relative-rotamers-evidence.json',
              'docs/research/ncac-new-system-evidence.json'):
        inputs[p]=sha(Path(p).read_bytes())
    hashes={p:sha(Path(p).read_bytes()) for p in ('scripts/audit_mean_comparability.py',
        'scripts/research_rosetta_repack.py','scripts/summarize_rosetta_fixed_scoring.py')}
    write(output,dict(schema='dvbfixer.mean-comparability-plan.v1',cases=cases,inputs=inputs,source_hashes=hashes,
        enrollment='ALL retained mean pools and corresponding parent-relative feasible pools; seven inherited nonexecutions retained',
        scope='saved retrospective readback, no new scores or candidate selection',production_enabled=False))


def audit(plan_path,output):
    raw=plan_path.read_bytes()
    record=json.loads(raw)
    for p,d in {**record['inputs'],**record['source_hashes']}.items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('bound saved state changed')
    reports={}
    snapshots=0
    for c in record['cases']:
        if c['status']!='readback_ready':
            reports[c['id']]=dict(status=c['status'],group=c['group'])
            continue
        mean,parent=c['mean'],c['parent']
        if Path(mean['source']).read_bytes()!=Path(parent['source']).read_bytes() or mean['generated']!=parent['generated'] or mean['target_key']!=parent['target_key']:
            raise ValueError('arms differ in source bytes/ownership/target')
        contexts=[]
        for view in (mean,parent):
            for cell in view['cells']:
                if not source_records_preserved(Path(view['source']).read_bytes(),Path(cell['pdb']).read_bytes(),view['generated']):
                    raise ValueError('original source records changed')
                s=json.loads(gzip.decompress(Path(cell['snapshot']).read_bytes()))
                contexts.append(context(s,view['generated'],view['target_key']))
                snapshots+=1
        equal={k:all(v[k]==contexts[0][k] for v in contexts) for k in contexts[0]}
        reports[c['id']]=dict(status='audited',group=c['group'],mean_snapshots=len(mean['cells']),parent_snapshots=len(parent['cells']),
                            context_equal=equal,all_context_equal=all(equal.values()))
    # Accuracy files are reporting only and cannot alter pool membership/context.
    means=json.loads(Path('docs/research/dunbrack-means-evidence.json').read_bytes())
    aggregate={}
    for group in ('old_eleven','new_inputs_eleven'):
        names=[n for n,c in reports.items() if c['group']==group and c['status']=='audited']
        metric={}
        for arm in ('prior_accuracy','scored_accuracy'):
            metric[arm]={k:statistics.median(means['cases'][n][arm]['rounded'][k] for n in names) for k in
                         ('exact_name_rmsd_angstrom','generated_only_val_symmetry_rmsd_angstrom')}
        aggregate[group]=dict(planned=sum(c['group']==group for c in reports.values()),audited=len(names),
            cross_arm_context_equal=sum(reports[n]['all_context_equal'] for n in names),mean_accuracy_medians=metric)
    write(output,dict(schema='dvbfixer.mean-comparability-audit.v1',cases=reports,aggregate=aggregate,
        snapshots=snapshots,plan_sha256=sha(raw),new_scores_evaluated=False,new_selectors=False,
        weights_modified=False,table_data_included=False,production_changed=False))
    print(json.dumps(dict(snapshots=snapshots,aggregate=aggregate),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('plan')
    p.add_argument('--output',type=Path,required=True)
    p=sub.add_parser('audit')
    p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.command=='plan':
        plan(args.output)
    else:
        audit(args.plan,args.output)
