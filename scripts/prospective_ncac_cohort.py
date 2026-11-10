"""New-input enrollment for the frozen whole-mask NCAC relative-rotamer arm.

No substitutions, sidechain repair of other observations, atom filtering, angle
tuning or native-derived selection. Preserve refusals, shortfalls and failures.
"""
from __future__ import annotations

import argparse
import json
from contextlib import redirect_stderr, redirect_stdout
from dataclasses import asdict
from pathlib import Path
from unittest.mock import patch

from scripts import research_template_causality as tracer
from scripts.branched_template_baseline import preflight, template_proposal
from scripts.compare_branched_completion import measure
from scripts.prospective_independent_cohort import project_first, selected_targets
from scripts.research_ncac_rotamers import ANGLES, rotate_branches
from scripts.research_rosetta_joint import sha, verify
from scripts.research_rosetta_repack import identity, write_json
from scripts.research_template_anchors import anchor_data, make_anchor_overlay

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.pdb_component_reconstruction import _parse_atom_lines
from dvbfixer.protein_completion.branched_source import plan_branched_source
from dvbfixer.protein_completion.contracts import Addition
from dvbfixer.protein_completion.service import complete_sidechain
from dvbfixer.protein_completion.validation import sidechain_bonds

SYSTEMS={'1PGB':'protein G B1 domain','1BPI':'bovine pancreatic trypsin inhibitor',
         '1CSP':'Bacillus subtilis cold-shock protein'}
PARENT=Path('tmp/pdbfixer-research/ncac-rotamers-score-plan-01/plan.json')


def whole_mask(name):
    if name=='VAL':
        return ['CB','CG1','CG2']
    if name=='THR':
        return ['CB','OG1','CG2']
    raise ValueError('only VAL/THR enrolled')


def enroll(output):
    parent_raw=PARENT.read_bytes()
    parent=json.loads(parent_raw)
    hashes=dict(parent['source_hashes'])
    for p in ('scripts/prospective_ncac_cohort.py','scripts/prospective_independent_cohort.py'):
        hashes[p]=sha(Path(p).read_bytes())
    for p,d in parent['source_hashes'].items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('frozen algorithm changed')
    output.parent.mkdir(parents=True,exist_ok=True)
    write_json(output,dict(schema='dvbfixer.prospective-ncac-enrollment.v1',
        systems=[dict(id=n,description=d,url=f'https://files.rcsb.org/download/{n}.pdb') for n,d in SYSTEMS.items()],
        selection='first coordinate MODEL and first encountered ATOM chain; first two VAL and first two THR identities in file order, whole-sidechain masks only; retain shortfalls/refusals',
        observations='reuse unchanged project_first; retain H/altloc/occupancy and touching link endpoints; do not sanitize into admission',
        parent_score_plan=str(PARENT),parent_score_plan_sha256=sha(parent_raw),source_hashes=hashes,
        environment=preflight(),angles_degrees=ANGLES,max_masks=12,replacement_structures_allowed=False,
        arms='stock shared seeded PDBFixer Reference; unchanged own three-staggered proposal pool; NCAC-normal plus fixed 0/+120/-120 distal rotations',
        geometry='unchanged common gate for every proposal; all failures retained; no native tuning',
        scoring='unchanged fixed-heavy ref2015 runner; ALL and ONLY common-feasible own/NCAC candidates, separate complete pools; stable index tie; no native control in runner',
        native_policy='read native reference coordinates for accuracy only after all score choices commit; geometry/proposals use masked inputs only',
        qualification_scope='prospective input/algorithm enrollment, deliberately chosen systems; not randomized/blinded or independent from Rosetta training/templates',
        production_enabled=False))


def freeze(enrollment_path,inputs,output):
    raw=enrollment_path.read_bytes()
    enrollment=json.loads(raw)
    for p,d in enrollment['source_hashes'].items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('enrolled algorithm changed')
    output.mkdir(parents=True,exist_ok=False)
    cases,systems,artifacts=[],{},{}

    def save(name,data):
        (output/name).write_bytes(data)
        artifacts[name]=sha(data)

    for system in enrollment['systems']:
        n=system['id']
        original=(inputs/(n+'.pdb')).read_bytes()
        save(n+'-original.pdb',original)
        try:
            reference,chain=project_first(original)
        except ValueError as error:
            systems[n]=dict(status='projection_refused',reason=str(error))
            continue
        save(n+'-reference.pdb',reference)
        targets,counts=selected_targets(reference)
        systems[n]=dict(status='projected',chain=chain,selected_counts=counts,count_shortfall=counts!={'VAL':2,'THR':2})
        for name,target in targets:
            mask=whole_mask(name)
            identifier=f'{n}-{name}{target.sequence_number}{target.insertion_code}-'+ '-'.join(mask)
            source=b''.join(line for line in reference.splitlines(keepends=True) if not (
                line.startswith(b'ATOM  ') and (line[21:22].decode(),line[22:26].strip().decode(),line[26:27].strip().decode())
                ==(target.chain_id,target.sequence_number,target.insertion_code) and line[12:16].strip().decode() in mask))
            save(identifier+'.pdb',source)
            cases.append(dict(id=identifier,fixture=n,target=asdict(target),residue_name=name,mask=mask,
                              source=identifier+'.pdb',reference=n+'-reference.pdb'))
    (output/'enrollment.json').write_bytes(raw)
    write_json(output/'plan.json',dict(schema='dvbfixer.prospective-ncac-masks.v1',enrollment=enrollment,
        enrollment_sha256=sha(raw),cases=cases,systems=systems,artifacts=artifacts,source_hashes=enrollment['source_hashes']))


def score_case(name,arm,source,request,rows,output,artifacts):
    """Transport canonical/source-derived target graph, never duplicate scoring."""
    identifier=name+'-'+arm
    payload=rows[0][2]
    lookup={(a.chain_id,a.sequence_number,a.insertion_code,a.name):
            identity(dict(chain=a.chain_id,resid=int(a.sequence_number),icode=a.insertion_code or ' ',
                          resname=a.residue_name,atom=a.name))
            for a in _parse_atom_lines(payload.decode().splitlines(keepends=True))}

    def key(atom):
        c=atom.component
        return lookup[(c.chain_id,c.sequence_number,c.insertion_code,atom.name)]

    cells=[]
    for i,audit,pdb in rows:
        filename=f'{identifier}/own-{i}.pdb'
        path=output/filename
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_bytes(pdb)
        artifacts[filename]=sha(pdb)
        if audit.accepted:
            cells.append(dict(id=f'own-{i}',arm='own',index=i,pdb=filename,
                old_source_accepted=True,common_accepted=True,common_audit=asdict(audit)))
    filename=f'{identifier}/source.pdb'
    (output/filename).write_bytes(source)
    artifacts[filename]=sha(source)
    target=request.residue
    return dict(id=identifier,status='admitted',cells=cells,
        target_key=(target.chain_id,int(target.sequence_number),target.insertion_code or ' ',request.residue_name),
        generated=[key(ExactAtomRef(target,n)) for n in request.missing_atoms],
        expected_external=[(key(a),key(b)) for a,b in request.links],
        target_internal=[(key(ExactAtomRef(target,a)),key(ExactAtomRef(target,b))) for a,b in sidechain_bonds(request.residue_name)
                         if (target.chain_id,target.sequence_number,target.insertion_code,a) in lookup
                         and (target.chain_id,target.sequence_number,target.insertion_code,b) in lookup])


def run(plan_root,output):
    raw=(plan_root/'plan.json').read_bytes()
    plan=json.loads(raw)
    verify(plan,plan_root)
    enrollment=plan['enrollment']
    if preflight()!=enrollment['environment']:
        raise ValueError('scientific environment differs from enrollment')
    parent_raw=Path(enrollment['parent_score_plan']).read_bytes()
    if sha(parent_raw)!=enrollment['parent_score_plan_sha256']:
        raise ValueError('parent score settings changed')
    parent=json.loads(parent_raw)
    output.mkdir(parents=True,exist_ok=False)
    (output/'executed-plan.json').write_bytes(raw)
    score_root=output/'score-plan'
    score_root.mkdir()
    results,cases,artifacts={},[],{}
    for original in plan['cases']:
        name=original['id']
        folder=output/name
        folder.mkdir()
        source=(plan_root/original['source']).read_bytes()
        (folder/'input.pdb').write_bytes(source)
        target=ComponentInstanceRef(**original['target'])
        try:
            request=plan_branched_source(source,target).request
            if set(request.missing_atoms)!=set(original['mask']):
                raise ValueError('actual missing atoms differ from enrolled mask')
        except ValueError as error:
            results[name]=dict(status='shared_admission_refused',reason=str(error))
            write_json(folder/'result.json',results[name])
            continue
        own=complete_sidechain(request)
        write_json(folder/'own-attempts.json',asdict(own))
        own_rows=[]
        for i,attempt in enumerate(own.attempts):
            audit,pdb=measure(source,request,attempt.candidate.additions,anchored=False)
            own_rows.append((i,audit,pdb))
        if not own_rows:
            raise ValueError('expected proposal inventory, even for empty feasible pool')
        cases.append(score_case(name,'baseline',source,request,own_rows,score_root,artifacts))
        row=dict(status='admitted',ncac_status='not_executed')
        try:
            with (folder/'stock.log').open('x') as log,redirect_stdout(log),redirect_stderr(log):
                stock=template_proposal(source,request,retain_snapshot=lambda s,folder=folder:write_json(folder/f'stock-seed-{s.seed}.json',asdict(s)))
            write_json(folder/'stock-proposal.json',asdict(stock))
            audit,pdb=measure(source,request,stock.additions,anchored=False)
            (folder/'stock-proposal.pdb').write_bytes(pdb)
            row['stock']=dict(status='evaluated',common_audit=asdict(audit))
        except Exception as error:
            row['stock']=dict(status='execution_blocker',error_type=type(error).__name__,error=str(error))
        declaration=dict(original,input_path=str(folder/'input.pdb'))
        try:
            data=anchor_data(declaration,enrollment['environment'])
            calls=[]
            hook=make_anchor_overlay(data['names'],data['reference_points_nm'],data['template_points_nm'],
                                     'ncac',tracer.proper_overlay,calls)
            original_fit=tracer.proper_overlay
            with patch.object(tracer,'proper_overlay',hook):
                trace=tracer.run_arm(declaration,dict(overlay='proper',relaxation='normal'),folder/'ncac-trace')
            if tracer.proper_overlay is not original_fit or len(calls)!=len(trace['overlays']):
                raise ValueError('NCAC hook restoration/transport differs')
            write_json(folder/'anchor-policy.json',dict(calls=calls,anchor_data=data,research_hook_restored=True))
            returned=[e for e in trace['events'] if e['stage']=='returned-seed'][-1]
            points=returned['target_points_angstrom']
            rows,proposals=[],[]
            for i,angle in enumerate(ANGLES):
                proposed=rotate_branches(points,angle)
                added=tuple(Addition(ExactAtomRef(target,n),n[0],tuple(proposed[n])) for n in request.missing_atoms)
                audit,pdb=measure(source,request,added,anchored=False)
                rows.append((i,audit,pdb))
                proposals.append(dict(index=i,relative_angle_degrees=angle,points=proposed,common_audit=asdict(audit)))
            write_json(folder/'ncac-proposals.json',dict(candidates=proposals))
            cases.append(score_case(name,'ncac',source,request,rows,score_root,artifacts))
            row.update(ncac_status='evaluated',trace_sha256=sha((folder/'ncac-trace/evidence.json').read_bytes()))
        except Exception as error:
            row.update(ncac_status='execution_or_protocol_blocker',error_type=type(error).__name__,error=str(error))
        results[name]=row
        write_json(folder/'result.json',row)
    write_json(output/'evidence.json',dict(schema='dvbfixer.prospective-ncac-run.v1',results=results,
         plan_sha256=sha(raw),plan_root=str(plan_root),source_hashes=plan['source_hashes'],native_coordinates_read=False))
    settings={k:parent[k] for k in ('expected_rosetta_revision','score_function','seed','radii','init_flags',
        'hydrogen_protocol','extra_heavy_policy','virtual_atom_policy','topology_scope','selection_policy','native_policy')}
    write_json(score_root/'plan.json',dict(schema='dvbfixer.rosetta-fixed-scoring-plan.v2',cases=cases,
        artifacts=artifacts,source_hashes=plan['source_hashes'],**settings,prospectively_enrolled=True,public_capability=False))
    print(json.dumps(dict(planned=len(results),admitted=sum(r['status']=='admitted' for r in results.values()),
        ncac_evaluated=sum(r.get('ncac_status')=='evaluated' for r in results.values()),
        score_views=len(cases),feasible_score_inputs=sum(len(c['cells']) for c in cases)),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('enroll')
    p.add_argument('--output',type=Path,required=True)
    p=sub.add_parser('freeze')
    p.add_argument('--enrollment',type=Path,required=True)
    p.add_argument('--inputs',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p=sub.add_parser('run')
    p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.command=='enroll':
        enroll(args.output)
    elif args.command=='freeze':
        freeze(args.enrollment,args.inputs,args.output)
    else:
        run(args.plan,args.output)
