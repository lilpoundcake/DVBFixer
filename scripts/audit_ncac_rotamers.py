"""Independent candidate/score-snapshot audit, then native accuracy reporting.

Never generate or rescore. Verify rigid branch transport from saved parent
coordinates and score every feasible member exactly once before held-out readback.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import statistics
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scripts.audit_ncac_expansion import accuracy
from scripts.compare_branched_completion import measure
from scripts.research_rosetta_joint import sha
from scripts.research_rosetta_repack import identity, observed_atoms, write_json
from scripts.rosetta_fixed_contract import choose_scored_own, import_audit
from scripts.summarize_rosetta_fixed_scoring import source_records_preserved

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.pdb_component_reconstruction import _component_ref, _parse_atom_lines
from dvbfixer.protein_completion.branched_source import plan_branched_source
from dvbfixer.protein_completion.contracts import Addition


def check_rotation(parent, actual, angle):
    if set(parent)!=set(actual):
        raise ValueError('rotated atom set changed')
    for n in ('N','CA','C','O','CB'):
        if parent[n]!=actual[n]:
            raise ValueError('rotation moved fixed backbone/CB')
    cb,ca=np.asarray(parent['CB']),np.asarray(parent['CA'])
    u=(cb-ca)/np.linalg.norm(cb-ca)
    for n in set(parent)-{'N','CA','C','O','CB'}:
        v,w=np.asarray(parent[n])-cb,np.asarray(actual[n])-cb
        axial_v,axial_w=np.dot(v,u),np.dot(w,u)
        pv,pw=v-axial_v*u,w-axial_w*u
        if not math.isclose(axial_v,axial_w,abs_tol=1e-12) or not math.isclose(np.linalg.norm(pv),np.linalg.norm(pw),abs_tol=1e-12):
            raise ValueError('branch axial/radial geometry changed')
        measured=math.degrees(math.atan2(np.dot(u,np.cross(pv,pw)),np.dot(pv,pw)))
        if not math.isclose(measured,angle,abs_tol=1e-9):
            raise ValueError('relative branch angle differs')
    distal=list(set(parent)-{'N','CA','C','O','CB'})
    if not math.isclose(math.dist(parent[distal[0]],parent[distal[1]]),math.dist(actual[distal[0]],actual[distal[1]]),abs_tol=1e-12):
        raise ValueError('distal pair did not move rigidly together')


def audit(plan_root,run_root,output):
    plan_raw=(plan_root/'plan.json').read_bytes()
    plan=json.loads(plan_raw)
    protocol_raw=(plan_root/'protocol.json').read_bytes()
    protocol=json.loads(protocol_raw)
    geometry=json.loads((plan_root/'geometry.json').read_bytes())
    for name,digest in plan['artifacts'].items():
        if sha((plan_root/name).read_bytes())!=digest:
            raise ValueError('candidate artifact changed')
    for name,digest in {**protocol['inputs'],**plan['source_hashes']}.items():
        if sha(Path(name).read_bytes())!=digest:
            raise ValueError('frozen input/implementation changed')
    if geometry['protocol_sha256']!=sha(protocol_raw) or geometry['native_coordinates_read']:
        raise ValueError('geometry/protocol binding differs')
    raw=(run_root/'evidence.json').read_bytes()
    run=json.loads(raw)
    if run['plan_sha256']!=sha(plan_raw) or (run_root/'executed-plan.json').read_bytes()!=plan_raw:
        raise ValueError('executed score plan differs')
    declarations={c['id']:c for c in protocol['cases']}
    if set(run['results'])!=set(declarations) or set(geometry['cases'])!=set(declarations) or {c['id'] for c in plan['cases']}!=set(declarations):
        raise ValueError('case inventory differs')
    parent_root=Path(protocol['run_path'])
    parent_run=json.loads((parent_root/'evidence.json').read_bytes())
    ledger,replay={},{}
    snapshots=0
    for case in plan['cases']:
        name=case['id']
        original=declarations[name]
        source=(plan_root/name/'source.pdb').read_bytes()
        if sha(source)!=original['input_sha256']:
            raise ValueError('original masked source changed')
        request=plan_branched_source(source,ComponentInstanceRef(**original['target'])).request
        trace_raw=(parent_root/(name+'-ncac')/'evidence.json').read_bytes()
        if sha(trace_raw)!=parent_run['results'][name+'-ncac']['trace_sha256']:
            raise ValueError('saved parent NCAC trace changed')
        parent=[e for e in json.loads(trace_raw)['events'] if e['stage']=='returned-seed'][-1]['target_points_angstrom']
        proposals=geometry['cases'][name]['candidates']
        if [c['relative_angle_degrees'] for c in proposals]!=[0.,120.,-120.]:
            raise ValueError('complete three-angle pool required')
        audited=[]
        feasible=[]
        for i,c in enumerate(proposals):
            if c['index']!=i or c['id']!=f'own-{i}' or c['arm']!='own':
                raise ValueError('candidate identity/order differs')
            check_rotation(parent,c['points'],c['relative_angle_degrees'])
            added=tuple(Addition(ExactAtomRef(request.residue,n),n[0],tuple(c['points'][n])) for n in request.missing_atoms)
            measured,pdb=measure(source,request,added,anchored=False)
            if json.loads(json.dumps(asdict(measured)))!=c['common_audit'] or pdb!=(plan_root/c['pdb']).read_bytes():
                raise ValueError('independent geometry/materialization differs')
            if c['common_accepted']!=measured.accepted or c['old_source_accepted']!=measured.accepted:
                raise ValueError('candidate feasibility differs')
            if not source_records_preserved(source,pdb,case['generated']):
                raise ValueError('candidate changed source records')
            if measured.accepted:
                feasible.append(c['id'])
            audited.append(dict(index=i,relative_angle_degrees=c['relative_angle_degrees'],common_accepted=measured.accepted))
        if [c['id'] for c in case['cells']]!=feasible:
            raise ValueError('runner plan omitted feasible or included rejected candidate')
        result=run['results'][name]
        scored=result['cells']
        if [c['id'] for c in scored]!=feasible:
            raise ValueError('scored pool omitted feasible candidate')
        extras,virtuals,context_h=[],[],[]
        for c,declaration in zip(scored,case['cells'],strict=True):
            if any(c[k]!=v for k,v in declaration.items()):
                raise ValueError('scored declaration changed')
            folder=run_root/name/c['id']
            if json.loads((folder/'score.json').read_bytes())!=c:
                raise ValueError('score record changed')
            snapshot_raw=(folder/'pose-snapshot.json.gz').read_bytes()
            if sha(snapshot_raw)!=c['snapshot_sha256']:
                raise ValueError('score snapshot changed')
            s=json.loads(gzip.decompress(snapshot_raw))
            expected=observed_atoms((plan_root/c['pdb']).read_bytes())
            outgoing={tuple(a if a[-1]=='C' else b)[:4] for a,b in case['expected_external'] if {a[-1],b[-1]}=={'C','N'}}
            terminal={(r['chain'],r['resid'],r['icode']) for r in s['residues'] if 'CtermProteinFull' in r['name']}
            expected_ids={identity(a) for a in expected}
            allowed=[identity(a) for a in s['heavy_atoms'] if a['atom']=='OXT' and identity(a)[:3] in terminal
                     and identity(a)[:4] not in outgoing and identity(a) not in expected_ids]
            ia=import_audit(expected,s['heavy_atoms'],case['expected_external'],s['external'],case['target_internal'],s['target_internal'],allowed)
            if json.loads(json.dumps(ia))!=c['import_audit'] or not ia['accepted'] or c['score_status']!='scored':
                raise ValueError('independent fixed-heavy/topology import blocker')
            if any(not v['virtual_by_rosetta_type'] for v in s['virtual_atoms']) or {identity(v) for v in s['virtual_atoms']} & {identity(a) for a in s['heavy_atoms']}:
                raise ValueError('physical/virtual inventory differs')
            if not c['virtual_coordinates_unchanged_by_scoring']:
                raise ValueError('scoring moved virtual coordinates')
            if not math.isclose(sum(t['raw']*t['weight'] for t in c['score_terms'].values()),c['total_score'],rel_tol=1e-10,abs_tol=1e-8):
                raise ValueError('score term total differs')
            extras.append(sorted((identity(a),a['xyz']) for a in s['heavy_atoms'] if identity(a) in allowed))
            virtuals.append(sorted((identity(a),a['xyz']) for a in s['virtual_atoms']))
            context_h.append(sorted(((h['chain'],h['resid'],h['icode'],h['atom']),h['xyz']) for h in s['hydrogens'] if (h['chain'],h['resid'],h['icode'])!=tuple(case['target_key'][:3])))
            snapshots+=1
        context_equal=all(not x or all(v==x[0] for v in x) for x in (extras,virtuals,context_h))
        selection=choose_scored_own(scored)
        if selection!=result['selection'] or selection!=json.loads((run_root/name/'selection-before-native.json').read_bytes()):
            raise ValueError('complete feasible selection replay differs')
        ledger[name]=dict(candidates=audited,selection=selection,context_coordinates_identical=context_equal,
                          feasible_count=len(feasible),ncac_previously_executed=original['ncac_previously_executed'])
        replay[name]=(request,source,proposals)
    selection_path=output.with_name(output.stem+'-selection.json')
    if output.exists() or selection_path.exists():
        raise FileExistsError('fresh audit and selection destinations required')
    write_json(selection_path,dict(cases=ledger,native_coordinates_read=False,score_plan_sha256=sha(plan_raw),
                                  run_sha256=sha(raw),validated_snapshots=snapshots,audit_script_sha256=sha(Path(__file__).read_bytes())))
    # All case selections and independent audits are committed before native PDBs.
    expansion=json.loads(Path(protocol['expansion_plan']).read_bytes())
    mask_root=Path(expansion['mask_root'])
    masks=json.loads((mask_root/'plan.json').read_bytes())
    previous_raw=Path(expansion['prior_audit']).read_bytes()
    if sha(previous_raw)!=expansion['inputs'][expansion['prior_audit']]:
        raise ValueError('prior own-scored evidence changed')
    previous=json.loads(previous_raw)
    parent_path=next(p for p in protocol['inputs'] if Path(p).name=='ncac-whole-mask-evidence.json')
    parent_evidence=json.loads(Path(parent_path).read_bytes())
    reports={}
    for name,(request,source,proposals) in replay.items():
        original=declarations[name]
        ref=(mask_root/original['reference']).read_bytes()
        if sha(ref)!=masks['artifacts'][original['reference']]:
            raise ValueError('native reporting reference changed')
        truth={a.name:a.position for a in _parse_atom_lines(ref.decode().splitlines(keepends=True)) if _component_ref(a)==request.residue}
        candidates=[]
        for c in proposals:
            added=tuple(Addition(ExactAtomRef(request.residue,n),n[0],tuple(c['points'][n])) for n in request.missing_atoms)
            candidates.append(dict(index=c['index'],common_accepted=c['common_accepted'],
                accuracy=accuracy(request,added,(plan_root/c['pdb']).read_bytes(),truth)))
        selected=ledger[name]['selection']['selected_index']
        own=previous['cases'][name]['own_scored']
        if not previous['cases'][name]['own_scored']['common_accepted']:
            raise ValueError('prior own comparator not common-feasible')
        zero_rmsd=candidates[0]['accuracy']['rounded']['exact_name_rmsd_angstrom']
        if zero_rmsd!=parent_evidence['cases'][name]['methods']['ncac']['accuracy']['rounded']['exact_name_rmsd_angstrom']:
            raise ValueError('zero-angle accuracy does not repeat unchanged parent')
        reports[name]=dict(ledger[name],candidates=candidates,selected_index=selected,
                           zero_rmsd=zero_rmsd,
                           prior_own_scored_rmsd=own['rounded_exact_name_rmsd_angstrom'],
                           selected_rmsd=None if selected is None else candidates[selected]['accuracy']['rounded']['exact_name_rmsd_angstrom'])
    aggregate={}
    for label,subset in [('all',list(reports.values())),('seen_four',[r for r in reports.values() if r['ncac_previously_executed']]),
                         ('new_seven',[r for r in reports.values() if not r['ncac_previously_executed']])]:
        matched=[r for r in subset if r['selected_index'] is not None and r['context_coordinates_identical']]
        aggregate[label]=dict(cases=len(subset),matched=len(matched),changed=sum(r['selected_index']!=0 for r in matched),
                improved=sum(r['selected_rmsd']<r['zero_rmsd']-1e-12 for r in matched),
                worsened=sum(r['selected_rmsd']>r['zero_rmsd']+1e-12 for r in matched),
                prior_own_scored_median=statistics.median(r['prior_own_scored_rmsd'] for r in matched) if matched else None,
                zero_median=statistics.median(r['zero_rmsd'] for r in matched) if matched else None,
                selected_median=statistics.median(r['selected_rmsd'] for r in matched) if matched else None)
    write_json(output,dict(schema='dvbfixer.ncac-relative-rotamer-audit.v1',cases=reports,aggregate=aggregate,
        selection_ledger_name=selection_path.name,selection_ledger_sha256=sha(selection_path.read_bytes()),
        validated_snapshots=snapshots,score_plan_sha256=sha(plan_raw),run_sha256=sha(raw),
        audit_script_sha256=sha(Path(__file__).read_bytes()),production_changed=False,scientific_qualification=False))
    print(json.dumps(aggregate,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    audit(args.plan,args.run,args.output)
