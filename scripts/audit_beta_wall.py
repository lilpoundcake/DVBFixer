"""Saved-only independent stereo/geometry/ownership audit for private beta stress.

Abstract actual-frame heavy-geometry views only; never restore original positions
or claim H geometry/general minimizer convergence/native reconstruction accuracy.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
from dataclasses import asdict
from pathlib import Path

from scripts.audit_legacy_start_controls import coordinate_delta, load_atoms
from scripts.audit_legacy_stereo_trace import determinant, state_report
from scripts.research_dunbrack_prior import sha, write
from scripts.research_legacy_comparator import verify

from dvbfixer.diagnose.geometry import signed_dihedral_degrees
from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.protein_completion.contracts import (
    Addition,
    Candidate,
    CompletionRequest,
    Observation,
)
from dvbfixer.protein_completion.validation import _audit_candidate


def refs(rows):
    return [ExactAtomRef(ComponentInstanceRef(1,a['key'][0],a['key'][1],a['key'][2]),a['key'][-1]) for a in rows]


def geometry(rows,bonds,target,rounded=False):
    material=[dict(a,xyz=[float(f'{x:.3f}') for x in a['xyz']]) if rounded else a for a in rows]
    identities=refs(material)
    residue=ComponentInstanceRef(**target)
    missing=('CB','OG1','CG2')
    observations=[]
    additions=[]
    for a,ref in zip(material,identities,strict=True):
        if a['element'] in {'H','D',None}:
            continue
        if ref.component==residue and ref.name in missing:
            additions.append(Addition(ref,a['element'],tuple(a['xyz'])))
        else:
            observations.append(Observation(ref,a['element'],tuple(a['xyz'])))
    heavy={o.atom for o in observations}|{a.atom for a in additions}
    links=tuple((identities[a],identities[b]) for a,b in bonds if identities[a] in heavy and identities[b] in heavy
        and identities[a].component!=identities[b].component)
    request=CompletionRequest(sha(json.dumps(material,sort_keys=True,allow_nan=False).encode()),residue,'THR',
        tuple(observations),missing,links)
    points={o.atom.name:o.xyz for o in observations if o.atom.component==residue}
    points.update({a.atom.name:a.xyz for a in additions})
    chi=signed_dihedral_degrees(*(points[n] for n in ('N','CA','CB','OG1')))
    raw=_audit_candidate(request,Candidate(request.digest,tuple(additions),chi),measured_chi=True,pdb_roundtrip=rounded)
    stereo=state_report(material,target)
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
    target_cip=stereo['target']['target_cip']
    accepted=raw.accepted and not global_failures and target_cip==dict(CA='S',CB='R')
    return dict(accepted=accepted,target_geometry=asdict(raw),stereo=stereo,
        global_stereo_failures=global_failures,rounded_all_coordinates=rounded,
        scope='actual-frame abstract heavy-only view, no original-observation or H-geometry guarantee')


def force_audit(source,before_path,after_path,arm,protocol):
    import openmm

    original=openmm.XmlSerializer.deserialize(source)
    before=openmm.XmlSerializer.deserialize(before_path.read_text())
    after=openmm.XmlSerializer.deserialize(after_path.read_text())
    count=original.getNumForces()
    extra=int(arm in ('wall','drive_wall'))+int(arm in ('drive','drive_wall'))
    if before.getNumForces()!=count+extra or after.getNumForces()!=count+extra:
        raise ValueError('private beta-force inventory differs')
    before_base=openmm.XmlSerializer.deserialize(openmm.XmlSerializer.serialize(before))
    for i in reversed(range(count,before_base.getNumForces())):
        before_base.removeForce(i)
    if openmm.XmlSerializer.serialize(before_base)!=openmm.XmlSerializer.serialize(original):
        raise ValueError('original full physical system changed before minimization')
    for i in range(count):
        force=original.getForce(i)
        if openmm.XmlSerializer.serialize(force)!=openmm.XmlSerializer.serialize(before.getForce(i)):
            raise ValueError('original physical force modified in phase1')
        if isinstance(force,openmm.CustomExternalForce):
            if force.getEnergyFunction()!='0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)':
                raise ValueError('unexpected original external force')
            for j in range(force.getNumParticles()):
                index,parameters=force.getParticleParameters(j)
                force.setParticleParameters(j,index,[parameters[0]/10,*parameters[1:]])
        if openmm.XmlSerializer.serialize(force)!=openmm.XmlSerializer.serialize(after.getForce(i)):
            raise ValueError('phase2 changed more than original restraint strengths')
    for i in range(count,before.getNumForces()):
        force=before.getForce(i)
        if openmm.XmlSerializer.serialize(force)!=openmm.XmlSerializer.serialize(after.getForce(i)):
            raise ValueError('beta stress force changed between phases')
        params={force.getGlobalParameterName(j):force.getGlobalParameterDefaultValue(j) for j in range(force.getNumGlobalParameters())}
        if force.getName()=='PRIVATE_BETA_WALL':
            expected=dict(bscale=protocol['scale_nm3'],bwall_k=protocol['wall_k_kj_mol'],
                          bfloor=protocol['wall_floor_a3']/1000/protocol['scale_nm3'])
            if params!=expected or force.getForceGroup()!=30:
                raise ValueError('beta wall constants differ from frozen protocol')
        elif force.getName()=='PRIVATE_BETA_DRIVE':
            if params!=dict(bscale=protocol['scale_nm3'],bdrive_k=protocol['drive_k_kj_mol']) or force.getForceGroup()!=31:
                raise ValueError('beta drive constants differ from frozen protocol')
        else:
            raise ValueError('unrecognized research force')
    if any(before.getParticleMass(i)!=original.getParticleMass(i) for i in range(original.getNumParticles())):
        raise ValueError('particle mass changed in private experiment')
    after_base=openmm.XmlSerializer.deserialize(openmm.XmlSerializer.serialize(after))
    for i in reversed(range(count,after_base.getNumForces())):
        after_base.removeForce(i)
    if openmm.XmlSerializer.serialize(after_base)!=openmm.XmlSerializer.serialize(original):
        raise ValueError('phase2 physical system differs beyond prescribed restraints')
    return dict(original_forces_preserved=True,phase2_only_restraints_reduced=True,research_forces=extra)


def audit(plan_path,root,output):
    raw=plan_path.read_bytes()
    plan=json.loads(raw)
    verify(plan)
    run_raw=(root/'evidence.json').read_bytes()
    run=json.loads(run_raw)
    if run['plan_sha256']!=sha(raw) or (root/'executed-plan.json').read_bytes()!=raw or run['coordinates_restored']:
        raise ValueError('private beta-wall execution binding differs')
    expected=[(c['id'],c['seed'],arm) for c in plan['cases'] for arm in plan['protocol']['arms']]
    if [(r['id'],r['seed'],r['arm']) for r in run['cases']]!=expected:
        raise ValueError('complete beta stress enrollment differs')
    cases={(c['id'],c['seed']):c for c in plan['cases']}
    reports=[]
    for row in run['cases']:
        case=cases[(row['id'],row['seed'])]
        folder=root/row['id']/f"seed-{row['seed']}"/row['arm']
        for name,digest in row['artifacts'].items():
            if sha((folder/name).read_bytes())!=digest:
                raise ValueError('retained private beta evidence changed')
        if json.loads((folder/'result.json').read_bytes())!=row:
            raise ValueError('retained beta cell result differs')
        report=dict(id=row['id'],seed=row['seed'],arm=row['arm'],status=row['status'])
        if row['status']=='completed':
            source=Path(case['source'])
            prefix=case['prefix']
            template=load_atoms(source/f'{prefix}-atoms.json.gz')
            bonds=json.loads((source/f'{prefix}-topology.json').read_bytes())['bonds']
            initial=load_atoms(folder/'initial.json.gz')
            delta=coordinate_delta(template,initial)
            if not delta['heavy']['exact'] or not delta['hydrogen']['exact']:
                raise ValueError('private beta arm did not start at exact saved coordinates')
            phases=[]
            for stage in ('phase1','phase2'):
                saved=json.loads(gzip.decompress((folder/f'{stage}.json.gz').read_bytes()))
                if not math.isfinite(saved['total_energy_kj_mol']):
                    raise ValueError('nonfinite private beta endpoint energy')
                rows=saved['atoms']
                movement=coordinate_delta(initial,rows)
                try:
                    measured=geometry(rows,bonds,case['target'])
                    rounded=geometry(rows,bonds,case['target'],rounded=True)
                    phases.append(dict(stage=stage,raw=measured,rounded=rounded,movement=movement,
                        accepted=measured['accepted'] and rounded['accepted']))
                except ValueError as error:
                    phases.append(dict(stage=stage,accepted=False,movement=movement,
                        geometry_audit_blocker=dict(type=type(error).__name__,reason=str(error))))
            forces=force_audit((source/f'{prefix}-system.xml').read_text(),folder/'system-before.xml',folder/'system-after.xml',row['arm'],plan['protocol'])
            report.update(phases=phases,force_audit=forces,initial_coordinates_exact=True,
                final_geometry_accepted=phases[-1]['accepted'],final_output_policy='eligible_for_research_only' if phases[-1]['accepted'] else 'refused_geometry')
        elif row['status']!=case['admission'] and row['status']!='execution_blocker':
            raise ValueError('undeclared admission/refusal outcome')
        reports.append(report)
    summaries={}
    for arm in plan['protocol']['arms']:
        rows=[r for r in reports if r['arm']==arm]
        completed=[r for r in rows if r['status']=='completed']
        summaries[arm]=dict(enrolled=len(rows),executed=len(completed),
            initial_refusals=sum(r['status'].startswith('refused') for r in rows),
            not_applicable=sum(r['status'].startswith('not_applicable') for r in rows),
            execution_blockers=sum(r['status']=='execution_blocker' for r in rows),
            final_geometry_passes=sum(r['final_geometry_accepted'] for r in completed),
            final_target_CB_R=sum(r['phases'][-1].get('raw',{}).get('stereo',{}).get('target',{}).get('target_cip',{}).get('CB')=='R' for r in completed))
    write(output,dict(schema='dvbfixer.private-beta-wall-audit.v1',cases=reports,summary=summaries,
        plan_sha256=sha(raw),run_sha256=sha(run_raw),audit_script_sha256=sha(Path(__file__).read_bytes()),
        native_coordinates_read=False,native_accuracy_reported=False,production_changed=False,profile_selected=False,
        original_failure_causally_explained=False,production_qualified=False,
        interpretation='artificial-drive stress and finite beta wall, all observations may move; final geometry rejection is separate from CB sign'))
    print(json.dumps(summaries,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    audit(args.plan,args.run,args.output)
