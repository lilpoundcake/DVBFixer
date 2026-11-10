"""Saved-only full legacy audit in its actual frame, then separate native report.

Heavy-only remasking is a diagnostic view of the real output, never an exported
reconstruction or observation restoration. Preserve raw H-containing outputs.
"""
from __future__ import annotations

import argparse
import json
import math
import statistics
from dataclasses import asdict
from pathlib import Path

from scripts.compare_branched_completion import measure
from scripts.research_dunbrack_prior import sha, write
from scripts.research_legacy_comparator import verify
from scripts.unfiltered_branched_comparison import reporting_scores

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.ffutils.dat import DatRecord
from dvbfixer.ffutils.variants import text_rename_variants_to_parent
from dvbfixer.pdb_component_reconstruction import _component_ref, _parse_atom_lines
from dvbfixer.protein_completion.branched_source import plan_branched_source
from dvbfixer.protein_completion.contracts import Addition


def atoms(raw):
    result={}
    for atom in _parse_atom_lines(raw.decode().splitlines(keepends=True)):
        if atom.element in {'H','D'}:
            continue
        key=(_component_ref(atom),atom.name)
        if key in result:
            raise ValueError('duplicate heavy identity in actual output')
        result[key]=atom
    return result


def motion(source,output):
    before,after=atoms(source),atoms(output)
    missing=before.keys()-after.keys()
    distances=[math.dist(a.position,after[key].position) for key,a in before.items() if key in after]
    return dict(original_heavy_count=len(before),matched_heavy_count=len(distances),
        missing_identities=[repr(k) for k in sorted(missing,key=repr)],
        additional_heavy_count=len(after.keys()-before.keys()),
        moved_at_pdb_precision=sum(d>0 for d in distances),
        rms_displacement_angstrom=math.sqrt(sum(d*d for d in distances)/len(distances)) if distances else None,
        maximum_displacement_angstrom=max(distances,default=None),aligned=False,observations_restored=False)


def diagnostic_projection(raw,target,mask):
    """Remove H and only declared target atoms; keep actual heavy coordinates."""
    parsed=_parse_atom_lines(raw.decode().splitlines(keepends=True))
    removed={a.serial for a in parsed if a.element in {'H','D'} or (_component_ref(a)==target and a.name in mask)}
    lines=[]
    for line in raw.decode().splitlines(keepends=True):
        if line.startswith(('ATOM  ','HETATM','ANISOU')) and int(line[6:11]) in removed:
            continue
        if line.startswith('CONECT'):
            serials=[int(line[i:i+5]) for i in range(6,len(line.rstrip()),5) if line[i:i+5].strip()]
            if any(s in removed for s in serials):
                continue
        lines.append(line)
    projected=''.join(lines).encode()
    actual=atoms(raw)
    for key,atom in atoms(projected).items():
        if actual[key].position!=atom.position:
            raise ValueError('diagnostic projection changed real heavy coordinates')
    return projected


def audit(plan_path,run_root,output):
    raw=plan_path.read_bytes()
    plan=json.loads(raw)
    verify(plan)
    if (run_root/'executed-plan.json').read_bytes()!=raw:
        raise ValueError('legacy executed plan differs')
    run_raw=(run_root/'evidence.json').read_bytes()
    run=json.loads(run_raw)
    if run['plan_sha256']!=sha(raw) or run['native_coordinates_read'] or run['observations_restored']:
        raise ValueError('legacy run scope differs')
    if [r['id'] for r in run['cases']]!=[c['id'] for c in plan['cases']]:
        raise ValueError('legacy enrolled inventory differs')
    output.mkdir(parents=True,exist_ok=False)
    results=[]
    for case,executed in zip(plan['cases'],run['cases'],strict=True):
        folder=run_root/case['id']
        for path,digest in executed['artifacts'].items():
            if sha((folder/path).read_bytes())!=digest:
                raise ValueError('raw legacy output/log changed')
        source=Path(case['source']).read_bytes()
        if source!=(folder/'input.pdb').read_bytes():
            raise ValueError('legacy input changed')
        record=dict(id=case['id'],execution_status=executed['status'],stages={})
        target=ComponentInstanceRef(**case['target'])
        for stage,filename in (('prepare','prepared.pdb'),('minimize','minimized.pdb')):
            path=folder/filename
            if not path.exists():
                record['stages'][stage]=dict(status='no_output')
                continue
            actual=path.read_bytes()
            # Shared variant normalization changes names only on a private copy.
            private=output/(case['id']+'-'+stage+'-variant-input.pdb')
            private.write_bytes(actual)
            normalized_path,variants=text_rename_variants_to_parent(str(private))
            normalized=Path(normalized_path).read_bytes()
            if normalized_path!=str(private):
                Path(normalized_path).unlink()
            derived=diagnostic_projection(normalized,target,case['mask'])
            derived_path=output/(case['id']+'-'+stage+'-heavy-remasked.pdb')
            derived_path.write_bytes(derived)
            result=dict(status='audited',raw_sha256=sha(actual),original_observation_motion=motion(source,actual),
                normalized_variant_residues=len(variants),diagnostic_projection_sha256=sha(derived),
                geometry_scope='heavy-only actual relaxed frame, target remasked; no original-coordinate restore; H geometry not covered')
            try:
                request=plan_branched_source(derived,target).request
                full=atoms(normalized)
                added=tuple(Addition(ExactAtomRef(target,n),full[(target,n)].element,full[(target,n)].position) for n in case['mask'])
                common,_=measure(derived,request,added,anchored=False)
                result['actual_frame_common_audit']=asdict(common)
            except ValueError as error:
                result['geometry_audit_blocker']=dict(type=type(error).__name__,reason=str(error))
            record['stages'][stage]=result
        dat=folder/'prepared.dat'
        if dat.exists():
            record['dat_added_atom_count']=len(DatRecord.load(dat).added_atoms)
        results.append(record)
    verify(plan)
    write(output/'selection-before-native.json',dict(schema='dvbfixer.legacy-comparator-audit.v1',cases=results,
        plan_sha256=sha(raw),run_sha256=sha(run_raw),audit_script_sha256=sha(Path(__file__).read_bytes()),
        native_coordinates_read=False,observations_restored=False,profile_selected=False,
        production_changed=False,scientific_qualification=False))
    print(json.dumps(dict(cases=len(results),execution_completed=sum(r['execution_status']=='completed' for r in results),
        final_geometry_passes=sum(r['stages'].get('minimize',{}).get('actual_frame_common_audit',{}).get('accepted',False) for r in results)),indent=2))


def report(plan_path,run_root,audit_root,output):
    plan=json.loads(plan_path.read_bytes())
    verify(plan)
    ledger_raw=(audit_root/'selection-before-native.json').read_bytes()
    ledger=json.loads(ledger_raw)
    if ledger['plan_sha256']!=sha(plan_path.read_bytes()) or ledger['audit_script_sha256']!=sha(Path(__file__).read_bytes()):
        raise ValueError('legacy pre-native audit binding differs')
    masks=json.loads((Path(plan['mask_root'])/'plan.json').read_bytes())
    results=[]
    for case,row in zip(plan['cases'],ledger['cases'],strict=True):
        if case['id']!=row['id']:
            raise ValueError('legacy report inventory differs')
        native_path=Path(plan['mask_root'])/case['reference']
        native=native_path.read_bytes()
        if sha(native)!=masks['artifacts'][case['reference']]:
            raise ValueError('legacy native reporting reference changed')
        target=ComponentInstanceRef(**case['target'])
        source=Path(case['source']).read_bytes()
        request=plan_branched_source(source,target).request
        truth={name:a.position for (component,name),a in atoms(native).items() if component==target}
        stages={}
        for stage,filename in (('prepare','prepared.pdb'),('minimize','minimized.pdb')):
            retained=row['stages'][stage]
            path=run_root/case['id']/filename
            if not path.exists():
                continue
            actual=path.read_bytes()
            if sha(actual)!=retained['raw_sha256']:
                raise ValueError('legacy actual output changed after audit')
            full=atoms(actual)
            added=tuple(Addition(ExactAtomRef(target,n),full[(target,n)].element,full[(target,n)].position) for n in case['mask'])
            stages[stage]=dict(original_frame_accuracy=reporting_scores(request,added,truth),
                observation_motion=retained['original_observation_motion'],
                actual_frame_common_audit=retained.get('actual_frame_common_audit'),
                geometry_audit_blocker=retained.get('geometry_audit_blocker'))
        results.append(dict(id=case['id'],stages=stages,execution_status=row['execution_status']))
    aggregate={}
    for stage in ('prepare','minimize'):
        rows=[r['stages'][stage] for r in results if stage in r['stages']]
        aggregate[stage]=dict(enrolled=len(results),reported=len(rows),
            original_frame_target_rmsd_median_angstrom=statistics.median(r['original_frame_accuracy']['exact_name_rmsd_angstrom'] for r in rows) if rows else None,
            actual_frame_geometry_passes=sum(bool(r['actual_frame_common_audit'] and r['actual_frame_common_audit']['accepted']) for r in rows),
            maximum_observed_motion_angstrom=max((r['observation_motion']['maximum_displacement_angstrom'] for r in rows),default=None))
    write(output,dict(schema='dvbfixer.legacy-comparator-report.v1',cases=results,aggregate=aggregate,
        plan_sha256=sha(plan_path.read_bytes()),audit_ledger_sha256=sha(ledger_raw),
        native_reporting_after_geometry_audit=True,observations_restored=False,
        profile_selected=False,production_changed=False,scientific_qualification=False))
    print(json.dumps(aggregate,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    for name in ('audit','report'):
        p=sub.add_parser(name)
        p.add_argument('--plan',type=Path,required=True)
        p.add_argument('--run',type=Path,required=True)
        p.add_argument('--output',type=Path,required=True)
        if name=='report':
            p.add_argument('--audit',type=Path,required=True)
    args=parser.parse_args()
    if args.command=='audit':
        audit(args.plan,args.run,args.output)
    else:
        report(args.plan,args.run,args.audit,args.output)
