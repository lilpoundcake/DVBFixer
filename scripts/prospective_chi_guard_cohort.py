"""Pre-download enrollment and transport for a private new-input comparison.

Reuse frozen reconstruction/scoring/chi kernels; no production capability, new
angles/weights, native-based selection or replacement of refused inputs.
"""
from __future__ import annotations

import argparse
import json
import math
from datetime import UTC, datetime
from pathlib import Path
from unittest.mock import patch
from urllib.request import urlopen

from scripts import prospective_ncac_cohort as preparation
from scripts.research_dunbrack_prior import sha, write

SYSTEMS={'1ENH':'engrailed homeodomain','1CRN':'crambin','1QYS':'Top7 designed protein'}
KERNELS=('scripts/research_dunbrack_prior.py','scripts/research_dunbrack_prior_weights.py',
    'scripts/research_dunbrack_means.py','scripts/research_generated_chi_relax.py',
    'scripts/research_chi_stopping_controls.py','scripts/research_chi_endpoint_guard.py',
    'scripts/audit_generated_chi_relax.py','scripts/audit_ncac_expansion.py')


def choose_complete(rows,guarded=False):
    """Selection policy frozen BEFORE download; no native information accepted."""
    if any(r['status']!='evaluated' for r in rows):
        return dict(status='incomplete_pool',selected_index=None)
    if guarded and any(r['guard_decision']!='endpoint_gradient_met' for r in rows):
        return dict(status='incomplete_guard_pool',selected_index=None)
    if any(not math.isfinite(r['score']) for r in rows):
        raise ValueError('finite completed-pool scores required')
    feasible=[r for r in rows if r['common_accepted']]
    selected=min(feasible,key=lambda r:(r['score'],r['index']))['index'] if feasible else None
    return dict(status='selected' if selected is not None else 'empty_geometry_pool',selected_index=selected)


def verify(enrollment):
    for p,d in enrollment['source_hashes'].items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('pre-download frozen algorithm/selection implementation changed')


def enroll(output):
    output.mkdir(parents=True,exist_ok=False)
    base_path=output/'preparation-enrollment.json'
    with patch.object(preparation,'SYSTEMS',SYSTEMS):
        preparation.enroll(base_path)
    if preparation.SYSTEMS==SYSTEMS:
        raise ValueError('new-system preparation hook not restored')
    base_raw=base_path.read_bytes()
    record=json.loads(base_raw)
    for filename in (*KERNELS,'scripts/prospective_chi_guard_cohort.py'):
        record['source_hashes'][filename]=sha(Path(filename).read_bytes())
    record['schema']='dvbfixer.prospective-chi-guard-cohort.v1'
    record['preparation_enrollment_sha256']=sha(base_raw)
    record['max_masks']=12
    record['arms']=dict(stock='shared seeded PDBFixer Reference, unchanged; no reflection/refinement, NOT full legacy',
        own='unchanged synthetic three-staggered pool, complete common-feasible-pool fixed ref2015 score',
        mean_prior='NCAC normal geometry with installed-library three well means, max raw family weight',
        mean_fixed='same complete common-feasible mean pool, unchanged fixed-heavy ref2015 import score',
        mean_c='same mean starts, original absolute1e-6/gmax1e-4/max200 chi control; raw refined-H score, final common gate',
        mean_guard='same mean starts/inner C, actual endpoint guard at most THREE TOTAL cycles; final common gate; ANY exhaustion/execution blocker blocks whole guard pool')
    record['selection_contract']='choose_complete, frozen in this module; geometry rejections are retained completed evaluations; execution/budget shortfalls block complete pools; stable index ties; no cross-arm best-profile choice'
    record['native_policy']='source projection/masking uses only frozen identity/file-order rule; scientific proposals consume masked input only; commit ALL arm choices before native accuracy reporting; never use reference coordinates to generate/rank/refine'
    record['library_contract']='same installed PyRosetta revision/binary as retained research, complete three-well API samples at source phi/psi; raw probability() accessor weights not a normalized probability; periodic controls; no published table parameters'
    record['numerical_contract']='no new tolerance/restart/angle/weight grid; compare ALL declared arms, choose none by known native errors; no convergence claim'
    record['scope']='new deliberately selected source inputs, first model/first ATOM chain only; NOT randomized/blinded/temporal/training-independent; Top7 is Rosetta-designed, NOT evidence independent of Rosetta design/training'
    record['source_retention']='original download bytes/URL/time/SHA256 and every projected/masked source retained inside project; no source normalization into admission; no source substitutions'
    record['execution_gates']='transport/plans must bind this enrollment and unchanged kernels; current source preparation callable, library/mean/control adapters and independent all-arm selection/readback required before scientific accuracy qualification'
    record['no_structure_downloaded_by_enrollment']=True
    verify(record)
    write(output/'enrollment.json',record)
    print(json.dumps(dict(systems=list(SYSTEMS),max_masks=12,structures_downloaded=False,
        comparison_arms=list(record['arms']),production_enabled=False),indent=2))


def fetch(enrollment_path,output):
    raw=enrollment_path.read_bytes()
    record=json.loads(raw)
    verify(record)
    output.mkdir(parents=True,exist_ok=False)
    (output/'enrollment.json').write_bytes(raw)
    provenance=[]
    for system in record['systems']:
        row=dict(id=system['id'],url=system['url'],retrieved_utc=datetime.now(UTC).isoformat())
        try:
            with urlopen(system['url'],timeout=60) as response:
                data=response.read(20*1024*1024+1)
                row.update(status_code=response.status,resolved_url=response.geturl())
            if len(data)>20*1024*1024:
                raise ValueError('source exceeds frozen20MiB download limit')
            path=output/(system['id']+'.pdb')
            with path.open('xb') as stream:
                stream.write(data)
            row.update(bytes=len(data),sha256=sha(data))
            if row['status_code']!=200 or not any(line.startswith(b'ATOM  ') for line in data.splitlines()):
                raise ValueError('download is not a coordinate PDB source')
            row['status']='downloaded'
        except Exception as error:
            row.update(status='download_blocker',error_type=type(error).__name__,error=str(error))
        provenance.append(row)
        write(output/(system['id']+'-download.json'),row)
    write(output/'provenance.json',dict(enrollment_sha256=sha(raw),systems=provenance,replacements=False))
    print(json.dumps(provenance,indent=2))
    if any(r['status']!='downloaded' for r in provenance):
        raise RuntimeError('source acquisition incomplete; retain failures, no replacement')


def freeze(enrollment_path,inputs,output):
    raw=enrollment_path.read_bytes()
    record=json.loads(raw)
    verify(record)
    provenance=json.loads((inputs/'provenance.json').read_bytes())
    if provenance['enrollment_sha256']!=sha(raw) or [r['id'] for r in provenance['systems']]!=[r['id'] for r in record['systems']]:
        raise ValueError('source acquisition enrollment/inventory differs')
    for row in provenance['systems']:
        if row['status']!='downloaded' or sha((inputs/(row['id']+'.pdb')).read_bytes())!=row['sha256']:
            raise ValueError('source acquisition missing/changed')
    preparation.freeze(enrollment_path,inputs,output)
    (output/'download-provenance.json').write_bytes((inputs/'provenance.json').read_bytes())


def prepare(plan_root,output):
    record=json.loads((plan_root/'plan.json').read_bytes())
    if record['enrollment']['schema']!='dvbfixer.prospective-chi-guard-cohort.v1':
        raise ValueError('new-input guarded comparison enrollment required')
    verify(record['enrollment'])
    preparation.run(plan_root,output)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('enroll')
    p.add_argument('--output',type=Path,required=True)
    p=sub.add_parser('fetch')
    p.add_argument('--enrollment',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p=sub.add_parser('freeze')
    p.add_argument('--enrollment',type=Path,required=True)
    p.add_argument('--inputs',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p=sub.add_parser('prepare')
    p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.command=='enroll':
        enroll(args.output)
    elif args.command=='fetch':
        fetch(args.enrollment,args.output)
    elif args.command=='freeze':
        freeze(args.enrollment,args.inputs,args.output)
    else:
        prepare(args.plan,args.output)
