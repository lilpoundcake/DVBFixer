"""Independent enrollment/freeze layer around unchanged reconstruction/scoring.

New protein systems are specified before input inspection or proposal execution.
No alternative structure, altloc choice or occupancy rewrite on admission failure.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict
from pathlib import Path

from scripts.unfiltered_branched_comparison import json_bytes, masks, select_targets

SYSTEMS = {'1CRN': 'crambin', '1AKI': 'hen egg-white lysozyme', '1BTA': 'barstar'}
POLICY = Path('tmp/pdbfixer-research/rosetta-fixed-plan-03/plan.json')


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def implementation_hashes() -> dict:
    paths = [Path(__file__).resolve(), *sorted(Path('src/dvbfixer/protein_completion').glob('*.py')),
             *[Path('scripts')/p for p in ('unfiltered_branched_comparison.py', 'compare_branched_completion.py',
                                          'branched_template_baseline.py', 'branched_stereo_oracle.py',
                                          'research_rosetta_fixed_scoring.py', 'research_rosetta_repack.py',
                                          'rosetta_fixed_contract.py', 'summarize_rosetta_fixed_scoring.py',
                                          'audit_unfiltered_reference.py')]]
    root = Path('.').resolve()
    return {str(p.resolve().relative_to(root)): sha(p.read_bytes()) for p in paths}


def enroll(output: Path) -> None:
    policy = json.loads(POLICY.read_bytes())
    for name, digest in policy['source_hashes'].items():
        if sha(Path(name).read_bytes()) != digest:
            raise ValueError('chosen scoring implementation changed since retrospective pilot')
    record = dict(schema='dvbfixer.prospective-branched-enrollment.v1',
                  systems=[dict(id=k, description=v, url=f'https://files.rcsb.org/download/{k}.pdb') for k, v in SYSTEMS.items()],
                  selection='first two observed VAL and first two observed THR identities in chain-A file order; four masks each',
                  coordinate_model='first coordinate MODEL in file, explicitly projected; unwrapped single model if no MODEL records',
                  projection='every chain-A ATOM byte in that model, retaining altloc/occupancy and chain TER; retain LINK/SSBOND/CONECT touching selected atoms, including excluded endpoints for fail-closed admission',
                  expected_masks=48, replacement_structures_allowed=False, outcome_filtering=False,
                  score_policy_plan_sha256=sha(POLICY.read_bytes()),
                  score_selection_policy=policy['selection_policy'],
                  hydrogen_protocol=policy['hydrogen_protocol'], init_flags=policy['init_flags'],
                  expected_rosetta_revision=policy['expected_rosetta_revision'],
                  source_hashes=implementation_hashes(),
                  qualification_scope='prospective follow-up of an already chosen scoring rule on three deliberately selected protein systems; not randomized or blinded; correlated masks',
                  production_enabled=False)
    with output.open('xb') as stream:
        stream.write(json_bytes(record))


def project(original: bytes) -> bytes:
    """Keep first-model chain A and relevant original bond annotations verbatim."""
    lines = original.decode('ascii').splitlines(keepends=True)
    has_models = any(l.startswith('MODEL ') for l in lines)
    active = not has_models
    first_seen, last_selected = False, False
    coordinate_lines = []
    selected_serials = set()
    for line in lines:
        if line.startswith('MODEL '):
            active = not first_seen
            first_seen = True
            continue
        if line.startswith('ENDMDL'):
            active = False
            last_selected = False
            continue
        if not active:
            continue
        if line.startswith(('ATOM  ', 'HETATM')):
            last_selected = line.startswith('ATOM  ') and line[21] == 'A'
            if last_selected:
                coordinate_lines.append(line)
                selected_serials.add(int(line[6:11]))
        elif line[:6].strip() == 'TER' and last_selected:
            coordinate_lines.append(line)
            last_selected = False
    if not selected_serials:
        raise ValueError('predeclared first-model chain A absent; do not substitute a chain/model')
    metadata, conect = [], []
    in_model, model_count = False, 0
    for line in lines:
        if line.startswith('MODEL '):
            in_model = True
            model_count += 1
            continue
        if line.startswith('ENDMDL'):
            in_model = False
            continue
        if in_model and model_count > 1:
            continue
        record = line[:6].strip()
        if ((record == 'SSBOND' and (line[15:16] == 'A' or line[29:30] == 'A'))
                or (record == 'LINK' and (line[21:22] == 'A' or line[51:52] == 'A'))):
            metadata.append(line)
        elif record == 'CONECT':
            fields = line[6:].rstrip('\n')
            values = [int(fields[i:i+5]) for i in range(0, len(fields), 5) if fields[i:i+5].strip()]
            if selected_serials.intersection(values):
                conect.append(line)
    if not coordinate_lines[-1].startswith('TER'):
        coordinate_lines.append('TER\n')
    return (''.join(metadata+coordinate_lines+conect)+'END\n').encode('ascii')


def amend(enrollment: Path, output: Path) -> None:
    """Explicit post-import-blocker amendment; never overwrite original enrollment."""
    raw = enrollment.read_bytes()
    record = json.loads(raw)
    current = implementation_hashes()
    for name, digest in record['source_hashes'].items():
        if name != 'scripts/prospective_branched_cohort.py' and current.get(name) != digest:
            raise ValueError('only projection serialization may change in this amendment')
    record.update(source_hashes=current,
                  parent_enrollment_sha256=sha(raw),
                  protocol_amendment=dict(reason='OpenMM cannot read CONECT before the coordinate model exists',
                                          change='place retained CONECT after coordinates; no atom or link changes',
                                          reconstruction_results_already_observed=True,
                                          new_cohort_scoring_not_started=True,
                                          cohort_and_scoring_rule_unchanged=True))
    with output.open('xb') as stream:
        stream.write(json_bytes(record))


def freeze(enrollment: Path, inputs: Path, output: Path) -> None:
    raw = enrollment.read_bytes()
    record = json.loads(raw)
    if record['source_hashes'] != implementation_hashes():
        raise ValueError('implementation changed after enrollment')
    payloads, cases = {}, []
    for system in record['systems']:
        label = system['id']
        original = (inputs/(label+'.pdb')).read_bytes()
        reference = project(original)
        payloads[label+'-original.pdb'], payloads[label+'-reference.pdb'] = original, reference
        lines = reference.decode('ascii').splitlines(keepends=True)
        for name, target in select_targets(lines):
            for mask in masks(name):
                identifier = f"{label}-{name}{target.sequence_number}{target.insertion_code}-{'-'.join(mask)}"
                source = ''.join(l for l in lines if not (l.startswith('ATOM  ') and l[21] == target.chain_id
                                 and l[22:26].strip() == target.sequence_number and l[26].strip() == target.insertion_code
                                 and l[12:16].strip() in mask)).encode('ascii')
                payloads[identifier+'.pdb'] = source
                cases.append(dict(id=identifier, fixture=label, target=asdict(target), residue_name=name,
                                  mask=mask, mode='anchored' if len(mask) == 1 else 'paired',
                                  source=identifier+'.pdb', reference=label+'-reference.pdb'))
    if len(cases) != record['expected_masks']:
        raise ValueError('full predeclared enrollment required; retain unavailable structure, do not replace')
    plan = dict(schema='dvbfixer.unfiltered-branched-plan.v1', case_count=len(cases), cases=cases,
                artifacts={n: sha(data) for n, data in payloads.items()}, source_hashes=record['source_hashes'],
                prospective_enrollment_sha256=sha(raw), enrollment=record,
                common_gate='unchanged measured-chi/raw/rounded/stereo gate',
                no_own_success_filter=True, production_enabled=False)
    output.mkdir(parents=True, exist_ok=False)
    for name, data in payloads.items():
        (output/name).write_bytes(data)
    (output/'enrollment.json').write_bytes(raw)
    (output/'plan.json').write_bytes(json_bytes(plan))
    print(f'Frozen all {len(cases)} masks before reconstruction/scoring')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('enroll')
    p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('freeze')
    p.add_argument('--enrollment', type=Path, required=True)
    p.add_argument('--inputs', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('amend-conect-order')
    p.add_argument('--enrollment', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'enroll':
        enroll(args.output)
    elif args.command == 'freeze':
        freeze(args.enrollment, args.inputs, args.output)
    else:
        amend(args.enrollment, args.output)


if __name__ == '__main__':
    main()
