"""New-system, pre-download enrollment around unchanged three-chi/template tools.

First-model/first-ATOM-chain selection is identity-only. Preserve unsupported
chemistry, H, altloc and occupancy, every refusal and all external link endpoints.
No replacement system, repair, score retuning or native-based candidate selection.
"""
from __future__ import annotations

import argparse
import json
import platform
from contextlib import redirect_stderr, redirect_stdout
from dataclasses import asdict
from pathlib import Path

from scripts.branched_template_baseline import preflight, template_proposal
from scripts.compare_branched_completion import measure
from scripts.prospective_branched_cohort import implementation_hashes
from scripts.research_rosetta_fixed_scoring import plan as build_score_plan
from scripts.research_rosetta_joint import POLICY, sha, verify
from scripts.research_rosetta_repack import write_json
from scripts.unfiltered_branched_comparison import json_bytes, masks

from dvbfixer.domain.structure_identity import ComponentInstanceRef
from dvbfixer.protein_completion.anchored_branch import complete_anchored_branch
from dvbfixer.protein_completion.anchored_source import plan_anchored_source
from dvbfixer.protein_completion.branched_source import plan_branched_source
from dvbfixer.protein_completion.service import complete_sidechain

SYSTEMS = {'2CI2': 'chymotrypsin inhibitor 2', '1R69': '434 repressor amino-terminal domain',
           '1CTF': 'ribosomal protein L7/L12 carboxy-terminal domain'}


def project_first(original: bytes) -> tuple[bytes, str]:
    lines = original.decode('ascii').splitlines(keepends=True)
    has_models = any(line.startswith('MODEL ') for line in lines)
    active, seen, chain, last = not has_models, False, None, False
    coordinates, serials = [], set()
    for line in lines:
        if line.startswith('MODEL '):
            active = not seen
            seen = True
            continue
        if line.startswith('ENDMDL'):
            active, last = False, False
            continue
        if not active:
            continue
        if line.startswith(('ATOM  ', 'HETATM')):
            if line.startswith('ATOM  ') and chain is None:
                chain = line[21]
            last = line.startswith('ATOM  ') and line[21] == chain
            if last:
                coordinates.append(line)
                serials.add(int(line[6:11]))
        elif line[:6].strip() == 'TER' and last:
            coordinates.append(line)
            last = False
    if chain is None or not coordinates:
        raise ValueError('first coordinate model has no ATOM chain; no substitution')
    metadata, conect, in_model, model_number = [], [], False, 0
    for line in lines:
        if line.startswith('MODEL '):
            in_model, model_number = True, model_number+1
            continue
        if line.startswith('ENDMDL'):
            in_model = False
            continue
        if in_model and model_number > 1:
            continue
        record = line[:6].strip()
        if ((record == 'LINK' and chain in (line[21:22], line[51:52]))
                or (record == 'SSBOND' and chain in (line[15:16], line[29:30]))):
            metadata.append(line)
        elif record == 'CONECT':
            fields = line[6:].rstrip('\n')
            endpoints = [int(fields[i:i+5]) for i in range(0, len(fields), 5) if fields[i:i+5].strip()]
            if serials.intersection(endpoints):
                conect.append(line)
    if not coordinates[-1].startswith('TER'):
        coordinates.append('TER\n')
    return (''.join(metadata+coordinates+conect)+'END\n').encode(), chain


def selected_targets(reference: bytes):
    seen, selected, counts = set(), [], {'VAL': 0, 'THR': 0}
    for line in reference.decode('ascii').splitlines():
        if not line.startswith('ATOM  ') or line[17:20] not in counts:
            continue
        name = line[17:20]
        ref = ComponentInstanceRef(1, line[21], line[22:26].strip(), line[26].strip())
        if ref not in seen and counts[name] < 2:
            selected.append((name, ref))
            counts[name] += 1
        seen.add(ref)
    return selected, counts


def enroll(output: Path):
    policy_raw = POLICY.read_bytes()
    policy = json.loads(policy_raw)
    for filename, digest in policy['source_hashes'].items():
        if sha(Path(filename).read_bytes()) != digest:
            raise ValueError('previously chosen scoring policy changed')
    hashes = implementation_hashes()
    hashes['scripts/prospective_independent_cohort.py'] = sha(Path(__file__).read_bytes())
    for filename in ('scripts/prospective_scoring_followup.py', 'src/dvbfixer/ffutils/geometry.py',
                     'src/dvbfixer/ffutils/dat.py', 'src/dvbfixer/pdb_component_reconstruction.py',
                     'src/dvbfixer/domain/structure_identity.py', 'src/dvbfixer/diagnose/geometry.py'):
        hashes[filename] = sha(Path(filename).read_bytes())
    record = dict(schema='dvbfixer.independent-system-enrollment.v1',
                  systems=[dict(id=k, description=v, url=f'https://files.rcsb.org/download/{k}.pdb') for k, v in SYSTEMS.items()],
                  selection='first coordinate MODEL, first encountered ATOM chain; up to first two VAL and first two THR identities, all four masks; retain count shortfalls',
                  observations='preserve chain ATOM bytes, H/altloc/occupancy, TER and LINK/SSBOND/CONECT touching selected atoms; excluded endpoints retained for fail-closed admission',
                  max_cases=48, expected_complete_system_cases=16, replacement_structures_allowed=False,
                  own='unchanged three-staggered proposals or exact anchored closure; clearance baseline and unchanged ref2015 complete-feasible-pool selector',
                  template='unchanged shared seeded PDBFixer helper on Reference; snapshots retained; no reflection or refinement; not full legacy',
                  common_gate='unchanged measured-chi/raw/rounded/stereo/sterics, not native-tuned or a calibrated gold standard',
                  score_policy_plan_sha256=sha(policy_raw), score_selection_policy=policy['selection_policy'],
                  hydrogen_protocol=policy['hydrogen_protocol'], init_flags=policy['init_flags'],
                  expected_rosetta_revision=policy['expected_rosetta_revision'], source_hashes=hashes,
                  qualification_scope='new deliberately selected systems/folds for this research arm; prospective input/policy enrollment, not randomized, blinded, temporal holdout or independence from Rosetta training/templates',
                  no_truth_in_proposals_or_selection=True, production_enabled=False)
    with output.open('xb') as stream:
        stream.write(json_bytes(record))


def freeze(enrollment_path: Path, inputs: Path, output: Path):
    raw = enrollment_path.read_bytes()
    enrollment = json.loads(raw)
    for filename, digest in enrollment['source_hashes'].items():
        if sha(Path(filename).read_bytes()) != digest:
            raise ValueError('enrolled implementation changed')
    output.mkdir(parents=True, exist_ok=False)
    artifacts, cases, systems = {}, [], {}

    def save(filename, data):
        (output/filename).write_bytes(data)
        artifacts[filename] = sha(data)

    for system in enrollment['systems']:
        identifier = system['id']
        original = (inputs/(identifier+'.pdb')).read_bytes()
        save(identifier+'-original.pdb', original)
        try:
            reference, chain = project_first(original)
        except ValueError as error:
            systems[identifier] = dict(status='projection_refused', reason=str(error), original_sha256=sha(original))
            continue
        save(identifier+'-reference.pdb', reference)
        targets, counts = selected_targets(reference)
        systems[identifier] = dict(status='projected', chain=chain, selected_counts=counts,
                                  count_shortfall=counts != {'VAL': 2, 'THR': 2}, original_sha256=sha(original))
        for name, target in targets:
            for mask in masks(name):
                case_id = f'{identifier}-{name}{target.sequence_number}{target.insertion_code}-'+ '-'.join(mask)
                source = b''.join(line for line in reference.splitlines(keepends=True) if not (
                    line.startswith(b'ATOM  ') and (line[21:22].decode(), line[22:26].strip().decode(), line[26:27].strip().decode())
                    == (target.chain_id, target.sequence_number, target.insertion_code) and line[12:16].strip().decode() in mask))
                filename = case_id+'.pdb'
                save(filename, source)
                cases.append(dict(id=case_id, fixture=identifier, target=asdict(target), residue_name=name, mask=mask,
                                  mode='anchored' if len(mask) == 1 else 'paired', source=filename, reference=identifier+'-reference.pdb'))
    record = dict(schema='dvbfixer.independent-system-mask-plan.v1', cases=cases, case_count=len(cases), systems=systems,
                  artifacts=artifacts, source_hashes=enrollment['source_hashes'], prospective_enrollment_sha256=sha(raw),
                  enrollment=enrollment, production_enabled=False)
    (output/'enrollment.json').write_bytes(raw)
    write_json(output/'plan.json', record)


def run(plan_root: Path, output: Path):
    raw = (plan_root/'plan.json').read_bytes()
    plan = json.loads(raw)
    verify(plan, plan_root)
    environment = preflight()  # A broken environment is not a scientific failure.
    output.mkdir(parents=True, exist_ok=False)
    (output/'executed-plan.json').write_bytes(raw)
    results = {}
    for case in plan['cases']:
        root = output/case['id']
        root.mkdir()
        source = (plan_root/case['source']).read_bytes()
        (root/'input.pdb').write_bytes(source)
        target = ComponentInstanceRef(**case['target'])
        anchored = case['mode'] == 'anchored'
        try:
            source_plan = (plan_anchored_source if anchored else plan_branched_source)(source, target)
            if set(source_plan.request.missing_atoms) != set(case['mask']):
                raise ValueError('declared deletion mask differs from actual missing canonical atoms')
        except ValueError as error:
            row = dict(status='shared_admission_refused', reason=str(error), own_executed=False, template_executed=False)
        else:
            request = source_plan.request
            (root/'source-plan.json').write_bytes(json_bytes(asdict(source_plan)))
            own = (complete_anchored_branch if anchored else complete_sidechain)(request)
            (root/'own-attempts.json').write_bytes(json_bytes(asdict(own)))
            row = dict(status='admitted', request_digest=request.digest, own_executed=True, template_executed=False,
                       own=dict(status=own.status, reason=own.reason))
            if own.selected is not None:
                audit, pdb = measure(source, request, own.selected.additions, anchored=anchored)
                (root/'own-proposal.pdb').write_bytes(pdb)
                row['own']['common_audit'] = asdict(audit)

            def retain(snapshot, root=root):
                (root/f'template-seed-{snapshot.seed}.json').write_bytes(json_bytes(asdict(snapshot)))

            try:
                with (root/'template.log').open('x') as log, redirect_stdout(log), redirect_stderr(log):
                    template = template_proposal(source, request, retain_snapshot=retain)
                (root/'template-proposal.json').write_bytes(json_bytes(asdict(template)))
                audit, pdb = measure(source, request, template.additions, anchored=anchored)
                (root/'template-proposal.pdb').write_bytes(pdb)
                row['template'] = dict(status='evaluated', common_audit=asdict(audit), helper_attempts=len(template.snapshots))
                row['template_executed'] = True
            except Exception as error:
                row['template'] = dict(status='execution_or_protocol_blocker', error_type=type(error).__name__, error=str(error))
                (root/'execution-error.json').write_bytes(json_bytes(row['template']))
        (root/'comparison.json').write_bytes(json_bytes(row))
        row['retained_files'] = {p.name: sha(p.read_bytes()) for p in sorted(root.iterdir())}
        results[case['id']] = row
    write_json(output/'evidence.json', dict(schema='dvbfixer.independent-system-comparison.v1', results=results,
                                           plan_root=str(plan_root.resolve()), retained_directory=str(output.resolve()),
                                           plan_sha256=sha(raw), environment=environment, python=platform.python_version()))


def prepare_scores(enrollment_path: Path, comparison_path: Path, output: Path):
    from scripts.prospective_scoring_followup import check_policy

    enrollment_raw = enrollment_path.read_bytes()
    enrollment = json.loads(enrollment_raw)
    comparison = json.loads(comparison_path.read_bytes())
    mask_raw = (Path(comparison['plan_root'])/'plan.json').read_bytes()
    masks_plan = json.loads(mask_raw)
    if sha(mask_raw) != comparison['plan_sha256'] or masks_plan['prospective_enrollment_sha256'] != sha(enrollment_raw):
        raise ValueError('prospective enrollment/comparison binding differs')
    blockers = [k for k, r in comparison['results'].items()
                if r['status'] == 'admitted' and r['template']['status'] == 'execution_or_protocol_blocker']
    if blockers:
        raise ValueError(f'retained template execution/protocol blockers prevent unchanged scoring adapter: {blockers}')
    output.mkdir(parents=True, exist_ok=False)
    base = output/'archived-builder-plan'
    build_score_plan(comparison_path, base)
    record = json.loads((base/'plan.json').read_bytes())
    check_policy(enrollment, record)
    for filename, digest in record['artifacts'].items():
        data = (base/filename).read_bytes()
        if sha(data) != digest:
            raise ValueError('score builder changed an artifact')
        path = output/filename
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
    record.update(retrospective=False, prospectively_enrolled=True, prospective_enrollment_sha256=sha(enrollment_raw),
                  archived_builder_plan_sha256=sha((base/'plan.json').read_bytes()),
                  provenance_adapter_sha256=sha(Path(__file__).read_bytes()), qualification_scope=enrollment['qualification_scope'])
    (output/'enrollment.json').write_bytes(enrollment_raw)
    write_json(output/'plan.json', record)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('enroll')
    p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('freeze')
    p.add_argument('--enrollment', type=Path, required=True)
    p.add_argument('--inputs', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('run')
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('prepare-scores')
    p.add_argument('--enrollment', type=Path, required=True)
    p.add_argument('--comparison', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'enroll':
        enroll(args.output)
    elif args.command == 'freeze':
        freeze(args.enrollment, args.inputs, args.output)
    elif args.command == 'run':
        run(args.plan, args.output)
    else:
        prepare_scores(args.enrollment, args.comparison, args.output)
