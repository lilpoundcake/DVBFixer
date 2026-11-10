"""Bind a prospective cohort to the already chosen, unchanged scoring protocol.

The legacy builder's retrospective label is retained in an archived base plan;
the final, not-yet-executed plan gets explicit enrollment provenance. No score
weights, eligibility, candidate generation or selection rule may change here.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

from scripts.research_rosetta_repack import write_json

POLICY = Path('tmp/pdbfixer-research/rosetta-fixed-plan-03/plan.json')


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def check_policy(enrollment: dict, plan: dict) -> None:
    policy_raw = POLICY.read_bytes()
    if sha(policy_raw) != enrollment['score_policy_plan_sha256']:
        raise ValueError('original chosen policy plan changed')
    policy = json.loads(policy_raw)
    for key in ('score_function', 'seed', 'radii', 'init_flags', 'hydrogen_protocol', 'extra_heavy_policy',
                'virtual_atom_policy', 'selection_policy', 'native_policy', 'topology_scope', 'expected_rosetta_revision'):
        if plan[key] != policy[key]:
            raise ValueError('previously chosen scoring policy changed')
    for enrolled, actual in [('score_selection_policy', 'selection_policy'), ('hydrogen_protocol', 'hydrogen_protocol'),
                              ('init_flags', 'init_flags'), ('expected_rosetta_revision', 'expected_rosetta_revision')]:
        if enrollment[enrolled] != plan[actual]:
            raise ValueError('previously chosen scoring policy changed')
    for path, digest in enrollment['source_hashes'].items():
        if sha(Path(path).read_bytes()) != digest:
            raise ValueError('implementation changed since prospective enrollment')
    for path, digest in plan['source_hashes'].items():
        if enrollment['source_hashes'].get(path) != digest:
            raise ValueError('scoring source not bound by enrollment')


def prepare(enrollment_path: Path, comparison_path: Path, output: Path) -> None:
    from scripts.research_rosetta_fixed_scoring import plan as build_plan

    enrollment_raw = enrollment_path.read_bytes()
    enrollment = json.loads(enrollment_raw)
    comparison = json.loads(comparison_path.read_bytes())
    masks_raw = (Path(comparison['plan_root'])/'plan.json').read_bytes()
    masks = json.loads(masks_raw)
    if sha(masks_raw) != comparison['plan_sha256'] or masks['prospective_enrollment_sha256'] != sha(enrollment_raw):
        raise ValueError('cohort is not bound to this prospective enrollment')
    output.mkdir(parents=True, exist_ok=False)
    base = output/'archived-builder-plan'
    build_plan(comparison_path, base)
    plan_raw = (base/'plan.json').read_bytes()
    record = json.loads(plan_raw)
    check_policy(enrollment, record)
    # All payloads remain byte-identical. Only the final plan's provenance changes.
    for name, digest in record['artifacts'].items():
        data = (base/name).read_bytes()
        if sha(data) != digest:
            raise ValueError('builder artifact changed')
        destination = output/name
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(data)
    record.update(retrospective=False, prospectively_enrolled=True,
                  prospective_enrollment_sha256=sha(enrollment_raw),
                  archived_builder_plan_sha256=sha(plan_raw),
                  provenance_adapter_sha256=sha(Path(__file__).read_bytes()),
                  qualification_scope=enrollment['qualification_scope'])
    (output/'enrollment.json').write_bytes(enrollment_raw)
    write_json(output/'plan.json', record)


def report_changed_choices(plan: dict, root: Path, changed: list) -> list:
    """Existing generated-only VAL symmetry reporting; no new selection policy."""
    from scripts.unfiltered_branched_comparison import reporting_scores

    from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
    from dvbfixer.pdb_component_reconstruction import _component_ref, _parse_atom_lines
    from dvbfixer.protein_completion.anchored_source import plan_anchored_source
    from dvbfixer.protein_completion.branched_source import plan_branched_source
    from dvbfixer.protein_completion.contracts import Addition

    cases = {c['id']: c for c in plan['cases']}
    reports = []
    for row in changed:
        case = cases[row['id']]
        target = ComponentInstanceRef(**case['target'])
        planner = plan_anchored_source if case['mode'] == 'anchored' else plan_branched_source
        request = planner((root/case['id']/'source.pdb').read_bytes(), target).request
        truth = {a.name: a.position for a in _parse_atom_lines((root/case['reference_reporting_only']).read_text().splitlines(keepends=True))
                 if _component_ref(a) == target}
        values = {}
        for label, index in [('previous', row['previous']), ('selected', row['selected'])]:
            cell = next(c for c in case['cells'] if c['arm'] == 'own' and c['index'] == index)
            points = {a.name: a.position for a in _parse_atom_lines((root/cell['pdb']).read_text().splitlines(keepends=True))
                      if _component_ref(a) == target}
            additions = tuple(Addition(ExactAtomRef(target, n), n[0], points[n]) for n in request.missing_atoms)
            values[label] = reporting_scores(request, additions, truth)
        reports.append(dict(id=case['id'], reporting_only=True, **values))
    return reports


def summarize(enrollment_path: Path, masks: Path, comparison: Path, scores: Path, run: Path, output: Path, *,
              native_audit: Path, amendment_audit: Path, own_replay_audit: Path, parent_enrollment: Path) -> None:
    enrollment_raw = enrollment_path.read_bytes()
    enrollment = json.loads(enrollment_raw)
    mask_plan = json.loads((masks/'plan.json').read_bytes())
    score_plan = json.loads((scores/'plan.json').read_bytes())
    if mask_plan['prospective_enrollment_sha256'] != sha(enrollment_raw) or score_plan['prospective_enrollment_sha256'] != sha(enrollment_raw):
        raise ValueError('enrollment binding changed')
    check_policy(enrollment, score_plan)
    raw = json.loads((run/'evidence.json').read_bytes())
    if raw['plan_sha256'] != sha((scores/'plan.json').read_bytes()):
        raise ValueError('executed scoring plan differs')
    # Read the frozen helper's independent readback, without changing its policy.
    from scripts.summarize_rosetta_fixed_scoring import summarize as readback

    inner = output.with_name(output.stem+'-readback.json')
    readback(scores, run, inner)
    result = json.loads(inner.read_bytes())
    reconstructed = json.loads(comparison.read_bytes())
    native = json.loads(native_audit.read_bytes())
    if native['comparison_sha256'] != sha(comparison.read_bytes()):
        raise ValueError('native audit belongs to another reconstruction record')
    parent_raw = parent_enrollment.read_bytes()
    if enrollment['parent_enrollment_sha256'] != sha(parent_raw):
        raise ValueError('amendment parent changed')
    amendment = json.loads(amendment_audit.read_bytes())
    if not amendment['case_declarations_identical']:
        raise ValueError('amendment changed enrollment')
    for name, row in amendment['artifacts'].items():
        if sha((masks/name).read_bytes()) != row['new_sha256'] or not row['nonconect_bytes_identical'] or not row['conect_records_identical']:
            raise ValueError('amendment payload audit differs')
    own_replay = json.loads(own_replay_audit.read_bytes())
    if not own_replay['source_statuses_identical'] or any(not r['all_generated_coordinates_and_audits_identical'] for r in own_replay['cases'].values()):
        raise ValueError('own proposals changed after annotation-order correction')
    by_system = {}
    for system in enrollment['systems']:
        names = [c['id'] for c in mask_plan['cases'] if c['fixture'] == system['id']]
        rows = [result['cases'][n] for n in names]
        by_system[system['id']] = dict(
            description=system['description'], planned=len(names),
            source_statuses=dict(Counter(r['status'] for r in rows)),
            refusal_reasons=dict(Counter(r['reason'] for r in rows if r['status'] == 'shared_admission_refused')),
            old_own_passes=sum(reconstructed['results'][n].get('own', {}).get('common_audit', {}).get('accepted', False) for n in names),
            template_passes=sum(reconstructed['results'][n].get('template', {}).get('common_audit', {}).get('accepted', False) for n in names),
            scored_selection_statuses=dict(Counter(r['selection']['status'] for r in rows if r['status'] == 'evaluated')))
    write_json(output, dict(schema='dvbfixer.prospective-scoring-followup.v1',
                            prospectively_enrolled=True, policy_retuned=False,
                            enrollment=enrollment, enrollment_sha256=sha(enrollment_raw),
                            mask_plan_sha256=sha((masks/'plan.json').read_bytes()),
                            score_plan_sha256=sha((scores/'plan.json').read_bytes()),
                            comparison_sha256=sha(comparison.read_bytes()),
                            independent_readback_path=str(inner), independent_readback_sha256=sha(inner.read_bytes()),
                            provenance_adapter_sha256=sha(Path(__file__).read_bytes()),
                            preparation_adapter_sha256=score_plan['provenance_adapter_sha256'],
                            parent_enrollment=json.loads(parent_raw),
                            protocol_amended=True, amendment_audit=amendment,
                            own_replay_audit=own_replay,
                            native_audit_sha256=sha(native_audit.read_bytes()),
                            native_control_count=len(native['results']),
                            native_control_passes=sum(r['reference_common_audit']['accepted'] for r in native['results'].values()),
                            per_system=by_system, aggregate=result['aggregate'],
                            changed_choice_reporting=report_changed_choices(score_plan, scores, result['aggregate']['selections_changed']),
                            scientific_qualification=False,
                            label_scope='prospective enrollment/policy follow-up, not a randomized or blinded qualified accuracy benchmark; legacy readback benchmark flag remains false'))
    print(json.dumps(by_system, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('prepare')
    p.add_argument('--enrollment', type=Path, required=True)
    p.add_argument('--comparison', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('summarize')
    for name in ('enrollment', 'masks', 'comparison', 'scores', 'run', 'output',
                 'native-audit', 'amendment-audit', 'own-replay-audit', 'parent-enrollment'):
        p.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'prepare':
        prepare(args.enrollment, args.comparison, args.output)
    else:
        summarize(args.enrollment, args.masks, args.comparison, args.scores, args.run, args.output,
                  native_audit=args.native_audit, amendment_audit=args.amendment_audit,
                  own_replay_audit=args.own_replay_audit, parent_enrollment=args.parent_enrollment)


if __name__ == '__main__':
    main()
