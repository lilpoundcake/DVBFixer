"""Saved-proposal independent-system readback, without reconstructing or rescoring.

Retain admission/geometry/protocol failures and exact source observations. Native
coordinates are reporting-only, after frozen score-selection replay.
"""
from __future__ import annotations

import argparse
import json
import math
import statistics
from collections import Counter
from dataclasses import asdict
from pathlib import Path

from scripts.audit_rosetta_joint import norm, restore_candidate
from scripts.compare_branched_completion import measure
from scripts.prospective_independent_cohort import project_first, selected_targets
from scripts.prospective_scoring_followup import check_policy
from scripts.research_rosetta_joint import sha, verify
from scripts.research_rosetta_repack import write_json
from scripts.unfiltered_branched_comparison import masks, reporting_scores

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.pdb_component_reconstruction import _component_ref, _parse_atom_lines
from dvbfixer.protein_completion.anchored_branch import audit_anchored_branch
from dvbfixer.protein_completion.anchored_source import plan_anchored_source
from dvbfixer.protein_completion.branched_source import plan_branched_source
from dvbfixer.protein_completion.contracts import Addition
from dvbfixer.protein_completion.validation import audit_candidate


def additions(raw):
    return tuple(Addition(ExactAtomRef(ComponentInstanceRef(**a['atom']['component']), a['atom']['name']),
                          a['element'], tuple(a['xyz'])) for a in raw)


def verify_template(request, proposal, root):
    observed = {o.atom: o for o in request.observations}
    expected = {ExactAtomRef(request.residue, n) for n in request.missing_atoms}
    snapshots = proposal['snapshots']
    if not snapshots or [s['seed'] for s in snapshots] != list(range(1, len(snapshots)+1)) or len(snapshots) > 5:
        raise ValueError('template seed/snapshot inventory differs')
    for snapshot in snapshots:
        if snapshot != json.loads((root/f"template-seed-{snapshot['seed']}.json").read_bytes()):
            raise ValueError('template snapshot payload differs')
        by_atom = {}
        for raw_atom, element, xyz in snapshot['atoms']:
            atom = ExactAtomRef(ComponentInstanceRef(**raw_atom['component']), raw_atom['name'])
            if atom in by_atom:
                raise ValueError('duplicate template atom identity')
            by_atom[atom] = element, tuple(xyz)
        if not observed.keys() <= by_atom.keys() or not expected <= by_atom.keys():
            raise ValueError('template omitted observations/requested atoms')
        movement = max(math.dist(o.xyz, by_atom[atom][1]) for atom, o in observed.items())
        if movement != snapshot['observed_max_movement_angstrom'] or movement > 1e-6:
            raise ValueError('template moved observations')
    actual = additions(proposal['additions'])
    if {a.atom for a in actual} != expected or len(actual) != len(expected):
        raise ValueError('template generated ownership differs')
    if any((a.element, a.xyz) != by_atom[a.atom] for a in actual):
        raise ValueError('template additions differ from last snapshot')
    extras = set(by_atom)-set(observed)-expected
    declared_extras = {ExactAtomRef(ComponentInstanceRef(**a['component']), a['name']) for a in proposal['temporary_extra_atoms']}
    if extras != declared_extras or any(a.name != 'OXT' for a in extras):
        raise ValueError('template rebuilt undeclared nonterminal atoms')
    return actual


def audit(masks_root, comparison_root, score_plan_root, score_root, score_readback, output):
    mask_raw = (masks_root/'plan.json').read_bytes()
    mask_plan = json.loads(mask_raw)
    verify(mask_plan, masks_root)
    enrollment_raw = (masks_root/'enrollment.json').read_bytes()
    if sha(enrollment_raw) != mask_plan['prospective_enrollment_sha256']:
        raise ValueError('input enrollment binding differs')
    enrollment = json.loads(enrollment_raw)
    comparison_raw = (comparison_root/'evidence.json').read_bytes()
    comparison = json.loads(comparison_raw)
    score_plan_raw = (score_plan_root/'plan.json').read_bytes()
    score_plan = json.loads(score_plan_raw)
    check_policy(enrollment, score_plan)
    verify(score_plan, score_plan_root)
    if (comparison['plan_sha256'] != sha(mask_raw) or (comparison_root/'executed-plan.json').read_bytes() != mask_raw
            or score_plan['prospective_enrollment_sha256'] != sha(enrollment_raw)
            or score_plan['prior_comparison_sha256'] != sha(comparison_raw)):
        raise ValueError('comparison/scoring enrollment chain differs')
    readback_raw = score_readback.read_bytes()
    readback = json.loads(readback_raw)
    score_raw = (score_root/'evidence.json').read_bytes()
    if (readback['plan_sha256'] != sha(score_plan_raw) or readback['raw_sha256'] != sha(score_raw)
            or readback['summarizer_sha256'] != sha(Path('scripts/summarize_rosetta_fixed_scoring.py').read_bytes())):
        raise ValueError('independent scoring readback differs')
    expected_ids, systems = [], {}
    for system in enrollment['systems']:
        identifier = system['id']
        original = (masks_root/(identifier+'-original.pdb')).read_bytes()
        tracked = Path(f'tests/fixtures/protein_context/independent_followup/{identifier}.pdb')
        if tracked.read_bytes() != original:
            raise ValueError('tracked original differs from frozen original')
        reference, chain = project_first(original)
        if reference != (masks_root/(identifier+'-reference.pdb')).read_bytes():
            raise ValueError('source projection no longer reproduces')
        targets, counts = selected_targets(reference)
        selected_count = 0
        for name, target in targets:
            for mask in masks(name):
                case_id = f'{identifier}-{name}{target.sequence_number}{target.insertion_code}-'+'-'.join(mask)
                declaration = next(c for c in mask_plan['cases'] if c['id'] == case_id)
                if declaration['target'] != asdict(target) or declaration['mask'] != list(mask):
                    raise ValueError('declared target/mask differs from input-identity selection')
                expected = b''.join(line for line in reference.splitlines(keepends=True) if not (
                    line.startswith(b'ATOM  ') and (line[21:22].decode(), line[22:26].strip().decode(), line[26:27].strip().decode())
                    == (target.chain_id, target.sequence_number, target.insertion_code) and line[12:16].strip().decode() in mask))
                if expected != (masks_root/declaration['source']).read_bytes():
                    raise ValueError('mask altered unselected source bytes')
                expected_ids.append(case_id)
                selected_count += 1
        systems[identifier] = dict(description=system['description'], chain=chain, selected_counts=counts,
                                  declared_cases=selected_count, original_sha256=sha(original))
    if expected_ids != [c['id'] for c in mask_plan['cases']] or set(expected_ids) != set(comparison['results']) or set(expected_ids) != set(readback['cases']):
        raise ValueError('case omitted or added to cohort')
    score_declarations = {c['id']: c for c in score_plan['cases']}
    cases, blockers = {}, []
    for case in mask_plan['cases']:
        name = case['id']
        root = comparison_root/name
        row = comparison['results'][name]
        for filename, digest in row['retained_files'].items():
            if sha((root/filename).read_bytes()) != digest:
                raise ValueError('comparison artifact digest differs')
        source = (root/'input.pdb').read_bytes()
        if source != (masks_root/case['source']).read_bytes() or source != (score_plan_root/name/'source.pdb').read_bytes():
            raise ValueError('arms did not receive identical masked bytes')
        target = ComponentInstanceRef(**case['target'])
        anchored = case['mode'] == 'anchored'
        try:
            source_plan = (plan_anchored_source if anchored else plan_branched_source)(source, target)
            if set(source_plan.request.missing_atoms) != set(case['mask']):
                raise ValueError('declared deletion mask differs from actual missing canonical atoms')
        except ValueError as error:
            if row['status'] != 'shared_admission_refused' or row['reason'] != str(error):
                raise ValueError('admission refusal replay differs')
            cases[name] = dict(status=row['status'], reason=row['reason'])
            continue
        request = source_plan.request
        if norm(asdict(source_plan)) != json.loads((root/'source-plan.json').read_bytes()):
            raise ValueError('source-bound request replay differs')
        own = json.loads((root/'own-attempts.json').read_bytes())
        own_candidates, passes = [], []
        for i, attempt in enumerate(own['attempts']):
            candidate = restore_candidate(attempt['candidate'])
            validator = audit_anchored_branch if anchored else audit_candidate
            audit_result = validator(request, candidate)
            if norm(asdict(audit_result)) != attempt['audit']:
                raise ValueError('saved own-candidate audit differs')
            own_candidates.append(candidate)
            if audit_result.accepted:
                passes.append(i)
        selected_clearance = max(passes, key=lambda i: min(0., own['attempts'][i]['audit']['minimum_nonbonded_clearance_angstrom'] or 0.)) if passes else None
        if selected_clearance != own['selected_index']:
            raise ValueError('original clearance winner differs')
        score_case = score_declarations[name]
        selected = readback['cases'][name]['selection']['selected_index']
        scored_cells = {c['id']: c for c in readback['cases'][name]['cells']}
        # Native/reference coordinates enter only after all proposal audits and
        # the unchanged independent score selector have already been replayed.
        truth = {a.name: a.position for a in _parse_atom_lines((masks_root/case['reference']).read_text().splitlines(keepends=True))
                 if _component_ref(a) == target}
        candidate_reports = []
        for i, candidate in enumerate(own_candidates):
            measured, pdb = measure(source, request, candidate.additions, anchored=anchored)
            cell = next(c for c in score_case['cells'] if c['arm'] == 'own' and c['index'] == i)
            if norm(asdict(measured)) != cell['common_audit'] or (score_plan_root/cell['pdb']).read_bytes() != pdb:
                raise ValueError('own rounded/common audit differs')
            candidate_reports.append(dict(index=i, source_accepted=own['attempts'][i]['audit']['accepted'],
                                          common_accepted=measured.accepted,
                                          rounded_exact_name_rmsd_angstrom=scored_cells[f'own-{i}']['reporting_only_rounded_exact_name_rmsd_angstrom'],
                                          scores=reporting_scores(request, candidate.additions, truth)))
        template_raw = json.loads((root/'template-proposal.json').read_bytes())
        template = verify_template(request, template_raw, root)
        template_audit, template_pdb = measure(source, request, template, anchored=anchored)
        template_cell = next(c for c in score_case['cells'] if c['arm'] == 'template')
        if (template_pdb != (root/'template-proposal.pdb').read_bytes()
                or template_pdb != (score_plan_root/template_cell['pdb']).read_bytes()
                or norm(asdict(template_audit)) != template_cell['common_audit']):
            raise ValueError('template geometry/source readback differs')
        native = tuple(Addition(ExactAtomRef(target, n), n[0], truth[n]) for n in request.missing_atoms)
        native_audit, _ = measure(source, request, native, anchored=anchored)
        clear_report = next((r for r in candidate_reports if r['index'] == selected_clearance), None)
        chosen = next((r for r in candidate_reports if r['index'] == selected), None)
        cases[name] = dict(status='evaluated', fixture=case['fixture'], mode=case['mode'],
                           clearance_index=selected_clearance, score_index=selected,
                           own_clearance=clear_report, own_scored=chosen, own_candidates=candidate_reports,
                           template=dict(common_accepted=template_audit.accepted, common_audit=asdict(template_audit),
                                         rounded_exact_name_rmsd_angstrom=scored_cells['template-0']['reporting_only_rounded_exact_name_rmsd_angstrom'],
                                         scores=reporting_scores(request, template, truth), retained_seed_states=len(template_raw['snapshots'])),
                           native_common_audit=asdict(native_audit), scientific_qualification=False)
    for identifier, system in systems.items():
        subset = [r for r in cases.values() if r.get('fixture') == identifier]
        system.update(own_clearance_passes=sum(r['own_clearance'] is not None and r['own_clearance']['common_accepted'] for r in subset),
                      own_scored_passes=sum(r['own_scored'] is not None and r['own_scored']['common_accepted'] for r in subset),
                      template_passes=sum(r['template']['common_accepted'] for r in subset),
                      native_passes=sum(r['native_common_audit']['accepted'] for r in subset))
    mutual = [r for r in cases.values() if r.get('status') == 'evaluated' and r['own_scored'] is not None
              and r['own_scored']['common_accepted'] and r['template']['common_accepted']]
    changed = [dict(id=n, previous=r['clearance_index'], selected=r['score_index'],
                    previous_scores=r['own_clearance']['scores'], selected_scores=r['own_scored']['scores'])
               for n, r in cases.items() if r.get('status') == 'evaluated' and r['own_clearance'] is not None and r['own_scored'] is not None
               and r['clearance_index'] != r['score_index']]
    stereo_rejections = []
    for name, row in cases.items():
        if row.get('status') != 'evaluated':
            continue
        measured = row['template']['common_audit']
        if 'non_L_or_nearplanar_ca' not in measured['raw']['findings']:
            continue
        volume = measured['raw']['ca_volume_angstrom3']
        ca_cip = measured['rdkit_target_cip'].get('CA')
        classification = ('confirmed_D_CA_by_negative_volume_and_R_CIP' if volume < -.5 and ca_cip == 'R'
                          else 'small_volume_CA_CIP_unassigned' if abs(volume) <= .5 and ca_cip is None
                          else 'other_stereo_disagreement')
        stereo_rejections.append(dict(id=name, classification=classification, raw_ca_volume_angstrom3=volume,
                                      rounded_target_cip=measured['rdkit_target_cip'],
                                      retained_seed_states=row['template']['retained_seed_states'],
                                      common_accepted=measured['accepted'], not_full_legacy_pipeline=True))
    aggregate = dict(declared_cases=len(cases), statuses=dict(Counter(r['status'] for r in cases.values())),
                     scoring=readback['aggregate'], per_system=systems, changed_choices=changed, mutual_pass_cases=len(mutual),
                     mutual_rounded_exact_name_median_own=statistics.median(r['own_scored']['rounded_exact_name_rmsd_angstrom'] for r in mutual) if mutual else None,
                     mutual_rounded_exact_name_median_template=statistics.median(r['template']['rounded_exact_name_rmsd_angstrom'] for r in mutual) if mutual else None,
                     template_ca_stereo_rejections=stereo_rejections,
                     retained_template_seed_states=sum(r['template']['retained_seed_states'] for r in cases.values() if r.get('status') == 'evaluated'),
                     protocol_blockers=blockers)
    write_json(output, dict(schema='dvbfixer.independent-system-audit.v1', aggregate=aggregate, cases=cases,
                            enrollment=enrollment, enrollment_sha256=sha(enrollment_raw), mask_plan_sha256=sha(mask_raw),
                            comparison_sha256=sha(comparison_raw), scoring_plan_sha256=sha(score_plan_raw),
                            scoring_readback_path=str(score_readback), scoring_readback_sha256=sha(readback_raw),
                            audit_script_sha256=sha(Path(__file__).read_bytes()), prospective_enrollment=True,
                            production_enabled=False, scientific_qualification=False,
                            scope=enrollment['qualification_scope']))
    print(json.dumps(aggregate, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for argument in ('masks', 'comparison', 'score-plan', 'scores', 'score-readback', 'output'):
        parser.add_argument('--'+argument, type=Path, required=True)
    args = parser.parse_args()
    audit(args.masks, args.comparison, args.score_plan, args.scores, args.score_readback, args.output)
