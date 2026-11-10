"""Supplement frozen joint readback with independent input-only enrollment audit.

Never generates coordinates or scores; rejected combinations remain retained.
Native coordinate checks are reporting-only, not proposal/selection inputs.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
from collections import Counter
from dataclasses import asdict, replace
from pathlib import Path

from scripts.audit_rosetta_joint import norm, replay
from scripts.prospective_branched_cohort import project
from scripts.research_rosetta_joint import mask_classes, sha, verify
from scripts.research_rosetta_repack import write_json

from dvbfixer.diagnose.geometry import signed_dihedral_degrees
from dvbfixer.domain.structure_identity import ExactAtomRef
from dvbfixer.pdb_component_reconstruction import _component_ref, _parse_atom_lines
from dvbfixer.protein_completion.branched import BRANCHES
from dvbfixer.protein_completion.contracts import Addition, Candidate, Observation
from dvbfixer.protein_completion.mixed_cohort import _audit, _joint, materialize_mixed_cohort
from dvbfixer.protein_completion.validation import _audit_candidate


def measured_native_audits(plan, candidates: tuple, *, rounded: bool) -> list[dict]:
    """Existing comparator geometry at measured chi; never generator admission.

    Context includes the other native additions. Only the grid/metadata requirement
    is replaced by measured chi, just as in the existing template/native comparator;
    sterics, bond/angle/stereo thresholds and observations are unchanged.
    """
    audits = []
    for i, (target, candidate) in enumerate(zip(plan.targets, candidates, strict=True)):
        request = target.request
        others = tuple(Observation(a.atom, a.element, a.xyz)
                       for j, c in enumerate(candidates) if j != i for a in c.additions)
        contextual = replace(request, observations=request.observations+others)
        points = {o.atom.name: o.xyz for o in request.observations if o.atom.component == request.residue}
        points.update({a.atom.name: a.xyz for a in candidate.additions})
        first = BRANCHES[request.residue_name].first
        chi = signed_dihedral_degrees(*(points[n] for n in ('N', 'CA', 'CB', first)))
        measured = replace(candidate, request_digest=contextual.digest, torsion_degrees=chi)
        audits.append(asdict(_audit_candidate(contextual, measured, anchored_branch=target.mode == 'anchored',
                                             measured_chi=True, pdb_roundtrip=rounded)))
    return norm(audits)


def audit_enrollment(root: Path) -> dict:
    raw = (root/'enrollment.json').read_bytes()
    record = json.loads(raw)
    verify(record, root)
    expected_policy = dict(metric='observed CA distance only; no sidechain coordinates',
                           max_ca_distance_angstrom=8., min_file_order_separation=3,
                           per_system_pair_limit=4, tie_break='file-order index pair',
                           pair_reuse_allowed=True, outcome_filtering=False, replace_empty_pools=False)
    if record['selection'] != expected_policy:
        raise ValueError('selector policy differs from declared contact protocol')
    declared = []
    for system in ('1CRN', '1AKI'):
        original = Path(f'tests/fixtures/protein_context/scoring_followup/{system}.pdb').read_bytes()
        if sha(original) != record['original_source_hashes'][system]:
            raise ValueError('original fixture changed')
        reference = (root/(system+'-reference.pdb')).read_bytes()
        if project(original) != reference:
            raise ValueError('projection changed input records')
        atoms = _parse_atom_lines(reference.decode().splitlines(keepends=True))
        grouped = {}
        for atom in atoms:
            grouped.setdefault(_component_ref(atom), []).append(atom)
        order = list(grouped)
        ranking = []
        for i, j in itertools.combinations(range(len(order)), 2):
            if j-i < 3 or order[i].chain_id != order[j].chain_id:
                continue
            left, right = grouped[order[i]], grouped[order[j]]
            if left[0].residue_name not in ('VAL', 'THR') or right[0].residue_name not in ('VAL', 'THR'):
                continue
            ca_left = [a.position for a in left if a.name == 'CA']
            ca_right = [a.position for a in right if a.name == 'CA']
            if len(ca_left) != 1 or len(ca_right) != 1:
                raise ValueError('ambiguous selector CA geometry')
            distance = math.dist(ca_left[0], ca_right[0])
            if distance <= 8.:
                ranking.append(dict(ca_distance_angstrom=distance, file_order_indices=[i, j],
                                    targets=[dict(residue_name=rows[0].residue_name, residue=asdict(ref))
                                             for rows, ref in ((left, order[i]), (right, order[j]))]))
        ranking.sort(key=lambda r: (r['ca_distance_angstrom'], r['file_order_indices']))
        if ranking != record['input_geometry_rankings'][system]:
            raise ValueError('independent CA-only ranking differs')
        for rank, pair in enumerate(ranking[:4]):
            for label, masks in mask_classes(tuple(t['residue_name'] for t in pair['targets'])):
                identifier = f'{system}-contact{rank}-{label}'
                matches = [c for c in record['cases'] if c['id'] == identifier]
                if len(matches) != 1:
                    raise ValueError('enrolled case omitted or duplicated')
                actual = matches[0]
                targets = [dict(t, mask=list(mask)) for t, mask in zip(pair['targets'], masks, strict=True)]
                if actual['targets'] != targets or actual['selection_geometry'] != pair or actual['selection_rank'] != rank:
                    raise ValueError('pair/mask declaration differs from input-only ranking')
                excluded = {(t['residue']['chain_id'], t['residue']['sequence_number'], t['residue']['insertion_code'], atom)
                            for t in targets for atom in t['mask']}
                expected = b''.join(line for line in reference.splitlines(keepends=True) if not (
                    line.startswith(b'ATOM  ') and (line[21:22].decode(), line[22:26].strip().decode(),
                                                   line[26:27].strip().decode(), line[12:16].strip().decode()) in excluded))
                if expected != (root/actual['source']).read_bytes():
                    raise ValueError('masked source does not preserve unselected records')
                declared.append(identifier)
    if declared != [c['id'] for c in record['cases']]:
        raise ValueError('case list omitted, reordered or added an enrolled case')
    return dict(enrollment_sha256=sha(raw), selection=expected_policy,
                eligible_pairs={s: len(r) for s, r in record['input_geometry_rankings'].items()},
                audited_cases=len(declared), original_and_masked_records_exact=True)


def audit(enrollment_root: Path, plan_root: Path, joint_audit: Path, output: Path) -> None:
    enrollment = audit_enrollment(enrollment_root)
    plan_raw = (plan_root/'plan.json').read_bytes()
    plan = json.loads(plan_raw)
    verify(plan, plan_root)
    raw = joint_audit.read_bytes()
    joint = json.loads(raw)
    if (joint['plan_sha256'] != sha(plan_raw) or plan['enrollment_sha256'] != enrollment['enrollment_sha256']
            or joint['audit_script_sha256'] != sha(Path('scripts/audit_rosetta_joint.py').read_bytes())
            or sha(Path(joint['raw_score_path']).read_bytes()) != joint['raw_score_sha256']):
        raise ValueError('joint readback binding differs')
    native_controls, rejected, failed_attempts = {}, [], []
    for case in plan['cases']:
        if case['status'] != 'planned':
            continue
        root = plan_root/case['id']
        source = (root/'source.pdb').read_bytes()
        saved = json.loads((root/'reconstruction.json').read_bytes())
        bound, result = replay(source, case, saved)
        reference = _parse_atom_lines((root/'reference.pdb').read_text().splitlines(keepends=True))
        native = []
        for target, individual in zip(bound.targets, result.targets, strict=True):
            request = target.request
            points = {a.name: a.position for a in reference if _component_ref(a) == request.residue}
            native.append(Candidate(request.digest, tuple(Addition(ExactAtomRef(request.residue, n), n[0], points[n])
                                                         for n in request.missing_atoms), None))
            failed_attempts.append(dict(case=case['id'], residue=asdict(request.residue),
                                        feasible_count=sum(a.audit.accepted for a in individual.attempts),
                                        attempts=[dict(index=i, accepted=a.audit.accepted, findings=a.audit.findings)
                                                  for i, a in enumerate(individual.attempts)]))
        raw_native = _joint(bound, tuple(native))
        measured_raw = measured_native_audits(bound, tuple(native), rounded=False)
        measured_rounded = measured_native_audits(bound, tuple(native), rounded=True)
        try:
            materialize_mixed_cohort(source, bound, tuple(native))
        except ValueError as error:
            generator_ok, generator_refusal = False, str(error)
        else:
            generator_ok, generator_refusal = True, None
        native_controls[case['id']] = dict(
            generator_policy_audits=norm([asdict(a) for a in raw_native]),
            generator_materialization_accepted=generator_ok, generator_refusal=generator_refusal,
            measured_geometry_raw=measured_raw, measured_geometry_rounded=measured_rounded,
            measured_geometry_accepted=all(a['accepted'] for a in measured_raw+measured_rounded),
            reporting_only=True, measured_chi_is_not_generator_admission=True)
        for cell, attempt in zip([c for c in case['cells'] if c['arm'] == 'own'], result.attempts, strict=True):
            if attempt.accepted:
                continue
            chosen = [t.attempts[i].candidate for t, i in zip(result.targets, attempt.indices, strict=True)]
            # Each proposal passes its original source-only audit; only adding
            # the other proposal can change this joint result.
            if not all(_audit(t, c).accepted for t, c in zip(bound.targets, chosen, strict=True)):
                raise ValueError('claimed joint-only rejection already fails source audit')
            distances = [dict(left=asdict(a.atom), right=asdict(b.atom), distance_angstrom=math.dist(a.xyz, b.xyz))
                         for a in chosen[0].additions for b in chosen[1].additions]
            rejected.append(dict(case=case['id'], index=cell['index'], indices=list(attempt.indices),
                                 source_only_audits_pass=True, joint_audits=norm([asdict(a) for a in attempt.audits]),
                                 cross_generated_distances=distances))
    findings = Counter(f.split(':')[0] for t in failed_attempts for a in t['attempts'] if not a['accepted'] for f in a['findings'])
    write_json(output, dict(schema='dvbfixer.exact-contact-audit.v1', enrollment_audit=enrollment,
                            joint_audit=joint, joint_audit_sha256=sha(raw),
                            native_geometry_controls=native_controls, target_pool_diagnostics=failed_attempts,
                            source_rejection_finding_counts=dict(findings), joint_only_rejections=rejected,
                            audit_script_sha256=sha(Path(__file__).read_bytes()),
                            public_capability=False, coupled_packing_advantage_demonstrated=False,
                            scope='input-backbone-selected cases on previously studied systems; correlated masks; no production qualification'))
    print(json.dumps(dict(enrollment=enrollment, aggregate=joint['aggregate'],
                          native_measured_geometry_passes=sum(n['measured_geometry_accepted'] for n in native_controls.values()),
                          generator_policy_native_materializations=sum(n['generator_materialization_accepted'] for n in native_controls.values()),
                          verified_joint_only_rejections=len(rejected), source_rejection_findings=dict(findings)), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--enrollment', type=Path, required=True)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--joint-audit', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    audit(args.enrollment, args.plan, args.joint_audit, args.output)
