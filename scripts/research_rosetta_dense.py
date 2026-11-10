"""Separate private nine-chi study on every frozen close-contact case.

Only the explicitly declared sampling policy changes. Existing production and
three-staggered audits are untouched. Measured-chi geometry is supplemented by
strict new-grid membership; anchored observations still admit only one closure.
"""
from __future__ import annotations

import argparse
import copy
import json
from dataclasses import asdict, replace
from pathlib import Path

from scripts.audit_rosetta_joint import norm
from scripts.branched_stereo_oracle import rdkit_labels
from scripts.research_rosetta_joint import sha, verify
from scripts.research_rosetta_repack import write_json

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.pdb_component_reconstruction import _component_ref, _parse_atom_lines
from dvbfixer.protein_completion.anchored_branch import (
    audit_anchored_branch,
    complete_anchored_branch,
)
from dvbfixer.protein_completion.branched import BRANCHES, TETRAHEDRAL_ANGLE
from dvbfixer.protein_completion.branched_source import _thr_beta_checks
from dvbfixer.protein_completion.cohort import MaterializedCohort, _enumerate_joint
from dvbfixer.protein_completion.contracts import (
    Addition,
    Attempt,
    Audit,
    Candidate,
    CompletionResult,
    Observation,
)
from dvbfixer.protein_completion.geometry import extend, l_cb
from dvbfixer.protein_completion.mixed_cohort import plan_mixed_cohort
from dvbfixer.protein_completion.source import (
    _global_ca,
    _insert_candidates,
    _layout,
    _readback_additions,
)
from dvbfixer.protein_completion.validation import _audit_candidate

CHIS = (-60., 60., 180., -80., -40., 40., 80., 160., -160.)
MAX_COMBINATIONS = 81
PARENT = Path('tmp/pdbfixer-research/contact-plan-01')


def audit_target(target, candidate, *, rounded=False):
    if target.mode == 'anchored':
        return audit_anchored_branch(target.request, candidate, pdb_roundtrip=rounded)
    if (isinstance(candidate.torsion_degrees, bool) or candidate.torsion_degrees not in CHIS
            or candidate.rotamer_evidence is not None):
        return Audit(False, ('dense_grid_or_authority_mismatch',))
    return _audit_candidate(target.request, candidate, measured_chi=True, pdb_roundtrip=rounded)


def joint_audits(plan, candidates, *, rounded=False):
    if not isinstance(candidates, tuple) or len(candidates) != len(plan.targets):
        raise ValueError('exact immutable target tuple required')
    audits = []
    for i, (target, candidate) in enumerate(zip(plan.targets, candidates, strict=True)):
        others = tuple(Observation(a.atom, a.element, a.xyz) for j, c in enumerate(candidates) if j != i for a in c.additions)
        request = replace(target.request, observations=target.request.observations+others)
        audits.append(audit_target(replace(target, request=request), replace(candidate, request_digest=request.digest), rounded=rounded))
    return tuple(audits)


def propose(target):
    if target.mode == 'anchored':
        return complete_anchored_branch(target.request)
    request = target.request
    points = {o.atom.name: o.xyz for o in request.observations if o.atom.component == request.residue}
    if 'CB' not in points:
        points['CB'] = l_cb(points['N'], points['CA'], points['C'])
    branch = BRANCHES[request.residue_name]
    attempts = []
    for chi in CHIS:
        proposal = dict(points)
        proposal[branch.first] = extend(points['N'], points['CA'], points['CB'], branch.first_length, TETRAHEDRAL_ANGLE, chi)
        proposal[branch.second] = extend(points['N'], points['CA'], points['CB'], branch.second_length, TETRAHEDRAL_ANGLE, chi+branch.second_chi_offset)
        additions = tuple(Addition(ExactAtomRef(request.residue, n), n[0], proposal[n]) for n in sorted(request.missing_atoms))
        candidate = Candidate(request.digest, additions, chi)
        attempts.append(Attempt(candidate, audit_target(target, candidate)))
    passing = [i for i, a in enumerate(attempts) if a.audit.accepted]
    selected = max(passing, key=lambda i: min(0., attempts[i].audit.minimum_nonbonded_clearance_angstrom or 0.)) if passing else None
    return CompletionResult('completed' if passing else 'rejected', request.digest, 'private_nine_chi_control', tuple(attempts), selected,
                            revision='private-nine-chi-v1', sampling_policy='three original chi plus +/-20 degrees; no empirical probabilities')


def materialize(source, plan, candidates):
    if not all(audit_target(t, c).accepted for t, c in zip(plan.targets, candidates, strict=True)):
        raise ValueError('source-only geometry rejected')
    if not all(a.accepted for a in joint_audits(plan, candidates)):
        raise ValueError('raw joint geometry rejected')
    requests = tuple(t.request for t in plan.targets)
    pdb = _insert_candidates(source, requests, candidates)
    added = _readback_additions(source, pdb, tuple(a for c in candidates for a in c.additions))
    rounded = tuple(replace(c, additions=tuple(a for a in added if a.atom.component == t.request.residue))
                    for t, c in zip(plan.targets, candidates, strict=True))
    if not all(a.accepted for a in joint_audits(plan, rounded, rounded=True)):
        raise ValueError('rounded joint geometry rejected')
    _, atoms, _ = _layout(pdb)
    _global_ca(atoms, allow_missing_cb=None)
    _thr_beta_checks(atoms, pending_target=None)
    parsed = _parse_atom_lines(pdb.decode().splitlines(keepends=True))
    for request in requests:
        points = {a.name: a.position for a in parsed if _component_ref(a) == request.residue}
        expected = {'CA': 'S', 'CB': 'R'} if request.residue_name == 'THR' else {'CA': 'S'}
        if rdkit_labels(points, request.residue_name) != expected:
            raise ValueError('independent target CIP mismatch')
    result = MaterializedCohort(pdb, sha(source), sha(pdb), plan.digest, added,
                                tuple((r.residue, r.residue_name) for r in requests))
    dat = result.dat_record()
    dat.description = 'Private nine-chi research; generated heavy atoms only; not production materialization'
    return pdb, dat


def enroll(output):
    parent_raw = (PARENT/'plan.json').read_bytes()
    parent = json.loads(parent_raw)
    verify(parent, PARENT)
    hashes = dict(parent['source_hashes'])
    hashes['scripts/research_rosetta_dense.py'] = sha(Path(__file__).read_bytes())
    hashes['scripts/rosetta_dense_diagnostics.py'] = sha(Path('scripts/rosetta_dense_diagnostics.py').read_bytes())
    hashes['scripts/compare_branched_completion.py'] = sha(Path('scripts/compare_branched_completion.py').read_bytes())
    record = dict(schema='dvbfixer.private-nine-chi-enrollment.v1', parent_plan_path=str(PARENT),
                  parent_plan_sha256=sha(parent_raw), source_hashes=hashes, artifacts={},
                  case_ids=[c['id'] for c in parent['cases']], chis=CHIS, max_combinations=MAX_COMBINATIONS,
                  scope='all previously frozen pairs/masks, no outcome filtering; retrospective systems, prospective new-grid execution',
                  anchored='same exact observed-partner algebraic closure; no new degrees of freedom',
                  geometry='existing measured-chi comparator plus exact declared-grid gate; raw/rounded/stereo/sterics unchanged',
                  score_policy=parent['enrollment']['score_policy'],
                  selection_policy='same full-pose minimum over COMPLETE source/raw/rounded-feasible set; any missing feasible score blocks selection',
                  comparisons=['coarse first-three subset within same unchanged contexts/scores', 'dense independent clearance',
                               'conditional factorized score E(i,j0)+E(i0,j)-E(i0,j0), first source-feasible indices; retain infeasible result'],
                  coupling='complete source-feasible Cartesian full-score finite difference; geometry-rejected cells remain diagnostic only',
                  public_capability=False)
    output.mkdir(parents=True, exist_ok=False)
    write_json(output/'enrollment.json', record)


def prepare(enrollment_root, output):
    raw = (enrollment_root/'enrollment.json').read_bytes()
    enrollment = json.loads(raw)
    verify(enrollment, enrollment_root)
    parent_root = Path(enrollment['parent_plan_path'])
    parent_raw = (parent_root/'plan.json').read_bytes()
    if sha(parent_raw) != enrollment['parent_plan_sha256'] or enrollment['chis'] != list(CHIS):
        raise ValueError('parent or declared grid changed')
    parent = json.loads(parent_raw)
    verify(parent, parent_root)
    output.mkdir(parents=True, exist_ok=False)
    cases, artifacts = [], {}

    def save(name, data):
        path = output/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        artifacts[name] = sha(data)

    for old in parent['cases']:
        case = copy.deepcopy(old)
        name = case['id']
        for filename in ('source.pdb', 'reference.pdb', 'native-control.pdb'):
            save(name+'/'+filename, (parent_root/name/filename).read_bytes())
        source = (output/name/'source.pdb').read_bytes()
        try:
            plan = plan_mixed_cohort(source, tuple(ComponentInstanceRef(**t['residue']) for t in case['targets']))
        except ValueError as error:
            cases.append(dict(case, status='admission_refused', reason=str(error)))
            continue
        targets = tuple(propose(t) for t in plan.targets)
        # The first three paired candidates and anchored closure must reproduce
        # saved old proposals/audits exactly; this comparison does not select.
        old_result = json.loads((parent_root/name/'reconstruction.json').read_bytes())
        for new, previous in zip(targets, old_result['targets'], strict=True):
            if norm([asdict(a) for a in new.attempts[:len(previous['attempts'])]]) != previous['attempts']:
                raise ValueError('coarse subset no longer reproduces old candidates/audits')
        result = _enumerate_joint(targets, lambda c, plan=plan: joint_audits(plan, c), plan.digest,
                                  max_combinations=MAX_COMBINATIONS, completed_reason='private_nine_chi_joint')
        save(name+'/source-plan.json', (json.dumps(asdict(plan), indent=2)+'\n').encode())
        save(name+'/reconstruction.json', (json.dumps(asdict(result), indent=2)+'\n').encode())
        cells = []
        for i, attempt in enumerate(result.attempts):
            candidates = tuple(t.attempts[j].candidate for t, j in zip(targets, attempt.indices, strict=True))
            pdb = _insert_candidates(source, tuple(t.request for t in plan.targets), candidates)
            eligible, reason = False, 'raw_joint_geometry_rejected'
            if attempt.accepted:
                try:
                    checked, dat = materialize(source, plan, candidates)
                except ValueError as error:
                    reason = str(error)
                else:
                    if checked != pdb:
                        raise ValueError('materialization differs from insertion')
                    path = output/name/f'joint-{i}.dat'
                    dat.save(path, verbose=False)
                    artifacts[str(path.relative_to(output))] = sha(path.read_bytes())
                    eligible, reason = True, None
            filename = f'{name}/joint-{i}.pdb'
            save(filename, pdb)
            cells.append(dict(id=f'joint-{i}', arm='own', index=i, indices=list(attempt.indices), pdb=filename,
                              old_source_accepted=attempt.accepted, common_accepted=eligible,
                              raw_joint_audits=norm([asdict(a) for a in attempt.audits]), rounded_refusal=reason))
        cells.extend(c for c in old['cells'] if c['arm'] == 'native_control')
        independent = [r.selected_index for r in targets]
        case.update(status='planned', cells=cells,
                    individual_attempt_counts=[len(r.attempts) for r in targets],
                    individual_feasible_counts=[sum(a.audit.accepted for a in r.attempts) for r in targets],
                    independent_indices=independent,
                    independent_joint_index=next((i for i, a in enumerate(result.attempts) if list(a.indices) == independent), None),
                    raw_clearance_joint_index=result.selected_index)
        cases.append(case)
    if [c['id'] for c in cases] != enrollment['case_ids']:
        raise ValueError('frozen case list changed')
    (output/'executed-enrollment.json').write_bytes(raw)
    write_json(output/'plan.json', dict(schema='dvbfixer.private-nine-chi-score-plan.v1', enrollment=enrollment,
                                      enrollment_sha256=sha(raw), source_hashes=enrollment['source_hashes'],
                                      artifacts=artifacts, cases=cases))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('enroll')
    p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('prepare')
    p.add_argument('--enrollment', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'enroll':
        enroll(args.output)
    else:
        prepare(args.enrollment, args.output)
