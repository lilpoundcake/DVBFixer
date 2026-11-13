"""Frozen two-state offline local-policy experiment; no new native execution."""
from __future__ import annotations

import argparse
import random
import time
from dataclasses import asdict, replace
from pathlib import Path

from scripts.audit_post_h_1vcu import target_geometry
from scripts.audit_post_h_local_noop import audit
from scripts.audit_thr2_native_capture import read
from scripts.post_hydrogen_local_noop import propose
from scripts.post_hydrogen_public_api import ownership
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify
from scripts.research_legacy_start_controls import compressed
from scripts.research_post_h_1vcu import check as check_original
from scripts.research_post_h_after_repair import input_from_table

ROOT = Path(__file__).resolve().parents[1]


def check(plan):
    verify(plan)
    parent = read(Path(plan['parent_plan']))
    verify(parent)
    check_original(read(Path(parent['parent_plan'])))


def freeze(output):
    parent_path = ROOT/'docs/research/post-hydrogen-after-repair-plan.json'
    parent = read(parent_path)
    verify(parent)
    check_original(read(Path(parent['parent_plan'])))
    evidence_path = ROOT/'docs/research/post-hydrogen-after-repair-evidence.json'
    evidence = read(evidence_path)
    if evidence['native_gate_both_admitted'] or [c['status'] for c in evidence['cases']] != ['proposed', 'refused']:
        raise ValueError('original failed strict gate must remain preserved')
    private = ROOT/'tmp/pdbfixer-research/post-hydrogen-after-repair-01'
    cases, inputs = [], {str(p):sha(p.read_bytes()) for p in
        (parent_path, evidence_path, ROOT/'docs/plans/post-hydrogen-local-noop.md')}
    for i, old in enumerate(parent['cases']):
        path = private/f'case-{i:02d}-repair.json.gz'
        digest = sha(path.read_bytes())
        if digest != evidence['artifacts'][path.name]:
            raise ValueError('retained repaired state digest differs')
        inputs[str(path)] = digest
        inputs[old['journal']] = sha(Path(old['journal']).read_bytes())
        cases.append(dict(index=i, arm=old['arm'], repair=str(path), journal=old['journal']))
    sources = [ROOT/'scripts'/name for name in ('post_hydrogen_local_noop.py',
        'audit_post_h_local_noop.py', 'research_post_h_local_noop.py')]
    sources.append(ROOT/'tests/test_post_hydrogen_local_noop.py')
    write(output, dict(schema='dvbfixer.private-post-h-local-noop-plan.v1', parent_plan=str(parent_path),
        inputs=inputs, source_hashes={str(p):sha(p.read_bytes()) for p in sources}, cases=cases,
        policy=dict(maximum_seconds=120, policy_calls=2, audits=2, new_contexts=0,
                    new_repairs=0, new_h_samples=0, new_minimizations=0,
                    native_continuation_requires_separate_freeze=True), production_changed=False))


def run(plan_path, root, output):
    import numpy as np

    plan = read(plan_path)
    check(plan)
    root.mkdir(parents=True, exist_ok=False)
    (root/'executed-plan.json').write_bytes(plan_path.read_bytes())
    rows, start = [], time.monotonic()
    for case in plan['cases']:
        row = dict(index=case['index'], arm=case['arm'], status='not-executed-budget')
        py_rng, np_rng = random.getstate(), np.random.get_state()
        if time.monotonic()-start <= plan['policy']['maximum_seconds']:
            try:
                saved, journal = read(Path(case['repair'])), read(Path(case['journal']))
                raw, before = input_from_table(saved['raw']), input_from_table(journal['before'])
                if saved['raw'] != journal['raw'] or saved['repaired_nm'] is None:
                    raise ValueError('retained raw/repair binding mismatch')
                owned, _ = ownership(before, raw, frozenset(), journal['reference_requested'])
                value = replace(owned, positions_nm=tuple(tuple(p) for p in saved['repaired_nm']))
                if any(p != q for a, p, q in zip(value.atoms, raw.positions_nm, value.positions_nm, strict=True)
                       if a.element not in ('H', 'D')):
                    raise ValueError('retained repair moved heavy atoms')
                proposal = propose(value)
                # Retain proposal before independent audit; an audit refusal must not lose it.
                compressed(root/f"case-{case['index']:02d}-proposal.json.gz", dict(
                    input=dict(atoms=[asdict(a) for a in value.atoms], positions_nm=value.positions_nm,
                        bonds=value.bonds, new_hydrogens=[asdict(a) for a in sorted(value.new_hydrogens)]),
                    proposal=dict(result=asdict(proposal.result), warnings=[asdict(w) for w in proposal.warnings],
                        original_new_hydrogens=[asdict(a) for a in sorted(proposal.original_new_hydrogens)],
                        intervention_hydrogens=[asdict(a) for a in sorted(proposal.intervention_hydrogens)])))
                row.update(status='proposal-awaiting-audit', warnings=[asdict(w) for w in proposal.warnings],
                    target_before=target_geometry(saved['raw']['atoms'], value.positions_nm),
                    target_after=target_geometry(saved['raw']['atoms'], proposal.result.positions_nm))
                independent = audit(value, proposal)
                row.update(status='proposed-and-audited', audit=independent)
            except Exception as error:
                row.update(status='refused-or-audit-error', error_type=type(error).__name__, error=str(error))
        current_np = np.random.get_state()
        row['python_rng_unchanged'] = random.getstate() == py_rng
        row['numpy_rng_unchanged'] = (np_rng[0] == current_np[0] and np.array_equal(np_rng[1], current_np[1])
                                      and np_rng[2:] == current_np[2:])
        rows.append(row)
        write(root/f"case-{case['index']:02d}-result.json", row)
    check(plan)
    write(output, dict(schema='dvbfixer.private-post-h-local-noop-evidence.v1',
        plan_sha256=sha(plan_path.read_bytes()), cases=rows,
        both_states_proposed_and_audited=all(r['status'] == 'proposed-and-audited'
            and r['python_rng_unchanged'] and r['numpy_rng_unchanged'] for r in rows) and len(rows) == 2,
        private_artifacts={p.name:sha(p.read_bytes()) for p in sorted(root.glob('*')) if p.is_file()},
        new_contexts=0, new_repairs=0, new_h_samples=0, new_minimizations=0,
        old_strict_gate_changed=False, production_ready=False,
        limits='Two retained-state offline proposals only, not native/live qualification. Compatible out-of-range bonds remain defects; old strict refusal and incomplete gradient attestation remain unchanged.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('plan', 'run'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--plan', type=Path)
    parser.add_argument('--root', type=Path)
    args = parser.parse_args()
    if args.command == 'plan':
        freeze(args.output.resolve())
    else:
        run(args.plan.resolve(), args.root.resolve(), args.output.resolve())
