"""Independent NumPy saved-return geometry audit; no proposal or native sampler."""
from __future__ import annotations

import argparse
from dataclasses import asdict
from pathlib import Path

from scripts.audit_thr2_native_capture import read
from scripts.research_dunbrack_prior import sha, write
from scripts.research_post_hydrogen_guard import check, load_return

from dvbfixer.ffutils.variants import VARIANT_TO_PARENT
from dvbfixer.model.diffusion.protein_atoms import CANONICAL_HEAVY_ATOMS


def audit(plan_path, evidence_path, root, output):
    import numpy as np

    plan, evidence = read(plan_path), read(evidence_path)
    check(plan)
    if evidence['plan_sha256'] != sha(plan_path.read_bytes()) or (root/'executed-plan.json').read_bytes() != plan_path.read_bytes():
        raise ValueError('static plan binding mismatch')
    if len(evidence['cases']) != len(plan['cases']):
        raise ValueError('incomplete static enrollment')
    for name, digest in evidence['artifacts'].items():
        if sha((root/name).read_bytes()) != digest:
            raise ValueError('proposal artifact changed')
    reports = []
    for index, (case, result) in enumerate(zip(plan['cases'], evidence['cases'], strict=True)):
        if result['index'] != index:
            raise ValueError('case order mismatch')
        report = dict(index=index, label=case['label'], mode=case['mode'], kind=case['kind'],
                      status=result['status'], requested_ff=case['requested_ff'])
        if result['status'] == 'proposed':
            value = load_return(case)  # Shared source admission, not the proposal rule.
            saved = read(root/result['proposal_file'])
            if (saved['atoms'] != [asdict(a) for a in value.atoms]
                    or saved['initial_nm'] != [list(p) for p in value.positions_nm]
                    or saved['bonds'] != [list(b) for b in value.bonds]
                    or saved['new_indices'] != [i for i, a in enumerate(value.atoms) if a.identity in value.new_hydrogens]):
                raise ValueError('proposal input not the admitted native return')
            before, after = np.asarray(saved['initial_nm']), np.asarray(saved['proposed_nm'])
            moved = np.where(np.any(before != after, axis=1))[0].tolist()
            if moved != saved['moved_indices'] or not set(moved) <= set(saved['new_indices']):
                raise ValueError('observation/moved mask mismatch')
            adjacent = [[] for _ in value.atoms]
            for a, b in value.bonds:
                adjacent[a].append(b)
                adjacent[b].append(a)
            supported, expected_moved, metrics = [], [], []
            for h in saved['new_indices']:
                parent, = adjacent[h]
                center = value.atoms[parent]
                inventory = CANONICAL_HEAVY_ATOMS.get(VARIANT_TO_PARENT.get(center.residue_name, center.residue_name), ())
                heavy = [i for i in adjacent[parent] if value.atoms[i].element not in ('H', 'D')]
                hs = [i for i in adjacent[parent] if value.atoms[i].element in ('H', 'D')]
                if (center.element != 'C' or len(heavy) != 3 or hs != [h] or center.identity.name not in inventory
                        or any(value.atoms[i].identity.component != center.identity.component for i in [h]+heavy)
                        or any(value.atoms[i].identity.name not in inventory for i in heavy)):
                    continue
                supported.append(h)
                a, b = before[heavy]-before[parent], before[heavy]-before[h]
                hv, fv = float(np.linalg.det(a)), float(np.linalg.det(b))
                radius = float(np.linalg.norm(before[h]-before[parent]))
                lengths = np.linalg.norm(a, axis=1)
                if np.any(lengths < 0.10) or np.any(lengths > 0.20) or not 0.05 <= radius <= 0.15:
                    raise ValueError('accepted outside-distance center')
                if abs(hv/np.prod(lengths)) <= 1e-4 or abs(fv/np.prod(np.linalg.norm(b, axis=1))) <= 1e-4:
                    raise ValueError('accepted degenerate center')
                if hv*fv < 0:
                    expected_moved.append(h)
                    raw = np.sum(before[parent]-before[heavy], axis=0)
                    wanted = before[parent]+radius*raw/np.linalg.norm(raw)
                    c = after[heavy]-after[h]
                    cv = float(np.linalg.det(c))
                    if (np.linalg.norm(after[h]-wanted) > 1e-14 or hv*cv <= 0
                            or abs(cv/np.prod(np.linalg.norm(c, axis=1))) <= 1e-4
                            or abs(float(np.linalg.norm(after[h]-after[parent]))-radius) > 1e-14):
                        raise ValueError('proposal formula/radius/sign mismatch')
                    others = [i for i in range(len(before)) if i not in (h, parent)]
                    if np.min(np.linalg.norm(after[others]-after[h], axis=1)) < 0.05:
                        raise ValueError('proposal gross overlap')
                    metrics.append(dict(atom=asdict(value.atoms[h].identity), radius_nm=radius,
                        heavy_volume_angstrom3=1000*hv, before_four_volume_angstrom3=1000*fv,
                        after_four_volume_angstrom3=1000*cv))
            if (expected_moved != moved or len(supported) != result['supported']
                    or len(saved['new_indices'])-len(supported) != result['unsupported']
                    or [asdict(value.atoms[i].identity) for i in moved] != result['moved']):
                raise ValueError('selectivity/enrollment mismatch')
            report.update(observations_exact=True, geometry_passed=True, supported=len(supported),
                          moved=len(moved), moved_geometry=metrics)
        else:
            report['error'] = result.get('error')
        reports.append(report)
    write(output, dict(schema='dvbfixer.private-post-h-static-audit.v1',
        plan_sha256=sha(plan_path.read_bytes()), evidence_sha256=sha(evidence_path.read_bytes()),
        script_sha256=sha(Path(__file__).read_bytes()), cases=reports,
        contexts=0, minimizations=0, production_changed=False,
        limits='Independent NumPy proposal geometry/preservation; shared source admission. Refusals remain outcomes, not passes. No dynamic claim.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('plan', 'evidence', 'root', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan.resolve(), args.evidence.resolve(), args.root.resolve(), args.output.resolve())
