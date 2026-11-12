"""Independent NumPy geometry/preservation readback; never proposes or minimizes."""
from __future__ import annotations

import argparse
from dataclasses import asdict
from pathlib import Path

from scripts.audit_thr2_native_capture import read
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify


def audit(plan_path, evidence_path, root, output):
    import numpy as np
    import openmm
    from openmm.unit import nanometer
    from scripts.qualify_hydrogen_initialization_guard import admit

    from dvbfixer.ffutils.variants import VARIANT_TO_PARENT
    from dvbfixer.model.diffusion.protein_atoms import CANONICAL_HEAVY_ATOMS

    plan, evidence = read(plan_path), read(evidence_path)
    verify(plan)
    if evidence['plan_sha256'] != sha(plan_path.read_bytes()) or (root/'executed-plan.json').read_bytes() != plan_path.read_bytes():
        raise ValueError('plan binding mismatch')
    if len(evidence['cases']) != len(plan['cases']):
        raise ValueError('incomplete enrollment')
    for name, digest in evidence['artifacts'].items():
        if sha((root/name).read_bytes()) != digest:
            raise ValueError('proposal artifact changed')
    records = []
    for expected, case in zip(plan['cases'], evidence['cases'], strict=True):
        if any(case[k] != expected[k] for k in ('kind', 'label', 'slots', 'h_calls', 'requested_ff')):
            raise ValueError('case enrollment changed')
        if len(case['results']) != len(expected['slots']):
            raise ValueError('missing slot outcome')
        for result in case['results']:
            report = dict(index=case['index'], slot=result['slot'], status=result['status'])
            if result['status'] == 'proposed':
                saved = read(root/result['proposal_file'])
                slot = expected['slots'][result['slot']]
                source = Path(expected['folder'])
                prefix = slot['prefix']
                source_positions = openmm.XmlSerializer.deserialize((source/f'{prefix}-state.xml').read_text()).getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
                value = admit(read(source/f'{prefix}-atoms.json.gz')['atoms'],
                    read(source/f"h-{slot['h_call']}-before.json.gz"), source_positions,
                    read(source/slot['bonds'])['bonds'],
                    openmm.XmlSerializer.deserialize((source/f'{prefix}-system.xml').read_text()))
                if (saved['atoms'] != [asdict(a) for a in value.atoms]
                        or saved['initial_nm'] != source_positions
                        or saved['bonds'] != [list(b) for b in value.bonds]
                        or saved['new_indices'] != [i for i, a in enumerate(value.atoms) if a.identity in value.new_hydrogens]):
                    raise ValueError('proposal input differs from hash-bound native source')
                before, after = np.asarray(saved['initial_nm']), np.asarray(saved['proposed_nm'])
                moved = np.where(np.any(before != after, axis=1))[0].tolist()
                if moved != saved['moved_indices'] or not set(moved) <= set(saved['new_indices']):
                    raise ValueError('observation/moved mask mismatch')
                adjacent = [[] for _ in saved['atoms']]
                for a, b in saved['bonds']:
                    adjacent[a].append(b)
                    adjacent[b].append(a)
                independent_supported, independent_moved = [], []
                for h in saved['new_indices']:
                    parent, = adjacent[h]
                    center = value.atoms[parent]
                    inventory = CANONICAL_HEAVY_ATOMS.get(VARIANT_TO_PARENT.get(center.residue_name, center.residue_name), ())
                    heavy = [i for i in adjacent[parent] if value.atoms[i].element not in ('H', 'D')]
                    hs = [i for i in adjacent[parent] if value.atoms[i].element in ('H', 'D')]
                    if (center.element != 'C' or len(heavy) != 3 or hs != [h]
                            or center.identity.name not in inventory
                            or any(value.atoms[i].identity.component != center.identity.component for i in [h]+heavy)
                            or any(value.atoms[i].identity.name not in inventory for i in heavy)):
                        continue
                    independent_supported.append(h)
                    heavy_vectors, four_vectors = before[heavy]-before[parent], before[heavy]-before[h]
                    distances = np.linalg.norm(heavy_vectors, axis=1)
                    radius = float(np.linalg.norm(before[h]-before[parent]))
                    if np.any(distances < 0.10) or np.any(distances > 0.20) or not 0.05 <= radius <= 0.15:
                        raise ValueError('accepted outside-distance center')
                    hv = float(np.linalg.det(heavy_vectors))
                    fv = float(np.linalg.det(four_vectors))
                    if (abs(hv/np.prod(np.linalg.norm(heavy_vectors, axis=1))) <= 1e-4
                            or abs(fv/np.prod(np.linalg.norm(four_vectors, axis=1))) <= 1e-4):
                        raise ValueError('accepted degenerate center')
                    if hv*fv < 0:
                        independent_moved.append(h)
                if (independent_moved != moved or len(independent_supported) != result['supported']
                        or len(saved['new_indices'])-len(independent_supported) != result['unsupported']):
                    raise ValueError('independent enrollment/selectivity mismatch')
                for h in moved:
                    parent, = adjacent[h]
                    heavy = [i for i in adjacent[parent] if i != h]
                    if len(heavy) != 3:
                        raise ValueError('not a three-heavy center')
                    hv = float(np.linalg.det(before[heavy]-before[parent]))
                    fv = float(np.linalg.det(before[heavy]-before[h]))
                    cv = float(np.linalg.det(after[heavy]-after[h]))
                    if hv*fv >= 0 or hv*cv <= 0:
                        raise ValueError('sign conflict not selectively resolved')
                    if abs(cv/np.prod(np.linalg.norm(after[heavy]-after[h], axis=1))) <= 1e-4:
                        raise ValueError('degenerate proposed tetrahedron')
                    old_radius = float(np.linalg.norm(before[h]-before[parent]))
                    radius = float(np.linalg.norm(after[h]-after[parent]))
                    raw = np.sum(before[parent]-before[heavy], axis=0)
                    wanted = before[parent]+old_radius*raw/np.linalg.norm(raw)
                    if abs(old_radius-radius) > 1e-14 or np.linalg.norm(wanted-after[h]) > 1e-14:
                        raise ValueError('radius/formula mismatch')
                    others = [i for i in range(len(before)) if i not in (h, parent)]
                    if np.min(np.linalg.norm(after[others]-after[h], axis=1)) < 0.05:
                        raise ValueError('gross overlap')
                identities = [saved['atoms'][i]['identity'] for i in moved]
                if identities != result['moved']:
                    raise ValueError('public moved identity mismatch')
                report.update(observations_exact=True, independent_geometry_passed=True, moved=len(moved))
            records.append(report)
    write(output, dict(schema='dvbfixer.private-initializer-static-audit.v1',
        plan_sha256=sha(plan_path.read_bytes()), evidence_sha256=sha(evidence_path.read_bytes()),
        script_sha256=sha(Path(__file__).read_bytes()), cases=len(plan['cases']), slots=records,
        contexts=0, minimizations=0, production_changed=False,
        limits='Independent moved-coordinate formula/sign/radius and observation checks only; refusals are retained, not successes or dynamic qualification.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('plan', 'evidence', 'root', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan.resolve(), args.evidence.resolve(), args.root.resolve(), args.output.resolve())
