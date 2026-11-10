"""Input-backbone-only enrollment of close two-target research controls.

Reuse the frozen exact runner without editing it or its historical evidence.
No admission, candidate construction, scoring or held-out sidechain coordinates
are consulted while ranking pairs. Empty/degenerate pools are never replaced.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
from dataclasses import asdict
from pathlib import Path

from scripts.prospective_branched_cohort import project
from scripts.research_rosetta_joint import POLICY, SYSTEMS, mask_classes, sha, source_hashes
from scripts.research_rosetta_repack import write_json

from dvbfixer.domain.structure_identity import ComponentInstanceRef

PAIR_LIMIT = 4
MAX_CA_DISTANCE = 8.0
MIN_FILE_ORDER_SEPARATION = 3


def ranked_pairs(reference: bytes) -> list[dict]:
    """Use CA only, with file-order ties; retain the whole eligible ranking."""
    residues = {}
    for line in reference.decode('ascii').splitlines():
        if not line.startswith('ATOM  '):
            continue
        key = (line[21], line[22:26].strip(), line[26].strip())
        row = residues.setdefault(key, dict(name=line[17:20], order=len(residues), ca=[]))
        if line[12:16].strip() == 'CA':
            row['ca'].append(tuple(float(line[i:i+8]) for i in (30, 38, 46)))
    candidates = []
    for (a, left), (b, right) in itertools.combinations(residues.items(), 2):
        if a[0] != b[0] or left['name'] not in ('VAL', 'THR') or right['name'] not in ('VAL', 'THR'):
            continue
        if right['order']-left['order'] < MIN_FILE_ORDER_SEPARATION:
            continue
        # Fail on ambiguous selector geometry, rather than using occupancy or an
        # outcome-based altloc choice. Scientific admission remains downstream.
        if len(left['ca']) != 1 or len(right['ca']) != 1:
            raise ValueError('selector requires one observed CA per eligible identity')
        distance = math.dist(left['ca'][0], right['ca'][0])
        if not math.isfinite(distance):
            raise ValueError('nonfinite selector geometry')
        if distance <= MAX_CA_DISTANCE:
            candidates.append(dict(ca_distance_angstrom=distance,
                                   file_order_indices=[left['order'], right['order']],
                                   targets=[dict(residue_name=r['name'], residue=asdict(ComponentInstanceRef(1, *k)))
                                            for k, r in ((a, left), (b, right))]))
    return sorted(candidates, key=lambda p: (p['ca_distance_angstrom'], p['file_order_indices']))


def enroll(output: Path) -> None:
    policy_raw = POLICY.read_bytes()
    policy = json.loads(policy_raw)
    for filename, digest in policy['source_hashes'].items():
        if sha(Path(filename).read_bytes()) != digest:
            raise ValueError('frozen score policy implementation changed')
    payloads, cases, rankings, originals = {}, [], {}, {}
    for system in SYSTEMS:
        original = Path(f'tests/fixtures/protein_context/scoring_followup/{system}.pdb').read_bytes()
        originals[system] = sha(original)
        reference = project(original)
        payloads[system+'-reference.pdb'] = reference
        rankings[system] = ranked_pairs(reference)
        for rank, pair in enumerate(rankings[system][:PAIR_LIMIT]):
            for label, masks in mask_classes(tuple(t['residue_name'] for t in pair['targets'])):
                identifier = f'{system}-contact{rank}-{label}'
                removed = {(t['residue']['chain_id'], t['residue']['sequence_number'], t['residue']['insertion_code']): set(mask)
                           for t, mask in zip(pair['targets'], masks, strict=True)}
                source = b''.join(line for line in reference.splitlines(keepends=True) if not (
                    line.startswith(b'ATOM  ') and line[12:16].strip().decode() in removed.get(
                        (line[21:22].decode(), line[22:26].strip().decode(), line[26:27].strip().decode()), set())))
                filename = identifier+'/source.pdb'
                payloads[filename] = source
                cases.append(dict(id=identifier, system=system, mask_class=label,
                                  targets=[dict(t, mask=mask) for t, mask in zip(pair['targets'], masks, strict=True)],
                                  source=filename, reference=system+'-reference.pdb',
                                  selection_rank=rank, selection_geometry=pair))
    hashes = source_hashes()
    hashes['scripts/research_rosetta_contacts.py'] = sha(Path(__file__).read_bytes())
    hashes['scripts/audit_rosetta_joint.py'] = sha(Path('scripts/audit_rosetta_joint.py').read_bytes())
    output.mkdir(parents=True, exist_ok=False)
    for filename, data in payloads.items():
        path = output/filename
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
    write_json(output/'enrollment.json', dict(
        schema='dvbfixer.exact-joint-contact-enrollment.v1', cases=cases,
        artifacts={name: sha(data) for name, data in payloads.items()}, source_hashes=hashes,
        original_source_hashes=originals, single_score_policy_sha256=sha(policy_raw), score_policy=policy,
        selection=dict(metric='observed CA distance only; no sidechain coordinates',
                       max_ca_distance_angstrom=MAX_CA_DISTANCE,
                       min_file_order_separation=MIN_FILE_ORDER_SEPARATION,
                       per_system_pair_limit=PAIR_LIMIT, tie_break='file-order index pair',
                       pair_reuse_allowed=True, outcome_filtering=False, replace_empty_pools=False),
        input_geometry_rankings=rankings, max_targets=2, max_combinations=9,
        candidate_gate='unchanged individual/raw-joint/rounded/global stereo guards; no rescue',
        selection_policy='frozen complete feasible-pool full-pose ref2015 minimum; native reporting only',
        coupling_diagnostic='frozen complete Cartesian finite difference; first feasible indices, not native',
        scope='deliberately close backbone pairs in previously studied systems; correlated masks; not independent-family qualification',
        public_capability=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    enroll(parser.parse_args().output)
