"""Exact original H-label repair then object-index mapping to full-System start.

No whole minimization, intervention-arm result reading or alias guessing.
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

from scripts.audit_methyl_diversity_capture import bijection
from scripts.research_dunbrack_prior import sha, write
from scripts.research_h_stage_transport import post_h_repair
from scripts.research_ha_minimizer_path import verify


def mapped_start(rows, bonds, points, mapping):
    bijection(mapping, len(rows))
    repaired, count = post_h_repair(rows, bonds, points)
    result = [None]*len(rows)
    for original, full in enumerate(mapping):
        result[full] = repaired[original]
    return result, count


def preflight(plan_path, capture_root, output):
    import openmm
    from openmm.unit import nanometer

    raw = plan_path.read_bytes()
    record = json.loads(raw)
    verify(record)
    run_raw = (capture_root/'evidence.json').read_bytes()
    run = json.loads(run_raw)
    if run['plan_sha256'] != sha(raw):
        raise ValueError('capture plan differs')
    results = []
    for cell in run['cases']:
        row = {k: cell[k] for k in ('id', 'seed', 'status')}
        if cell['status'] != 'captured':
            results.append(row)
            continue
        source = capture_root/cell['id']/f"seed-{cell['seed']}"
        for name, digest in cell['artifacts'].items():
            if sha((source/name).read_bytes()) != digest:
                raise ValueError('capture changed')
        topology = json.loads((source/'H-topology.json').read_bytes())
        whole = json.loads((source/'whole-topology.json').read_bytes())
        path = json.loads(gzip.decompress((source/'H-path.json.gz').read_bytes()))
        points, repairs = mapped_start(topology['atoms'], topology['bonds'], path['final_nm'], whole['H_to_whole_indices'])
        state = openmm.XmlSerializer.deserialize((source/'whole-state.xml').read_text())
        expected = [list(p) for p in state.getPositions().value_in_unit(nanometer)]
        row.update(exact_original_nm_start=points == expected, repairs=repairs)
        results.append(row)
    write(output, dict(schema='dvbfixer.graph-methyl-transport-preflight.v1', cases=results,
        exact_starts=sum(c.get('exact_original_nm_start', False) for c in results),
        captured=sum(c['status'] == 'captured' for c in results), plan_sha256=sha(raw),
        run_sha256=sha(run_raw), script_sha256=sha(Path(__file__).read_bytes()),
        transport_helper_sha256=sha(Path('scripts/research_h_stage_transport.py').read_bytes()),
        production_changed=False, whole_minimization_executed=False))
    if any(not c.get('exact_original_nm_start') for c in results if c['status'] == 'captured'):
        raise ValueError('original graph H-to-whole transport not exact; do not continue downstream')
    print(json.dumps(dict(captured=sum(c['status'] == 'captured' for c in results), exact_starts=sum(c.get('exact_original_nm_start', False) for c in results)), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--capture-plan', type=Path, required=True)
    parser.add_argument('--capture-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    preflight(args.capture_plan, args.capture_root, args.output)
