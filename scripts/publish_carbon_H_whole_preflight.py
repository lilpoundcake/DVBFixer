"""Publish metadata-only whole-plan transport admission; no coordinates/Context."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify


def publish(plan_path, output):
    raw = plan_path.read_bytes()
    plan = json.loads(raw)
    verify(plan)
    cases = []
    for c in plan['cases']:
        row = {k:c[k] for k in ('id', 'seed', 'conditioning', 'context_kind', 'status')}
        if c['status'] == 'admitted':
            row['whole_preflight'] = c['whole_preflight']
        cases.append(row)
    if (len(cases) != 40 or sum(c['status'] == 'refused_source' for c in cases) != 6
            or sum(c.get('whole_preflight', {}).get('exact_original_nm_start', False) for c in cases) != 34):
        raise ValueError('complete40/exact34 whole transport required')
    write(output, dict(schema='dvbfixer.carbon-H-whole-preflight.v1', cases=cases,
        policy=plan['policy'], plan_sha256=sha(raw), source_hashes=plan['source_hashes'],
        publication_script_sha256=sha(Path(__file__).read_bytes()),
        whole_minimization_executed=False, production_changed=False, scientific_qualification=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    publish(args.plan, args.output)
