"""Explicit source-only inventory amendment after diagnostic cherry-picks.

Retain the original inventory and offline template probe unchanged. Refuse any
PDBFixer call/import/count change; this amendment only records the new hash of
the read-only diagnose consumer, not new scientific qualification.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

from scripts.inventory_pdbfixer_usage import inventory


def amend(root, output):
    historical_path = root/'docs/research/pdbfixer-usage-inventory.json'
    historical_raw = historical_path.read_bytes()
    historical = json.loads(historical_raw)
    current = inventory(root)
    old = {r['path']: r for r in historical['files']}
    new = {r['path']: r for r in current['files']}
    if old.keys() != new.keys() or historical['totals'] != current['totals']:
        raise ValueError('consumer inventory/counts changed beyond diagnostic-only amendment')
    changed = []
    for name in old:
        if old[name]['imports'] != new[name]['imports'] or old[name]['uses'] != new[name]['uses']:
            raise ValueError('PDBFixer imports/calls changed beyond source-hash amendment')
        if old[name]['source_sha256'] != new[name]['source_sha256']:
            changed.append(dict(path=name, before_source_sha256=old[name]['source_sha256'],
                                after_source_sha256=new[name]['source_sha256']))
    if [r['path'] for r in changed] != ['src/dvbfixer/diagnose/pipeline.py']:
        raise ValueError('unexpected source hash changes for this explicit diagnostic transfer')
    result = dict(schema='dvbfixer.pdbfixer-usage-source-amendment.v1',
        historical_inventory=dict(path='docs/research/pdbfixer-usage-inventory.json',
            sha256=hashlib.sha256(historical_raw).hexdigest()),
        reason='diagnostic-only 039e1ba/9e71372 cherry-picks as 471844a/8eeeada; unchanged PDBFixer consumers',
        source_git_head=subprocess.check_output(['git', '-C', str(root), 'rev-parse', 'HEAD'], text=True).strip(),
        generator_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        files=current['files'], totals=current['totals'], changed_source_bindings=changed,
        historical_template_probe_rerun=False, scientific_results_changed=False,
        network_allowed=False, coordinates_generated=False)
    with output.open('x') as stream:
        json.dump(result, stream, sort_keys=True, indent=2)
        stream.write('\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = amend(Path(__file__).resolve().parents[1], args.output)
    print(json.dumps(dict(totals=result['totals'], changed_source_bindings=result['changed_source_bindings']), indent=2))
