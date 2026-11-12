"""Append current source bindings without rewriting historical OH inventories."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from scripts.inventory_pdbfixer_usage import inventory
from scripts.research_dunbrack_prior import sha, write

ROOT = Path(__file__).resolve().parents[1]
PARENT = Path('docs/research/pdbfixer-usage-inventory-hydrogen-postflight.json')


def generate(output):
    previous = json.loads((ROOT/PARENT).read_bytes())
    current = inventory(ROOT)
    if current['totals'] != previous['totals'] or len(current['files']) != len(previous['files']):
        raise ValueError('PDBFixer use inventory changed')
    changes = []
    for old, new in zip(previous['files'], current['files'], strict=True):
        if old['path'] != new['path']:
            raise ValueError('PDBFixer source paths changed')
        for key in ('imports', 'uses'):
            a = [{k:v for k, v in row.items() if k != 'line'} for row in old[key]]
            b = [{k:v for k, v in row.items() if k != 'line'} for row in new[key]]
            if a != b:
                raise ValueError('PDBFixer import/use semantics changed')
        if old != new:
            if old['path'] != 'src/dvbfixer/prepare/pipeline.py':
                raise ValueError('Unexpected production source binding changed')
            changes.append(dict(path=new['path'], before_sha256=old['source_sha256'],
                                after_sha256=new['source_sha256']))
    write(output, dict(current, historical_inventory=dict(path=str(PARENT), sha256=sha((ROOT/PARENT).read_bytes())),
        changed_source_bindings=changes, generator_source_sha256=sha(Path(__file__).read_bytes()),
        historical_template_probe_rerun=False, historical_scientific_results_changed=False,
        production_glycam_fallback_changed=True))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    generate(parser.parse_args().output)
