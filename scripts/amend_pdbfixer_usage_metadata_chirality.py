"""Register pure metadata-L factoring; keep historical PDBFixer inventories intact."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from scripts.inventory_pdbfixer_usage import inventory

ROOT = Path(__file__).resolve().parents[1]
PARENT = Path('docs/research/pdbfixer-usage-inventory-comparison.json')
CHANGED = 'src/dvbfixer/chemical_reconstruction/template_comparison.py'


def semantics(record):
    return {key: [{k: v for k, v in row.items() if k != 'line'} for row in record[key]] for key in ('imports', 'uses')}


def generate(output: Path):
    previous = json.loads((ROOT/PARENT).read_bytes())
    current = inventory(ROOT)
    if current['totals'] != previous['totals'] or len(current['files']) != len(previous['files']):
        raise ValueError('unexpected PDBFixer use delta')
    changes = []
    for before, after in zip(previous['files'], current['files'], strict=True):
        if before['path'] != after['path'] or semantics(before) != semantics(after):
            raise ValueError('unexpected import/call semantics change')
        if before != after:
            changes.append({'path': after['path'], 'old_source_sha256': before['source_sha256'],
                            'new_source_sha256': after['source_sha256']})
    if [c['path'] for c in changes] != [CHANGED]:
        raise ValueError('expected only private comparator metadata-L factoring')
    result = dict(current, historical_inventory={'path': str(PARENT), 'sha256': hashlib.sha256((ROOT/PARENT).read_bytes()).hexdigest()},
        changed_source_bindings=changes, generator_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        reason='shared metadata-only assert_all_l adapter; no PDBFixer constructor/rebuild/method semantics changed',
        historical_template_probe_rerun=False, historical_scientific_results_changed=False,
        public_backend_or_dependency_replacement=False)
    with output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    generate(parser.parse_args().output)
