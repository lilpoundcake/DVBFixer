"""Register the private WP2P baseline without rewriting historical use inventories."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from scripts.inventory_pdbfixer_usage import inventory

ROOT = Path(__file__).resolve().parents[1]
PARENT = Path('docs/research/pdbfixer-usage-inventory-glycam-fallback.json')
ADDED = 'src/dvbfixer/chemical_reconstruction/template_comparison.py'


def generate(output: Path):
    previous = json.loads((ROOT/PARENT).read_bytes())
    current = inventory(ROOT)
    old_files = [f for f in current['files'] if f['path'] != ADDED]
    added = [f for f in current['files'] if f['path'] == ADDED]
    if old_files != previous['files'] or len(added) != 1:
        raise ValueError('unexpected historical-source or dependency change')
    delta = {kind: current['totals'][kind]-previous['totals'][kind] for kind in current['totals']}
    if delta != {'constructor': 1, 'seeded_rebuild': 1, 'method': 1}:
        raise ValueError('unexpected private comparator use semantics')
    def sha(data):
        return hashlib.sha256(data).hexdigest()
    result = dict(current, historical_inventory={'path': str(PARENT), 'sha256': sha((ROOT/PARENT).read_bytes())},
        added_private_sources=[ADDED], added_totals=delta, generator_source_sha256=sha(Path(__file__).read_bytes()),
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
