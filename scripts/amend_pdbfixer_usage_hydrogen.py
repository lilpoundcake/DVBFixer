"""Current AST source amendment; never overwrite historical scientific inventories."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from scripts.inventory_pdbfixer_usage import inventory

ROOT = Path(__file__).resolve().parents[1]
OLD = Path('docs/research/pdbfixer-usage-inventory-stereo-diagnostics.json')


def semantic(records):
    return [{k:v for k, v in row.items() if k != 'line'} for row in records]


def generate(output):
    raw = (ROOT/OLD).read_bytes()
    historical = json.loads(raw)
    current = inventory(ROOT)
    if current['totals'] != historical['totals'] or len(current['files']) != len(historical['files']):
        raise ValueError('unexpected PDBFixer use inventory change, not a source-only amendment')
    allowed = {'src/dvbfixer/ffutils/geometry.py', 'src/dvbfixer/prepare/pipeline.py', 'src/dvbfixer/diagnose/structural.py'}
    changes = []
    for before, after in zip(historical['files'], current['files'], strict=True):
        if before['path'] != after['path'] or any(semantic(before[key]) != semantic(after[key]) for key in ('imports', 'uses')):
            raise ValueError('PDBFixer import/use semantics changed; explicit audit required')
        if before != after:
            if after['path'] not in allowed:
                raise ValueError('unexpected changed production source binding')
            changes.append(dict(path=after['path'], before_sha256=before['source_sha256'], after_sha256=after['source_sha256']))
    result = dict(current, historical_inventory=dict(path=str(OLD), sha256=hashlib.sha256(raw).hexdigest()),
        changed_source_bindings=changes, generator_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        historical_template_probe_rerun=False, historical_scientific_results_changed=False,
        production_hydrogen_postflight_changed=True,
        limits='AST inventory only; old research inventories and template evidence remain immutable')
    with output.open('x') as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write('\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    generate(parser.parse_args().output)
