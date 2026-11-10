"""Read-only H/whole boundary lineage audit; no graph-arm sampling."""
from __future__ import annotations

import argparse
import gzip
import json
import math
from collections import Counter
from pathlib import Path

from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify


def bijection(mapping, size):
    if len(mapping) != size or any(type(i) is not int for i in mapping) or sorted(mapping) != list(range(size)):
        raise ValueError('captured lineage is not a complete index bijection')


def validate_graph(rows, bonds, ownership):
    adjacent = [set() for _ in rows]
    for a, b in bonds:
        if a == b or not 0 <= a < len(rows) or not 0 <= b < len(rows) or b in adjacent[a]:
            raise ValueError('invalid capture graph')
        adjacent[a].add(b)
        adjacent[b].add(a)
    owned = set(ownership['new_h']) & set(ownership['mobile'])
    expected = []
    for i, row in enumerate(rows):
        hs = [j for j in adjacent[i] if rows[j]['element'] in {'H', 'D'}]
        others = adjacent[i]-set(hs)
        if row['element'] == 'C' and len(hs) == 3 and len(others) == 1 and set(hs) <= owned:
            expected.append(i)
    groups = ownership['graph']['groups']
    if [g['parent'] for g in groups] != sorted(expected, key=lambda i: tuple(rows[i]['key'])):
        raise ValueError('graph methyl enumeration/ordering differs')
    for group in groups:
        p, axis, frame = group['parent'], group['axis'], group['frame']
        if set(group['hydrogens']) != adjacent[p]-{axis} or rows[axis]['element'] in {'H', 'D'}:
            raise ValueError('graph methyl H/axis inventory differs')
        candidates = [i for i in adjacent[axis]-{p} if rows[i]['element'] not in {'H', 'D'}]
        candidates.sort(key=lambda i: (rows[i]['key'][:4] != rows[p]['key'][:4], tuple(rows[i]['key'])))
        if not candidates or frame != candidates[0] or group['parent_key'] != rows[p]['key'] or group['frame_key'] != rows[frame]['key']:
            raise ValueError('graph methyl frame/identity differs')
    return len(groups)


def audit(plan_path, root, output):
    import openmm
    from openmm.unit import dalton, nanometer

    raw = plan_path.read_bytes()
    plan = json.loads(raw)
    verify(plan)
    run_raw = (root/'evidence.json').read_bytes()
    run = json.loads(run_raw)
    expected = [(c['id'], s) for c in plan['cases'] for s in plan['policy']['seeds']]
    if (root/'executed-plan.json').read_bytes() != raw or run['plan_sha256'] != sha(raw) or [(c['id'], c['seed']) for c in run['cases']] != expected:
        raise ValueError('capture full enrollment/binding differs')
    results = []
    for cell in run['cases']:
        folder = root/cell['id']/f"seed-{cell['seed']}"
        for name, digest in cell['artifacts'].items():
            if sha((folder/name).read_bytes()) != digest:
                raise ValueError('capture artifact changed')
        row = {k: cell[k] for k in ('id', 'seed', 'status', 'returncode')}
        if cell['status'] != 'captured':
            results.append(row)
            continue
        event = json.loads((folder/'capture.json').read_bytes())
        if event['plan_sha256'] != sha(raw) or event['whole_minimization_executed'] or not event['hooks_restored'] or event['H_calls'] != 1:
            raise ValueError('capture observer/boundary contract differs')
        topology = json.loads((folder/'H-topology.json').read_bytes())
        ownership = json.loads((folder/'H-ownership.json').read_bytes())
        path = json.loads(gzip.decompress((folder/'H-path.json.gz').read_bytes()))
        rows, bonds = topology['atoms'], topology['bonds']
        if rows != ownership['rows'] or bonds != ownership['bonds']:
            raise ValueError('H topology/ownership inventories differ')
        system = openmm.XmlSerializer.deserialize((folder/'H-system.xml').read_text())
        if system.getNumParticles() != len(rows) or system.getNumConstraints() or system.usesPeriodicBoundaryConditions() or any(system.isVirtualSite(i) for i in range(len(rows))):
            raise ValueError('unsupported captured H System')
        mobile = [i for i in range(len(rows)) if system.getParticleMass(i).value_in_unit(dalton) > 0]
        if mobile != ownership['mobile'] or any(rows[i]['element'] not in {'H', 'D'} for i in mobile):
            raise ValueError('H System mobile/ownership differs')
        before, after = path['initial_nm'], path['final_nm']
        if len(before) != len(rows) or len(after) != len(rows) or any(len(p) != 3 or not all(math.isfinite(x) for x in p) for p in before+after):
            raise ValueError('captured H coordinate inventory differs')
        heavy = [i for i, r in enumerate(rows) if r['element'] not in {'H', 'D'}]
        if any(before[i] != after[i] for i in heavy):
            raise ValueError('native H minimizer moved heavy atoms')
        input_rows = json.loads(gzip.decompress((folder/'H-input.json.gz').read_bytes()))['atoms']
        keys = {tuple(r['key']) for r in input_rows}
        new_h = [i for i, r in enumerate(rows) if r['element'] in {'H', 'D'} and tuple(r['key']) not in keys]
        if new_h != ownership['new_h']:
            raise ValueError('new-H ownership differs from captured pre-H observations')
        group_count = validate_graph(rows, bonds, ownership)
        if group_count != event['groups']:
            raise ValueError('group summary differs')
        returned = json.loads(gzip.decompress((folder/'H-return.json.gz').read_bytes()))['atoms']
        post_map = ownership['post_H_lineage']
        bijection(post_map, len(rows))
        for i, j in enumerate(post_map):
            if returned[j]['key'] != ownership['post_H_labels'][j]['key'] or returned[j]['element'] != rows[i]['element']:
                raise ValueError('post-H label lineage differs')
            if any(not math.isclose(x*10, y, abs_tol=1e-12, rel_tol=0) for x, y in zip(after[i], returned[j]['xyz'], strict=True)):
                raise ValueError('H endpoint does not match native addHydrogens return')
        whole = json.loads((folder/'whole-topology.json').read_bytes())
        whole_path = json.loads(gzip.decompress((folder/'whole-start.json.gz').read_bytes()))
        mapping = whole['H_to_whole_indices']
        bijection(mapping, len(rows))
        mapped_bonds = {tuple(sorted((mapping[a], mapping[b]))) for a, b in bonds}
        if mapped_bonds != {tuple(sorted(b)) for b in whole['bonds']}:
            raise ValueError('native H-to-whole bond lineage differs')
        state = openmm.XmlSerializer.deserialize((folder/'whole-state.xml').read_text())
        whole_nm = [list(p) for p in state.getPositions().value_in_unit(nanometer)]
        if whole_nm != whole_path['positions_nm']:
            raise ValueError('whole State/nm snapshot differs')
        label_changes = []
        for i, j in enumerate(mapping):
            full = whole['atoms'][j]
            if full['element'] != rows[i]['element'] or full['key'][:3] != rows[i]['key'][:3]:
                raise ValueError('native object lineage changed residue identity/element')
            if full['key'] != rows[i]['key']:
                label_changes.append(dict(index=i, H_key=rows[i]['key'], whole_key=full['key']))
            if i in heavy and before[i] != whole_nm[j]:
                raise ValueError('post-H bridge moved a heavy observation')
        row.update(groups=group_count, excluded_groups=len(ownership['graph']['exclusions']), atoms=len(rows),
            new_H=len(new_h), heavy=len(heavy), exact_heavy_path=True, complete_bond_lineage=True,
            label_changes=label_changes, H_system_sha256=sha((folder/'H-system.xml').read_bytes()))
        results.append(row)
    summary = dict(enrolled=len(results), statuses=dict(Counter(r['status'] for r in results)),
        captured_groups=sum(r.get('groups', 0) for r in results))
    write(output, dict(schema='dvbfixer.methyl-diversity-capture-audit.v1', cases=results, summary=summary,
        plan_sha256=sha(raw), run_sha256=sha(run_raw), auditor_sha256=sha(Path(__file__).read_bytes()),
        amendment=plan.get('amendment'), production_changed=False, graph_H_arms_evaluated=False,
        whole_minimization_executed=False, scientific_qualification=False))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan, args.run, args.output)
