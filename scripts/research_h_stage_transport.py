"""Research-only faithful post-H repair bridge, no production orchestration.

Reconstruct the bound atom/bond topology without renaming identities, apply the
existing production repair once, require exact original whole-system start.
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

from scripts.audit_legacy_start_controls import coordinate_delta, load_atoms
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify


def topology_from_rows(rows, bonds):
    from openmm.app import Topology, element

    topology = Topology()
    chains, residues, atoms = {}, {}, []
    seen = set()
    for row in rows:
        key = tuple(row['key'])
        if key in seen or len(key) != 5:
            raise ValueError('ambiguous transport atom identity')
        seen.add(key)
        chain_id, resid, icode, resname, name = key
        chain = chains.get(chain_id)
        if chain is None:
            chain = topology.addChain(chain_id)
            chains[chain_id] = chain
        residue_key = key[:4]
        residue = residues.get(residue_key)
        if residue is None:
            residue = topology.addResidue(resname, chain, resid, insertionCode=icode)
            residues[residue_key] = residue
        atom = topology.addAtom(name, element.Element.getBySymbol(row['element']), residue)
        if atom.index != len(atoms):
            raise ValueError('transport atom indices changed during topology reconstruction')
        atoms.append(atom)
    seen_bonds = set()
    for a, b in bonds:
        if min(a, b) < 0 or max(a, b) >= len(atoms) or a == b or tuple(sorted((a, b))) in seen_bonds:
            raise ValueError('invalid or duplicate transport bond')
        seen_bonds.add(tuple(sorted((a, b))))
        topology.addBond(atoms[a], atoms[b])
    if [a.index for a in topology.atoms()] != list(range(len(rows))):
        raise ValueError('transport topology iteration would reorder captured atoms')
    return topology


def post_h_repair(rows, bonds, points_nm):
    from openmm import Vec3
    from openmm.unit import nanometer

    from dvbfixer.ffutils.geometry import repair_misplaced_hydrogens

    if len(rows) != len(points_nm):
        raise ValueError('transport coordinate inventory differs')
    topology = topology_from_rows(rows, bonds)
    positions = [Vec3(*p) for p in points_nm]*nanometer
    repaired = repair_misplaced_hydrogens(topology, positions)
    result = [list(p) for p in positions.value_in_unit(nanometer)]
    if any(result[i] != list(p) for i, p in enumerate(points_nm) if rows[i]['element'] not in {'H', 'D'}):
        raise ValueError('existing post-H repair moved a heavy atom')
    return result, repaired


def preflight(plan_path, output):
    import openmm
    from openmm.unit import nanometer

    raw = plan_path.read_bytes()
    plan = json.loads(raw)
    verify(plan)
    results = []
    for case in plan['cases']:
        source = Path(case['source'])
        rows = load_atoms(source/'h-internal-before.json.gz')
        path = json.loads(gzip.decompress(Path(case['path']).read_bytes()))
        bonds = json.loads((source/'topology.json').read_bytes())['bonds']
        points, repaired = post_h_repair(rows, bonds, path['final_nm'])
        after = [dict(r, xyz=[10*x for x in p]) for r, p in zip(rows, points, strict=True)]
        delta = coordinate_delta(load_atoms(source/'phase-1-before.json.gz'), after)
        original_state = openmm.XmlSerializer.deserialize((source/'initial-state.xml').read_text())
        original_nm = [list(p) for p in original_state.getPositions().value_in_unit(nanometer)]
        nm_exact = points == original_nm
        result = dict(id=case['id'], seed=case['seed'], source_arm=case['source_arm'],
            original_H_repairs=repaired, whole_system_start_difference=delta,
            nm_start_exact=nm_exact, exact=nm_exact and all(delta[k]['exact'] for k in ('heavy', 'hydrogen')))
        results.append(result)
    report = dict(schema='dvbfixer.h-stage-transport-preflight.v1', cases=results,
        exact_original_starts=sum(r['exact'] for r in results), enrolled=len(results),
        plan_sha256=sha(raw), script_sha256=sha(Path(__file__).read_bytes()),
        production_geometry_sha256=sha(Path('src/dvbfixer/ffutils/geometry.py').read_bytes()),
        production_changed=False, whole_system_minimization=False)
    write(output, report)
    print(json.dumps(dict(enrolled=len(results), exact_original_starts=report['exact_original_starts']), indent=2))
    if not all(r['exact'] for r in results):
        raise ValueError('post-H transport failed exact original-start preflight; do not proceed to downstream sampling')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    preflight(args.plan, args.output)
