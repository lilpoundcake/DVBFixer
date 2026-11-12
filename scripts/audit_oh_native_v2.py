"""Offline-only ordered chain-run transport amendment for the frozen OH auditor."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from unittest.mock import patch

from scripts import audit_oh_native as frozen
from scripts.research_dunbrack_prior import sha, write

ROOT = Path(__file__).resolve().parents[1]
AMENDMENT = ROOT/'docs/plans/oh-native-audit-amendment.md'


def restore_ordered(saved):
    from openmm import Vec3, unit
    from openmm.app import Topology, element
    from scripts.research_legacy_start_controls import snapshot_atoms

    topology = Topology()
    atoms, seen_residues = [], set()
    previous_chain, previous_residue, chain, residue = None, None, None, None
    for row in saved['atoms']:
        chain_id, resid, icode, resname, name = row['key']
        key = tuple(row['key'][:4])
        if chain is None or chain_id != previous_chain:
            chain = topology.addChain(chain_id)
            previous_chain = chain_id
        if key != previous_residue:
            if key in seen_residues:
                raise ValueError('ambiguous noncontiguous repeated residue identity')
            seen_residues.add(key)
            residue = topology.addResidue(resname, chain, resid, icode)
            previous_residue = key
        atoms.append(topology.addAtom(name, element.Element.getBySymbol(row['element']), residue))
    for a, b in saved['bonds']:
        topology.addBond(atoms[a], atoms[b])
    if len(atoms) != len(saved['positions_nm']):
        raise ValueError('saved atom/position count differs')
    positions = [Vec3(*p) for p in saved['positions_nm']]*unit.nanometer
    if [r['key'] for r in snapshot_atoms(topology, positions)] != [r['key'] for r in saved['atoms']]:
        raise ValueError('ordered identity transport differs')
    if [[a.index, b.index] for a, b in topology.bonds()] != saved['bonds']:
        raise ValueError('ordered explicit bond transport differs')
    return topology, positions


def audit(plan, run, raw_output, output):
    from scripts.research_ha_minimizer_path import verify

    record = json.loads(plan.read_bytes())
    verify(record)  # Still validates the unmodified original auditor and observer.
    transport = []
    for case in sorted(run.glob('case-*')):
        for path in sorted(case.glob('repair-*-before.json.gz')):
            saved = frozen.read(path)
            top, _ = restore_ordered(saved)
            transport.append(dict(case=case.name, artifact=path.name,
                snapshot_sha256=sha(path.read_bytes()), atoms=top.getNumAtoms(),
                chain_id_runs=[chain.id for chain in top.chains()], residues=top.getNumResidues()))
    with patch.object(frozen, 'restore', restore_ordered):
        frozen.audit(plan, run, raw_output)
    result = json.loads(raw_output.read_bytes())
    result['schema'] = 'dvbfixer.oh-native-audit.v2'
    result['audit_transport'] = dict(adapter_sha256=sha(Path(__file__).read_bytes()),
        amendment_sha256=sha(AMENDMENT.read_bytes()), frozen_scalar_audit_sha256=sha(raw_output.read_bytes()),
        native_rerun=False, snapshot_coordinates_changed=False, inventory=transport,
        limitation='Original chain-occurrence indices were not serialized; ordered ID runs preserve unique contiguous residue ownership and explicit bonds, not all native chain-object boundaries.')
    write(output, result)
    verify(record)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('plan', 'run', 'raw-output', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan.resolve(), args.run.resolve(), args.raw_output.resolve(), args.output.resolve())
