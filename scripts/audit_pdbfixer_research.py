"""Independent readback of research PDBs, observations, chirality and geometry."""
from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np
from openmm import app, unit

from dvbfixer.diagnose.chemistry import (
    check_backbone_bond_angles,
    check_bond_lengths,
    check_peptide_omegas,
)


def identity(atom):
    r = atom.residue
    return r.chain.id, r.id, r.insertionCode, atom.name


def audit(directory):
    pdb = app.PDBFile(str(directory/'returned.pdb'))
    source = app.PDBFile(str(directory/'input.pdb'))
    xyz = np.asarray(pdb.positions.value_in_unit(unit.nanometer))
    original = {identity(a): np.asarray(source.positions[a.index].value_in_unit(unit.nanometer)) for a in source.topology.atoms()}
    output = {identity(a): xyz[a.index] for a in pdb.topology.atoms()}
    missing = sorted(set(original)-set(output))
    if missing:
        raise ValueError(f'observed atoms lost: {missing}')
    movement = max(float(np.linalg.norm(output[k]-p)) for k,p in original.items())
    centres = []
    for residue in pdb.topology.residues():
        atoms = {a.name: xyz[a.index] for a in residue.atoms()}
        if residue.name == 'GLY' or not {'N', 'CA', 'C', 'CB'} <= atoms.keys():
            continue
        n, ca, c, cb = (atoms[k] for k in ('N', 'CA', 'C', 'CB'))
        volume = float(np.linalg.det(np.stack([n-ca, c-ca, cb-ca])))
        if volume <= 1e-6:
            centres.append([residue.chain.id, residue.id, residue.name, volume])
    findings = []
    for check in (check_bond_lengths, check_backbone_bond_angles, check_peptide_omegas):
        findings.extend(asdict(f) for f in check(pdb.topology, pdb.positions))
    return {'max_observed_movement_nm': movement,
            'nonpositive_or_nearplanar_centres': centres, 'geometry_findings': findings,
            'returned_pdb_sha256': hashlib.sha256((directory/'returned.pdb').read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--sidechains', type=Path, required=True)
    parser.add_argument('--pair', type=Path, required=True)
    parser.add_argument('--pair-source', type=Path)
    parser.add_argument('--provenance', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    evidence = json.loads((args.run/'evidence.json').read_text())
    evidence['sidechains'] = json.loads((args.sidechains/'evidence.json').read_text())
    evidence['hydrogen_pair_ablation'] = json.loads((args.pair/'evidence.json').read_text())
    pair_source = args.pair_source or args.run
    evidence['hydrogen_pair_source'] = {'directory': str(pair_source), 'controls': {}}
    for seed in (1, 7, 19, 42, 101):
        name = f'h-seed{seed}-iter50-strip0-ff1'
        hashes = {}
        for artifact in ('stage-0.npy', 'system.xml', 'atoms.json'):
            data = (pair_source/name/artifact).read_bytes()
            if data != (args.run/name/artifact).read_bytes():
                raise ValueError(f'pair proposal differs from recorded placement: {name}/{artifact}')
            hashes[artifact] = hashlib.sha256(data).hexdigest()
        evidence['hydrogen_pair_source']['controls'][name] = hashes
    if args.provenance:
        evidence['provenance'] = json.loads(args.provenance.read_text())
    evidence['audit'] = {}
    for root, group in ((args.run, evidence['heavy']), (args.sidechains, evidence['sidechains']['cohort'])):
        for name, result in group.items():
            if result['outcome'] != 'returned':
                raise ValueError(f'unexecuted scientific outcome: {name}: {result["outcome"]}')
            evidence['audit'][name] = audit(root/name)
            result['raw_result_sha256'] = hashlib.sha256((root/name/'result.json').read_bytes()).hexdigest()
    # Full callbacks stay local; portable record retains the first negative one
    # plus phase boundaries, with exact raw result digests for later comparison.
    for group in (evidence['heavy'], evidence['sidechains']['cohort']):
        for name, result in group.items():
            callbacks = result.pop('minimizer_iterations')
            result['callback_count'] = len(callbacks)
            result['first_negative_callback'] = next((r for r in callbacks if any(v[4] < -1e-6 for v in r['volumes_nm3'])), None)
            for event in result['events']:
                event['volumes_nm3'] = [v for v in event['volumes_nm3'] if v[1] in ('7','10','11','12','13','20','27','28','46')]
    evidence['scope'] = 'Private causal controls; no backend enablement or broad qualification'
    evidence['raw_directories'] = [str(args.run), str(args.sidechains), str(args.pair)]
    evidence['audit_script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.write_text(json.dumps(evidence, indent=2, ensure_ascii=False, allow_nan=False))
    print(json.dumps({k: {'movement_nm': v['max_observed_movement_nm'],
                           'non_L': v['nonpositive_or_nearplanar_centres'],
                           'geometry_findings': len(v['geometry_findings'])} for k,v in evidence['audit'].items()}, indent=2))


if __name__ == '__main__':
    main()
