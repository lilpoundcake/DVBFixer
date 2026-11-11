#!/usr/bin/env python3
"""Visualization-only coupled-group ensemble; never an accepted structure bundle.

MODEL indices select an entire joint proposal for gaps 2/3, not independent loop
donors. Noncontiguous source context gets explicit TER boundaries. No generated
coordinates from the other eight gaps enter this cropped visualization.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from dvbfixer.model.diffusion.contract import AtomIdentity

SEEDS = (7, 42, 101, 102, 103, 104)


def atom_records(text):
    result = {}
    for line in text.splitlines():
        if not line.startswith(('ATOM  ', 'HETATM')):
            continue
        atom = AtomIdentity(line[21], line[22:26].strip(), line[26].strip(), line[12:16].strip())
        if atom in result:
            raise ValueError('duplicate atom identity in ensemble input')
        result[atom] = line
    return result


def model_lines(axis, records, *, model_index):
    """Coordinate-preserving PDB export; explicit breaks between source spans."""
    output = [f'MODEL     {model_index:4d}\n']
    previous = None
    for serial, atom in enumerate(axis, 1):
        if atom not in records:
            raise ValueError('ensemble input omits a requested context atom')
        residue = (atom.chain, atom.residue_number, atom.insertion_code)
        if previous is not None and residue != previous:
            if residue[0] != previous[0] or int(residue[1]) != int(previous[1])+1 or residue[2] or previous[2]:
                output.append('TER\n')
        line = records[atom]
        output.append(line[:6]+f'{serial:5d}'+line[11:]+'\n')
        previous = residue
    output.extend(('TER\n', 'ENDMDL\n'))
    return ''.join(output)


def run(root, output):
    panel = root / '1he8-rfd3-context-development-panel-attempt02'
    completion = json.loads((panel / 'completed.json').read_text())
    if completion['executed'] != 5:
        raise ValueError('all newly declared context calls must finish first')
    context = root / '1he8-rfd3-pair-context-input/manifest.json'
    manifest = json.loads(context.read_text())
    source_bytes = (root / '1he8-visual-compatible-selection/source.pdb').read_bytes()
    if hashlib.sha256(source_bytes).hexdigest() != manifest['source_sha256']:
        raise ValueError('ensemble source differs from declared conditioning')
    original = atom_records(source_bytes.decode())
    axis = tuple(AtomIdentity(**a['atom']) for a in manifest['axis'])
    fixed = {AtomIdentity(**a['atom']) for a in manifest['axis'] if not a['generated']}
    if len(set(axis)) != len(axis) or len(fixed) != manifest['observed_atoms'] or len(axis)-len(fixed) != manifest['generated_atoms']:
        raise ValueError('context ownership/atom inventory differs')
    models = []
    custody = []
    for index, seed in enumerate(SEEDS, 1):
        path = (root / '1he8-rfd3-context-assembly-recovered/joint-context/raw.pdb' if seed == 7 else
                panel / f'seed-{seed}/raw-assembly.pdb')
        data = path.read_bytes()
        records = atom_records(data.decode())
        if any(records[a] != original[a] for a in fixed):
            raise ValueError('context observation changed before visualization')
        model = model_lines(axis, records, model_index=index)
        readback = atom_records(model)
        if set(readback) != set(axis) or any(readback[a][30:54] != records[a][30:54] for a in axis):
            raise ValueError('ensemble export changed coordinates/identities')
        models.append(model)
        custody.append({'model': index, 'seed': seed, 'raw_assembly': str(path),
                        'raw_sha256': hashlib.sha256(data).hexdigest(), 'atoms': len(axis)})
    output.mkdir(exist_ok=False)
    pdb = ('REMARK 999 VISUALIZATION ONLY: UNREVIEWED RAW RFD3 COUPLED GAP2/GAP3 CONTROLS\n'
           'REMARK 999 NOT SCIENTIFIC ACCEPTANCE; SAME MODEL MUST OWN BOTH GAPS\n'+''.join(models)+'END\n')
    target = output / 'joint-gap2-gap3_raw_all-seeds.pdb'
    target.write_text(pdb)
    report = {'scope': 'visualization-only raw coupled groups plus fixed observed context; no other eight generated donors',
        'model_to_seed': {str(i): s for i, s in enumerate(SEEDS, 1)}, 'joint_gap_indices': [2,3],
        'source_sha256': manifest['source_sha256'], 'context_manifest_sha256': hashlib.sha256(context.read_bytes()).hexdigest(),
        'ensemble_sha256': hashlib.sha256(target.read_bytes()).hexdigest(), 'models': custody,
        'source_atom_coordinates_unchanged': True, 'model_group_atomic': True,
        'noncontiguous_source_spans_separated_by_TER': True, 'original_native_chain_topology_not_inferred': True,
        'seed7': 'original retained run with explicit hash-bound source-label recovery; other seeds from width-patched clone',
        'user_visual_approved': False, 'geometry_accepted': False, 'qualification': False}
    (output / 'manifest.json').write_text(json.dumps(report, indent=2)+'\n')
    (output / 'README.txt').write_text(
        'Unreviewed raw proposals, NOT accepted structures.\n'
        'MODEL 1/2/3/4/5/6 = seeds 7/42/101/102/103/104.\n'
        'Select both loops A436-456 and A489-497 from the SAME MODEL; they were generated jointly.\n'
        'Fixed source context is included in the original frame; omitted source spans have explicit TER.\n'
        'The other eight generated gaps are excluded. Original cc94 visual approvals do not apply.\n')
    print(target)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    run(args.root.absolute(), args.output.absolute())
