"""Private canonical CB feasibility control and PDBFixer sidechain cohort.

The CB construction uses fixed N/CA/C, bond length/angles and an explicit positive
signed-volume branch. It is not a full residue builder or a production feature.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from research_pdbfixer import heavy_case


def construct_l_cb(n, ca, c, length=0.153, angle_n=110.5, angle_c=110.1):
    n, ca, c = (np.asarray(p, dtype=float) for p in (n, ca, c))
    if not all(np.isfinite(p).all() and p.shape == (3,) for p in (n, ca, c)):
        raise ValueError('finite 3D anchors required')
    u, v = n-ca, c-ca
    if min(np.linalg.norm(u), np.linalg.norm(v)) < 1e-10:
        raise ValueError('coincident anchors')
    u, v = u/np.linalg.norm(u), v/np.linalg.norm(v)
    normal = np.cross(u, v)
    if np.linalg.norm(normal) < 1e-6:
        raise ValueError('degenerate anchors')
    normal /= np.linalg.norm(normal)
    cos_uv = np.dot(u, v)
    coefficients = np.linalg.solve([[1, cos_uv], [cos_uv, 1]],
                                   length*np.cos(np.radians([angle_n, angle_c])))
    parallel = coefficients[0]*u+coefficients[1]*v
    height2 = length**2-np.dot(parallel, parallel)
    if height2 <= 1e-12:
        raise ValueError('incompatible fixed-anchor angles: no nonplanar solution')
    return ca+parallel+math.sqrt(height2)*normal


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[1]
    original = (root/'tests/fixtures/protein_context/1UBQ.pdb').read_text()
    subset = [l for l in original.splitlines(True) if l.startswith('ATOM') and l[21] == 'A' and int(l[22:26]) <= 72]
    results = {'cohort': {}, 'cb_controls': {}}
    # Declared before execution: complete sidechain deletion at six sites;
    # original held-out CB goes to scoring only, never reconstruction.
    for resid in (7, 13, 20, 27, 28, 46):
        text = ''.join(l for l in subset if int(l[22:26]) != resid or l[12:16].strip() in ('N', 'CA', 'C', 'O'))+'TER\nEND\n'
        for proper, guard in ((False, False), (True, False), (True, True)):
            name = f'res{resid}-proper{int(proper)}-guard{int(guard)}'
            results['cohort'][name] = heavy_case(text, args.output/name, proper, guard)
        name = f'res{resid}-backbone-proper-guard'
        results['cohort'][name] = heavy_case(text, args.output/name, True, True, True)
        atoms = {l[12:16].strip(): np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])/10
                 for l in subset if int(l[22:26]) == resid}
        cb = construct_l_cb(atoms['N'], atoms['CA'], atoms['C'])
        results['cb_controls'][str(resid)] = {
            'cb_nm': cb.tolist(),
            'volume_nm3': float(np.dot(np.cross(atoms['N']-atoms['CA'],atoms['C']-atoms['CA']),cb-atoms['CA'])),
            'bond_nm': float(np.linalg.norm(cb-atoms['CA'])),
            'held_out_cb_error_angstrom': float(np.linalg.norm(cb-atoms['CB'])*10)}
    (args.output/'evidence.json').write_text(json.dumps(results, indent=2, allow_nan=False))
    print(json.dumps({k: {'outcome': v['outcome'], 'final_d': v['final_d']} for k,v in results['cohort'].items()}, indent=2))
    print(json.dumps(results['cb_controls'], indent=2))


if __name__ == '__main__':
    main()
