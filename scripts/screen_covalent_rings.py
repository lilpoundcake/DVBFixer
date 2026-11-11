#!/usr/bin/env python3
"""Experimental canonical-residue small-ring diagnostics from explicit topology.

RDKit supplies a graph cycle basis only; no proximity bond perception, aromatic
assignment or force-field chemistry is reimplemented. This is not exhaustive
macrocycle/disulfide topology, a trajectory proof or a public acceptance gate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict
from pathlib import Path

from dvbfixer.model.diffusion.anchor_acceptance import _structure
from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest
from dvbfixer.model.diffusion.covalent_ring_probe import probe_covalent_ring
from dvbfixer.model.diffusion.gap_portfolio import _coordinates
from dvbfixer.model.diffusion.protein_atoms import CANONICAL_HEAVY_ATOMS


def canonical_cycle(cycle):
    values = tuple(cycle)
    if len(set(values)) != len(values):
        raise ValueError('graph cycle repeats an atom')
    return min(direction[i:]+direction[:i] for direction in (values, values[::-1]) for i in range(len(values)))


def screen(text, *, generated=frozenset(), max_rings=512):
    from rdkit import Chem, rdBase

    structure = _structure(text, 'small-ring-screen.pdb')
    xyz, _ = _coordinates(text)
    def identity(atom):
        return AtomIdentity(atom.residue.chain.id, atom.residue.id, atom.residue.insertionCode.strip(), atom.name)
    canonical = {identity(atom): atom for atom in structure.pdb.topology.atoms()
                 if atom.residue.name in CANONICAL_HEAVY_ATOMS and
                 (atom.name in CANONICAL_HEAVY_ATOMS[atom.residue.name] or atom.name == 'OXT')}
    atoms = tuple(sorted(canonical))
    if len(atoms) > 50_000 or not set(atoms) <= xyz.keys():
        raise ValueError('canonical topology/coordinates differ or exceed atom budget')
    bonds_set = set()
    excluded_crosslinks = 0
    for a,b in structure.pdb.topology.bonds():
        if identity(a) not in canonical or identity(b) not in canonical:
            continue
        if a.residue != b.residue:
            peptide = ({a.name,b.name} == {'C','N'} and a.residue.chain == b.residue.chain
                       and abs(a.residue.index-b.residue.index) == 1)
            if not peptide:
                excluded_crosslinks += 1
                continue
        bonds_set.add(frozenset((identity(a),identity(b))))
    bonds = frozenset(bonds_set)
    if len(bonds) > 50_000:
        raise ValueError('small-ring topology edge budget exceeded')
    graph = Chem.RWMol()
    index = {atom: i for i, atom in enumerate(atoms)}
    for atom in atoms:
        element = canonical[atom].element
        if element is None:
            raise ValueError('canonical atom has no element')
        graph.AddAtom(Chem.Atom(element.atomic_number))
    for edge in sorted(bonds, key=lambda pair: tuple(sorted(pair))):
        a,b = sorted(edge)
        graph.AddBond(index[a], index[b], Chem.BondType.SINGLE)
    # Bond orders are intentionally irrelevant here. Do not sanitize/reperceive
    # chemistry: this is cycle discovery on the already supplied covalent graph.
    cycles = set()
    excluded_cycles = 0
    for values in Chem.GetSymmSSSR(graph.GetMol()):
        cycle = tuple(atoms[int(i)] for i in values)
        residues = {(a.chain, a.residue_number, a.insertion_code) for a in cycle}
        if not 3 <= len(cycle) <= 8 or len(residues) != 1:
            excluded_cycles += 1
            continue
        cycles.add(canonical_cycle(cycle))
    if len(cycles) > max_rings:
        raise ValueError('small-ring cycle budget exceeded')
    findings = [asdict(probe_covalent_ring(cycle, {a:tuple(v) for a,v in xyz.items()}, bonds, generated))
                for cycle in sorted(cycles)]
    return {'method': 'experimental-canonical-topology-cycle-basis-ring-screen-v1',
        'rdkit_version': rdBase.rdkitVersion, 'canonical_atoms': len(atoms), 'topology_bonds': len(bonds),
        'excluded_topology_atoms': structure.pdb.topology.getNumAtoms()-len(atoms),
        'excluded_nonpeptide_interresidue_links': excluded_crosslinks,
        'topology_provenance': 'OpenMM PDBFile canonical-residue topology; intramolecular cycles and same-chain file-neighbor peptides only; not independent chemical validation of inherited topology',
        'excluded_residue_names': sorted({r.name for r in structure.pdb.topology.residues()
                                          if r.name not in CANONICAL_HEAVY_ATOMS}),
        'cycles_considered': len(cycles), 'cycles_outside_scope': excluded_cycles,
        'undefined_surfaces': sum(r['undefined_reason'] is not None for r in findings),
        'piercings': sum(len(r['piercings']) for r in findings), 'rings': findings,
        'limitations': 'cycle basis only, canonical intra-residue rings3-8 only; noncanonical rings, macrocycles, disulfide topology and unsaved trajectories unqualified; no acceptance gate',
        'qualified': False, 'publication_gate': False}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('candidate', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--request', type=Path)
    args = parser.parse_args()
    data = args.candidate.read_bytes()
    request = DiffusionRequest.from_json(args.request.read_text()) if args.request else None
    result = screen(data.decode(), generated=frozenset(request.generated_atoms) if request else frozenset())
    result['candidate_sha256'] = hashlib.sha256(data).hexdigest()
    result['ownership_mask_request_sha256'] = hashlib.sha256(args.request.read_bytes()).hexdigest() if args.request else None
    with args.output.open('x') as stream:
        stream.write(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('cycles_considered','undefined_surfaces','piercings')}))
