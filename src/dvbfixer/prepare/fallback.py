"""Conservative connectivity recovery for existing GLYCAM H before fallback.

The template force field here supplies graph metadata only. It is never used
to parameterize or minimize the forcefield-free hydrogen-placement attempt.
"""
from __future__ import annotations

from typing import Any

from dvbfixer.residue_registry import GLYCAM_PROTEIN_RESNAMES, is_glycam_sugar


class GlycamFallbackTopologyError(RuntimeError):
    """An existing GLYCAM H cannot be bound unambiguously without guessing."""


def restore_glycam_fallback_bonds(topology: Any, positions: Any, *, verbose: bool = False) -> int:
    """Restore missing intra-residue template bonds; leave atoms/positions fixed.

    Only GLYCAM residues containing an unbonded existing H are selected. Build
    an isolated metadata view through the canonical GLYCAM wrapper, validate
    every affected residue, then commit additions atomically. Existing links
    (especially glycosylation/disulfide/peptide links) are never removed or
    inferred anew. Unknown ligand residues do not enter this route.
    """
    from openmm.app import Topology

    from dvbfixer.ffutils import build_glycam_system

    atoms = list(topology.atoms())
    if len(atoms) != len(positions):
        raise GlycamFallbackTopologyError('Fallback atom/position counts differ')
    adjacent: dict[int, list[Any]] = {atom.index: [] for atom in atoms}
    existing = set()
    for a, b in topology.bonds():
        adjacent[a.index].append(b)
        adjacent[b.index].append(a)
        existing.add(frozenset((a.index, b.index)))
    affected = [res for res in topology.residues()
        if (is_glycam_sugar(res.name) or res.name in GLYCAM_PROTEIN_RESNAMES)
        and any(atom.element is not None and atom.element.symbol == 'H' and not adjacent[atom.index]
                for atom in res.atoms())]
    if not affected:
        return 0

    def refuse(res: Any, detail: str) -> None:
        raise GlycamFallbackTopologyError(
            f'{res.chain.id}/{res.name}{res.id}{res.insertionCode.strip()}: {detail}')

    # One private chain per selected residue prevents fabricated peptide links
    # between non-neighboring selected glycoprotein residues. Never infer
    # distance-based external links in this metadata view.
    view = Topology()
    view_to_original = {}
    atom_maps = {}
    for res in affected:
        atom_map = {}
        copied = view.addResidue(res.name, view.addChain(res.chain.id), res.id, res.insertionCode)
        for atom in res.atoms():
            if atom.name in atom_map:
                refuse(res, f'duplicate atom name {atom.name}')
            atom_map[atom.name] = atom
            clone = view.addAtom(atom.name, atom.element, copied, atom.id, atom.formalCharge)
            view_to_original[clone.index] = atom
        atom_maps[res] = atom_map
    # GLYCAM's NLN/OLS/OLT templates reference AMBER protein atom types, so
    # load the normal dependency set even though this view supplies only bonds.
    try:
        ff = build_glycam_system(['amber14-all.xml', 'amber14/GLYCAM_06j-1.xml',
                                 'amber14/tip3pfb.xml'], view, verbose=verbose)
    except Exception as error:
        raise GlycamFallbackTopologyError('Cannot load GLYCAM fallback bond templates') from error
    view_links = {frozenset((view_to_original[a.index].index, view_to_original[b.index].index))
                  for a, b in view.bonds() if a.residue is b.residue}
    additions = []
    for res in affected:
        atom_map = atom_maps[res]
        template = ff._templates.get(res.name)
        if template is None:
            refuse(res, 'no exact GLYCAM template')
        template_atoms = {atom.name: atom for atom in template.atoms}
        if len(template_atoms) != len(template.atoms):
            refuse(res, 'ambiguous template atom names')
        for name, atom in atom_map.items():
            expected = template_atoms.get(name)
            if expected is None or atom.element != expected.element:
                refuse(res, f'unknown or incompatible template atom {name}')
        parents: dict[str, list[str]] = {}
        for i, j in template.bonds:
            a, b = template.atoms[i], template.atoms[j]
            if a.element is not None and a.element.symbol == 'H':
                a, b = b, a
            if b.element is not None and b.element.symbol == 'H':
                parents.setdefault(b.name, []).append(a.name)
            if a.name in atom_map and b.name in atom_map:
                original_a, original_b = atom_map[a.name], atom_map[b.name]
                pair = frozenset((original_a.index, original_b.index))
                if pair not in view_links:
                    refuse(res, 'GLYCAM wrapper did not establish required template bonds')
                if pair not in existing:
                    additions.append((original_a, original_b))
        for atom in res.atoms():
            if atom.element is None or atom.element.symbol != 'H':
                continue
            candidates = parents.get(atom.name, [])
            if len(candidates) != 1 or candidates[0] not in atom_map:
                refuse(res, f'missing or ambiguous parent for {atom.name}')
            parent = atom_map[candidates[0]]
            if parent.element is None or parent.element.symbol == 'H':
                refuse(res, f'invalid heavy parent for {atom.name}')
            if adjacent[atom.index] and adjacent[atom.index] != [parent]:
                refuse(res, f'incompatible existing bonds for {atom.name}')
    for a, b in additions:
        topology.addBond(a, b)
    if additions:
        print(f'  [prep] restored {len(additions)} GLYCAM intra-residue template bonds '
              'before forcefield-free H fallback (coordinates unchanged)')
    return len(additions)
