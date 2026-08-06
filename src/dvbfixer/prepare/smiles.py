"""Optional SMILES-guided hydrogen addition for isolated ligand residues."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any


class SmilesPreparationError(ValueError):
    """A supplied ligand SMILES cannot be applied safely to the PDB graph."""


def parse_smiles_mappings(values: Iterable[str]) -> dict[str, str]:
    """Parse repeatable ``RESNAME=SMILES`` CLI values and validate SMILES.

    RDKit remains optional for the normal prepare path.  It is imported only
    when at least one mapping is supplied.
    """
    values = list(values)
    if not values:
        return {}
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise SmilesPreparationError(
            "--smiles requires RDKit; install the rdkit package or omit --smiles"
        ) from exc

    mappings: dict[str, str] = {}
    for raw in values:
        if "=" not in raw:
            raise SmilesPreparationError(
                f"invalid --smiles value {raw!r}; expected RESNAME=SMILES"
            )
        resname, smiles = raw.split("=", 1)
        resname = resname.strip()
        smiles = smiles.strip()
        if not resname or not smiles:
            raise SmilesPreparationError(
                f"invalid --smiles value {raw!r}; residue name and SMILES are required"
            )
        if any(ch.isspace() for ch in resname):
            raise SmilesPreparationError(
                f"invalid residue name {resname!r} in --smiles mapping"
            )
        if resname in mappings:
            raise SmilesPreparationError(
                f"duplicate --smiles mapping for residue name {resname!r}"
            )
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise SmilesPreparationError(
                f"invalid SMILES for residue {resname}: {smiles!r}"
            )
        if len(Chem.GetMolFrags(mol)) != 1:
            raise SmilesPreparationError(
                f"SMILES for residue {resname} must describe one connected molecule"
            )
        mappings[resname] = smiles
    return mappings


def _heavy_molecule(smiles: str) -> Any:
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:  # Already checked by parse_smiles_mappings; defensive API guard.
        raise SmilesPreparationError(f"invalid SMILES: {smiles!r}")
    mol = Chem.RemoveHs(mol)
    if not mol.GetNumAtoms():
        raise SmilesPreparationError("SMILES must contain at least one heavy atom")
    return mol


def _graph_mappings(pdb_atoms: list[Any], pdb_edges: set[tuple[int, int]], mol: Any) -> list[dict[int, int]]:
    """Return target-index -> SMILES-index graph isomorphisms.

    PDB bond orders are deliberately ignored: supplying SMILES exists to
    replace those geometry-derived guesses, while connectivity and elements
    provide the only safe atom correspondence available without atom names in
    SMILES.
    """
    n = len(pdb_atoms)
    if mol.GetNumAtoms() != n:
        return []

    pdb_adj = [set() for _ in range(n)]
    for a, b in pdb_edges:
        pdb_adj[a].add(b)
        pdb_adj[b].add(a)
    smi_adj = [set() for _ in range(n)]
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        smi_adj[a].add(b)
        smi_adj[b].add(a)

    atomic_numbers = [atom.element.atomic_number for atom in pdb_atoms]
    candidates = {
        ti: [si for si, atom in enumerate(mol.GetAtoms())
             if atom.GetAtomicNum() == atomic_numbers[ti]
             and len(smi_adj[si]) == len(pdb_adj[ti])]
        for ti in range(n)
    }
    if any(not options for options in candidates.values()):
        return []

    order = sorted(range(n), key=lambda i: (len(candidates[i]), -len(pdb_adj[i]), i))
    assignments: list[dict[int, int]] = []
    current: dict[int, int] = {}
    used: set[int] = set()

    def visit(depth: int) -> None:
        if len(assignments) > 4096:
            raise SmilesPreparationError(
                "ligand graph has too many symmetric atom mappings to validate safely"
            )
        if depth == n:
            assignments.append(dict(current))
            return
        ti = order[depth]
        for si in candidates[ti]:
            if si in used:
                continue
            compatible = True
            for other_t, other_s in current.items():
                if ((other_t in pdb_adj[ti]) != (other_s in smi_adj[si])):
                    compatible = False
                    break
            if not compatible:
                continue
            current[ti] = si
            used.add(si)
            visit(depth + 1)
            used.remove(si)
            del current[ti]

    visit(0)
    return assignments


def _chemical_signature(mapping: dict[int, int], mol: Any) -> tuple[Any, ...]:
    inverse = {si: ti for ti, si in mapping.items()}
    atom_signatures = []
    for ti in sorted(mapping):
        atom = mol.GetAtomWithIdx(mapping[ti])
        atom_signatures.append((
            ti, atom.GetFormalCharge(), atom.GetIsAromatic(), atom.GetTotalNumHs(),
        ))
    bond_signatures = []
    for bond in mol.GetBonds():
        a = inverse[bond.GetBeginAtomIdx()]
        b = inverse[bond.GetEndAtomIdx()]
        bond_signatures.append((
            min(a, b), max(a, b), float(bond.GetBondTypeAsDouble()), bond.GetIsAromatic(),
        ))
    return tuple(atom_signatures), tuple(sorted(bond_signatures))


def _bond_metadata(rdkit_bond: Any) -> tuple[Any, float]:
    from openmm.app.topology import Aromatic, Double, Single, Triple

    order = float(rdkit_bond.GetBondTypeAsDouble())
    if rdkit_bond.GetIsAromatic():
        return Aromatic, 1.5
    if order == 2.0:
        return Double, 2.0
    if order == 3.0:
        return Triple, 3.0
    return Single, 1.0


def add_hydrogens_from_smiles(topology: Any, positions: Any,
                              smiles_by_resname: dict[str, str],
                              verbose: bool = False) -> tuple[Any, Any]:
    """Replace H on mapped isolated residues using authoritative SMILES.

    Heavy atoms, their PDB names, residue identities, ordering, and coordinates
    are retained.  Every matching residue instance must map safely or the
    operation fails before a new topology is constructed.
    """
    if not smiles_by_resname:
        return topology, positions

    from rdkit import Chem
    from rdkit.Geometry import Point3D

    from dvbfixer.ffutils import PROTEIN_RESIDUES, SOLVENT_IONS, is_glycam_residue

    forbidden = PROTEIN_RESIDUES | SOLVENT_IONS
    for resname in smiles_by_resname:
        if resname in forbidden or is_glycam_residue(resname):
            raise SmilesPreparationError(
                f"--smiles target {resname!r} is not an isolated non-GLYCAM ligand residue"
            )

    residues = list(topology.residues())
    target_residues = [r for r in residues if r.name in smiles_by_resname]
    present = {r.name for r in target_residues}
    missing = sorted(set(smiles_by_resname) - present)
    if missing:
        raise SmilesPreparationError(
            "no ligand residue found for --smiles target(s): " + ", ".join(missing)
        )

    all_bonds = list(topology.bonds())
    generated: dict[int, dict[str, Any]] = {}
    for residue in target_residues:
        heavy = [a for a in residue.atoms() if a.element.symbol != "H"]
        local_index = {atom.index: i for i, atom in enumerate(heavy)}
        edges: set[tuple[int, int]] = set()
        for bond in all_bonds:
            a, b = bond[0], bond[1]
            if a.element.symbol == "H" or b.element.symbol == "H":
                continue
            a_here = a.residue is residue
            b_here = b.residue is residue
            if a_here != b_here:
                raise SmilesPreparationError(
                    f"mapped ligand {residue.chain.id}/{residue.name}{residue.id} "
                    "has an external heavy-atom bond; version 1 supports isolated residues only"
                )
            if a_here and b_here:
                i, j = local_index[a.index], local_index[b.index]
                edges.add((min(i, j), max(i, j)))

        mol = _heavy_molecule(smiles_by_resname[residue.name])
        mappings = _graph_mappings(heavy, edges, mol)
        if not mappings:
            elements = "".join(a.element.symbol for a in heavy)
            raise SmilesPreparationError(
                f"PDB graph for {residue.chain.id}/{residue.name}{residue.id} "
                f"does not match its SMILES (PDB heavy elements: {elements})"
            )
        signatures = {_chemical_signature(mapping, mol) for mapping in mappings}
        if len(signatures) != 1:
            raise SmilesPreparationError(
                f"PDB graph for {residue.chain.id}/{residue.name}{residue.id} "
                "has multiple chemically distinct SMILES atom mappings"
            )
        mapping = mappings[0]
        inverse = {si: ti for ti, si in mapping.items()}

        conformer = Chem.Conformer(mol.GetNumAtoms())
        for target_i, smiles_i in mapping.items():
            p = positions[heavy[target_i].index]
            conformer.SetAtomPosition(smiles_i, Point3D(
                float(p.x) * 10.0, float(p.y) * 10.0, float(p.z) * 10.0,
            ))
        mol.RemoveAllConformers()
        mol.AddConformer(conformer, assignId=True)
        molh = Chem.AddHs(mol, addCoords=True)
        confh = molh.GetConformer()

        h_by_parent: dict[int, list[tuple[float, float, float]]] = {}
        for atom in molh.GetAtoms():
            if atom.GetAtomicNum() != 1:
                continue
            neighbors = list(atom.GetNeighbors())
            if len(neighbors) != 1 or neighbors[0].GetIdx() not in inverse:
                raise SmilesPreparationError(
                    f"RDKit produced an unassignable H for residue {residue.name}{residue.id}"
                )
            parent_target = inverse[neighbors[0].GetIdx()]
            hp = confh.GetAtomPosition(atom.GetIdx())
            xyz = (hp.x / 10.0, hp.y / 10.0, hp.z / 10.0)
            if all(abs(v) < 1e-8 for v in xyz):
                parent_pos = positions[heavy[parent_target].index]
                xyz = (float(parent_pos.x) + 0.10,
                       float(parent_pos.y) + 0.05,
                       float(parent_pos.z) + 0.05)
            h_by_parent.setdefault(heavy[parent_target].index, []).append(xyz)

        bond_specs = []
        for bond in mol.GetBonds():
            ta = inverse[bond.GetBeginAtomIdx()]
            tb = inverse[bond.GetEndAtomIdx()]
            bond_specs.append((heavy[ta].index, heavy[tb].index, *_bond_metadata(bond)))
        generated[residue.index] = {
            "hydrogens": h_by_parent,
            "bonds": bond_specs,
        }

    from openmm import Vec3
    from openmm.app import Topology, element
    from openmm.app.topology import Single
    from openmm.unit import Quantity, nanometer

    new_top = Topology()
    box_vectors = topology.getPeriodicBoxVectors()
    if box_vectors is not None:
        new_top.setPeriodicBoxVectors(box_vectors)
    new_positions = []
    old_to_new: dict[int, Any] = {}
    new_h_bonds: list[tuple[Any, Any]] = []
    added = 0

    for chain in topology.chains():
        new_chain = new_top.addChain(chain.id)
        for residue in chain.residues():
            new_residue = new_top.addResidue(
                residue.name, new_chain, residue.id, residue.insertionCode,
            )
            managed = generated.get(residue.index)
            used_names = {a.name for a in residue.atoms() if a.element.symbol != "H"}
            for atom in residue.atoms():
                if managed is not None and atom.element.symbol == "H":
                    continue
                new_atom = new_top.addAtom(atom.name, atom.element, new_residue)
                old_to_new[atom.index] = new_atom
                p = positions[atom.index]
                new_positions.append(Vec3(float(p.x), float(p.y), float(p.z)))
                if managed is None:
                    continue
                h_positions = managed["hydrogens"].get(atom.index, [])
                if atom.name.startswith(("O", "N", "S")):
                    base = f"H{atom.name}"
                elif len(atom.name) > 1:
                    base = f"H{atom.name[1:]}"
                else:
                    base = "H"
                for ordinal, h_pos in enumerate(h_positions, 1):
                    candidate = base if ordinal == 1 else f"{base}{ordinal}"
                    suffix = ordinal
                    while candidate in used_names:
                        suffix += 1
                        candidate = f"{base}{suffix}"
                    used_names.add(candidate)
                    h_atom = new_top.addAtom(candidate, element.hydrogen, new_residue)
                    new_positions.append(Vec3(*h_pos))
                    new_h_bonds.append((new_atom, h_atom))
                    added += 1

    for bond in all_bonds:
        a, b = bond[0], bond[1]
        if a.residue.index in generated and b.residue is a.residue:
            continue
        na, nb = old_to_new.get(a.index), old_to_new.get(b.index)
        if na is not None and nb is not None:
            new_top.addBond(na, nb, type=bond.type, order=bond.order)
    for residue_index, specs in generated.items():
        del residue_index
        for a_idx, b_idx, bond_type, order in specs["bonds"]:
            new_top.addBond(old_to_new[a_idx], old_to_new[b_idx],
                            type=bond_type, order=order)
    for parent, hydrogen in new_h_bonds:
        new_top.addBond(parent, hydrogen, type=Single, order=1.0)

    if verbose:
        print(f"  SMILES added {added} H atoms to {len(target_residues)} ligand residue(s)")
    return new_top, Quantity(new_positions, nanometer)
