"""Independent research-only RDKit CIP labels; never a coordinate generator.

Use a local free-amino-acid graph to order CA/CB ligands. The auxiliary terminal
oxygen is solely a graph/conformer placeholder; no oracle coordinate enters a
proposal. Atom maps transport coordinates only and are removed before CIP.
"""
from __future__ import annotations

from collections.abc import Mapping

from dvbfixer.protein_completion.contracts import Point


def rdkit_labels(points: Mapping[str, Point], residue_name: str) -> dict[str, str]:
    from rdkit import Chem
    from rdkit.Chem import rdCIPLabeler

    if residue_name not in {"VAL", "THR"}:
        raise ValueError("branched stereo oracle supports VAL/THR only")
    first = "[OH:7]" if residue_name == "THR" else "[CH3:7]"
    molecule = Chem.MolFromSmiles(f"[NH2:1][CH:2]([C:3](=[O:4])[OH:5])[CH:6]({first})[CH3:8]")
    if molecule is None:
        raise ValueError("independent CIP graph unavailable")
    names = {1: "N", 2: "CA", 3: "C", 4: "O", 5: "OXT", 6: "CB",
             7: "OG1" if residue_name == "THR" else "CG1", 8: "CG2"}
    xyz = dict(points)
    xyz["OXT"] = tuple(2*points["C"][i]-points["O"][i] for i in range(3))
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    conformer.Set3D(True)
    atom_names = {}
    for atom in molecule.GetAtoms():
        name = names[atom.GetAtomMapNum()]
        atom_names[atom.GetIdx()] = name
        conformer.SetAtomPosition(atom.GetIdx(), xyz[name])
        atom.SetAtomMapNum(0)  # Never let coordinate-index labels alter ligand ranking.
    molecule.AddConformer(conformer)
    Chem.AssignAtomChiralTagsFromStructure(molecule, replaceExistingTags=True)
    rdCIPLabeler.AssignCIPLabels(molecule)
    return {atom_names[a.GetIdx()]: a.GetProp("_CIPCode") for a in molecule.GetAtoms() if a.HasProp("_CIPCode")}
