"""Fail-closed authoritative fixed-ligand context for diffusion sampling."""

from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from types import SimpleNamespace

from dvbfixer.model.diffusion.contract import (
    AtomIdentity,
    HeterogenContext,
    ResidueIdentity,
)
from dvbfixer.prepare.smiles import (
    SmilesPreparationError,
    _chemical_signature,
    _graph_mappings,
    _heavy_molecule,
)


class HeterogenContextError(ValueError):
    """A retained heterogen cannot be represented with authoritative chemistry."""


def build_smiles_heterogen_contexts(
    pdb_text: str,
    smiles_by_resname: dict[str, str],
) -> tuple[HeterogenContext, ...]:
    """Map isolated HETATM residues to explicit authoritative SMILES graphs."""
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise HeterogenContextError(
            "heterogen-conditioned diffusion requires RDKit"
        ) from exc

    serial_atoms: dict[int, tuple[ResidueIdentity, AtomIdentity, str]] = {}
    residue_atoms: dict[ResidueIdentity, list[tuple[int, AtomIdentity, str]]] = defaultdict(list)
    residue_names: dict[ResidueIdentity, str] = {}
    conect: set[tuple[int, int]] = set()
    for line_number, line in enumerate(pdb_text.splitlines(), 1):
        record = line[:6].strip()
        if record == "HETATM":
            if len(line) < 78:
                raise HeterogenContextError(
                    f"truncated HETATM record on line {line_number}"
                )
            if line[16].strip():
                raise HeterogenContextError(
                    "heterogen-conditioned diffusion rejects alternate locations"
                )
            try:
                serial = int(line[6:11])
                occupancy = float(line[54:60])
            except ValueError as exc:
                raise HeterogenContextError(
                    f"malformed HETATM record on line {line_number}"
                ) from exc
            if occupancy != 1.0:
                raise HeterogenContextError(
                    "heterogen-conditioned diffusion requires occupancy 1.00"
                )
            element = line[76:78].strip().upper()
            if not element or element == "H":
                raise HeterogenContextError(
                    "heterogen context must contain explicit heavy elements and no hydrogens"
                )
            residue = ResidueIdentity(line[21], line[22:26].strip(), line[26].strip())
            identity = AtomIdentity(
                residue.chain,
                residue.residue_number,
                residue.insertion_code,
                line[12:16].strip(),
            )
            resname = line[17:20].strip()
            if serial in serial_atoms:
                raise HeterogenContextError("heterogen atom serials must be unique")
            serial_atoms[serial] = (residue, identity, element)
            residue_atoms[residue].append((serial, identity, element))
            previous = residue_names.setdefault(residue, resname)
            if previous != resname:
                raise HeterogenContextError("heterogen residue name is inconsistent")
        elif record == "CONECT":
            try:
                serials = [
                    int(line[index : index + 5])
                    for index in range(6, len(line), 5)
                    if line[index : index + 5].strip()
                ]
            except ValueError as exc:
                raise HeterogenContextError("malformed heterogen CONECT record") from exc
            if serials:
                conect.update(
                    tuple(sorted((serials[0], other)))
                    for other in serials[1:]
                    if other != serials[0]
                )

    if not residue_atoms:
        if smiles_by_resname:
            raise HeterogenContextError(
                "--diffusion-heterogen-smiles was supplied but input has no HETATM residue"
            )
        return ()
    missing_authority = sorted(set(residue_names.values()) - set(smiles_by_resname))
    if missing_authority:
        raise HeterogenContextError(
            "retained heterogens require authoritative SMILES mappings: "
            + ", ".join(missing_authority)
        )

    contexts: list[HeterogenContext] = []
    periodic_table = Chem.GetPeriodicTable()
    for residue, entries in residue_atoms.items():
        resname = residue_names[residue]
        mol = _heavy_molecule(smiles_by_resname[resname])
        serial_to_index = {serial: index for index, (serial, _identity, _element) in enumerate(entries)}
        internal_edges: set[tuple[int, int]] = set()
        entry_serials = set(serial_to_index)
        for first, second in conect:
            if bool(first in entry_serials) != bool(second in entry_serials):
                raise HeterogenContextError(
                    "covalently attached heterogens are unsupported by the fixed-context slice"
                )
            if first in entry_serials:
                internal_edges.add(tuple(sorted((
                    serial_to_index[first],
                    serial_to_index[second],
                ))))
        pdb_atoms = [
            SimpleNamespace(
                element=SimpleNamespace(
                    atomic_number=periodic_table.GetAtomicNumber(element.title())
                )
            )
            for _serial, _identity, element in entries
        ]
        try:
            mappings = _graph_mappings(pdb_atoms, internal_edges, mol)
            signatures = {_chemical_signature(mapping, mol) for mapping in mappings}
        except SmilesPreparationError as exc:
            raise HeterogenContextError(str(exc)) from exc
        if not mappings:
            raise HeterogenContextError(
                f"SMILES graph is incompatible with retained residue {resname}"
            )
        if len(signatures) != 1:
            raise HeterogenContextError(
                f"SMILES mapping is chemically ambiguous for retained residue {resname}"
            )
        canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        identities = tuple(identity for _serial, identity, _element in entries)
        digest_payload = {
            "canonical_smiles": canonical,
            "atoms": [
                [atom.chain, atom.residue_number, atom.insertion_code, atom.atom_name]
                for atom in identities
            ],
        }
        graph_sha256 = hashlib.sha256(
            json.dumps(
                digest_payload,
                ensure_ascii=False,
                sort_keys=True,
                separators=(",", ":"),
            ).encode("utf-8")
        ).hexdigest()
        contexts.append(HeterogenContext(
            residue=residue,
            residue_name=resname,
            atoms=identities,
            canonical_smiles=canonical,
            authoritative_graph_sha256=graph_sha256,
        ))
    return tuple(contexts)
