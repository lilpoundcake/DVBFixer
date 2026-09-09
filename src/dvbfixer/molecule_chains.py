"""PDB adapter for assigning unique chain IDs to molecular components."""

from __future__ import annotations

import re
from collections import defaultdict
from urllib.parse import unquote

from dvbfixer.domain.structure_identity import allocate_chain_ids

_ASSEMBLY_CHAINS = re.compile(r"CHAINS?:\s*(.+)$")
_EXCLUDED_COMPONENTS = {"HOH", "WAT", "TIP", "TIP3", "SOL", "NA", "CL", "K", "MG", "CA", "ZN"}


def _residue_key(line: str) -> tuple[str, str, str, str]:
    return (line[21:22], line[22:26], line[26:27], line[17:20].strip())


def _metadata_chain_ids(lines: list[str]) -> set[str]:
    reserved: set[str] = set()
    for line in lines:
        if not line.startswith("REMARK 350"):
            continue
        match = _ASSEMBLY_CHAINS.search(line)
        if match:
            reserved.update(
                token.strip() for token in match.group(1).split(",")
                if len(token.strip()) == 1
            )
    return reserved


def assign_unique_molecule_chains(
    lines: list[str], alphabet: list[str], *, source_lines: list[str] | None = None
) -> tuple[list[str], int]:
    """Give each connected non-polymer component a unique PDB chain identifier.

    Heterogen-to-heterogen CONECT records join residues into one molecule. Links
    to polymer atoms are deliberately ignored for ownership. Waters and ions
    have already been removed by the split pipeline and are never reassigned.
    """
    hetero_residues: dict[tuple[str, str, str, str], list[int]] = defaultdict(list)
    polymer_chains: set[str] = set()
    for line in lines:
        if line.startswith("ATOM  "):
            polymer_chains.add(line[21:22])
        elif line.startswith("HETATM"):
            try:
                serial = int(line[6:11])
            except ValueError:
                continue
            key = _residue_key(line)
            hetero_residues[key].append(serial)

    if not hetero_residues:
        return lines, 0

    parent = {key: key for key in hetero_residues}

    def find(key):  # noqa: ANN001, ANN202
        while parent[key] != key:
            parent[key] = parent[parent[key]]
            key = parent[key]
        return key

    def union(left, right):  # noqa: ANN001, ANN202
        left, right = find(left), find(right)
        if left != right:
            parent[right] = left

    # mmCIF normalization emits authoritative label_asym component identities
    # in residue traversal order. Join residues carrying the same label.
    cif_labels = [
        unquote(line.split(maxsplit=4)[4].strip())
        for line in lines if line.startswith("REMARK 999 DVBFIXER CIF_COMPONENT ")
    ]
    ordered_residues = list(hetero_residues)
    if cif_labels and len(cif_labels) == len(ordered_residues):
        first_by_label: dict[str, tuple[str, str, str, str]] = {}
        for key, label in zip(ordered_residues, cif_labels):
            if label in first_by_label:
                union(first_by_label[label], key)
            else:
                first_by_label[label] = key

    # The legacy split renderer preserves CONECT records byte-for-byte while
    # reserializing coordinates. Consequently those serials cannot safely be
    # interpreted after rendering. PDB residue boundaries are therefore the
    # conservative component boundary here. The mmCIF adapter can supply its
    # authoritative label_asym grouping before this stage in a future writer.

    groups: dict[tuple[str, str, str, str], list[tuple[str, str, str, str]]] = defaultdict(list)
    for key in hetero_residues:
        groups[find(key)].append(key)
    components = sorted(groups.values(), key=lambda group: min(hetero_residues[key][0] for key in group))
    reserved = polymer_chains | _metadata_chain_ids(lines)
    assigned = allocate_chain_ids(len(components), alphabet=alphabet, reserved=reserved)
    replacement = {key: chain for group, chain in zip(components, assigned) for key in group}
    source_target: dict[
        tuple[str, str, str, str], tuple[str, str, str, str]
    ] = {}
    if source_lines is not None:
        source_residues = list(dict.fromkeys(
            _residue_key(line) for line in source_lines
            if line.startswith("HETATM") and line[17:20].strip() not in _EXCLUDED_COMPONENTS
        ))
        if len(source_residues) == len(ordered_residues):
            source_target = {
                source: rendered
                for source, rendered in zip(source_residues, ordered_residues)
            }
    rewritten: list[str] = []
    inserted_remark = False
    for line in lines:
        new_line = line
        if line.startswith("HETATM"):
            new_line = line[:21] + replacement[_residue_key(line)] + line[22:]
        elif line.startswith("ANISOU"):
            key = _residue_key(line)
            target = source_target.get(key, key)
            chain = replacement.get(target)
            if chain:
                new_line = line[:21] + chain + target[1] + target[2] + line[27:]
        elif line.startswith("HET   "):
            key = _residue_key(line)
            target = source_target.get(key, key)
            chain = replacement.get(target)
            if chain:
                new_line = line[:21] + chain + target[1] + target[2] + line[27:]
        elif line.startswith("LINK  "):
            chars = list(line)
            for chain_pos, resid_slice, icode_pos, resname_slice in (
                (21, slice(22, 26), 26, slice(17, 20)),
                (51, slice(52, 56), 56, slice(47, 50)),
            ):
                try:
                    key = (line[chain_pos], line[resid_slice], line[icode_pos], line[resname_slice].strip())
                except IndexError:
                    continue
                target = source_target.get(key, key)
                chain = replacement.get(target)
                if chain:
                    chars[chain_pos] = chain
                    chars[resid_slice] = target[1]
                    chars[icode_pos] = target[2]
            new_line = "".join(chars)
        if not inserted_remark and line.startswith(("ATOM  ", "HETATM", "MODEL ")):
            rewritten.append(
                f"REMARK 999 DVBFIXER UNIQUE_MOLECULE_CHAINS {len(components)} COMPONENTS\n"
            )
            inserted_remark = True
        rewritten.append(new_line)
    return rewritten, len(components)
