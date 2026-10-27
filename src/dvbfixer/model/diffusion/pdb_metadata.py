"""Reconcile dependent PDB metadata with the selected target and actual records.

Coordinates, atom identity, occupancy, ANISOU and valid covalent links are not
rewritten. Source experimental metadata remains recoverable from the digest-bound
input artifact; absent-chain/range annotations must not describe the new model.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

from dvbfixer.biological_assembly import filter_biological_assembly_records
from dvbfixer.model.diffusion.contract import DiffusionContractError, DiffusionRequest
from dvbfixer.model.diffusion.scope import _CANONICAL_ONE_TO_THREE

ResidueKey = tuple[str, str, str]
# chain, resSeq start/stop, insertion-code column, residue-name start/stop
ResidueField = tuple[int, int, int, int, int, int]
_FIELDS: dict[str, tuple[ResidueField, ...]] = {
    "HELIX": ((19, 21, 25, 25, 15, 18), (31, 33, 37, 37, 27, 30)),
    "SHEET": ((21, 22, 26, 26, 17, 20), (32, 33, 37, 37, 28, 31)),
    "TURN": ((19, 20, 24, 24, 15, 18), (30, 31, 35, 35, 26, 29)),
    "CISPEP": ((15, 17, 21, 21, 11, 14), (29, 31, 35, 35, 25, 28)),
    "SSBOND": ((15, 17, 21, 21, 11, 14), (29, 31, 35, 35, 25, 28)),
    "LINK": ((21, 22, 26, 26, 17, 20), (51, 52, 56, 56, 47, 50)),
    "MODRES": ((16, 18, 22, 22, 12, 15),),
    "SEQADV": ((16, 18, 22, 22, 12, 15),),
    "HET": ((12, 13, 17, 17, 7, 10),),
    "DBREF": ((12, 14, 18, 18, 0, 0), (12, 20, 24, 24, 0, 0)),
    "DBREF1": ((12, 14, 18, 18, 0, 0), (12, 20, 24, 24, 0, 0)),
}
_SHEET_REGISTRATION: tuple[ResidueField, ...] = ((49, 50, 54, 54, 45, 48), (64, 65, 69, 69, 60, 63))
_SITE_FIELDS: tuple[ResidueField, ...] = ((22, 23, 27, 27, 18, 21), (33, 34, 38, 38, 29, 32),
                                        (44, 45, 49, 49, 40, 43), (55, 56, 60, 60, 51, 54))


@dataclass(frozen=True)
class MetadataNormalization:
    text: str
    changes: tuple[str, ...]


def reconcile_pdb_metadata(
    text: str, request: DiffusionRequest, *, add_missing_seqres: bool = True, annotate: bool = True,
) -> MetadataNormalization:
    """Normalize only sequence/identity-dependent records, not arbitrary headers.

    Validation uses add_missing_seqres=False: absent optional annotations are
    allowed, but supplied annotations must be consistent. Production materializers
    always emit the authoritative request sequence, including cropped targets.
    """
    lines = text.splitlines(keepends=True)
    targets = {target.chain: tuple(_CANONICAL_ONE_TO_THREE[letter] for letter in target.sequence)
               for target in request.target_sequences}
    residues: dict[ResidueKey, str] = {}
    atoms: dict[ResidueKey, set[str]] = defaultdict(set)
    for line in lines:
        if line[:6].strip() in {"ATOM", "HETATM"} and len(line) >= 27:
            key = (line[21], line[22:26].strip(), line[26].strip())
            residues[key] = line[17:20].strip()
            atoms[key].add(line[12:16].strip())
    chains = {key[0] for key in residues}
    changes: list[str] = []

    def valid(line: str, fields: tuple[ResidueField, ...]) -> bool:
        for chain, start, stop, code, name_start, name_stop in fields:
            if len(line.rstrip("\r\n")) < stop:
                return False
            key = (line[chain], line[start:stop].strip(), line[code:code + 1].strip())
            if key not in residues or (name_stop and line[name_start:name_stop].strip() != residues[key]):
                return False
        return True

    seqres: dict[str, list[str]] = defaultdict(list)
    for line in lines:
        if line.startswith("SEQRES"):
            seqres[line[11:12]].append(line)
    replacement_seqres: list[str] = []
    for chain in dict.fromkeys((*seqres, *(targets if add_missing_seqres else ()))):
        old = seqres.get(chain, [])
        if chain not in chains:
            changes.extend(["SEQRES"] * len(old))
            continue
        if chain not in targets:
            replacement_seqres.extend(old)
            continue
        expected = targets[chain]
        try:
            consistent = (tuple(name for line in old for name in line[19:].split()) == expected and
                          all(int(line[13:17]) == len(expected) and int(line[7:10]) == i
                              for i, line in enumerate(old, 1)))
        except ValueError:
            consistent = False
        if consistent:
            replacement_seqres.extend(old)
        else:
            if len(expected) > 9999:
                raise DiffusionContractError("target SEQRES exceeds PDB residue-count capacity")
            replacement_seqres.extend(f"SEQRES {i // 13 + 1:3d} {chain} {len(expected):4d}  " + " ".join(expected[i:i + 13]) + "\n"
                                      for i in range(0, len(expected), 13))
            changes.append(f"SEQRES:{chain}")

    assembly = filter_biological_assembly_records(lines, chains)
    original_assembly = [line for line in lines if line.startswith("REMARK 350")]
    if assembly != original_assembly:
        changes.append("REMARK350")

    # Sheets/sites are relational, multi-record annotations. Drop incomplete
    # groups rather than inventing strand registration or an incomplete site.
    groups: dict[tuple[str, str], list[str]] = defaultdict(list)
    for line in lines:
        record = line[:6].strip()
        if record in {"SHEET", "SITE"}:
            groups[record, line[11:14]].append(line)
    valid_groups = set()
    for (record, identifier), records in groups.items():
        fields = _FIELDS[record] if record == "SHEET" else _SITE_FIELDS
        count = 0
        good = True
        for line in records:
            used = fields if record == "SHEET" else tuple(field for field in fields if line[field[4]:field[5]].strip())
            good = good and bool(used) and valid(line, used)
            count += len(used) if record == "SITE" else 1
            if record == "SHEET" and line[38:40].strip() not in {"", "0"}:
                good = good and valid(line, _SHEET_REGISTRATION)
        try:
            good = good and all(int(line[14:16] if record == "SHEET" else line[15:17]) == count for line in records)
        except ValueError:
            good = False
        if good:
            valid_groups.add((record, identifier))

    dbref2_keys = {(line[7:11], line[12:13]) for line in lines if line.startswith("DBREF2")}
    valid_dbref1 = {(line[7:11], line[12:13]) for line in lines if line.startswith("DBREF1") and valid(line, _FIELDS["DBREF1"])} & dbref2_keys
    heterogens = {line[17:20].strip() for line in lines if line.startswith("HETATM")}
    output: list[str] = []
    sequence_written = False
    assembly_written = False
    last_coordinate: str | None = None
    for original in lines:
        line = original
        record = line[:6].strip()
        if record == "SEQRES":
            if not sequence_written:
                output.extend(replacement_seqres)
                sequence_written = True
            continue
        if not seqres and not sequence_written and record in {"MODRES", "HET", "HETNAM", "FORMUL", "HELIX", "SHEET", "SSBOND", "CRYST1", "MODEL", "ATOM", "HETATM", "END"}:
            output.extend(replacement_seqres)
            sequence_written = True
        if line.startswith("REMARK 350"):
            if not assembly_written:
                output.extend(assembly)
                assembly_written = True
            continue
        keep = True
        if record in {"SHEET", "SITE"}:
            keep = (record, line[11:14]) in valid_groups
        elif record in _FIELDS:
            keep = valid(line, _FIELDS[record])
            if record == "DBREF1":
                keep = keep and (line[7:11], line[12:13]) in valid_dbref1
            if record == "LINK" and keep:
                keep = all(line[start:stop].strip() in atoms[(line[chain], line[number:number + 4].strip(), line[code:code + 1].strip())]
                           for start, stop, chain, number, code in ((12, 16, 21, 22, 26), (42, 46, 51, 52, 56)))
        elif record == "DBREF2":
            keep = (line[7:11], line[12:13]) in valid_dbref1
        elif record in {"HETNAM", "HETSYN", "FORMUL"}:
            keep = line[12:15].strip() in heterogens if record == "FORMUL" else line[11:14].strip() in heterogens
        elif line.startswith(("REMARK 465", "REMARK 470")):
            # Table prose is retained; interpret only residue data rows.
            is_missing_residue = line.startswith("REMARK 465")
            start, stop, code = (21, 26, 26) if is_missing_residue else (20, 24, 24)
            number = line[start:stop].strip()
            if number.lstrip("-").isdigit() and line[15:18].strip():
                key = (line[19:20], number, line[code:code + 1].strip())
                if is_missing_residue:
                    keep = key not in residues and key[0] in chains and key[0] not in targets
                else:
                    missing = [name for name in line[25:].split() if name not in atoms.get(key, set())]
                    keep = key in residues and bool(missing)
                    if keep and missing != line[25:].split():
                        line = line[:25] + " " + " ".join(f"{name:>4}" for name in missing) + "\n"
        elif record == "TER" and last_coordinate is not None and line[17:20].strip():
            body = line.rstrip("\r\n").ljust(27)
            if body[17:27] != last_coordinate[17:27]:
                line = body[:17] + last_coordinate[17:27] + body[27:] + ("\r\n" if line.endswith("\r\n") else "\n")
        if not keep:
            changes.append(record if record != "REMARK" else line[:10].strip())
            continue
        if line != original:
            changes.append(record)
        output.append(line)
        if record in {"ATOM", "HETATM"}:
            last_coordinate = line
        elif record in {"TER", "ENDMDL", "MODEL"}:
            last_coordinate = None
    if not sequence_written:
        output[0:0] = replacement_seqres
    if changes and annotate:
        note = f"REMARK 999 DVBFIXER METADATA_SYNC {len(changes)} dependent record update(s)\n"
        first_atom = next((i for i, line in enumerate(output) if line[:6].strip() in {"MODEL", "ATOM", "HETATM"}), len(output))
        output.insert(first_atom, note)
    return MetadataNormalization("".join(output), tuple(changes))
