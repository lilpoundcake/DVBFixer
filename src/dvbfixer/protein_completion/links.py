"""Strict PDB link metadata transport; no chemistry inference or symmetry expansion."""
from __future__ import annotations

import math
from dataclasses import dataclass

from dvbfixer.domain.structure_identity import ExactAtomRef
from dvbfixer.pdb_component_reconstruction import PdbAtom, _component_ref


@dataclass(frozen=True)
class DeclaredLink:
    record: str
    left: ExactAtomRef
    right: ExactAtomRef
    reported_distance: float | None


def parse_declared_links(lines: list[str], atoms: list[PdbAtom]) -> tuple[DeclaredLink, ...]:
    """Resolve both endpoints exactly, including chain case and insertion code.

    No missing endpoint is invented. Blank/1555 symmetry means the deposited
    frame; any other crystallographic operation requires separate support.
    """
    index = {(a.chain_id, a.sequence_number, a.insertion_code, a.name): a for a in atoms}
    if len(index) != len(atoms):
        raise ValueError("ambiguous declared-link endpoints")
    links = []
    for line in lines:
        record = line[:6].strip()
        if record not in {"LINK", "SSBOND"}:
            continue
        if len(line) > 6 and line[6] not in {" ", "\n"}:
            raise ValueError("malformed LINK/SSBOND prefix")
        if len(line.rstrip("\n")) < 72:
            raise ValueError("truncated LINK/SSBOND record")
        if any(line[start:stop].strip() not in {"", "1555"} for start, stop in ((59, 65), (66, 72))):
            raise ValueError("nonidentity link symmetry requires crystal-frame support")

        def endpoint(chain: str, number: str, code: str, name: str, residue_name: str) -> ExactAtomRef:
            atom = index.get((chain, number.strip(), code.strip(), name.strip()))
            if atom is None or atom.residue_name != residue_name.strip() or atom.alternate_location:
                raise ValueError("missing/mismatched LINK/SSBOND endpoint")
            return ExactAtomRef(_component_ref(atom), atom.name)

        if record == "LINK":
            if line[16].strip() or line[46].strip():
                raise ValueError("LINK altloc requires explicit selection")
            left = endpoint(line[21], line[22:26], line[26], line[12:16], line[17:20])
            right = endpoint(line[51], line[52:56], line[56], line[42:46], line[47:50])
        else:
            try:
                serial = int(line[7:10])
            except ValueError as exc:
                raise ValueError("invalid SSBOND serial") from exc
            if not 1 <= serial <= 999 or line[11:14] != "CYS" or line[25:28] != "CYS":
                raise ValueError("SSBOND requires canonical CYS endpoints and serial")
            left = endpoint(line[15], line[17:21], line[21], "SG", "CYS")
            right = endpoint(line[29], line[31:35], line[35], "SG", "CYS")
        if left == right:
            raise ValueError("self LINK/SSBOND")
        field = line[73:78].strip()
        try:
            distance = float(field) if field else None
        except ValueError as exc:
            raise ValueError("invalid LINK/SSBOND distance") from exc
        if distance is not None and (not math.isfinite(distance) or not 0 < distance < 10):
            raise ValueError("invalid LINK/SSBOND distance")
        links.append(DeclaredLink(record, left, right, distance))
    return tuple(links)
