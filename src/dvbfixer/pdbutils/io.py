"""Line-level PDB read/write helpers.

Split out of ``pdbutils.py`` in Phase 1.4 of the revision plan. These
are the small utilities that touch PDB text directly — CONECT serial
remap, atom-serial map, insert-before-END — used by prepare / pull /
transplant when merging CONECT from a differently-serialised source.

The inference engine (``infer_conect_records`` + the OpenBabel wiring
+ domain overrides) lives in :mod:`dvbfixer.pdbutils.inference`. All
public names are re-exported from :mod:`dvbfixer.pdbutils` so callers
using ``from dvbfixer.pdbutils import X`` keep working.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path


def update_conect(line: str, serial_map: Mapping[int, int]) -> str | None:
    """CONECT: remap atom serial numbers. Each serial is 5 chars starting at col 6.

    Returns None if any referenced serial isn't a real atom in this file
    (not in `serial_map`) — some deposited PDBs carry stale/dangling
    CONECT records left over from atoms that no longer exist (confirmed
    on a real structure: an orphaned "CONECT 5111 5112 5113"-style
    record with no matching ATOM/HETATM line at all). Silently passing
    such an old serial number through UNCHANGED (the previous fallback,
    ``serial_map.get(old_serial, old_serial)``) is actively harmful once
    renumbered: that stale number can coincide with the NEW serial
    dvbfixer assigns to a real, unrelated atom, fabricating a spurious
    bond between atoms that were never actually connected. Dropping the
    whole record is safer than keeping a partially-wrong one.
    """
    parts = []
    parts.append(line[:6])  # "CONECT"
    i = 6
    while i + 5 <= len(line.rstrip('\n')):
        s = line[i:i + 5].strip()
        if s and s.isdigit():
            old_serial = int(s)
            new_serial = serial_map.get(old_serial)
            if new_serial is None:
                return None
            parts.append(f"{new_serial:5d}")
        else:
            parts.append(line[i:i + 5])
        i += 5
    return ''.join(parts) + '\n'


def build_serial_map(pdb_path: str | Path) -> dict[tuple[str, str, str], int]:
    """Build ``(chain, resid, atomname) -> serial`` from a PDB file.

    Used with :func:`remap_conect_records` to bridge CONECT records
    between two PDBs whose atoms overlap by identity but not by serial
    numbering (typical after concatenating an acceptor + a graft).
    """
    serial_map: dict[tuple[str, str, str], int] = {}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                serial = int(line[6:11])
                chain = line[21]
                resid = line[22:26].strip()
                atomname = line[12:16].strip()
                serial_map[(chain, resid, atomname)] = serial
    return serial_map


def remap_conect_records(
    input_path: str | Path,
    new_serial_map: dict[tuple[str, str, str], int],
) -> list[str]:
    """Read CONECT from ``input_path``, remap each serial to
    ``new_serial_map[(chain, resid, atomname)]``.

    Skips CONECT lines that reference atoms not present in the target
    map (typical when a residue was dropped). Returns the remapped
    CONECT lines ready to be written into the target file.
    """
    with open(input_path) as f:
        lines = f.readlines()

    old_serial_to_key: dict[int, tuple[str, str, str]] = {}
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            serial = int(line[6:11])
            chain = line[21]
            resid = line[22:26].strip()
            atomname = line[12:16].strip()
            old_serial_to_key[serial] = (chain, resid, atomname)

    result: list[str] = []
    for line in lines:
        if not line.startswith("CONECT"):
            continue
        serials: list[int] = []
        s = line[6:]
        while len(s) >= 5:
            chunk = s[:5].strip()
            if chunk:
                serials.append(int(chunk))
            s = s[5:]
        if len(serials) < 2:
            continue

        new_serials: list[int] = []
        all_mapped = True
        for old_s in serials:
            key = old_serial_to_key.get(old_s)
            if key and key in new_serial_map:
                new_serials.append(new_serial_map[key])
            else:
                all_mapped = False
                break

        if all_mapped and len(new_serials) >= 2:
            conect = f"CONECT{new_serials[0]:5d}"
            for ns in new_serials[1:]:
                conect += f"{ns:5d}"
            result.append(conect.ljust(80) + "\n")

    return result


def append_before_end(output_path: str | Path, extra_lines: list[str]) -> None:
    """Insert ``extra_lines`` immediately before the ``END`` record.

    No-op if ``extra_lines`` is empty. Used to append CONECT / TER
    blocks after the atom section without disturbing the ``END``
    sentinel.
    """
    if not extra_lines:
        return
    with open(output_path) as f:
        lines = f.readlines()
    with open(output_path, "w") as f:
        for line in lines:
            if line.startswith("END"):
                f.writelines(extra_lines)
            f.write(line)
