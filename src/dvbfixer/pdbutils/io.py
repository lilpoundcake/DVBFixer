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

from pathlib import Path


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
