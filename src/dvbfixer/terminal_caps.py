"""Build neutral ACE/NME peptide terminal caps in PDB files."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from dvbfixer.tempfiles import make_temp_path

PROTEIN_RESNAMES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "ASH",
    "GLH", "CYX", "CYM", "LYN", "LSN", "MSE", "NLN", "OLS", "OLT",
}
CAP_RESNAMES = {"ACE", "NME", "NHE", "NH2"}


@dataclass
class _Residue:
    chain: str
    resid: int
    icode: str
    resname: str
    line_indices: list[int]
    atoms: dict[str, tuple[int, np.ndarray]]


# Ideal ff19SB fragments.  The three anchor atoms define a local peptide
# backbone frame; only cap heavy atoms are copied into the target structure.
_N_ANCHORS = np.array(((3.555, 3.970, 0.0), (4.853, 4.614, 0.0),
                       (4.713, 6.129, 0.0)))
_ACE_HEAVY = {
    "CH3": np.array((2.000, 2.090, 0.0)),
    "C": np.array((3.427, 2.641, 0.0)),
    "O": np.array((4.391, 1.877, 0.0)),
}
_C_ANCHORS = np.array(((5.846, 6.835, 0.0), (5.846, 8.284, 0.0),
                       (7.273, 8.814, 0.0)))
_NME_HEAVY = {
    "N": np.array((7.417, 10.142, 0.0)),
    "C": np.array((8.722, 10.770, 0.0)),
}


def _transform(reference: np.ndarray, target: np.ndarray,
               points: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Rigidly align a reference backbone frame to target coordinates."""
    ref_center = reference.mean(axis=0)
    target_center = target.mean(axis=0)
    covariance = (reference - ref_center).T @ (target - target_center)
    u, _s, vt = np.linalg.svd(covariance)
    rotation = vt.T @ u.T
    if np.linalg.det(rotation) < 0:
        vt[-1, :] *= -1
        rotation = vt.T @ u.T
    return {
        name: (coord - ref_center) @ rotation.T + target_center
        for name, coord in points.items()
    }


def _parse_residues(lines: list[str]) -> list[_Residue]:
    residues: list[_Residue] = []
    last_key = None
    segment = 0
    for index, line in enumerate(lines):
        if line.startswith("TER"):
            segment += 1
            last_key = None
            continue
        if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 54:
            continue
        try:
            resid = int(line[22:26])
            xyz = np.array((float(line[30:38]), float(line[38:46]),
                            float(line[46:54])))
            serial = int(line[6:11])
        except ValueError:
            continue
        key = (segment, line[21], resid, line[26], line[17:20].strip())
        if key != last_key:
            residues.append(_Residue(
                chain=line[21], resid=resid, icode=line[26].strip(),
                resname=line[17:20].strip(), line_indices=[], atoms={},
            ))
            last_key = key
        residue = residues[-1]
        residue.line_indices.append(index)
        residue.atoms[line[12:16].strip()] = (serial, xyz)
    return residues


def _atom_line(serial: int, atom: str, resname: str, chain: str, resid: int,
               xyz: np.ndarray, element: str) -> str:
    # PDB atom names beginning with a letter occupy columns 14-16 for
    # one-character elements.  This also lets OpenMM apply pdbNames.xml
    # aliases consistently (ACE/NME H1/H2/H3 on the later H-addition pass).
    atom_field = f" {atom:<3s}" if len(atom) < 4 else f"{atom:<4s}"
    return (f"ATOM  {serial:5d} {atom_field}{' ':1s}{resname:>3s} {chain}"
            f"{resid:4d}    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
            f"  1.00  0.00          {element:>2s}  \n")


def _allocate_resid(used: set[int], start: int, direction: int) -> int:
    candidate = start + direction
    while -999 <= candidate <= 9999 and candidate in used:
        candidate += direction
    if not -999 <= candidate <= 9999:
        raise ValueError("no unused PDB residue number is available for a terminal cap")
    used.add(candidate)
    return candidate


def add_terminal_caps_to_pdb(
    input_path: str | Path,
    *,
    chain_ids: list[str] | None = None,
    verbose: bool = False,
) -> Path:
    """Return a temporary PDB with ACE/NME added to selected protein chains.

    ``chain_ids=None`` caps every protein chain.  ``_`` selects a blank chain.
    Existing ACE/NME caps are retained, making repeated application idempotent.
    """
    source = Path(input_path)
    lines = source.read_text().splitlines(keepends=True)
    residues = _parse_residues(lines)
    requested = None if chain_ids is None else {
        " " if chain == "_" else chain for chain in chain_ids
    }

    # Protein segments are contiguous residue runs separated by TER records.
    segments: list[list[_Residue]] = []
    current: list[_Residue] = []
    previous_end = -2
    for residue in residues:
        if current and (residue.chain != current[-1].chain
                        or residue.line_indices[0] != previous_end + 1):
            segments.append(current)
            current = []
        current.append(residue)
        previous_end = residue.line_indices[-1]
    if current:
        segments.append(current)

    candidates = [segment for segment in segments
                  if any(r.resname in PROTEIN_RESNAMES for r in segment)]
    available = {r.chain for segment in candidates for r in segment
                 if r.resname in PROTEIN_RESNAMES}
    if requested is not None:
        missing = requested - available
        if missing:
            labels = ", ".join("_" if c == " " else c for c in sorted(missing))
            print(f"WARNING: --cap-chain selection not found among protein chains: "
                  f"{labels}; continuing without those chain(s)")
            requested -= missing

    max_serial = max((serial for r in residues for serial, _ in r.atoms.values()),
                     default=0)
    next_serial = max_serial + 1
    used_by_chain: dict[str, set[int]] = {}
    for residue in residues:
        used_by_chain.setdefault(residue.chain, set()).add(residue.resid)

    insert_before: dict[int, list[str]] = {}
    insert_after: dict[int, list[str]] = {}
    remove_indices: set[int] = set()
    new_conects: list[str] = []
    added = 0

    for segment in candidates:
        protein = [r for r in segment if r.resname in PROTEIN_RESNAMES]
        if not protein:
            continue
        chain = protein[0].chain
        if requested is not None and chain not in requested:
            continue
        first, last = protein[0], protein[-1]
        first_pos = segment.index(first)
        last_pos = segment.index(last)
        prev_res = segment[first_pos - 1] if first_pos else None
        next_res = segment[last_pos + 1] if last_pos + 1 < len(segment) else None

        if prev_res is not None and prev_res.resname in CAP_RESNAMES:
            if prev_res.resname != "ACE":
                raise ValueError(
                    f"chain {chain or '_'} already has unsupported N-cap {prev_res.resname}"
                )
        else:
            missing_atoms = {"N", "CA", "C"} - set(first.atoms)
            if missing_atoms:
                raise ValueError(
                    f"cannot cap chain {chain or '_'}: first residue {first.resname}{first.resid} "
                    f"is missing backbone anchor(s) {', '.join(sorted(missing_atoms))}"
                )
            target = np.array([first.atoms[name][1] for name in ("N", "CA", "C")])
            coords = _transform(_N_ANCHORS, target, _ACE_HEAVY)
            cap_resid = _allocate_resid(used_by_chain[chain], first.resid, -1)
            serials: dict[str, int] = {}
            cap_lines = []
            for name, element in (("CH3", "C"), ("C", "C"), ("O", "O")):
                serials[name] = next_serial
                cap_lines.append(_atom_line(next_serial, name, "ACE", chain,
                                            cap_resid, coords[name], element))
                next_serial += 1
            insert_before.setdefault(first.line_indices[0], []).extend(cap_lines)
            n_serial = first.atoms["N"][0]
            new_conects.extend([
                f"CONECT{serials['CH3']:5d}{serials['C']:5d}\n",
                f"CONECT{serials['C']:5d}{serials['CH3']:5d}{serials['O']:5d}{n_serial:5d}\n",
                f"CONECT{serials['O']:5d}{serials['C']:5d}\n",
                f"CONECT{n_serial:5d}{serials['C']:5d}\n",
            ])
            added += 1

        if next_res is not None and next_res.resname in CAP_RESNAMES:
            if next_res.resname != "NME":
                raise ValueError(
                    f"chain {chain or '_'} already has unsupported C-cap {next_res.resname}"
                )
        else:
            missing_atoms = {"N", "CA", "C"} - set(last.atoms)
            if missing_atoms:
                raise ValueError(
                    f"cannot cap chain {chain or '_'}: last residue {last.resname}{last.resid} "
                    f"is missing backbone anchor(s) {', '.join(sorted(missing_atoms))}"
                )
            target = np.array([last.atoms[name][1] for name in ("N", "CA", "C")])
            coords = _transform(_C_ANCHORS, target, _NME_HEAVY)
            cap_resid = _allocate_resid(used_by_chain[chain], last.resid, 1)
            serials = {}
            cap_lines = []
            for name, element in (("N", "N"), ("C", "C")):
                serials[name] = next_serial
                cap_lines.append(_atom_line(next_serial, name, "NME", chain,
                                            cap_resid, coords[name], element))
                next_serial += 1
            insert_after.setdefault(last.line_indices[-1], []).extend(cap_lines)
            c_serial = last.atoms["C"][0]
            new_conects.extend([
                f"CONECT{c_serial:5d}{serials['N']:5d}\n",
                f"CONECT{serials['N']:5d}{c_serial:5d}{serials['C']:5d}\n",
                f"CONECT{serials['C']:5d}{serials['N']:5d}\n",
            ])
            if "OXT" in last.atoms:
                remove_indices.update(
                    line_idx for line_idx in last.line_indices
                    if lines[line_idx][12:16].strip() == "OXT"
                )
                # Existing CONECT records referencing OXT are filtered below.
            added += 1

    if not added:
        return source

    out: list[str] = []
    for index, line in enumerate(lines):
        if line.startswith("CONECT"):
            # Regenerate only cap connectivity; retain existing records unless
            # they reference an atom removed above (currently terminal OXT).
            try:
                refs = [int(line[i:i + 5]) for i in range(6, len(line), 5)
                        if line[i:i + 5].strip()]
            except ValueError:
                refs = []
            removed_serials = {
                int(lines[i][6:11]) for i in remove_indices
                if lines[i].startswith(("ATOM  ", "HETATM"))
            }
            if any(ref in removed_serials for ref in refs):
                continue
        if index in insert_before:
            out.extend(insert_before[index])
        if index not in remove_indices:
            out.append(line)
        if index in insert_after:
            out.extend(insert_after[index])

    end_index = next((i for i, line in enumerate(out) if line.startswith("END")), len(out))
    out[end_index:end_index] = new_conects
    output = make_temp_path(suffix=".pdb", prefix="dvbfixer_caps_")
    output.write_text("".join(out))
    if verbose:
        print(f"Added {added} terminal cap residue(s) (ACE/NME) → {output}")
    return output
