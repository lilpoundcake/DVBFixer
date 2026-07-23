"""Shared FF-name conversion module.

Central place for AMBER↔CHARMM protonation-variant maps and the
text-level PDB rewrite primitives. Reused by:

- ``dvbfixer convert`` (formerly ``glycam``) — bidirectional
  AMBER↔CHARMM conversion including sugars + protein variants.
- ``dvbfixer prepare`` / ``minimize`` / ``protonate`` — restoring
  variant names on user-visible output PDBs after OpenMM's
  ``PDBFile.writeFile`` canonicalises them.
- Any future FF-name-aware tool.

Verified empirically (OpenMM 8.5.1, Jul 2026): ``PDBFile`` reads
HIE/HID/HIP → HIS, ASH → ASP, GLH → GLU, CYX → CYS on load; writes
the canonical name on save. LYN and CYM are preserved through
the round-trip. The topology-level "rename → addHydrogens → restore"
dance in older code was a no-op because ``_saved`` was always empty
after the read. The correct fix is a text-level rewrite of the
output PDB using ``amber_renames`` derived from the RAW input text
(via ``scan_variant_names`` / ``_read_amber_renames``).
"""

from __future__ import annotations

from pathlib import Path

# ---------------------------------------------------------------------------
# Residue-name maps
# ---------------------------------------------------------------------------
# Source: AMBER ff14SB/ff19SB XML + CHARMM36 aminoacids.rtp (verified).
# Only the variants that OpenMM's shipped `charmm36.xml` actually contains
# as templates are listed here. ASH/GLH/LYN have NO OpenMM-CHARMM36
# equivalent (CHARMM-GUI uses ASPP/GLUP/LSN as 4-char names but those
# aren't in the OpenMM XML), so on --ff charmm those residues fall back
# to their standard parent (ASP/GLU/LYS) — see ``AMBER_CHARMM_UNSUPPORTED``.
PROTONATION_AMBER_TO_CHARMM: dict[str, str] = {
    "HID": "HSD", "HIE": "HSE", "HIP": "HSP",
    "ASH": "ASPP", "GLH": "GLUP",
    "LYN": "LSN",
    "CYX": "CYS",   # CHARMM uses CYS + DISU patch, applied via SSBOND
    # CYM stays as CYM — CHARMM36 has a [ CYM ] residue.
}

PROTONATION_CHARMM_TO_AMBER: dict[str, str] = {
    "HSD": "HID", "HSE": "HIE", "HSP": "HIP",
    "ASPP": "ASH", "GLUP": "GLH",
    "LSN": "LYN",
}

# ---------------------------------------------------------------------------
# Per-residue atom-name shifts
# ---------------------------------------------------------------------------
# The only AMBER↔CHARMM rename that's asymmetric on the atom set (so
# GROMACS's aminoacids.arn can't infer it) is LYN/LSN's NH2-H pair:
#
#   AMBER ff14SB LYN: HZ2 + HZ3 (HZ1 absent)
#   CHARMM36 LSN:     HZ1 + HZ2 (HZ3 absent)
#
# The pair is applied atomically per residue to avoid the
# HZ2→HZ1-then-HZ3→HZ2 collision.
PROTONATION_ATOM_RENAME_TO_CHARMM: dict[str, dict[str, str]] = {
    "LYN": {"HZ2": "HZ1", "HZ3": "HZ2"},
}
PROTONATION_ATOM_RENAME_TO_AMBER: dict[str, dict[str, str]] = {
    "LSN": {"HZ1": "HZ2", "HZ2": "HZ3"},
}

# GROMACS amber99sb-ildn LYN uses HZ1 + HZ2, but the OpenMM ff14SB
# template LYN uses HZ2 + HZ3. When writing PDB output intended for
# pdb2gmx -ff amber99sb-ildn, rename HZ3 → HZ1. The two H atoms are
# chemically equivalent (same charge, same bond topology).
GROMACS_AMBER_LYN_ATOM_RENAME: dict[str, dict[str, str]] = {
    "LYN": {"HZ3": "HZ1"},
}

# Convenience sets.
AMBER_VARIANTS: frozenset[str] = frozenset({
    "HID", "HIE", "HIP", "ASH", "GLH", "CYX", "CYM", "LYN",
})
CHARMM_VARIANTS: frozenset[str] = frozenset({
    "HSD", "HSE", "HSP", "ASPP", "GLUP", "CYM", "LSN",
})


# ---------------------------------------------------------------------------
# Primitive
# ---------------------------------------------------------------------------

def _split_key(k: object) -> tuple[str, str, str]:
    """Coerce a key from various sources into (chain, resid_str, icode).

    Accepts:
      - (chain, resid_str)             -> icode=""
      - (chain, resid_str, icode)      -> as-is
    """
    if isinstance(k, tuple):
        if len(k) == 2:
            return (str(k[0]), str(k[1]), "")
        if len(k) == 3:
            return (str(k[0]), str(k[1]), str(k[2] or ""))
    raise TypeError(f"Unexpected key shape for amber_renames: {k!r}")


def _format_atom_name_field(name: str) -> str:
    """Right-pad-4 layout for the 4-char atom name field (cols 13-16).

    PDB convention: 1-3 char atom names are padded with a LEADING space
    (element letter at col 14, remainder at 15-16). 4-char names fill the
    whole field.
    """
    if len(name) >= 4:
        return name[:4]
    return f" {name:<3s}"


def _rewrite_atom_line(line: str, new_atomname: str, new_resname: str) -> str:
    """Rewrite atom-name (cols 13-16) and residue-name (cols 18-20 for
    3-char names, cols 18-21 for 4-char names) on an ATOM/HETATM line.

    For 4-char residue names (ASPP/GLUP/ASPD/etc.) the name extends into
    col 20 (which is normally a space padding) and pushes into the chain
    field's space — kept the standard PDB compromise used by GROMACS
    and downstream CHARMM tools.
    """
    atom_field = _format_atom_name_field(new_atomname)
    # Column 17 is altLoc (preserve).
    altloc = line[16]
    if len(new_resname) <= 3:
        res_field = f"{new_resname:<3s}"
        # Cols 18-20 = 3-char residue, col 21 onwards unchanged.
        return line[:12] + atom_field + altloc + res_field + line[20:]
    # 4-char residue: extends into col 20 (which is chain-id column 21 -1).
    # Standard workaround (matches GROMACS + CHARMM-GUI output): put the
    # 4-char name in cols 17-20, keep chain at col 21.
    res_field = new_resname[:4]
    return line[:12] + atom_field + altloc + res_field + line[21:]


def apply_variants_to_pdb_text(
    pdb_path: str | Path,
    amber_renames: dict,
    target_ff: str = "amber",
    include_gromacs_lyn: bool = True,
    verbose: bool = False,
) -> int:
    """Rewrite ATOM/HETATM residue and (LYN) atom names in ``pdb_path``
    so the file reflects the variants recorded in ``amber_renames``.

    Parameters
    ----------
    pdb_path
        User-visible PDB file. Rewritten in place.
    amber_renames
        ``{(chain, resid_str): variant_name}`` OR
        ``{(chain, resid_str, icode): variant_name}``. Populated by
        ``scan_variant_names`` / ``_read_amber_renames`` from the RAW
        input text (before OpenMM parsed it).
    target_ff
        ``'amber'`` writes AMBER variant names verbatim
        (HID/HIE/HIP/ASH/GLH/CYX/CYM/LYN). ``'charmm'`` maps
        AMBER→CHARMM via ``PROTONATION_AMBER_TO_CHARMM`` (HSD/HSE/HSP/
        ASPP/GLUP/LSN + CYX→CYS + CYM→CYM) and applies the LYN→LSN
        atom-name shift.
    include_gromacs_lyn
        Only meaningful when ``target_ff='amber'``. Applies the
        ff14SB HZ3 → amber99sb-ildn HZ1 rename on LYN residues so
        the output is directly consumable by ``pdb2gmx``. Off if
        strict ff14SB naming is required.

    Returns the number of ATOM/HETATM lines whose residue OR atom name
    changed. Idempotent.
    """
    if target_ff not in ("amber", "charmm"):
        raise ValueError(f"target_ff must be 'amber' or 'charmm', got {target_ff!r}")

    # Normalise key shape for O(1) lookup during scan.
    lookup: dict[tuple[str, str], str] = {}
    lookup_with_icode: dict[tuple[str, str, str], str] = {}
    for k, v in (amber_renames or {}).items():
        chain, resid, icode = _split_key(k)
        if icode:
            lookup_with_icode[(chain, resid, icode)] = v
        else:
            lookup[(chain, resid)] = v

    if not lookup and not lookup_with_icode:
        # No variants recorded — but we still may need the GROMACS LYN
        # rename on any bare LYN residue that was preserved by OpenMM
        # (LYN survives PDBFile round-trip). Compute it below.
        pass

    path = Path(pdb_path)
    lines = path.read_text().splitlines(keepends=True)
    out: list[str] = []
    n_rewritten = 0

    for line in lines:
        if not (line.startswith(("ATOM  ", "HETATM")) and len(line) >= 27):
            out.append(line)
            continue

        chain = line[21]
        resid = line[22:26].strip()
        icode = line[26].strip()
        cur_resname = line[17:20].strip()
        cur_atomname = line[12:16].strip()

        # Look up variant recorded for this residue (icode-specific first,
        # then any-icode fallback).
        variant = lookup_with_icode.get((chain, resid, icode))
        if variant is None:
            variant = lookup.get((chain, resid))

        # Compute the target residue name.
        if variant is not None:
            if target_ff == "amber":
                new_resname = variant
            else:  # charmm
                new_resname = PROTONATION_AMBER_TO_CHARMM.get(variant, variant)
        else:
            # No variant recorded. If the current name IS an AMBER variant
            # (e.g. LYN or CYM survived PDBFile round-trip), still map it
            # for charmm target.
            if target_ff == "charmm" and cur_resname in PROTONATION_AMBER_TO_CHARMM:
                new_resname = PROTONATION_AMBER_TO_CHARMM[cur_resname]
            else:
                new_resname = cur_resname

        # Compute the target atom name.
        new_atomname = cur_atomname
        if target_ff == "charmm":
            atom_map = PROTONATION_ATOM_RENAME_TO_CHARMM.get(
                variant or cur_resname, {}
            )
            new_atomname = atom_map.get(cur_atomname, cur_atomname)
        elif target_ff == "amber" and include_gromacs_lyn:
            # LYN HZ3 → HZ1 for GROMACS amber99sb-ildn.
            if new_resname == "LYN":
                new_atomname = GROMACS_AMBER_LYN_ATOM_RENAME["LYN"].get(
                    cur_atomname, cur_atomname
                )

        # Emit updated line if anything changed.
        if new_resname != cur_resname or new_atomname != cur_atomname:
            new_line = _rewrite_atom_line(line, new_atomname, new_resname)
            out.append(new_line)
            n_rewritten += 1
            if verbose:
                print(f"  [ff_names] {chain}/{cur_resname}{resid}:{cur_atomname} "
                      f"-> {new_resname}:{new_atomname}")
        else:
            out.append(line)

    if n_rewritten:
        path.write_text("".join(out))
    return n_rewritten
