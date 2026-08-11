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

# Per-residue atom renames: OpenMM ff14SB (IUPAC) → GROMACS
# amber99sb-ildn. Only the "3" atom of each methylene needs a
# rename — HB3 → HB1 leaves HB2 alone, and the result {HB1, HB2}
# is exactly the GROMACS canonical (the two methylene H atoms are
# chemically equivalent, so pdb2gmx accepts either ordering).
# Naturally idempotent: once HB3 is gone, second pass finds nothing.
#
# Sub-agent survey (Jul 2026, verified against OpenMM XML +
# `FF/amber99sb-ildn-lipid21.ff/aminoacids.rtp`) confirmed this is
# the complete set of intra-residue differences the GROMACS
# `aminoacids.arn` file does NOT rewrite at pdb2gmx read time.
GROMACS_AMBER_ATOM_RENAMES: dict[str, dict[str, str]] = {
    # No methylene / α-branched.
    "ALA": {}, "VAL": {}, "THR": {},
    # β-methylene only.
    "SER": {"HB3": "HB1"},
    "CYS": {"HB3": "HB1"},
    "CYX": {"HB3": "HB1"},
    "CYM": {"HB3": "HB1"},
    "ASN": {"HB3": "HB1"},
    "ASP": {"HB3": "HB1"},
    "ASH": {"HB3": "HB1"},
    "PHE": {"HB3": "HB1"},
    "TYR": {"HB3": "HB1"},
    "TRP": {"HB3": "HB1"},
    "HIS": {"HB3": "HB1"},
    "HID": {"HB3": "HB1"},
    "HIE": {"HB3": "HB1"},
    "HIP": {"HB3": "HB1"},
    "LEU": {"HB3": "HB1"},
    # β + γ methylene.
    "GLU": {"HB3": "HB1", "HG3": "HG1"},
    "GLH": {"HB3": "HB1", "HG3": "HG1"},
    "GLN": {"HB3": "HB1", "HG3": "HG1"},
    "MET": {"HB3": "HB1", "HG3": "HG1"},
    # β + γ + δ methylene.
    "PRO": {"HB3": "HB1", "HG3": "HG1", "HD3": "HD1"},
    "HYP": {"HB3": "HB1", "HG3": "HG1", "HD3": "HD1"},
    "ARG": {"HB3": "HB1", "HG3": "HG1", "HD3": "HD1"},
    # β + γ + δ + ε methylene.
    "LYS": {"HB3": "HB1", "HG3": "HG1", "HD3": "HD1", "HE3": "HE1"},
    "LYN": {"HB3": "HB1", "HG3": "HG1", "HD3": "HD1", "HE3": "HE1",
            # LYN NZ: ff14SB HZ2+HZ3 → GROMACS HZ1+HZ2 (HZ3 → HZ1,
            # HZ2 stays; same trick as the methylenes).
            "HZ3": "HZ1"},
    # Special cases.
    "GLY": {"HA3": "HA1"},           # α-methylene on GLY
    "ILE": {"HG13": "HG11"},         # γ-methylene only
    # Cap residues. OpenMM writes PDB-style H1/H2/H3 (and C for the NME
    # methyl carbon), while the GROMACS AMBER RTP uses HH31/HH32/HH33 and CH3.
    "ACE": {"H1": "HH31", "H2": "HH32", "H3": "HH33"},
    "NME": {"H1": "HH31", "H2": "HH32", "H3": "HH33", "C": "CH3"},
}

# Kept for back-compat with 0.6.x callers that referenced it directly.
GROMACS_AMBER_LYN_ATOM_RENAME: dict[str, dict[str, str]] = {
    "LYN": {"HZ3": "HZ1"},
}

# Backbone amide H rename applied on every protein residue when the
# target FF is CHARMM. GROMACS charmm36 RTP uses HN for the amide
# hydrogen; ff14SB uses H.
CHARMM_UNIVERSAL_ATOM_RENAMES: dict[str, str] = {
    "H": "HN",
}

# OpenMM's cap templates use portable PDB atom names. GROMACS CHARMM36 RTP
# entries use the CHARMM methyl and amide-hydrogen names instead.
GROMACS_CHARMM_CAP_ATOM_RENAMES: dict[str, dict[str, str]] = {
    "ACE": {"H1": "HH31", "H2": "HH32", "H3": "HH33"},
    "NME": {
        "H1": "HH31", "H2": "HH32", "H3": "HH33",
        "C": "CH3", "H": "HN",
    },
}

# N-terminal ATOM rename: applied on the first residue of each
# protein chain when target is GROMACS amber (only if the raw H is
# still present — skip if Modeller already placed H1/H2/H3 for the
# charged NH3+ terminus).
_NTER_H_RENAME = {"H": "H1"}

# N-terminal proline (NPRO) backbone rename: OpenMM has H2/H3;
# GROMACS uses H1/H2. Single-source rename H3 → H1 leaves H2 alone,
# result {H1, H2} matches GROMACS (the two ring-N H atoms are
# chemically equivalent). Naturally idempotent.
_NPRO_BACKBONE_H_RENAME = {"H3": "H1"}

# C-terminal ATOM rename: applied on the last residue of each
# protein chain when target is GROMACS amber. Atomic per residue
# to avoid O→OC2-then-OXT→OC2 collision.
_CTER_O_RENAME = {"O": "OC2", "OXT": "OC1"}

# Nucleic acid renames (DNA + RNA) — verified against OpenMM
# amber14/DNA.OL15.xml, amber14/RNA.OL3.xml vs GROMACS dna.rtp /
# rna.rtp. Covers the confirmed methylene shifts on the sugar
# backbone. Rare cases (terminal 5'/3' OH atoms) not yet verified.
_NA_METHYLENE = {"H2'": "H2'1", "H2''": "H2'2",
                 "H5'": "H5'1", "H5''": "H5'2"}
GROMACS_AMBER_NA_ATOM_RENAMES: dict[str, dict[str, str]] = {
    "DA": dict(_NA_METHYLENE), "DC": dict(_NA_METHYLENE),
    "DG": dict(_NA_METHYLENE), "DT": dict(_NA_METHYLENE),
    "A": dict(_NA_METHYLENE), "C": dict(_NA_METHYLENE),
    "G": dict(_NA_METHYLENE), "U": dict(_NA_METHYLENE),
    # RNA-only: OpenMM HO'2 → GROMACS HO2' (element letter position).
    "A_RNA_EXTRA": {"HO'2": "HO2'"},   # merged below into A/C/G/U
    "C_RNA_EXTRA": {"HO'2": "HO2'"},
    "G_RNA_EXTRA": {"HO'2": "HO2'"},
    "U_RNA_EXTRA": {"HO'2": "HO2'"},
}
# Fold the RNA-extras into A/C/G/U so callers see one map per residue.
for _r in ("A", "C", "G", "U"):
    GROMACS_AMBER_NA_ATOM_RENAMES[_r].update(
        GROMACS_AMBER_NA_ATOM_RENAMES.pop(f"{_r}_RNA_EXTRA")
    )

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

# All residue names dvbfixer considers "protein" for the purposes of
# N/C-terminal detection. Includes standard AAs + AMBER variants + a
# few common non-standard names.
_PROTEIN_RESNAMES: frozenset[str] = frozenset({
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "HYP",
    # AMBER variants
    "HID", "HIE", "HIP", "ASH", "GLH", "CYX", "CYM", "LYN",
    # CHARMM variants
    "HSD", "HSE", "HSP", "ASPP", "GLUP", "LSN",
    # Non-standard
    "MSE", "SEC", "PYL", "SEP", "TPO", "PTR", "CSO", "CSD", "CME",
})
_TERMINAL_CHAIN_RESNAMES = _PROTEIN_RESNAMES | {"ACE", "NME", "NHE", "NH2"}


def _effective_amber_rename_map(
    resname: str,
    atom_set: set[str],
    is_nter: bool = False,
    is_cter: bool = False,
) -> dict[str, str]:
    """Compute the atom-name rename map to apply to a single residue
    when writing GROMACS amber99sb-ildn compatible output.

    All renames are single-source (HB3 → HB1, not the pair HB2 → HB1
    + HB3 → HB2), so repeat application is naturally idempotent: once
    the source atom is gone, the rename can't fire again.
    """
    m: dict[str, str] = {}

    # Standard per-residue renames from the lookup table (methylenes,
    # LYN NZ, ILE γ, GLY α, ACE/NME cap H).
    m.update(GROMACS_AMBER_ATOM_RENAMES.get(resname, {}))

    # Terminal renames.
    if is_nter and resname in _PROTEIN_RESNAMES:
        if resname in ("PRO", "HYP"):
            # NPRO backbone ring-N H3 → H1 (H2 stays). Only fires when
            # H3 is present.
            if "H3" in atom_set:
                m.update(_NPRO_BACKBONE_H_RENAME)
        else:
            # H → H1 only if H exists and H1 doesn't (Modeller may
            # already have placed H1/H2/H3 for a charged terminus).
            if "H" in atom_set and "H1" not in atom_set:
                m.update(_NTER_H_RENAME)
    if is_cter and resname in _PROTEIN_RESNAMES:
        # Only apply if both O and OXT are present (unshifted state).
        # Once renamed to OC1/OC2, second pass finds neither and skips.
        if "O" in atom_set and "OXT" in atom_set:
            m.update(_CTER_O_RENAME)

    return m


def _classify_terminal_residues(
    lines: list[str],
) -> dict[tuple[str, str, str], set[str]]:
    """Return {(chain, resid_str, icode): {'nter', 'cter'}} for the
    first / last protein residue of each chain in the input PDB text.

    Single-residue chains get BOTH 'nter' and 'cter' in the set —
    the terminal renames for each end apply independently.

    "Protein" = residue name in :data:`_PROTEIN_RESNAMES`. HETATMs and
    solvent are ignored so waters near a chain terminus don't get
    misclassified.
    """
    per_chain: dict[str, list[tuple[str, str]]] = {}
    seen: dict[str, set[tuple[str, str]]] = {}
    for line in lines:
        if not (line.startswith(("ATOM  ", "HETATM")) and len(line) >= 27):
            continue
        resname = line[17:20].strip()
        if resname not in _TERMINAL_CHAIN_RESNAMES:
            continue
        chain = line[21]
        resid = line[22:26].strip()
        icode = line[26].strip()
        key = (resid, icode)
        seen_set = seen.setdefault(chain, set())
        if key in seen_set:
            continue
        seen_set.add(key)
        per_chain.setdefault(chain, []).append(key)

    out: dict[tuple[str, str, str], set[str]] = {}
    for chain, residues in per_chain.items():
        if not residues:
            continue
        first_resid, first_icode = residues[0]
        last_resid, last_icode = residues[-1]
        out.setdefault((chain, first_resid, first_icode), set()).add("nter")
        out.setdefault((chain, last_resid, last_icode), set()).add("cter")
    return out


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
    include_gromacs_shifts: bool = True,
    verbose: bool = False,
    *,
    include_gromacs_lyn: bool | None = None,  # deprecated alias
) -> int:
    """Rewrite ATOM/HETATM residue and atom names in ``pdb_path`` so the
    file is directly GROMACS-canonical for ``target_ff``.

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
        (HID/HIE/HIP/ASH/GLH/CYX/CYM/LYN) plus applies GROMACS
        amber99sb-ildn atom-name shifts (methylene HB2/HB3 → HB1/HB2,
        GLY HA shift, ILE γ shift, LYN NZ shift, N-terminal H → H1,
        C-terminal O/OXT → OC1/OC2, ACE/NME cap H). ``'charmm'`` maps
        AMBER→CHARMM residue names + LYN→LSN atom shift +
        universal backbone H → HN.
    include_gromacs_shifts
        Enable all GROMACS-compatibility rewrites (per-residue methylene
        shifts, terminal renames, cap renames, backbone H→HN for CHARMM).
        Off if strict ff14SB naming is required.
    include_gromacs_lyn
        Deprecated alias for ``include_gromacs_shifts``. Emits no
        warning yet; will be removed in a future release.

    Returns the number of ATOM/HETATM lines whose residue OR atom name
    changed. Idempotent.
    """
    if target_ff not in ("amber", "charmm"):
        raise ValueError(f"target_ff must be 'amber' or 'charmm', got {target_ff!r}")

    if include_gromacs_lyn is not None:
        include_gromacs_shifts = include_gromacs_lyn

    # Normalise key shape for O(1) lookup during scan.
    lookup: dict[tuple[str, str], str] = {}
    lookup_with_icode: dict[tuple[str, str, str], str] = {}
    for k, v in (amber_renames or {}).items():
        chain, resid, icode = _split_key(k)
        if icode:
            lookup_with_icode[(chain, resid, icode)] = v
        else:
            lookup[(chain, resid)] = v

    path = Path(pdb_path)
    lines = path.read_text().splitlines(keepends=True)

    # Terminal-residue detection (needed for N/C-term renames).
    terminal_map = _classify_terminal_residues(lines) if include_gromacs_shifts else {}

    # Per-residue atom-name sets, used to decide whether a shift-style
    # rename (HB2/HB3 → HB1/HB2, etc.) should trigger. Shifts collide
    # on themselves under repeated application; gating on presence of
    # the "highest-numbered source" atom (HB3, HZ3, HA3, ...) makes
    # the pass idempotent.
    residue_atoms: dict[tuple[str, str, str], set[str]] = {}
    for line in lines:
        if not (line.startswith(("ATOM  ", "HETATM")) and len(line) >= 27):
            continue
        key = (line[21], line[22:26].strip(), line[26].strip())
        residue_atoms.setdefault(key, set()).add(line[12:16].strip())

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
        terminal_pos = terminal_map.get((chain, resid, icode), set())

        if target_ff == "charmm":
            # LYN → LSN NZ H shift, cap names, plus backbone H → HN.
            atom_map = dict(PROTONATION_ATOM_RENAME_TO_CHARMM.get(
                variant or cur_resname, {}
            ))
            if include_gromacs_shifts:
                atom_map.update(GROMACS_CHARMM_CAP_ATOM_RENAMES.get(
                    new_resname, {}
                ))
            if include_gromacs_shifts and cur_resname in _PROTEIN_RESNAMES:
                atom_map.setdefault(
                    "H", CHARMM_UNIVERSAL_ATOM_RENAMES["H"]
                )
            new_atomname = atom_map.get(cur_atomname, cur_atomname)
        elif target_ff == "amber" and include_gromacs_shifts:
            # Build the effective per-residue rename map, gated on the
            # actual atom set so shifts don't collide on repeat passes.
            atoms_here = residue_atoms.get((chain, resid, icode), set())
            combined = _effective_amber_rename_map(
                new_resname, atoms_here,
                is_nter=("nter" in terminal_pos),
                is_cter=("cter" in terminal_pos),
            )
            # Nucleic-acid renames (H2'/H2'' etc.) — source atoms are
            # unique so no shift-collision; safe to apply unconditionally.
            if cur_resname in GROMACS_AMBER_NA_ATOM_RENAMES:
                for src, tgt in GROMACS_AMBER_NA_ATOM_RENAMES[cur_resname].items():
                    combined.setdefault(src, tgt)
            new_atomname = combined.get(cur_atomname, cur_atomname)

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
