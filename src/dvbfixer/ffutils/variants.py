"""Shared handling of protonation-variant residue names.

Every OpenMM-using tool has to deal with the same core problem: AMBER and
CHARMM protonation variants (LYN, HIE/HSD, ASH/ASPP, CYX, ...) are NOT
recognised by OpenMM's ``PDBFile`` parser as amino acids, so when they
appear mid-chain they silently block peptide-bond inference from the
previous residue's C to their own N. Downstream ``addHydrogens`` then
fails with a confusing "missing external C atom" error on the residue
BEFORE the variant.

The fix is three-step:

1. Before OpenMM parses the file, rewrite variant residue names to their
   standard parents (LYS/HIS/ASP/GLU/CYS) in the raw PDB text and remember
   what was renamed. This is what :func:`text_rename_variants_to_parent`
   does — text-level, no OpenMM in the loop.
2. Call ``Modeller.addHydrogens(..., variants=[...])`` — the variants list
   is what tells OpenMM to place HZ3 for LYS-only, skip HE2 for HID, etc.
   Build it from the saved map via :func:`build_variants_list`.
3. After ``addHydrogens`` returns a new topology, restore the variant
   residue names on the new topology via
   :func:`restore_variants_post_addhydrogens`. The old topology's residue
   references are stale after the call.

Previously each of ``prepare``, ``minimize``, ``protonate``, and
``transplant`` re-implemented this dance with subtly different variant
tables and dict shapes. This module is the one source of truth.
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Any

# AMBER protonation variants — these appear in PDBs produced by dvbfixer's
# own ``protonate`` step, by AmberTools, or by users who annotated a
# residue with a non-default protonation state.
AMBER_VARIANT_TO_PARENT: dict[str, str] = {
    "HID": "HIS",
    "HIE": "HIS",
    "HIP": "HIS",
    "ASH": "ASP",
    "GLH": "GLU",
    "CYX": "CYS",
    "CYM": "CYS",
    "LYN": "LYS",
}

# CHARMM protonation variants — appear in CHARMM-GUI output.
CHARMM_VARIANT_TO_PARENT: dict[str, str] = {
    "HSD": "HIS",
    "HSE": "HIS",
    "HSP": "HIS",
    "ASPP": "ASP",
    "GLUP": "GLU",
    "LSN": "LYS",
}

# Combined AMBER + CHARMM lookup used at all text-level rename sites.
VARIANT_TO_PARENT: dict[str, str] = {
    **AMBER_VARIANT_TO_PARENT,
    **CHARMM_VARIANT_TO_PARENT,
}

# Just the names (not the mapping) — used by ``minimize`` when scanning a
# PDB for residues that need addHydrogens variants but the caller doesn't
# care about the parent.
AMBER_VARIANTS: frozenset[str] = frozenset(AMBER_VARIANT_TO_PARENT.keys())
CHARMM_VARIANTS: frozenset[str] = frozenset(CHARMM_VARIANT_TO_PARENT.keys())
ALL_VARIANTS: frozenset[str] = AMBER_VARIANTS | CHARMM_VARIANTS


ReskeyText = tuple[str, str, str]  # (chain, resseq, icode) as strings, whitespace-stripped
SavedMap = dict[ReskeyText, str]


def _reskey_from_atom_line(line: str) -> ReskeyText:
    chain = line[21] if len(line) > 21 else " "
    resseq_field = line[22:26]
    try:
        resseq = str(int(resseq_field))
    except ValueError:
        resseq = resseq_field.strip() or "0"
    icode = line[26].strip() if len(line) > 26 else ""
    return (chain, resseq, icode)


def scan_variant_names(pdb_path: str | Path) -> SavedMap:
    """Scan a PDB file for AMBER/CHARMM protonation variant residue names.

    Returns ``{(chain, resseq, icode): original_name}``. Used by callers
    (e.g. ``minimize``) that need the map to build an addHydrogens variants
    list WITHOUT rewriting the file — OpenMM's ``PDBFile`` reader will
    normalise the names on load, so this must be called on the raw PDB
    text first.
    """
    saved: SavedMap = {}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            resname = line[17:20].strip()
            if resname in ALL_VARIANTS:
                saved[_reskey_from_atom_line(line)] = resname
    return saved


def text_rename_variants_to_parent(
    pdb_path: str | Path,
    *,
    verbose: bool = False,
) -> tuple[str, SavedMap]:
    """Rewrite AMBER/CHARMM variant residue names to their standard parents.

    Operates on the raw PDB text. Returns ``(out_path, saved_map)`` where
    ``saved_map`` is ``{(chain, resseq, icode): original_variant}``. If no
    lines needed rewriting, returns ``(str(pdb_path), {})`` — the caller
    can use the original file directly.

    Needed so OpenMM's ``PDBFile`` parser can infer peptide bonds through
    mid-chain variants without failing on the adjacent residue.
    """
    pdb_path = str(pdb_path)
    with open(pdb_path) as f:
        lines = f.readlines()

    n_rewrite = 0
    saved: SavedMap = {}
    out_lines: list[str] = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 20:
            resname = line[17:20].strip()
            parent = VARIANT_TO_PARENT.get(resname)
            if parent is not None and parent != resname:
                saved[_reskey_from_atom_line(line)] = resname
                # 3-char left-aligned — preserves column alignment even when
                # the source name was 4 chars (ASPP/GLUP/LSN).
                line = line[:17] + f"{parent:<3s}" + line[20:]
                n_rewrite += 1
        out_lines.append(line)

    if n_rewrite == 0:
        return pdb_path, {}

    with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as tmp:
        tmp.writelines(out_lines)
        out_path = tmp.name
    if verbose:
        print(
            f"  [variants] pre-renamed {n_rewrite} variant line(s) to standard "
            f"parents so OpenMM can infer peptide bonds; original names "
            f"restored after addHydrogens"
        )
    return out_path, saved


def build_variants_list(
    topology: Any,
    saved: SavedMap,
) -> list[str | None] | None:
    """Build the ``variants=`` list for ``Modeller.addHydrogens``.

    Returns a list aligned with ``topology.residues()``: variant name for
    residues in ``saved``, ``None`` otherwise. If no residue in the
    topology matches an entry in ``saved``, returns ``None`` — OpenMM
    accepts ``variants=None`` as "use defaults".
    """
    if not saved:
        return None
    variants: list[str | None] = []
    hit = False
    for res in topology.residues():
        key = (res.chain.id, str(res.id), "")
        # ``res.id`` is a plain integer-as-string in OpenMM; iCode is not
        # preserved in ``res.id``. Match on (chain, resseq) with any icode
        # by trying the empty-icode key first, then falling back to any
        # icode variant.
        if key in saved:
            variants.append(saved[key])
            hit = True
            continue
        # Fallback: same chain + resseq, non-empty icode.
        fallback = next(
            (name for (c, r, i), name in saved.items() if c == res.chain.id and r == str(res.id)),
            None,
        )
        variants.append(fallback)
        if fallback is not None:
            hit = True
    return variants if hit else None


def rename_variants_to_parent_in_topology(topology: Any) -> dict[tuple[str, str], str]:
    """Rename variant residues in an OpenMM topology to standard parents.

    Used when the caller already has a topology (not a file) and wants to
    prep it for ``addHydrogens``. Returns ``{(chain_id, res_id): original}``
    keyed by OpenMM's own identifiers so the post-pass can look residues
    up in the NEW topology produced by ``addHydrogens``.
    """
    saved: dict[tuple[str, str], str] = {}
    for res in topology.residues():
        parent = VARIANT_TO_PARENT.get(res.name)
        if parent is not None and parent != res.name:
            saved[(res.chain.id, res.id)] = res.name
            res.name = parent
    return saved


def restore_variants_post_addhydrogens(
    topology: Any,
    saved: dict[tuple[str, str], str] | SavedMap,
) -> None:
    """Restore variant residue names on the topology produced by addHydrogens.

    Accepts either shape of ``saved``:
    - ``{(chain_id, res_id): name}`` from
      :func:`rename_variants_to_parent_in_topology`
    - ``{(chain, resseq, icode): name}`` from
      :func:`text_rename_variants_to_parent` /
      :func:`scan_variant_names`
    """
    for res in topology.residues():
        key_topo = (res.chain.id, res.id)
        if key_topo in saved:
            res.name = saved[key_topo]  # type: ignore[index]
            continue
        key_text = (res.chain.id, str(res.id), "")
        if key_text in saved:
            res.name = saved[key_text]  # type: ignore[index]
            continue
        # Match on chain+resseq ignoring icode as a last-resort fallback.
        # Callers routinely merge text-shape (3-tuple) keys into
        # topology-shape (2-tuple) dicts, so skip 2-tuples here — those
        # were already covered by the direct lookup on line 226.
        for key, name in saved.items():
            if not (isinstance(key, tuple) and len(key) == 3):
                continue
            c, r, _i = key
            if isinstance(r, str) and c == res.chain.id and r == str(res.id):
                res.name = name
                break


def fix_lyn_hz_naming(
    topology: Any,
    saved: dict[tuple[str, str], str] | SavedMap,
    extra_lyn_keys: dict[Any, str] | None = None,
) -> int:
    """Rename ``HZ1 -> HZ3`` on every LYN residue in ``topology``.

    OpenMM's ``hydrogens.xml`` (used by ``Modeller.addHydrogens``) gates
    HZ3 by ``variant="LYS"``, so calling addHydrogens with the LYN
    variant produces ``HZ1 + HZ2``. But the authoritative reference —
    the AMBER14 ff14SB LYN residue template in
    ``openmm/app/data/amber14/protein.ff14SB.xml`` — defines the two
    NZ hydrogens as ``HZ2 + HZ3`` (HZ1 absent). HZ1 and HZ3 are
    chemically equivalent (same charge, same bond topology) so we
    rename in place. See ``docs/known-issues.md`` for the full
    reasoning.

    Identifies LYN residues from ``saved`` (any shape supported by
    :func:`restore_variants_post_addhydrogens`) plus optional
    ``extra_lyn_keys`` (e.g. PROPKA-assigned LYN that wasn't in the input
    PDB). Returns the number of atoms renamed.
    """
    lyn_ids: set[tuple[str, str]] = set()

    for key, name in saved.items():
        if name != "LYN":
            continue
        if isinstance(key, tuple) and len(key) == 2:
            lyn_ids.add(key)
        elif isinstance(key, tuple) and len(key) == 3:
            chain, resseq, _icode = key
            lyn_ids.add((chain, str(resseq)))

    if extra_lyn_keys:
        for key, name in extra_lyn_keys.items():
            if name != "LYN":
                continue
            if isinstance(key, tuple) and len(key) >= 2:
                chain, resseq = key[0], key[1]
                lyn_ids.add((chain, str(resseq)))

    renamed = 0
    for res in topology.residues():
        if (res.chain.id, res.id) not in lyn_ids:
            continue
        for atom in res.atoms():
            if atom.name == "HZ1":
                atom.name = "HZ3"
                renamed += 1
                break
    return renamed
