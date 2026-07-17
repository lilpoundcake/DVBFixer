"""Mutation parsing and cleanup for ``dvbfixer prepare``.

Split out of the flat ``prepare.py`` in the Phase 2.3 follow-up work.

- :func:`parse_mutations` — turns a list of ``--mutate CHAIN:RESNUM:NEW_AA``
  specs into ``(substitutions_by_chain, variant_overrides, deletions)``.
  ``NEW_AA == "del"`` is a deletion; anything else is a substitution.
  Insertion codes are supported in ``RESNUM`` for deletions only —
  PDBFixer's ``applyMutations`` can't address iCode residues directly.
- :func:`apply_deletions_to_pdb_text` — the substantive part. Runs on
  raw PDB text BEFORE PDBFixer parses. For each deletion/substitution
  target: BFS from the sidechain anchor (ND2 / OG / OG1 / OH / SG / NZ /
  NH2 depending on residue type) into HETATM territory to pick up
  attached glycans, drop LINK records naming the affected residue,
  repair disulfide partners (CYX → CYS + drop HG). Substitution
  cleanup is filtered when the new AA's standard parent matches the
  old resname (e.g. CYS → CYX is just a protonation-variant rename).
- :func:`apply_mutations` — thin wrapper over ``fixer.applyMutations``
  for the substitution path.
"""

from __future__ import annotations

import sys
from pathlib import Path

from dvbfixer.ffutils.variants import VARIANT_TO_PARENT as _VARIANT_TO_STANDARD


def _split_resnum_icode(resnum_str):
    """Split a resnum string like '100A' into (100, 'A'). '100' → (100, ' ')."""
    import re
    m = re.fullmatch(r'(-?\d+)([A-Za-z]?)', resnum_str.strip())
    if not m:
        return None, None
    return int(m.group(1)), (m.group(2).upper() or ' ')


def parse_mutations(mutate_args):
    """Parse --mutate arguments into PDBFixer mutation format and deletion list.

    Input format: ['A:39:ALA', 'B:100:GLY', 'A:83:HIP', 'H:446:del', 'H:100A:del']
    Handles AMBER protonation variants (HIP, ASH, GLH, etc.) and deletions.

    Returns: (mutations_by_chain, variant_overrides, deletions)
      mutations_by_chain: dict chain_id -> [(resnum, standard_aa)]
      variant_overrides: dict (chain_id, resnum) -> variant_name
      deletions: list of (chain_id, resseq_int, icode_str) tuples
    """
    from collections import defaultdict
    mutations_by_chain = defaultdict(list)
    variant_overrides = {}
    deletions = []
    sub_keys = set()      # (chain, resseq_int, icode) — to detect conflict with del
    del_keys = set()
    for spec in mutate_args:
        parts = spec.split(':')
        if len(parts) != 3:
            print(f"Error: invalid --mutate format '{spec}' (expected "
                  f"CHAIN:RESNUM:NEW_AA or CHAIN:RESNUM:del)", file=sys.stderr)
            sys.exit(1)
        chain, resnum, new_aa = parts

        if new_aa.lower() == 'del':
            resseq, icode = _split_resnum_icode(resnum)
            if resseq is None:
                print(f"Error: invalid resnum '{resnum}' in deletion '{spec}'",
                      file=sys.stderr)
                sys.exit(1)
            key = (chain, resseq, icode)
            if key in sub_keys:
                print(f"Error: residue {chain}:{resnum} cannot be both "
                      f"substituted and deleted", file=sys.stderr)
                sys.exit(1)
            del_keys.add(key)
            deletions.append(key)
            continue

        new_aa_upper = new_aa.upper()
        # Disallow icodes in substitution targets — PDBFixer's applyMutations
        # expects bare integer resnums for the residue id.
        resseq_chk, icode_chk = _split_resnum_icode(resnum)
        if resseq_chk is not None and icode_chk != ' ':
            print(f"Error: substitution mutation '{spec}' uses an insertion "
                  f"code ('{resnum}') — PDBFixer can't address insertion-code "
                  f"residues directly. Renumber the chain first.", file=sys.stderr)
            sys.exit(1)
        if resseq_chk is not None:
            sub_key = (chain, resseq_chk, ' ')
            if sub_key in del_keys:
                print(f"Error: residue {chain}:{resnum} cannot be both "
                      f"substituted and deleted", file=sys.stderr)
                sys.exit(1)
            sub_keys.add(sub_key)

        # If it's a protonation variant, record it and use standard name for PDBFixer
        standard_aa = _VARIANT_TO_STANDARD.get(new_aa_upper, new_aa_upper)
        if new_aa_upper != standard_aa:
            variant_overrides[(chain, resnum)] = new_aa_upper

        mutations_by_chain[chain].append((resnum, standard_aa))
    return mutations_by_chain, variant_overrides, deletions


def apply_deletions_to_pdb_text(input_path, deletions, verbose=False,
                                  substitution_cleanups=None):
    """Apply residue cleanups to the PDB at the raw-text level.

    Two flavours of target:

    * **deletions** — fully remove the residue's atoms. The targeted residue
      itself disappears from the output. Triggered by `--mutate X:N:del`.
    * **substitution_cleanups** — keep the residue, but remove dependent
      atoms/records that won't survive the upcoming substitution. Triggered
      by `--mutate X:N:NEW_AA` when the new AA can't carry the old AA's
      sidechain-mediated bonds (e.g. ASN→ALA loses N-linked glycan;
      CYS→ALA loses the disulfide bridge). Pass as a list of
      `(chain, resseq, icode, old_resname, new_resname)` tuples; cleanup
      is only run when the parent residue names differ (so CYS→CYX, which
      is a protonation-variant rename, does NOT trigger cleanup and the
      SS bond is preserved).

    Both flavours trigger:

    * Glycan walk: BFS through CONECT from the residue's sidechain anchor
      (ASN ND2, SER OG, THR OG1, TYR OH, CYS SG, LYS NZ, ARG NH2) into
      attached HETATM territory — the whole glycan tree is removed.
    * Disulfide partner repair: CYX→CYS rename + drop HG (so addHydrogens
      regenerates it); SSBOND record dropped.
    * LINK record drop with warning; partner left as-is.

    Operates on raw text BEFORE the PDBFixer pipeline.

    Returns: (cleaned_path, removed_residues_meta) where cleaned_path is a
    temp file (the original path if no targets), and meta is a list of
    dicts ready for the .dat file.
    """
    import math
    import tempfile as _tf

    substitution_cleanups = substitution_cleanups or []
    if not deletions and not substitution_cleanups:
        return input_path, []

    with open(input_path) as f:
        lines = f.readlines()

    # --- Pass A: index residues, CONECT, SSBOND, LINK ---------------------
    # (chain, resseq, icode) → list of line indices for that residue's atoms
    res_lines = {}      # reskey → [line_idx, ...]
    res_resname = {}    # reskey → resname
    serial_to_reskey = {}  # int serial → reskey
    serial_to_atomname = {}
    serial_is_hetatm = {}
    serial_to_coord = {}
    chain_seq_order = {}  # chain_id → list of reskeys in file order

    for i, line in enumerate(lines):
        if not (line.startswith('ATOM  ') or line.startswith('HETATM')):
            continue
        chain = line[21]
        try:
            resseq = int(line[22:26].strip())
        except ValueError:
            continue
        icode = line[26] if len(line) > 26 else ' '
        resname = line[17:20].strip()
        atomname = line[12:16].strip()
        try:
            serial = int(line[6:11].strip())
        except ValueError:
            serial = None
        reskey = (chain, resseq, icode)

        if reskey not in res_lines:
            res_lines[reskey] = []
            res_resname[reskey] = resname
            chain_seq_order.setdefault(chain, []).append(reskey)
        res_lines[reskey].append(i)

        if serial is not None:
            serial_to_reskey[serial] = reskey
            serial_to_atomname[serial] = atomname
            serial_is_hetatm[serial] = line.startswith('HETATM')
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                serial_to_coord[serial] = (x, y, z)
            except ValueError:
                pass

    # CONECT graph
    conect_graph = {}
    conect_line_indices = {}  # int_serial → set of line indices (for filtering)
    conect_line_serials = {}  # line_idx → list of int serials
    for i, line in enumerate(lines):
        if not line.startswith('CONECT'):
            continue
        try:
            serials = [int(line[6 + j * 5: 6 + (j + 1) * 5].strip())
                       for j in range(min(11, (len(line.rstrip()) - 6) // 5))
                       if line[6 + j * 5: 6 + (j + 1) * 5].strip()]
        except ValueError:
            continue
        if not serials:
            continue
        source = serials[0]
        partners = serials[1:]
        for p in partners:
            conect_graph.setdefault(source, set()).add(p)
            conect_graph.setdefault(p, set()).add(source)
        for s in serials:
            conect_line_indices.setdefault(s, set()).add(i)
        conect_line_serials[i] = serials

    # SSBOND records — parse partner residues
    ssbond_line_indices = []
    ssbond_pairs = []  # (line_idx, reskey1, reskey2)
    for i, line in enumerate(lines):
        if not line.startswith('SSBOND'):
            continue
        try:
            ch1 = line[15]
            rs1 = int(line[17:21].strip())
            ic1 = line[21] if len(line) > 21 else ' '
            ch2 = line[29]
            rs2 = int(line[31:35].strip())
            ic2 = line[35] if len(line) > 35 else ' '
        except (ValueError, IndexError):
            continue
        rk1 = (ch1, rs1, ic1)
        rk2 = (ch2, rs2, ic2)
        ssbond_line_indices.append(i)
        ssbond_pairs.append((i, rk1, rk2))

    # LINK records — parse partner residues (cols are similar to SSBOND)
    link_pairs = []  # (line_idx, reskey1, reskey2)
    for i, line in enumerate(lines):
        if not line.startswith('LINK'):
            continue
        try:
            ch1 = line[21]
            rs1 = int(line[22:26].strip())
            ic1 = line[26] if len(line) > 26 else ' '
            ch2 = line[51]
            rs2 = int(line[52:56].strip())
            ic2 = line[56] if len(line) > 56 else ' '
        except (ValueError, IndexError):
            continue
        rk1 = (ch1, rs1, ic1)
        rk2 = (ch2, rs2, ic2)
        link_pairs.append((i, rk1, rk2))

    # --- Pass B: resolve each deletion target ---------------------------
    SIDECHAIN_ANCHORS = {
        'ASN': 'ND2', 'GLN': 'NE2',
        'SER': 'OG', 'THR': 'OG1', 'TYR': 'OH', 'HYP': 'OD1',
        'CYS': 'SG', 'CYX': 'SG', 'CYM': 'SG',
        'LYS': 'NZ', 'ARG': 'NH2',
    }

    deletion_set = set(deletions)
    atoms_to_remove = set()        # serial numbers
    residues_to_remove = set()     # reskeys (incl. attached glycans)
    cyx_partners_to_repair = set() # reskeys (rename CYX→CYS, drop HG)
    dropped_ssbond_lines = set()
    dropped_link_lines = set()
    removed_residues_meta = []

    # Unify deletion + substitution-cleanup targets. Each entry is
    # (reskey, kind, old_resname, new_resname). For substitutions, only
    # run cleanup when the new AA's parent residue name differs from the
    # old — so CYS→CYX (a protonation variant rename) does NOT trigger
    # cleanup and the SS bond is preserved.
    targets = []
    for dkey in deletions:
        targets.append((dkey, "delete", None, None))
    for entry in substitution_cleanups:
        ch, rs, ic, old_rn, new_rn = entry
        rkey = (ch, rs, ic)
        if old_rn == new_rn:
            continue  # CYS→CYX style rename — no cleanup
        targets.append((rkey, "substitution", old_rn, new_rn))

    for tkey, tkind, t_oldrn, t_newrn in targets:
        if tkey not in res_lines:
            chain, rs, ic = tkey
            ic_disp = '' if ic == ' ' else ic
            new_aa_disp = 'del' if tkind == 'delete' else (t_newrn or 'NEW_AA')
            print(f"Warning: residue {chain}:{rs}{ic_disp} not found in input — "
                  f"skipping --mutate {chain}:{rs}{ic_disp}:{new_aa_disp}",
                  file=sys.stderr)
            continue
        if tkind == "delete":
            dkey = tkey
            residues_to_remove.add(dkey)
            for ln_idx in res_lines[dkey]:
                line = lines[ln_idx]
                try:
                    s = int(line[6:11].strip())
                    atoms_to_remove.add(s)
                except ValueError:
                    pass
        else:
            # Substitution: keep the residue, only its dependent atoms are
            # cleaned up. dkey is reused as the local name for the target.
            dkey = tkey

        resname = res_resname[dkey]
        # --- Glycan walk (BFS) from the sidechain anchor ---
        glycan_resnames = []
        anchor_name = SIDECHAIN_ANCHORS.get(resname)
        if anchor_name:
            anchor_serial = None
            for ln_idx in res_lines[dkey]:
                if lines[ln_idx][12:16].strip() == anchor_name:
                    try:
                        anchor_serial = int(lines[ln_idx][6:11].strip())
                    except ValueError:
                        pass
                    break
            if anchor_serial is not None:
                # BFS over CONECT, crossing only into HETATM atoms of other residues
                visited = {anchor_serial}
                queue = [anchor_serial]
                while queue:
                    s = queue.pop()
                    for p in conect_graph.get(s, ()):
                        if p in visited:
                            continue
                        visited.add(p)
                        p_reskey = serial_to_reskey.get(p)
                        if p_reskey is None:
                            continue
                        # Only cross into HETATM residues != deleted residue
                        if not serial_is_hetatm.get(p, False):
                            continue
                        if p_reskey == dkey:
                            queue.append(p)  # same residue — keep walking
                            continue
                        if p_reskey in residues_to_remove:
                            queue.append(p)
                            continue
                        residues_to_remove.add(p_reskey)
                        glycan_resnames.append({
                            "chain": p_reskey[0],
                            "resid": str(p_reskey[1]),
                            "icode": p_reskey[2] if p_reskey[2] != ' ' else '',
                            "resname": res_resname.get(p_reskey, '?'),
                        })
                        for line_idx in res_lines.get(p_reskey, []):
                            try:
                                ss = int(lines[line_idx][6:11].strip())
                                atoms_to_remove.add(ss)
                                queue.append(ss)
                            except ValueError:
                                pass

        # --- Disulfide partner (CYS/CYX) ---
        disulfide_partner = None
        if resname in ('CYS', 'CYX'):
            for ln_idx_, rk1, rk2 in ssbond_pairs:
                if rk1 == dkey or rk2 == dkey:
                    dropped_ssbond_lines.add(ln_idx_)
                    partner = rk2 if rk1 == dkey else rk1
                    if partner != dkey and partner not in residues_to_remove:
                        cyx_partners_to_repair.add(partner)
                        disulfide_partner = partner
            # Also check CONECT SG→SG edges
            sg_serial = None
            for ln_idx in res_lines[dkey]:
                if lines[ln_idx][12:16].strip() == 'SG':
                    try:
                        sg_serial = int(lines[ln_idx][6:11].strip())
                    except ValueError:
                        pass
                    break
            if sg_serial is not None:
                for p in conect_graph.get(sg_serial, ()):
                    if serial_to_atomname.get(p) == 'SG':
                        p_reskey = serial_to_reskey.get(p)
                        if p_reskey and p_reskey != dkey and p_reskey not in residues_to_remove:
                            cyx_partners_to_repair.add(p_reskey)
                            if disulfide_partner is None:
                                disulfide_partner = p_reskey

        # --- LINK records mentioning this residue ---
        for ln_idx_, rk1, rk2 in link_pairs:
            if rk1 == dkey or rk2 == dkey:
                dropped_link_lines.add(ln_idx_)
                partner = rk2 if rk1 == dkey else rk1
                if partner not in residues_to_remove and partner != dkey:
                    print(f"Warning: LINK record between {dkey} and {partner} "
                          f"will be dropped — partner ({res_resname.get(partner, '?')}) "
                          f"left as-is. Inspect manually if it needs repair.",
                          file=sys.stderr)

        # --- Terminal vs internal classification (deletions only) ---
        if tkind == "delete":
            chain = dkey[0]
            order = chain_seq_order.get(chain, [])
            try:
                idx_in_chain = order.index(dkey)
            except ValueError:
                idx_in_chain = -1
            prev_reskey = None
            next_reskey = None
            # Walk backward to find prev that isn't being deleted
            for j in range(idx_in_chain - 1, -1, -1):
                if order[j] not in deletion_set:
                    prev_reskey = order[j]
                    break
            for j in range(idx_in_chain + 1, len(order)):
                if order[j] not in deletion_set:
                    next_reskey = order[j]
                    break

            gap_type = "internal"
            if prev_reskey is None and next_reskey is None:
                gap_type = "whole_chain"
            elif prev_reskey is None:
                gap_type = "terminal_N"
            elif next_reskey is None:
                gap_type = "terminal_C"
        else:
            # Substitution cleanup: residue stays — no gap classification.
            prev_reskey = None
            next_reskey = None
            gap_type = "substitution"

        # Compute the bridge distance (prev.C → next.N) for internal gaps
        gap_distance = None
        if gap_type == "internal":
            prev_c = None
            next_n = None
            for ln_idx in res_lines.get(prev_reskey, []):
                if lines[ln_idx][12:16].strip() == 'C':
                    try:
                        s = int(lines[ln_idx][6:11].strip())
                        prev_c = serial_to_coord.get(s)
                    except ValueError:
                        pass
                    break
            for ln_idx in res_lines.get(next_reskey, []):
                if lines[ln_idx][12:16].strip() == 'N':
                    try:
                        s = int(lines[ln_idx][6:11].strip())
                        next_n = serial_to_coord.get(s)
                    except ValueError:
                        pass
                    break
            if prev_c is not None and next_n is not None:
                gap_distance = math.sqrt(
                    sum((a - b) ** 2 for a, b in zip(prev_c, next_n))
                )

        chain_, rs_, ic_ = dkey
        meta = {
            "chain": chain_,
            "resid": str(rs_),
            "icode": ic_ if ic_ != ' ' else '',
            "resname": resname,
            "removed_atoms": len(res_lines[dkey]) if tkind == "delete" else 0,
            "gap_type": gap_type,
            "prev_residue": (f"{prev_reskey[0]}/{res_resname.get(prev_reskey, '?')}"
                             f"/{prev_reskey[1]}{prev_reskey[2].strip()}"
                             if prev_reskey else None),
            "next_residue": (f"{next_reskey[0]}/{res_resname.get(next_reskey, '?')}"
                             f"/{next_reskey[1]}{next_reskey[2].strip()}"
                             if next_reskey else None),
            "gap_distance_A": (round(gap_distance, 2)
                               if gap_distance is not None else None),
            "linked_glycan_residues": glycan_resnames,
            "disulfide_partner_repaired": (
                f"{disulfide_partner[0]}/{res_resname.get(disulfide_partner, '?')}"
                f"/{disulfide_partner[1]}{disulfide_partner[2].strip()}"
                if disulfide_partner else None),
        }
        if tkind == "substitution":
            meta["substituted_to"] = t_newrn
        removed_residues_meta.append(meta)
        if verbose:
            ic_disp = '' if ic_ == ' ' else ic_
            if tkind == "delete":
                print(f"  Delete {chain_}:{rs_}{ic_disp} ({resname}) — "
                      f"{gap_type}{', ' + str(len(glycan_resnames)) + ' glycan res' if glycan_resnames else ''}"
                      f"{', SS partner ' + str(disulfide_partner) + ' → CYS' if disulfide_partner else ''}")
            elif glycan_resnames or disulfide_partner:
                print(f"  Substitute {chain_}:{rs_}{ic_disp} ({resname}→{t_newrn}) — "
                      f"cleanup: "
                      f"{str(len(glycan_resnames)) + ' glycan res' if glycan_resnames else ''}"
                      f"{', SS partner ' + str(disulfide_partner) + ' → CYS' if disulfide_partner else ''}")
            if gap_distance is not None and gap_distance > 5.0:
                print(f"    WARNING: gap C(i-1)→N(i+1) distance = "
                      f"{gap_distance:.2f} Å — chain may need pulling/modeling")

    # --- Pass C: write the cleaned PDB --------------------------------------
    out_fd, out_path = _tf.mkstemp(suffix='_deldel.pdb')
    import os as _os
    _os.close(out_fd)

    def _rewrite_conect(line_idx):
        serials = conect_line_serials.get(line_idx, [])
        kept = [s for s in serials if s not in atoms_to_remove]
        if len(kept) < 2:
            return None
        # Pad each serial to 5 chars
        parts = ['CONECT']
        for s in kept:
            parts.append(f"{s:>5d}")
        return ''.join(parts) + '\n'

    with open(out_path, 'w') as f:
        for i, line in enumerate(lines):
            if line.startswith(('ATOM  ', 'HETATM')):
                try:
                    s = int(line[6:11].strip())
                except ValueError:
                    s = None
                if s is not None and s in atoms_to_remove:
                    continue
                # CYX partner repair: rename CYX → CYS on this line, drop HG
                chain = line[21]
                try:
                    rs = int(line[22:26].strip())
                except ValueError:
                    f.write(line)
                    continue
                ic = line[26] if len(line) > 26 else ' '
                reskey = (chain, rs, ic)
                if reskey in cyx_partners_to_repair:
                    atomname = line[12:16].strip()
                    if atomname == 'HG':
                        continue  # drop existing HG (will be re-added)
                    if line[17:20].strip() == 'CYX':
                        line = line[:17] + 'CYS' + line[20:]
                f.write(line)
            elif line.startswith('CONECT'):
                rewritten = _rewrite_conect(i)
                if rewritten is not None:
                    f.write(rewritten)
            elif line.startswith('SSBOND') and i in dropped_ssbond_lines:
                continue
            elif line.startswith('LINK') and i in dropped_link_lines:
                continue
            else:
                f.write(line)

    if verbose:
        print(f"  Deletion-cleaned PDB: {out_path}")
    return Path(out_path), removed_residues_meta


def apply_mutations(fixer, mutations_by_chain, verbose=False):
    """Apply point mutations using PDBFixer."""
    res_lookup = {}
    for res in fixer.topology.residues():
        res_lookup[(res.chain.id, res.id)] = res.name

    for chain_id, muts in mutations_by_chain.items():
        pdbfixer_muts = []
        for resnum, new_aa in muts:
            old_aa = res_lookup.get((chain_id, resnum))
            if old_aa is None:
                print(f"Warning: residue {chain_id}:{resnum} not found, skipping mutation")
                continue
            mut_str = f"{old_aa}-{resnum}-{new_aa}"
            pdbfixer_muts.append(mut_str)
            if verbose:
                print(f"  Mutation: {chain_id}:{old_aa}{resnum} -> {new_aa}")
        if pdbfixer_muts:
            fixer.applyMutations(pdbfixer_muts, chain_id)
