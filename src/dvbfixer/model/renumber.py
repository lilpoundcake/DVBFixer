"""Residue-number mapping for `dvbfixer model`.

Split out of the flat ``model.py`` in the Phase 2.2 follow-up work.
Contains the three-tier deterministic strategy for mapping Modeller's
continuous 1..N numbering back onto the user's original resseq /
insertion-code space.

Strategy (tried in order):

1. **K-finder** (:func:`_find_seqres_offset_by_resseq`) — the preferred
   path. Picks ``K`` = number of N-terminal SEQRES residues absent from
   ATOM, then assigns SEQRES position ``N`` → ``resseq = first_resseq +
   N - K``. Trusts input resseq jumps as authoritative gap positions,
   tolerates up to 10% letter mismatches (mutations), and handles
   N-term extras, internal gaps, and C-term extras uniformly.
2. **Needleman-Wunsch** (:func:`_align_atoms_to_seqres`) — semi-global
   fallback with affine gaps (open=-10, extend=-1, X-neutral, free end
   gaps on the SEQRES side). Used when K-finder can't converge because
   of excessive letter mismatches.
3. **Interpolated gap-fill** (:func:`_interpolate_gaps`) — walks
   Modeller's alignment mask and interpolates gap positions between
   flanking originals. Kept as a last resort for chains that defeat
   both deterministic paths (icode-bearing antibody Kabat numbering,
   etc.).

:func:`build_resnum_mapping` is the entry point; it picks the strategy
per chain and returns ``{(chain, model_resnum): (orig_resseq,
orig_icode)}`` for every atom the caller needs to renumber.
"""

from __future__ import annotations

from dvbfixer.model.cli import AA3TO1


def _get_original_resids_per_chain(original_lines, chain_order):
    """Get ordered list of unique (resSeq, iCode) per chain from original PDB."""
    chain_set = set(chain_order)
    result = {ch: [] for ch in chain_order}
    seen = {ch: set() for ch in chain_order}
    for line in original_lines:
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        ch = line[21]
        if ch not in chain_set:
            continue
        key = (int(line[22:26].strip()), line[26])
        if key not in seen[ch]:
            seen[ch].add(key)
            result[ch].append(key)
    return result


def _get_original_resletters_per_chain(original_lines, chain_order):
    """Per protein chain: ordered list of (resseq, icode, AA-letter) for ATOM records."""
    chain_set = set(chain_order)
    result = {ch: [] for ch in chain_order}
    seen = {ch: set() for ch in chain_order}
    for line in original_lines:
        if not line.startswith("ATOM  "):
            continue
        ch = line[21]
        if ch not in chain_set:
            continue
        resseq = int(line[22:26].strip())
        icode = line[26]
        key = (resseq, icode)
        if key not in seen[ch]:
            resname = line[17:20].strip()
            letter = AA3TO1.get(resname, 'X')
            seen[ch].add(key)
            result[ch].append((resseq, icode, letter))
    return result


def _find_seqres_offset_by_resseq(orig_residues, seqres_seq, tolerance=0.1):
    """Find offset K such that SEQRES[K + (resseq - first_resseq)] == atom letter
    for as many atoms as possible.

    K is the number of N-terminal SEQRES residues absent from ATOM (signal
    peptide etc.). This formula assumes input resseq is monotonic and
    contiguous within each stretch — i.e. the resseq jumps in the input
    correspond exactly to missing SEQRES positions. That makes it the
    "right" alignment for the user's perspective (gap-fill resseqs match
    the gap in the input).

    Picks the K with the highest match count. Accepts K when
    matches >= (1 - tolerance) * len(atoms), allowing a few mutations.

    Returns K (int >= 0) or None if no usable offset exists (icodes
    present, resseq exceeds SEQRES, or too many mismatches).
    """
    if not orig_residues:
        return 0
    if any(icode != ' ' for _, icode, _ in orig_residues):
        return None
    first_resseq = orig_residues[0][0]
    L = len(seqres_seq)
    max_offset = orig_residues[-1][0] - first_resseq
    if max_offset >= L:
        return None
    best_K = None
    best_matches = -1
    threshold = (1.0 - tolerance) * len(orig_residues)
    for K in range(L - max_offset):
        matches = 0
        for resseq, _, letter in orig_residues:
            pos = K + (resseq - first_resseq)
            if seqres_seq[pos] == letter or letter == 'X':
                matches += 1
        if matches > best_matches:
            best_matches = matches
            best_K = K
        if matches == len(orig_residues):
            return K  # perfect match — short-circuit
    if best_matches >= threshold:
        return best_K
    return None


def _align_atoms_to_seqres(orig_residues, seqres_seq):
    """Align ATOM letters to SEQRES letters via semi-global Needleman-Wunsch.

    The result is a list of length `len(orig_residues)` giving the SEQRES
    position (0-indexed) each ATOM residue maps to, or `None` for that entry
    if the residue could not be aligned. Returns `None` for the whole call
    when no usable alignment exists (e.g. ATOMs longer than SEQRES, no
    letter information).

    Semi-global = free end gaps on the SEQRES side, no skipping of ATOM
    letters. Affine gaps with a heavy gap-open penalty so consecutive
    gaps stay clumped (align2d's blosum scoring can split them, which is
    exactly the FcgRI failure mode this is bypassing).

    Robustness vs the previous `_find_seqres_offset`:
    - tolerates point mutations (mismatch score absorbed, alignment still
      placed in the only sensible spot)
    - tolerates `X` (unknown) residues with zero cost — they don't anchor
      or punish
    - handles N-terminal AND C-terminal SEQRES extras automatically via
      free end gaps
    - falls back gracefully when SEQRES is shorter than ATOMs
    """
    Q = len(orig_residues)
    R = len(seqres_seq)
    if Q == 0:
        return []
    if R == 0 or Q > R:
        return None

    query = [r[2] for r in orig_residues]

    MATCH, MISMATCH, X_NEUTRAL = 2, -1, 0
    GAP_OPEN, GAP_EXTEND = -10, -1
    NEG = float('-inf')

    # Three matrices for affine gaps:
    #   M[i][j]  = best score ending with query[i-1] aligned to seqres[j-1]
    #   IX[i][j] = best score ending with a SEQRES gap (skip seqres[j-1]; not used since gap-in-seqres = skip-query, disallowed)
    #   IY[i][j] = best score ending with a query gap (skip seqres[j-1])
    # Because we forbid skipping query letters, IX is dropped. Only IY is needed.
    M = [[NEG] * (R + 1) for _ in range(Q + 1)]
    IY = [[NEG] * (R + 1) for _ in range(Q + 1)]
    M[0][0] = 0
    # free end gaps on SEQRES side at the START: M[0][j] = 0 (any j)
    for j in range(R + 1):
        M[0][j] = 0
        IY[0][j] = NEG

    bt_M = [[None] * (R + 1) for _ in range(Q + 1)]
    bt_IY = [[None] * (R + 1) for _ in range(Q + 1)]

    for i in range(1, Q + 1):
        for j in range(1, R + 1):
            q, r = query[i - 1], seqres_seq[j - 1]
            if q == 'X' or r == 'X':
                sub = X_NEUTRAL
            elif q == r:
                sub = MATCH
            else:
                sub = MISMATCH
            # M[i][j]: match/mismatch from diag
            from_M = M[i - 1][j - 1] + sub
            from_IY = IY[i - 1][j - 1] + sub
            if from_M >= from_IY:
                M[i][j] = from_M
                bt_M[i][j] = 'M'
            else:
                M[i][j] = from_IY
                bt_M[i][j] = 'IY'
            # IY[i][j]: gap in query (skip seqres[j-1])
            open_iy = M[i][j - 1] + GAP_OPEN
            ext_iy = IY[i][j - 1] + GAP_EXTEND
            if open_iy >= ext_iy:
                IY[i][j] = open_iy
                bt_IY[i][j] = 'M'
            else:
                IY[i][j] = ext_iy
                bt_IY[i][j] = 'IY'

    # Free end gaps on SEQRES side at the END: best ending = max of M[Q][j] for j in 1..R
    best_j = 0
    best_score = NEG
    for j in range(1, R + 1):
        if M[Q][j] > best_score:
            best_score = M[Q][j]
            best_j = j
    if best_score == NEG:
        return None

    result = [None] * Q
    i, j = Q, best_j
    state = 'M'
    while i > 0 and j > 0:
        if state == 'M':
            result[i - 1] = j - 1
            prev = bt_M[i][j]
            i -= 1
            j -= 1
            state = prev
        else:  # IY
            prev = bt_IY[i][j]
            j -= 1
            state = prev
    if any(r is None for r in result):
        return None
    return result


def _interpolate_gaps(full_resids, region_len, renumber_from_1=False):
    """Fill `None` entries in `full_resids[:region_len]` by interpolating
    between flanking non-None entries.

    Handles:
    - internal gap (left and right templates present): place sequentially
      from `left + 1`. If there isn't enough room between `left` and
      `right` in the input's OWN numbering (the author never reserved
      resseqs for this gap — common for a genuinely-absent, unresolved
      loop/linker whose neighbors were numbered as if it didn't exist),
      shift every already-placed resid from this gap onward by the
      deficit instead of numbering backward into already-used resseqs
      (WARN emitted; see below).
    - N-terminal gap (no left): number backward from the first template.
      If that would produce non-positive resseqs (input's first template
      resseq ≤ gap_len), fall back to shifting the whole chain by
      ``1 - min_resseq`` so the first residue is 1 (WARN emitted).
    - C-terminal gap (no right within region): number forward from the
      last template.
    - all-None region: number sequentially from 1.

    If ``renumber_from_1=True``, the N-terminal branch unconditionally
    shifts to start at 1 regardless of whether the preserved numbering
    would fit.

    Operates in-place on `full_resids`. Entries beyond `region_len` are
    left alone (HETATM `.` slots are filled by the trailing sequential
    pass in `build_resnum_mapping`).
    """
    region_len = min(region_len, len(full_resids))
    i = 0
    n_term_shift_warned = False
    while i < region_len:
        if full_resids[i] is not None:
            i += 1
            continue
        gap_start = i
        while i < region_len and full_resids[i] is None:
            i += 1
        gap_end = i
        gap_len = gap_end - gap_start

        left = None
        for j in range(gap_start - 1, -1, -1):
            if full_resids[j] is not None:
                left = full_resids[j][0]
                break
        right = None
        for j in range(gap_end, region_len):
            if full_resids[j] is not None:
                right = full_resids[j][0]
                break

        if left is not None and right is not None:
            available = right - left - 1
            if available >= gap_len:
                for k in range(gap_len):
                    full_resids[gap_start + k] = (left + k + 1, ' ')
            else:
                # Not enough room in the input's own numbering for this
                # gap (e.g. an unresolved loop/linker whose flanking
                # residues were numbered back-to-back by the depositor,
                # as if the gap didn't exist). Numbering backward from
                # `right` would collide with resseqs already assigned to
                # `left` and everything before it — place the gap
                # sequentially from `left + 1` instead (matching the
                # "enough room" branch above) and shift every
                # already-placed resid from `right` onward forward by the
                # deficit, so nothing collides. This mirrors the
                # N-terminal branch's shift-the-rest-of-the-chain pattern
                # below.
                deficit = gap_len - available
                for k in range(gap_len):
                    full_resids[gap_start + k] = (left + k + 1, ' ')
                for j in range(gap_end, region_len):
                    if full_resids[j] is None:
                        continue
                    rs, ic = full_resids[j]
                    full_resids[j] = (rs + deficit, ic)
                print(f"  [renumber] internal gap of {gap_len} residue(s) "
                      f"has only {available} resseq(s) of room between "
                      f"original resSeq {left} and {right} — shifting "
                      f"every downstream residue in this chain by "
                      f"+{deficit} to avoid resSeq collisions. Original "
                      f"author numbering is not fully preserved past this "
                      f"point.")
        elif left is None and right is not None:
            candidates = [right - gap_len + k for k in range(gap_len)]
            if min(candidates) > 0 and not renumber_from_1:
                for k in range(gap_len):
                    full_resids[gap_start + k] = (candidates[k], ' ')
            else:
                # Shift the whole chain so the first residue is 1.
                # Compute shift so min(candidates) becomes 1.
                shift = 1 - min(candidates)
                if not n_term_shift_warned:
                    print(f"  [renumber] N-terminal fill would yield "
                          f"non-positive resseqs (min={min(candidates)}); "
                          f"shifting chain by +{shift} so first residue "
                          f"is 1. Pass --no-renumber-from-1 to preserve "
                          f"original numbering (currently ignored on "
                          f"non-positive fill).")
                    n_term_shift_warned = True
                for k in range(gap_len):
                    full_resids[gap_start + k] = (candidates[k] + shift, ' ')
                # Also shift EVERY already-placed resid in this region.
                for j in range(region_len):
                    if j >= gap_start and j < gap_end:
                        continue
                    if full_resids[j] is None:
                        continue
                    rs, ic = full_resids[j]
                    full_resids[j] = (rs + shift, ic)
        elif left is not None and right is None:
            for k in range(gap_len):
                full_resids[gap_start + k] = (left + k + 1, ' ')
        else:
            for k in range(gap_len):
                full_resids[gap_start + k] = (k + 1, ' ')


def build_resnum_mapping(per_chain_masks, all_chains, protein_chains, original_lines,
                          protein_seq_map=None, renumber_from_1=False):
    """Build mapping: (model_chain, model_resSeq) -> (original_resSeq, icode).

    Strategy (preferred): when input chain numbering is contiguous within the
    target sequence length and free of insertion codes, use a DETERMINISTIC
    resseq-based mapping: target-sequence position N (0-indexed) →
    resseq = first_resseq + N. Modeller produces one residue per target
    position, so this gives the correct resseq for both template and
    gap-filled positions without depending on align2d's gap placement
    (which can put a gap one position too early/late and then create cascade
    collisions in the dedup step, pushing modeled residues to the end of
    the chain).

    Fallback (when input has insertion codes, N-terminal SEQRES extras, or
    other non-contiguous numbering): walk align2d's mask one True position
    at a time, consuming `orig_rids` in order, and interpolate gap-filled
    positions between flanking originals (with N-term/C-term branching).

    Trailing positions beyond the protein region (HETATM `.` slots) are
    filled with sequential resseqs after the last protein residue.

    Only builds mapping for protein chains; non-protein chains are skipped.
    """
    protein_set = set(protein_chains)
    orig_resids = _get_original_resids_per_chain(original_lines, protein_chains)
    orig_letters = _get_original_resletters_per_chain(original_lines, protein_chains)
    mapping = {}
    protein_seq_map = protein_seq_map or {}

    for ci, chain in enumerate(all_chains):
        if ci >= len(per_chain_masks):
            continue
        mask = per_chain_masks[ci]

        if chain not in protein_set:
            # Non-protein chain — no renumbering needed
            continue

        orig_rids = orig_resids.get(chain, [])
        full_resids = [None] * len(mask)  # Each entry is (resSeq, iCode)

        seq_len = len(protein_seq_map.get(chain, ''))
        seqres_seq = protein_seq_map.get(chain, '')

        # Mapping strategy (tried in order):
        #
        # 1. resseq-based K-finder (preferred). Picks K = N-terminal-extras
        #    count and maps position N → resseq (first_resseq + N - K).
        #    Trusts the resseq jumps in the input PDB — when there's a
        #    jump from resseq 218 to 224 in chain A, that gap of 5 maps
        #    to SEQRES positions 198..202 (i.e. resseqs 219..223). This is
        #    what the user wants: gap-fill residues land inside the gap,
        #    not at the C-terminus.
        #
        # 2. Letter-only NW alignment (mutation-tolerant fallback). Used
        #    when K-finder bails (too many letter mismatches). NW handles
        #    point mutations but doesn't use resseq jumps, so it may place
        #    the gap one position off when the surrounding letters are
        #    ambiguous (e.g. repetitive sequence).
        #
        # 3. align2d's mask (last resort). Uses Modeller's alignment +
        #    flank interpolation. Susceptible to align2d's gap placement
        #    choices.
        seqres_index = None
        used = None
        if seq_len > 0 and orig_rids:
            K = _find_seqres_offset_by_resseq(orig_letters.get(chain, []), seqres_seq)
            if K is not None:
                used = 'K-finder'
                first_resseq = orig_rids[0][0]
                # Build seqres_index from K + (resseq - first_resseq)
                seqres_index = []
                for resseq, icode, _letter in orig_letters.get(chain, []):
                    if icode == ' ':
                        seqres_index.append(K + (resseq - first_resseq))
                    else:
                        seqres_index = None
                        break
            if seqres_index is None:
                # Try NW as fallback
                nw_result = _align_atoms_to_seqres(orig_letters.get(chain, []), seqres_seq)
                if nw_result is not None:
                    used = 'NW'
                    seqres_index = nw_result

        atom_only = orig_letters.get(chain, [])
        used_deterministic = (
            seqres_index is not None and len(seqres_index) == len(atom_only)
        )
        if used_deterministic:
            # Place each ATOM (protein-only) residue at its aligned SEQRES position.
            for k, (resseq, icode, _letter) in enumerate(atom_only):
                pos = seqres_index[k]
                if pos is None or pos < 0 or pos >= len(mask):
                    continue
                full_resids[pos] = (resseq, icode)
            _interpolate_gaps(full_resids, seq_len,
                              renumber_from_1=renumber_from_1)

            # Place HETATM residues (attached glycans, ligands) at the
            # trailing `.` slots beyond the protein region, preserving
            # their ORIGINAL resseqs (so an N-linked NAG keeps the resseq
            # it had in the input PDB rather than getting renumbered to
            # the next sequential integer after the last protein residue).
            #
            # BUT: "original resseq" here means whatever the upstream
            # standalone `renumber.py` step assigned earlier in the
            # pipeline — for a chain with no SEQRES records, that step
            # numbers a HETATM sequentially right after the chain's
            # ATOM-only residue count (confirmed: `renumber.py`'s
            # no-SEQRES branch does exactly this). If the true, FASTA-
            # complete sequence (`seq_len`, known here but NOT to that
            # earlier step) is LONGER than the ATOM-only count — i.e.
            # there's a real gap that `_interpolate_gaps` above just
            # filled — the newly gap-filled protein residues can
            # legitimately need exactly the resseq range the ligand was
            # naively given "one past the end". Confirmed on a real
            # structure (test/lipid/7x35_r_u.pdb): chain A had 267 ATOM
            # residues but a 278-residue FASTA sequence; a HETATM ligand
            # (PLM, palmitic acid) got resseq 268 from the upstream
            # no-SEQRES numbering, then `_interpolate_gaps` assigned the
            # first of the 11 gap-filled C-terminal residues that SAME
            # resseq 268 — a real protein residue and a ligand ending up
            # with an identical (chain, resseq), which silently
            # corrupted BOTH residues' atoms wherever anything downstream
            # looks up coordinates by (chain, resseq, atomname) without
            # also checking resname.
            protein_resseqs = {rs for rs, _ic in full_resids[:seq_len] if rs is not None}
            atom_keys = {(r, i) for r, i, _ in atom_only}
            hetatm_rids = [rid for rid in orig_rids if rid not in atom_keys]
            het_pos = seq_len
            next_safe_resseq = (max(protein_resseqs) if protein_resseqs else 0) + 1
            for het_rid in hetatm_rids:
                while het_pos < len(mask) and full_resids[het_pos] is not None:
                    het_pos += 1
                if het_pos >= len(mask):
                    break
                if het_rid[0] in protein_resseqs:
                    print(f"  WARNING: HETATM at chain {chain} original resseq "
                          f"{het_rid[0]}{het_rid[1].strip()} collides with a "
                          f"gap-filled protein residue at the same resseq — "
                          f"renumbering the HETATM to {next_safe_resseq} instead "
                          f"of preserving its original number.")
                    het_rid = (next_safe_resseq, ' ')
                    next_safe_resseq += 1
                full_resids[het_pos] = het_rid
                het_pos += 1
        else:
            # Fallback: align2d mask-based consumption + interpolation
            orig_idx = 0
            for i, is_template in enumerate(mask):
                if is_template and orig_idx < len(orig_rids):
                    full_resids[i] = orig_rids[orig_idx]
                    orig_idx += 1

            # Fill gap runs via the same interpolation logic the
            # deterministic path uses above — do not reimplement it here.
            # This used to be a second, hand-rolled copy of the same
            # gap-filling logic that had drifted out of sync (e.g. no
            # `renumber_from_1` handling, and the same collision bug in
            # its "not enough room" branch that `_interpolate_gaps` fixes).
            _interpolate_gaps(full_resids, len(mask),
                              renumber_from_1=renumber_from_1)

        # Fill any remaining None entries (mostly HETATM `.` slots beyond
        # the protein region) with sequential resseqs after the last
        # assigned residue.
        last_num = 0
        for i in range(len(full_resids)):
            if full_resids[i] is not None:
                last_num = full_resids[i][0]
            else:
                last_num += 1
                full_resids[i] = (last_num, ' ')

        # Deduplicate: only on the mask-based fallback path. The
        # deterministic path places residues at their natural SEQRES
        # positions with the original resseqs; running dedup would shift
        # HETATM resseqs (which can be numerically below protein resseqs)
        # to the wrong values.
        if not used_deterministic:
            for i in range(1, len(full_resids)):
                if full_resids[i][0] <= full_resids[i-1][0]:
                    full_resids[i] = (full_resids[i-1][0] + 1, ' ')

        # Modeller numbers all chains continuously: 1..N_total
        offset = sum(len(per_chain_masks[j]) for j in range(ci))
        for pos_in_chain in range(len(mask)):
            model_resnum = offset + pos_in_chain + 1
            mapping[(chain, model_resnum)] = full_resids[pos_in_chain]

    return mapping

