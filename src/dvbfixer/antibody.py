"""Antibody-aware numbering for `dvbfixer renumber`.

Supports five schemes: Kabat, Chothia, IMGT, Martin, EU.

For V-domains (variable region) the four V-domain schemes (Kabat / Chothia /
IMGT / Martin) are produced by ANARCI directly. EU's V-domain is identical to
Kabat's so it's handled the same way.

For C-domains (constant region) only EU has a fully-defined position table for
IgG (Edelman 1969 numbering). Kabat / Chothia / Martin don't define C-domain
positions — when one of those is requested, EU is used for C-domains and a
note is printed. IMGT C-domain numbering is a separate table; for now we
fall back to EU there too.

Reference sequences (human germline) are embedded inline below. They cover
human IgG1 heavy-chain constant region, human kappa light constant, and
human lambda light constant. For input chains matching a non-human / non-IgG1
construct, the Needleman-Wunsch alignment still works as long as the chain
shares >70% identity with the reference — it absorbs single-residue allotype
variation as substitutions (X-neutral) without losing the EU positions.

Partial chains (a structure missing the V-domain, or hinge-truncated Fc, or
just FR1 of the variable region) are handled automatically: the V-domain
detection and the C-domain alignments run independently and only emit
positions for residues that actually align.
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# EU-numbered reference sequences (human germline)
# ---------------------------------------------------------------------------
# Heavy-chain constant region (UniProt P01857, human IgG1).
# EU positions 118..447 (one entry per letter; total 330 residues).
# CH1: 118-215, hinge: 216-230, CH2: 231-340, CH3: 341-447.
IGG1_HEAVY_CONST_EU_START = 118
IGG1_HEAVY_CONST_SEQ = (
    "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLG"
    "TQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDP"
    "EVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQV"
    "YTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSC"
    "SVMHEALHNHYTQKSLSLSPGK"
)

# Kappa light constant (UniProt P01834, human IGKC).
# EU positions 108..214 (107 residues).
CK_EU_START = 108
CK_SEQ = (
    "RTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKA"
    "DYEKHKVYACEVTHQGLSSPVTKSFNRGEC"
)

# Lambda light constant (UniProt P0CG04, human IGLC2). EU positions 108..214 (107 residues).
CL_EU_START = 108
CL_SEQ = (
    "GQPKANPTVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADGSPVKAGVETTKPSKQSNNKYAASSYLSLTPE"
    "QWKSHRSYSCQVTHEGSTVEKTVAPTECS"
)


_ALL_SCHEMES = {"kabat", "chothia", "imgt", "martin", "eu", "aho"}
_ANARCI_V_SCHEMES = {
    # Map user-facing scheme → scheme name accepted by anarci.anarci()
    "kabat": "kabat",
    "chothia": "chothia",
    "imgt": "imgt",
    "martin": "martin",
    "aho": "aho",
    # EU's V-domain numbering matches Kabat's V-domain numbering.
    "eu": "kabat",
}


def _have_anarci() -> bool:
    try:
        import anarci  # noqa: F401
        return True
    except ImportError:
        return False


# ---------------------------------------------------------------------------
# Needleman-Wunsch alignment (lightweight; constant-region sized inputs only)
# ---------------------------------------------------------------------------

def _nw_align(query: str, ref: str) -> list[int] | None:
    """Semi-global NW alignment of `query` against `ref` (both AA strings).

    Returns a list of length `len(query)` where `result[i]` is the 0-based
    position in `ref` that query letter `i` aligns to (or `None` for that
    entry if the position couldn't be placed). Free end gaps on the REF
    side — the query can sit at any offset inside the reference.

    Returns `None` if the alignment is infeasible (`query` longer than
    `ref`, or every position scores below 0).
    """
    Q = len(query)
    R = len(ref)
    if Q == 0:
        return []
    if R == 0 or Q > R:
        return None

    MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND, X_NEUTRAL = 2, -1, -10, -1, 0
    NEG = float("-inf")
    M = [[NEG] * (R + 1) for _ in range(Q + 1)]
    IY = [[NEG] * (R + 1) for _ in range(Q + 1)]
    for j in range(R + 1):
        M[0][j] = 0
    bt_M = [[None] * (R + 1) for _ in range(Q + 1)]
    bt_IY = [[None] * (R + 1) for _ in range(Q + 1)]

    for i in range(1, Q + 1):
        for j in range(1, R + 1):
            q, r = query[i - 1], ref[j - 1]
            if q == "X" or r == "X":
                sub = X_NEUTRAL
            elif q == r:
                sub = MATCH
            else:
                sub = MISMATCH
            from_M = M[i - 1][j - 1] + sub
            from_IY = IY[i - 1][j - 1] + sub
            if from_M >= from_IY:
                M[i][j] = from_M
                bt_M[i][j] = "M"
            else:
                M[i][j] = from_IY
                bt_M[i][j] = "IY"
            open_iy = M[i][j - 1] + GAP_OPEN
            ext_iy = IY[i][j - 1] + GAP_EXTEND
            if open_iy >= ext_iy:
                IY[i][j] = open_iy
                bt_IY[i][j] = "M"
            else:
                IY[i][j] = ext_iy
                bt_IY[i][j] = "IY"

    best_j = 0
    best_score = NEG
    for j in range(1, R + 1):
        if M[Q][j] > best_score:
            best_score = M[Q][j]
            best_j = j
    if best_score == NEG or best_score <= 0:
        return None

    result: list[int | None] = [None] * Q
    i, j = Q, best_j
    state = "M"
    while i > 0 and j > 0:
        if state == "M":
            result[i - 1] = j - 1
            prev = bt_M[i][j]
            i -= 1
            j -= 1
            state = prev
        else:  # IY: skip ref position
            prev = bt_IY[i][j]
            j -= 1
            state = prev

    return result


# ---------------------------------------------------------------------------
# V-domain numbering via ANARCI
# ---------------------------------------------------------------------------

def _number_v_domain(seq: str, scheme: str):
    """Run ANARCI against `seq` and return a list of (input_idx, resseq, icode)
    for every residue ANARCI placed in the V-domain, or None when the chain
    has no detectable V-domain.

    `scheme` is the user-requested scheme name (kabat/chothia/imgt/martin/eu/aho).
    """
    if scheme not in _ANARCI_V_SCHEMES:
        return None
    if not _have_anarci():
        return None
    anarci_scheme = _ANARCI_V_SCHEMES[scheme]
    try:
        from anarci import anarci
    except ImportError:
        return None
    try:
        numbered, alignment_details, _ = anarci(
            [("query", seq)], scheme=anarci_scheme,
            allow={"H", "K", "L"},
        )
    except Exception:
        return None
    if not numbered or not numbered[0]:
        return None
    placements = []
    # ANARCI may return multiple hits (e.g. scFv with two V-domains).
    for hit in numbered[0]:
        nums, start, end = hit
        cursor = start
        for (resseq, icode), aa in nums:
            if aa == "-":
                continue
            if cursor > end:
                break
            placements.append((cursor, resseq, icode if icode != " " else " "))
            cursor += 1
    return placements


# ---------------------------------------------------------------------------
# C-domain numbering against bundled EU references
# ---------------------------------------------------------------------------

_C_REFERENCES = [
    ("heavy", IGG1_HEAVY_CONST_SEQ, IGG1_HEAVY_CONST_EU_START),
    ("kappa", CK_SEQ, CK_EU_START),
    ("lambda", CL_SEQ, CL_EU_START),
]


def _number_c_domain(seq: str):
    """Align `seq` against the three C-domain references (heavy IgG1, kappa,
    lambda) and pick the best-scoring alignment. Returns a list of
    (input_idx, resseq, icode) for residues that aligned, plus the
    reference name. Returns (None, None) if no alignment placed >=50% of
    residues.
    """
    best = None
    best_score = 0
    best_ref = None
    for ref_name, ref_seq, ref_start in _C_REFERENCES:
        idx = _nw_align(seq, ref_seq)
        if idx is None:
            continue
        # Score = number of placed positions where letters match
        placed = [(i, p) for i, p in enumerate(idx) if p is not None]
        matches = sum(1 for i, p in placed if seq[i] == ref_seq[p])
        if matches > best_score:
            best_score = matches
            best = (idx, ref_seq, ref_start)
            best_ref = ref_name

    if best is None or best_score < max(8, len(seq) // 2):
        return None, None

    idx, ref_seq, ref_start = best
    placements = []
    for i, p in enumerate(idx):
        if p is None:
            continue
        # EU position = ref_start + p (0-indexed p, ref_start is 1-indexed)
        placements.append((i, ref_start + p, " "))
    return placements, best_ref


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def number_chain(seq: str, scheme: str, verbose: bool = False):
    """Compute antibody numbering for one chain's amino-acid sequence.

    `seq` is a 1-letter-code string. `scheme` is one of the five user-facing
    scheme names (kabat / chothia / imgt / martin / eu). Returns:

      {
        "is_antibody": bool,
        "placements": [(input_idx, resseq, icode), ...],  # one entry per
                                                          # placed residue
        "segments": [(start_idx, end_idx, kind, scheme_used), ...],
        "warnings": [...],
      }

    `input_idx` is the 0-based index into `seq`. Residues not placed in any
    domain are absent from `placements`; the caller should fall back to a
    different strategy (e.g. SEQRES-based renumbering) for those positions.
    """
    scheme = scheme.lower()
    if scheme not in _ALL_SCHEMES:
        raise ValueError(f"Unknown antibody scheme '{scheme}'. "
                         f"Valid: {sorted(_ALL_SCHEMES)}")

    warnings: list[str] = []
    placements = []
    segments = []
    occupied = set()

    # V-domain pass
    v = _number_v_domain(seq, scheme)
    if v:
        placements.extend(v)
        for (i, _r, _ic) in v:
            occupied.add(i)
        v_idxs = [p[0] for p in v]
        segments.append((min(v_idxs), max(v_idxs), "V", scheme))
        if verbose:
            print(f"    V-domain: positions {min(v_idxs)}-{max(v_idxs)} ({scheme})")

    # C-domain pass — run on the residues NOT placed by ANARCI (a contiguous
    # post-V slice in the common case; the whole sequence for C-only chains).
    if not v:
        post_seq = seq
        post_offset = 0
    else:
        # Use everything after the last V-domain index
        v_end = max(p[0] for p in v) + 1
        post_seq = seq[v_end:]
        post_offset = v_end

    if post_seq:
        c_place, c_ref = _number_c_domain(post_seq)
        if c_place:
            # If the V-domain scheme (IMGT/Martin/Aho) extends past EU
            # position 117, EU's first CH1 position (118) would collide
            # with V-domain numbering. In that case, shift the EU
            # numbering forward so that the first C-domain residue is at
            # max(V) + 5 (preserving a small gap as a visual break).
            shift = 0
            if v:
                max_v_resseq = max(p[1] for p in v)
                first_c_resseq = c_place[0][1]
                if first_c_resseq <= max_v_resseq:
                    shift = (max_v_resseq + 5) - first_c_resseq
                    if verbose:
                        print(f"    C-domain ({c_ref}): EU position {first_c_resseq} "
                              f"collides with V-domain (last V={max_v_resseq}); "
                              f"shifting EU numbering by +{shift}")
                        warnings.append(
                            f"Mixed {scheme.upper()}-V / EU-C numbering: EU "
                            f"shifted by +{shift} to avoid collision with V-domain."
                        )
            elif scheme not in ("eu",) and verbose:
                print(f"    C-domain ({c_ref}): using EU numbering "
                      f"(scheme '{scheme}' is V-domain only)")
            for (i, resseq, icode) in c_place:
                global_i = i + post_offset
                if global_i in occupied:
                    continue
                placements.append((global_i, resseq + shift, icode))
                occupied.add(global_i)
            c_idxs = [c_place[0][0] + post_offset, c_place[-1][0] + post_offset]
            seg_scheme = "eu" if shift == 0 else f"eu+{shift}"
            segments.append((min(c_idxs), max(c_idxs), f"C/{c_ref}", seg_scheme))
            if verbose:
                print(f"    C-domain: positions {min(c_idxs)}-{max(c_idxs)} "
                      f"(EU, reference={c_ref})")

    is_antibody = len(placements) > 0
    if not is_antibody:
        warnings.append("no antibody domain detected — falling back to default")

    placements.sort()
    return {
        "is_antibody": is_antibody,
        "placements": placements,
        "segments": segments,
        "warnings": warnings,
    }


def number_chain_to_mapping(atom_residues, scheme, verbose=False):
    """Convenience wrapper: given the chain's ordered list of
    (resseq, icode, resname) tuples (as returned by `renumber.get_atom_residues`),
    produce the renumbering map `{(old_resseq, old_icode): (new_resseq, new_icode)}`.

    Residues that weren't placed by the antibody numbering keep their old
    keys absent from the map; the caller should merge with a fallback.
    """
    # Build 1-letter sequence
    from dvbfixer.renumber import _AA3TO1  # type: ignore[attr-defined]
    seq = "".join(_AA3TO1.get(rn, "X") for (_rs, _ic, rn) in atom_residues)
    info = number_chain(seq, scheme, verbose=verbose)
    mapping = {}
    for (idx, resseq, icode) in info["placements"]:
        old_rs, old_ic, _ = atom_residues[idx]
        mapping[(old_rs, old_ic)] = (resseq, icode)
    return mapping, info
