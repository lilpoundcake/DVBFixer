"""Shared protein-sequence alignment for renumbering and modeling."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class SequenceAlignment:
    """Placement of each observed residue at a zero-based reference index."""

    positions: tuple[int, ...]
    score: int
    ambiguous: bool
    substitutions: tuple[tuple[int, int, str, str], ...]
    internal_gaps: tuple[tuple[int, int], ...]

    @property
    def start(self) -> int:
        return self.positions[0]

    @property
    def end(self) -> int:
        """Exclusive end of the aligned reference interval."""
        return self.positions[-1] + 1


def align_observed_to_reference(
    observed: str,
    reference: str,
    *,
    match_score: int = 2,
    mismatch_score: int = -3,
    gap_open_score: int = -5,
    gap_extend_score: int = -1,
) -> SequenceAlignment | None:
    """Affine semi-global alignment of an observed chain to its reference.

    Every observed residue is consumed. Reference overhangs are free and
    skipped internal reference residues represent missing structural density.
    Equal-scoring placements resolve deterministically to the leftmost end.
    """
    observed = observed.upper()
    reference = reference.upper()
    q_len, r_len = len(observed), len(reference)
    if not observed or q_len > r_len:
        return None

    neg = -10**12
    m = [[neg] * (r_len + 1) for _ in range(q_len + 1)]
    gap = [[neg] * (r_len + 1) for _ in range(q_len + 1)]
    count_m = [[0] * (r_len + 1) for _ in range(q_len + 1)]
    count_gap = [[0] * (r_len + 1) for _ in range(q_len + 1)]
    back_m = [[""] * (r_len + 1) for _ in range(q_len + 1)]
    back_gap = [[""] * (r_len + 1) for _ in range(q_len + 1)]

    for j in range(r_len + 1):
        m[0][j] = 0  # free reference prefix
        count_m[0][j] = 1

    for i in range(1, q_len + 1):
        for j in range(1, r_len + 1):
            q, r = observed[i - 1], reference[j - 1]
            sub = (0 if q == "X" or r == "X" else
                   match_score if q == r else mismatch_score)

            from_m = m[i - 1][j - 1]
            from_gap = gap[i - 1][j - 1]
            best = max(from_m, from_gap)
            m[i][j] = best + sub
            back_m[i][j] = "M" if from_m >= from_gap else "G"
            paths = 0
            if from_m == best:
                paths += count_m[i - 1][j - 1]
            if from_gap == best:
                paths += count_gap[i - 1][j - 1]
            count_m[i][j] = min(2, paths)

            opened = m[i][j - 1] + gap_open_score
            extended = gap[i][j - 1] + gap_extend_score
            best_gap = max(opened, extended)
            gap[i][j] = best_gap
            back_gap[i][j] = "M" if opened >= extended else "G"
            paths = 0
            if opened == best_gap:
                paths += count_m[i][j - 1]
            if extended == best_gap:
                paths += count_gap[i][j - 1]
            count_gap[i][j] = min(2, paths)

    best_score = max(m[q_len][1:])  # free reference suffix
    best_ends = [j for j in range(1, r_len + 1) if m[q_len][j] == best_score]
    if not best_ends:
        return None
    best_j = best_ends[0]
    ambiguous = len(best_ends) > 1 or sum(
        count_m[q_len][j] for j in best_ends
    ) > 1

    positions = [-1] * q_len
    i, j, state = q_len, best_j, "M"
    while i > 0:
        if j <= 0:
            return None
        if state == "M":
            positions[i - 1] = j - 1
            state = back_m[i][j]
            i -= 1
            j -= 1
        else:
            state = back_gap[i][j]
            j -= 1
    if any(p < 0 for p in positions):
        return None

    substitutions = tuple(
        (i, pos, observed[i], reference[pos])
        for i, pos in enumerate(positions)
        if observed[i] != "X" and reference[pos] != "X"
        and observed[i] != reference[pos]
    )
    internal_gaps = tuple(
        (left + 1, right)
        for left, right in zip(positions, positions[1:])
        if right > left + 1
    )
    return SequenceAlignment(
        positions=tuple(positions),
        score=best_score,
        ambiguous=ambiguous,
        substitutions=substitutions,
        internal_gaps=internal_gaps,
    )


def format_alignment_diagnostic(chain: str, result: SequenceAlignment) -> str:
    """Return a compact summary using one-based reference coordinates."""
    gaps = ", ".join(
        f"{start + 1}-{end}" for start, end in result.internal_gaps
    ) or "none"
    suffix = ", ambiguous-best" if result.ambiguous else ""
    return (
        f"chain {chain}: reference {result.start + 1}-{result.end}, "
        f"score {result.score}, substitutions {len(result.substitutions)}, "
        f"internal gaps {gaps}{suffix}"
    )
