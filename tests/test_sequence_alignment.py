"""Regression coverage for shared FASTA/SEQRES alignment."""

from __future__ import annotations

from pathlib import Path

import pytest
from Bio import SeqIO

from dvbfixer import zbs
from dvbfixer.model.pipeline import get_atom_sequence, trim_terminal_gaps
from dvbfixer.renumber import main as renumber_main
from dvbfixer.sequence_alignment import align_observed_to_reference

_FIXTURES = Path(__file__).resolve().parent / "fixtures" / "numbering"


def _chain_c_reference() -> str:
    records = {
        record.id: str(record.seq)
        for record in SeqIO.parse(_FIXTURES / "8b01_renamed.fasta", "fasta")
    }
    return records["8B01_C"]


def _unique_resseqs(path: Path) -> list[int]:
    result = []
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM  "):
            continue
        resseq = int(line[22:26])
        if not result or result[-1] != resseq:
            result.append(resseq)
    return result


def test_8b01_chain_c_aligns_to_exact_contiguous_fasta_suffix() -> None:
    reference = _chain_c_reference()
    for name in ("8b01_a_b.pdb", "8b01_a_u.pdb"):
        lines = (_FIXTURES / name).read_text().splitlines(keepends=True)
        observed = get_atom_sequence(lines, "C")
        alignment = align_observed_to_reference(observed, reference)

        assert alignment is not None
        assert alignment.positions == tuple(range(403, 507))
        assert alignment.internal_gaps == ()
        assert not alignment.ambiguous


def test_no_terminal_trims_8b01_without_creating_internal_loop() -> None:
    reference = _chain_c_reference()
    for name in ("8b01_a_b.pdb", "8b01_a_u.pdb"):
        lines = (_FIXTURES / name).read_text().splitlines(keepends=True)
        observed = get_atom_sequence(lines, "C")
        trimmed = trim_terminal_gaps({"C": reference}, lines)

        assert trimmed == {"C": observed}
        assert len(trimmed["C"]) == 104


def test_no_terminal_keeps_real_internal_gap_while_trimming_ends() -> None:
    reference = "QQACDEFGHIKLMNPQRR"
    observed = "ACDEIKLMNPQ"  # reference FGH is unresolved internally
    alignment = align_observed_to_reference(observed, reference)

    assert alignment is not None
    assert (alignment.start, alignment.end) == (2, 16)
    assert alignment.internal_gaps == ((6, 9),)
    assert reference[alignment.start:alignment.end] == "ACDEFGHIKLMNPQ"


def test_locked_placement_preserves_internal_gap(tmp_path: Path) -> None:
    from dvbfixer.model.modeller_run import _fix_terminal_alignment, parse_alignment

    observed = "A" * 194 + "T"
    target = "A" * 194 + "G" * 10 + "T"
    raw = tmp_path / "raw.pir"
    fixed = tmp_path / "fixed.pir"
    raw.write_text(
        ">P1;template\nstructureX:template::::::::\n"
        f"{observed}*\n>P1;target\nsequence:target::::::::\n{target}*\n",
        encoding="utf-8",
    )

    _fix_terminal_alignment(
        raw,
        fixed,
        "template",
        {"C"},
        ["C"],
        [observed],
        [target],
        observed_target_indices_by_chain={"C": (*range(194), 204)},
    )

    assert parse_alignment(fixed) == "A" * 194 + "-" * 10 + "T"


def test_point_mutation_does_not_shift_downstream_residues() -> None:
    alignment = align_observed_to_reference("AGATVL", "AGSTVL")

    assert alignment is not None
    assert alignment.positions == tuple(range(6))
    assert alignment.substitutions == ((2, 2, "A", "S"),)


def test_ambiguous_repeat_uses_leftmost_best_placement() -> None:
    alignment = align_observed_to_reference("AAA", "AAAAA")

    assert alignment is not None
    assert alignment.positions == (0, 1, 2)
    assert alignment.ambiguous


def test_numbering_resolves_equal_sequence_scores_without_assuming_an_absolute_offset():
    observed = "WFPLAPCSTAALGKW"
    reference = "WFPLAPCSRSTSESTAALGKW"
    baseline = align_observed_to_reference(observed, reference)
    assert baseline is not None and baseline.ambiguous
    for offset in (-30, 100, 1000):
        ids = tuple((offset + i + (6 if i > 7 else 0), "") for i in range(len(observed)))
        result = align_observed_to_reference(observed, reference, residue_ids=ids)
        assert result is not None
        assert result.score == baseline.score
        assert result.positions[7] == 7  # retain the observed Cys–Ser peptide
        assert result.internal_gaps == ((8, 14),)
        assert reference[8:14] == "RSTSES"
        assert result.tie_breaker


def test_observed_peptide_geometry_overrides_a_previous_bad_renumbering():
    observed = "WFPLAPCSTAALGKW"
    reference = "WFPLAPCSRSTSESTAALGKW"
    # Wrong legacy numbering put the numerical jump before the observed SER.
    ids = tuple((100 + i + (6 if i >= 7 else 0), "") for i in range(len(observed)))
    continuity = tuple(False if i == 7 else True for i in range(len(observed) - 1))
    result = align_observed_to_reference(observed, reference, residue_ids=ids,
        peptide_continuity=continuity)
    assert result is not None and result.internal_gaps == ((8, 14),)
    assert result.positions[7] == 7


def test_evidence_never_overrides_a_better_sequence_score():
    baseline = align_observed_to_reference("ACDGHIK", "ACDEFGHIK")
    result = align_observed_to_reference("ACDGHIK", "ACDEFGHIK",
        peptide_continuity=(True,) * 6, residue_ids=tuple((i, "") for i in range(7)))
    assert result == baseline


def test_evidence_ranking_matches_exhaustive_small_alignments():
    import itertools
    import random

    rng = random.Random(7)
    for _ in range(100):
        size = rng.randint(2, 4)
        observed = "".join(rng.choice("AS") for _ in range(size))
        reference = "".join(rng.choice("AS") for _ in range(size + rng.randint(0, 3)))
        continuity = tuple(rng.choice((True, False, None)) for _ in range(size - 1))
        numbers = [10]
        for _boundary in range(size - 1):
            numbers.append(numbers[-1] + rng.choice((1, 3)))

        def rank(positions):
            score = sum(2 if observed[i] == reference[pos] else -3 for i, pos in enumerate(positions))
            peptide_splits = number_splits = 0
            for i, (left, right) in enumerate(zip(positions, positions[1:])):
                skipped = right - left - 1
                if skipped:
                    score += -5 - (skipped - 1)
                    peptide_splits += continuity[i] is True
                    number_splits += numbers[i + 1] == numbers[i] + 1
            return score, -peptide_splits, -number_splits

        candidates = list(itertools.combinations(range(len(reference)), size))
        best = max(map(rank, candidates))
        result = align_observed_to_reference(observed, reference,
            residue_ids=tuple((number, "") for number in numbers), peptide_continuity=continuity)
        assert result is not None
        assert rank(result.positions) == best
        assert result.score == best[0]
        assert result.ambiguous == (sum(rank(positions) == best for positions in candidates) > 1)


@pytest.mark.parametrize("kwargs", [
    {"residue_ids": ((1, ""),)}, {"peptide_continuity": (True,)},
])
def test_alignment_evidence_length_must_match_observed_residues(kwargs):
    with pytest.raises(ValueError, match="evidence must match"):
        align_observed_to_reference("ACD", "ACD", **kwargs)


def test_5j13_renumber_and_diffusion_preserve_the_observed_cys_ser_peptide(tmp_path):
    from dvbfixer.model.diffusion.contract import ResidueIdentity
    from dvbfixer.model.diffusion_cli import _target_sequences, build_cli_diffusion_request

    source = Path(__file__).parent / "fixtures/protein_context/5j13_selected.pdb"
    before = source.read_bytes()
    numbered = tmp_path / "numbered.pdb"
    renumber_main([str(source), "-o", str(numbered)])
    original_ser = next(line for line in before.decode().splitlines() if line.startswith("ATOM")
        and line[21] == "C" and line[22:26].strip() == "140" and line[12:16].strip() == "CA")
    mapped_ser = next(line for line in numbered.read_text().splitlines() if line.startswith("ATOM")
        and line[21] == "C" and line[30:54] == original_ser[30:54])
    assert mapped_ser[22:26].strip() == "168"
    numbered_text = numbered.read_text()
    misleading_numbering = "".join(
        line[:22] + " 174" + line[26:]
        if line.startswith("ATOM") and line[21] == "C" and line[22:26].strip() == "168" else line
        for line in numbered_text.splitlines(keepends=True))
    for name, text in (("original", before.decode()), ("renumbered", numbered_text), ("old-bad-numbering", misleading_numbering)):
        canonical = tmp_path / f"{name}-protein.pdb"
        canonical.write_text("".join(line for line in text.splitlines(keepends=True)
            if not line.startswith(("HETATM", "CONECT", "LINK  "))))
        request = build_cli_diffusion_request(canonical, _target_sequences(canonical, None), seeds=(7,),
            profile="protpardelle-1c-mps-local-v1")
        gaps = [gap for gap in request.gaps if gap.chain == "C"]
        gap = next(gap for gap in gaps if gap.target_interval.start == 168)
        assert gap.target_interval.stop == 174
        target = next(target.sequence for target in request.target_sequences if target.chain == "C")
        assert target[168:174] == "RSTSES"
        left_number = {"original": "140", "renumbered": "168", "old-bad-numbering": "174"}[name]
        assert gap.left_anchor == ResidueIdentity("C", left_number)
        assert gap.right_anchor == ResidueIdentity("C", "147" if name == "original" else "175")
        assert len(request.gaps) == 10
        assert sum(len(gap.generated_residues) for gap in request.gaps) == 124
    assert source.read_bytes() == before


def test_fasta_renumber_maps_both_8b01_inputs_to_full_sequence(tmp_path: Path) -> None:
    fasta = _FIXTURES / "8b01_renamed.fasta"
    for name in ("8b01_a_b.pdb", "8b01_a_u.pdb"):
        output = tmp_path / name
        renumber_main([
            str(_FIXTURES / name), "--fasta", str(fasta), "-o", str(output)
        ])
        assert _unique_resseqs(output) == list(range(404, 508))


def test_zbs_propagates_fasta_to_renumber(monkeypatch, tmp_path: Path) -> None:
    source = tmp_path / "input.pdb"
    source.write_text("END\n")
    fasta = tmp_path / "input.fasta"
    fasta.write_text(">chain_A\nA\n")
    captured = []

    def fake_renumber(argv):
        captured.extend(argv)
        output = Path(argv[argv.index("-o") + 1])
        output.write_text(source.read_text())

    monkeypatch.setattr("dvbfixer.renumber.main", fake_renumber)
    zbs.main([
        str(source), "--fasta", str(fasta), "--skip-model",
        "--skip-prepare", "--skip-minimize", "--no-align-to-input",
        "--no-postflight",
    ])

    assert captured[captured.index("--fasta") + 1] == str(fasta)
