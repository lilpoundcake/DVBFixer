"""Regression coverage for shared FASTA/SEQRES alignment."""

from __future__ import annotations

from pathlib import Path

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
