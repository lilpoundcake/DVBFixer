from pathlib import Path

from dvbfixer.model.pipeline import parse_fasta


def test_fasta_chain_ids_are_case_sensitive(tmp_path: Path) -> None:
    fasta = tmp_path / "two-case-distinct-chains.fasta"
    fasta.write_text(">chain_D\nAAAA\n>chain_d\nGGGG\n")
    assert parse_fasta(fasta) == {"D": "AAAA", "d": "GGGG"}
