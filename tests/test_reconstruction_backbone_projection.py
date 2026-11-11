from pathlib import Path

from scripts.screen_reconstruction_backbones import project

FIXTURE = Path(__file__).parent / "fixtures/protein_context/1UBQ.pdb"


def atoms():
    return [line for line in FIXTURE.read_text().splitlines(True) if line.startswith("ATOM  ") and int(line[22:26]) in (1, 2)]


def test_peptide_geometry_not_residue_number_adjacency_controls_projection():
    lines = [line[:22] + "  98" + line[26:] if int(line[22:26]) == 2 else line for line in atoms()]
    segments, _, metadata = project("".join(lines))
    assert len(segments) == 5  # two internal bonds per residue plus peptide C-N
    assert not metadata["unconnected_boundaries"]


def test_ter_is_a_break_even_when_coordinates_are_close():
    lines = atoms()
    index = next(i for i, line in enumerate(lines) if int(line[22:26]) == 2)
    segments, _, metadata = project("".join(lines[:index]) + "TER\n" + "".join(lines[index:]))
    assert len(segments) == 4
    assert metadata["unconnected_boundaries"][0]["connectivity"] is False


def test_chain_case_and_insertion_codes_survive_projection():
    lines = [line[:21] + ("D" if int(line[22:26]) == 1 else "d") + line[22:26] + "A" + line[27:] for line in atoms()]
    segments, _, _ = project("".join(lines))
    assert {s.first.chain for s in segments} == {"D", "d"}
    assert {s.first.insertion_code for s in segments} == {"A"}
    assert len(segments) == 4
