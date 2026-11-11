from pathlib import Path

import pytest
from scripts.prepare_mac_reconstruction_smoke import prepare as prepare_smoke
from scripts.prepare_rfd3_local_projection import prepare

from dvbfixer.model.diffusion_cli import build_cli_diffusion_request

FIXTURE = Path(__file__).parent / "fixtures/protein_context/1UBQ.pdb"


def request(tmp_path):
    inputs = tmp_path / "input"
    manifest = prepare_smoke(FIXTURE, inputs)
    source = "".join(line for line in (inputs / "source.pdb").read_text().splitlines(True)
                     if not (line.startswith("ATOM  ") and 20 <= int(line[22:26]) <= 21))
    path = tmp_path / "source.pdb"
    path.write_text(source)
    parsed = build_cli_diffusion_request(path, {"A": manifest["sequence"]}, seeds=(7,), profile="protpardelle-1c-cc94-mps-v1")
    request_path = tmp_path / "request.json"
    request_path.write_text(parsed.to_json())
    return request_path, path, parsed


def test_local_joint_projection_uses_established_sequence_and_zero_missing_coords(tmp_path):
    req, source, parsed = request(tmp_path)
    out = tmp_path / "projection"
    manifest = prepare(req, source, (0, 1), out)
    assert manifest["generated_atoms"] == len(parsed.generated_atoms)
    fixed = {line[21:27]+line[12:16]:line for line in source.read_text().splitlines() if line.startswith("ATOM  ")}
    rows = [line for line in (out / "conditioning.pdb").read_text().splitlines() if line.startswith("ATOM  ")]
    assert [int(line[6:11]) for line in rows] == list(range(1, len(rows)+1))
    for line in rows:
        number = int(line[22:26])
        if 10 <= number <= 13 or 20 <= number <= 21:
            assert [float(line[i:i+8]) for i in (30, 38, 46)] == [0., 0., 0.]
        else:
            assert line[30:54] == fixed[line[21:27]+line[12:16]][30:54]


def test_projection_never_silently_discards_an_unselected_internal_gap(tmp_path):
    req, source, _ = request(tmp_path)
    with pytest.raises(ValueError, match="unselected gap"):
        prepare(req, source, (0,), tmp_path / "projection")


def test_projection_binds_original_source_bytes(tmp_path):
    req, source, _ = request(tmp_path)
    source.write_text(source.read_text()+"REMARK altered\n")
    with pytest.raises(ValueError, match="digest"):
        prepare(req, source, (0, 1), tmp_path / "projection")
