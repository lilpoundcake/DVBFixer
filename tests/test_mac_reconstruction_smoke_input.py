import json
from pathlib import Path

import pytest
from scripts.prepare_mac_reconstruction_smoke import prepare

FIXTURE = Path(__file__).parent / "fixtures/protein_context/1UBQ.pdb"


def test_mac_smoke_hides_reference_and_preserves_sequence_and_observations(tmp_path):
    out = tmp_path / "input"
    report = prepare(FIXTURE, out)
    assert report["observed_atoms"] == 546 and report["generated_atoms"] == 28
    assert len(report["sequence"]) == 72
    spec = json.loads((out / "rfd3-input.json").read_text())["ubq72_gap10_13"]
    assert spec["select_fixed_atoms"]["A10-13"] == ""
    assert spec["select_unfixed_sequence"] is False
    fixed = {(line[21:27], line[12:16]): line for line in (out / "source.pdb").read_text().splitlines()
             if line.startswith("ATOM  ")}
    original = {(line[21:27], line[12:16]): line for line in FIXTURE.read_text().splitlines()
                if line.startswith("ATOM  ")}
    rows = [line for line in (out / "rfd3-sequence-dummy.pdb").read_text().splitlines() if line.startswith("ATOM  ")]
    assert [int(line[6:11]) for line in rows] == list(range(1, 575))
    for line in rows:
        key = line[21:27], line[12:16]
        if 10 <= int(line[22:26]) <= 13:
            assert key not in fixed
            assert [float(line[i:i + 8]) for i in (30, 38, 46)] == [0, 0, 0]
        else:
            assert line[30:54] == fixed[key][30:54] == original[key][30:54]


def test_mac_smoke_does_not_rewrite_existing_input(tmp_path):
    out = tmp_path / "input"
    prepare(FIXTURE, out)
    before = (out / "manifest.json").read_bytes()
    with pytest.raises(FileExistsError):
        prepare(FIXTURE, out)
    assert (out / "manifest.json").read_bytes() == before
