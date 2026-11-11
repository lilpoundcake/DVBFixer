import hashlib
import json
from pathlib import Path

import pytest
from scripts.prepare_reference_reconstruction_case import prepare

FIXTURE = Path(__file__).parent / "fixtures/protein_context/1UBQ.pdb"


def test_reference_case_is_new_and_does_not_feed_hidden_coordinates(tmp_path):
    before = hashlib.sha256(FIXTURE.read_bytes()).hexdigest()
    output = tmp_path / "case"
    result = prepare(FIXTURE, output, chain="A", first=1, last=72, hidden_first=10, hidden_last=13, case_id="control")
    assert result["source_admissible"] and result["reference_baseline_accepted"]
    assert result["generated_atoms"] == 28 and result["observed_atoms"] == 546
    assert hashlib.sha256(FIXTURE.read_bytes()).hexdigest() == before
    engine = (output / "rfd3/rfd3-input.json").read_text()
    assert "NOT-CONDITIONING" not in engine
    for line in (output / "rfd3/conditioning.pdb").read_text().splitlines():
        if line.startswith("ATOM  ") and 10 <= int(line[22:26]) <= 13:
            assert [float(line[i:i+8]) for i in (30,38,46)] == [0., 0., 0.]
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["files"]["source.pdb"] == hashlib.sha256((output / "source.pdb").read_bytes()).hexdigest()
    with pytest.raises(FileExistsError):
        prepare(FIXTURE, output, chain="A", first=1, last=72, hidden_first=10, hidden_last=13, case_id="control")


def test_reference_projection_never_normalizes_partial_occupancies(tmp_path):
    with pytest.raises(ValueError, match="partial occupancy"):
        prepare(FIXTURE, tmp_path / "case", chain="A", first=1, last=76, hidden_first=10, hidden_last=13, case_id="refused")
