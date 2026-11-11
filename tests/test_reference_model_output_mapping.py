import json
from pathlib import Path

import pytest
from scripts.audit_mac_reconstruction_smoke import audit
from scripts.prepare_reference_reconstruction_case import prepare
from scripts.remap_reference_model_output import run

FIXTURE = Path(__file__).parent / "fixtures/protein_context/1UBQ.pdb"


@pytest.mark.parametrize("model", ["rfd3","boltz"])
def test_native_ordinal_mapping_uses_source_numbering_not_reference_coordinates(tmp_path, model):
    case = tmp_path / "case"
    prepare(FIXTURE,case,chain="A",first=3,last=72,hidden_first=10,hidden_last=13,case_id="fake-output-control")
    rows = (case / "reference-projection-NOT-CONDITIONING.pdb").read_text().splitlines(True)
    prediction = tmp_path / "prediction.pdb"
    prediction.write_text("".join(line[:22]+f"{int(line[22:26])-2:4d}"+line[26:] if line.startswith("ATOM  ") else line for line in rows))
    (tmp_path / "prediction.json").write_text(json.dumps({"diffused_index_map":{f"A{i}":f"A{i-2}" for i in range(3,73)}}))
    mapped = run(case,prediction,tmp_path / "mapping",model=model)
    result = audit(case,mapped,tmp_path / "audit")
    assert result["heldout_generated_rmsd_angstrom"] < 1e-10
    assert result["assessments"]["source-restored"]["geometry_passed"]
    receipt = json.loads((tmp_path / "mapping/mapping.json").read_text())
    assert not receipt["reference_coordinates_read"] and not receipt["source_restored"]


def test_native_observed_drift_is_not_hidden_by_output_remapping(tmp_path):
    case = tmp_path / "case"
    prepare(FIXTURE,case,chain="A",first=3,last=72,hidden_first=10,hidden_last=13,case_id="fake-drift-control")
    rows = []
    for line in (case / "reference-projection-NOT-CONDITIONING.pdb").read_text().splitlines(True):
        if line.startswith("ATOM  "):
            number = int(line[22:26])
            if number == 20 and line[12:16].strip() == "CA":
                line = line[:30]+f"{float(line[30:38])+1:8.3f}"+line[38:]
            line = line[:22]+f"{number-2:4d}"+line[26:]
        rows.append(line)
    prediction = tmp_path / "prediction.pdb"
    prediction.write_text("".join(rows))
    run(case,prediction,tmp_path / "mapping",model="boltz")
    receipt = json.loads((tmp_path / "mapping/mapping.json").read_text())
    assert receipt["native_observed_max_angstrom"] > .5
    assert not receipt["source_restored"]
