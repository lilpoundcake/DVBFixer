"""Renumbering custody control with reference-derived fake predictions, NOT sampling."""
import gzip
import json
from pathlib import Path

import gemmi
import numpy as np
import pytest
from scripts.materialize_rfd3_local_controls import native_group
from scripts.prepare_mac_reconstruction_smoke import prepare as prepare_smoke
from scripts.prepare_rfd3_local_projection import prepare

from dvbfixer.model.diffusion_cli import build_cli_diffusion_request

FIXTURE = Path(__file__).parent / "fixtures/protein_context/1UBQ.pdb"


def test_local_native_identity_map_restores_non_one_based_source_numbers(tmp_path):
    original = tmp_path / "original"
    manifest = prepare_smoke(FIXTURE, original)
    source = tmp_path / "source.pdb"
    source.write_bytes((original / "source.pdb").read_bytes())
    request = build_cli_diffusion_request(source, {"A": manifest["sequence"]}, seeds=(7,), profile="protpardelle-1c-cc94-mps-v1")
    request_path = tmp_path / "request.json"
    request_path.write_text(request.to_json())
    inputs = tmp_path / "inputs"
    projection = prepare(request_path, source, (0,), inputs)
    residues = sorted({int(a["atom"]["residue_number"]) for a in projection["axis"]})
    numbering = {number:i+1 for i, number in enumerate(residues)}
    rows = [line for line in (original / "reference-projection-NOT-CONDITIONING.pdb").read_text().splitlines(True)
            if line.startswith("ATOM  ") and int(line[22:26]) in numbering]
    native = "".join(line[:22]+f"{numbering[int(line[22:26])]:4d}"+line[26:] for line in rows)
    predictions = tmp_path / "run/predictions"
    predictions.mkdir(parents=True)
    structure = gemmi.read_pdb_string(native)
    with gzip.open(predictions / "prediction.cif.gz", "wt") as stream:
        stream.write(structure.make_mmcif_document().as_string())
    (predictions / "prediction.json").write_text(json.dumps({"diffused_index_map": {f"A{k}":f"A{v}" for k,v in numbering.items()}}))
    coordinates, receipt = native_group(inputs, tmp_path / "run")
    assert set(coordinates) == set(request.generated_atoms)
    expected = {line[21:27]+line[12:16]:np.array([float(line[i:i+8]) for i in (30,38,46)]) for line in rows}
    for atom, xyz in coordinates.items():
        key = f"{atom.chain}{int(atom.residue_number):4d} "+f" {atom.atom_name:<3}"
        assert np.allclose(xyz, expected[key])
    assert receipt["native_observed_drift_max_angstrom"] < 1e-10
    broken = {f"A{k}":f"A{v}" for k,v in numbering.items()}
    broken.pop(f"A{residues[0]}")
    (predictions / "prediction.json").write_text(json.dumps({"diffused_index_map": broken}))
    with pytest.raises(ValueError, match="absent/truncated"):
        native_group(inputs, tmp_path / "run")
