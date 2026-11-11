import json

import pytest
from scripts.audit_mac_reconstruction_smoke import inventory
from scripts.mac_validation_inputs import prepare, transport
from scripts.mac_validation_policy import SEEDS, require_policy
from scripts.prepare_mac_comparison_expansion import prepare as prepare_case


@pytest.fixture
def validation_case(tmp_path):
    original, output = tmp_path / "original", tmp_path / "case"
    spec = {"id": "independent-test-control", "fixture": "tests/fixtures/protein_context/2CI2.pdb",
        "spans": [["I", 19, 83]], "hidden": [["I", 28, 31]]}
    prepare_case(spec, original, allow_observed_terminal_oxt=True)
    prepare(original, output)
    return output


def test_validation_preparation_preserves_observed_oxt_outside_native_canonical_axis(validation_case):
    case = validation_case
    receipt = json.loads((case / "transport.json").read_text())
    assert receipt["fixed_observed_terminal_passthrough"] == [{"chain": "I", "residue_number": "83", "atom_name": "OXT"}]
    assert any(a.atom_name == "OXT" for a in inventory(case / "source.pdb"))
    assert not any(a.atom_name == "OXT" for a in inventory(case / "rfd3/conditioning.pdb"))
    assert receipt["reference_coordinates_encoded"] is False


def test_unique_native_chain_ordinal_transport_retains_original_chain_case_and_oxt(validation_case, tmp_path):
    case = validation_case
    reference = (case / "reference-projection-NOT-CONDITIONING.pdb").read_text()
    lines = []
    for line in reference.splitlines(keepends=True):
        if not line.startswith("ATOM  ") or line[12:16].strip() == "OXT":
            continue
        lines.append(line[:21]+"B"+f"{int(line[22:26])-18:4d}"+line[26:])
    prediction = tmp_path / "synthetic-native.pdb"
    prediction.write_text("".join(lines)+"TER\nEND\n")
    output = tmp_path / "transport"
    report = transport(case, prediction, output, model="boltz")
    assert report["fixed_atoms_passed_through_only"] == 1
    assert report["custody"]["native_chain"] == "B"
    assert report["custody"]["original_chain_preserved_in_output"] == "I"
    assert {a.chain for a in inventory(output / "source-restored.pdb")} == {"I"}
    observed_oxt = next(line for line in (case / "source.pdb").read_text().splitlines() if line.startswith("ATOM  ") and line[12:16].strip() == "OXT")
    assert observed_oxt in (output / "source-restored.pdb").read_text().splitlines()


def test_reserved_seed_admission_requires_hash_and_exact_frozen_scope(tmp_path, monkeypatch):
    import hashlib
    policy = {"protocol": "private-fixed-source-validation-v1", "seeds": SEEDS, "models": ["rfd3"],
        "observed_motion": "none", "automatic_candidate_selection": False, "public_capability_enabled": False,
        "validation_seed_sampling_authorized": True}
    path = tmp_path / "policy.json"
    path.write_text(json.dumps(policy))
    monkeypatch.setenv("DVBFIXER_VALIDATION_POLICY", str(path))
    monkeypatch.setenv("DVBFIXER_VALIDATION_POLICY_SHA256", "0"*64)
    with pytest.raises(ValueError, match="digest"):
        require_policy("rfd3", SEEDS[0])
    monkeypatch.setenv("DVBFIXER_VALIDATION_POLICY_SHA256", hashlib.sha256(path.read_bytes()).hexdigest())
    assert require_policy("rfd3", SEEDS[0])["seed"] == SEEDS[0]
    with pytest.raises(ValueError, match="does not admit"):
        require_policy("rfd3", 7)
    with pytest.raises(ValueError, match="does not admit"):
        require_policy("protenix", SEEDS[0])
