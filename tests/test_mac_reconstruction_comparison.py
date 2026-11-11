from dataclasses import replace
from pathlib import Path

import pytest
from scripts.mac_reconstruction_comparison import (
    DEVELOPMENT_SEEDS,
    VALIDATION_SEEDS,
    modeller_alignment,
    normalize_modeller,
    target_axis,
)
from scripts.prepare_reference_reconstruction_case import prepare
from scripts.run_mac_modeller_proposal import run

from dvbfixer.model.diffusion.contract import DiffusionRequest

FIXTURE = Path(__file__).parent / "fixtures/protein_context/1UBQ.pdb"


@pytest.fixture
def case(tmp_path):
    output = tmp_path / "case"
    prepare(FIXTURE, output, chain="A", first=1, last=72, hidden_first=10, hidden_last=13, case_id="contract")
    return output, DiffusionRequest.from_json((output / "evaluation-mask.json").read_text())


def test_alignment_is_exact_and_reference_blind(case):
    output, request = case
    (output / "reference-projection-NOT-CONDITIONING.pdb").unlink()
    alignment = modeller_alignment(request)
    assert alignment.splitlines()[2].count("-") == 4
    assert alignment.splitlines()[2][9:13] == "----"
    assert alignment.splitlines()[5].rstrip("*") == request.target_sequences[0].sequence
    assert len(target_axis(request)["A"]) == 72
    assert not set(DEVELOPMENT_SEEDS) & set(VALIDATION_SEEDS)


def test_native_mapping_changes_only_identity_fields(case):
    output, request = case
    raw = (output / "reference-projection-NOT-CONDITIONING.pdb").read_text()
    native = "".join(line[:21]+"z"+line[22:]+"\n" for line in raw.splitlines() if line.startswith("ATOM  "))
    mapped, report = normalize_modeller(request, native)
    assert [line[30:] for line in mapped.splitlines() if line.startswith("ATOM  ")] == [
        line[30:] for line in native.splitlines() if line.startswith("ATOM  ")]
    assert all(line[21] == "A" for line in mapped.splitlines() if line.startswith("ATOM  "))
    assert report["source_restored"] is False and report["chirality_repaired"] is False


def test_native_mapping_refuses_missing_duplicate_or_wrong_atoms(case):
    output, request = case
    raw = (output / "reference-projection-NOT-CONDITIONING.pdb").read_text()
    atom = next(line for line in raw.splitlines(keepends=True) if line.startswith("ATOM  "))
    with pytest.raises(ValueError, match="missing requested"):
        normalize_modeller(request, raw.replace(atom, "", 1))
    with pytest.raises(ValueError, match="duplicate"):
        normalize_modeller(request, raw.replace(atom, atom+atom, 1))
    with pytest.raises(ValueError, match="target sequence"):
        normalize_modeller(request, raw.replace(atom, atom[:17]+"GLY"+atom[20:], 1))


def test_native_terminal_oxt_is_explicitly_excluded(case):
    output, request = case
    raw = (output / "reference-projection-NOT-CONDITIONING.pdb").read_text()
    atom = [line for line in raw.splitlines(keepends=True) if line.startswith("ATOM  ")][-1]
    extra = atom[:12]+" OXT"+atom[16:]
    _, report = normalize_modeller(request, raw.replace("TER\n", extra+"TER\n"))
    assert len(report["excluded_native_terminal_OXT"]) == 1


def test_axis_refuses_incomplete_placement(case):
    _, request = case
    placement = request.sequence_placements[0]
    incomplete = replace(placement, observed_target_indices=placement.observed_target_indices[1:],
        observed_residues=placement.observed_residues[1:])
    with pytest.raises(ValueError, match="every target ordinal"):
        target_axis(replace(request, sequence_placements=(incomplete,)))


def test_reserved_seed_refused_before_creating_outputs(tmp_path):
    with pytest.raises(ValueError, match="development seeds"):
        run(tmp_path / "missing-request", tmp_path / "missing-source", tmp_path / "output", 211)
    assert not (tmp_path / "output").exists()
