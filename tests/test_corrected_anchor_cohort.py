"""Provenance and shared-minimization controls for the frozen replay script."""

from pathlib import Path
from types import SimpleNamespace

import pytest

pytest.importorskip("openmm")
pytest.importorskip("pdbfixer")

from scripts import run_corrected_anchor_cohort as cohort

from dvbfixer.model.diffusion.anchor_acceptance import _structure
from dvbfixer.model.diffusion.anchor_policy import AnchorMovementLimits
from dvbfixer.model.diffusion_cli import _target_sequences, build_cli_diffusion_request

FIXTURE = Path(__file__).parent / "fixtures/protein_context/1UBQ.pdb"


def test_upstream_provenance_matches_observed_coordinates_and_b_values():
    original = FIXTURE.read_text()
    normalized = "".join(line for line in original.splitlines(keepends=True) if line.startswith("ATOM  "))
    evidence = cohort.verify_experimental_origin(normalized, original)
    assert evidence["method_source"] == "tracked-PDB-header"
    assert evidence["experimental_methods"] == ["X-RAY DIFFRACTION"]
    modified = "".join(line[:60] + "999.00" + line[66:] for line in normalized.splitlines(keepends=True))
    with pytest.raises(ValueError, match="do not match"):
        cohort.verify_experimental_origin(modified, original)


def test_stripped_headers_require_matching_external_experimental_metadata():
    original = "".join(line for line in FIXTURE.read_text().splitlines(keepends=True) if line.startswith("ATOM  "))
    with pytest.raises(ValueError, match="method unavailable"):
        cohort.verify_experimental_origin(original, original)
    metadata = {"entry": {"id": "1UBQ"}, "exptl": [{"method": "X-RAY DIFFRACTION"}]}
    evidence = cohort.verify_experimental_origin(original, original, entry_metadata=metadata, expected_entry="1UBQ")
    assert evidence["method_source"] == "RCSB-entry-metadata"
    with pytest.raises(ValueError, match="wrong entry identity"):
        cohort.verify_experimental_origin(original, original, entry_metadata=metadata, expected_entry="9CQD")
    metadata["exptl"] = [{"method": "ELECTRON MICROSCOPY"}]
    with pytest.raises(ValueError, match="method unavailable"):
        cohort.verify_experimental_origin(original, original, entry_metadata=metadata, expected_entry="1UBQ")


def test_cohort_does_not_treat_shared_stage_one_as_independent_minimizations(tmp_path, monkeypatch):
    original = FIXTURE.read_text()
    protein = "".join(line for line in original.splitlines(keepends=True) if not line.startswith(("HETATM", "CONECT")))
    source = "".join(line for line in protein.splitlines(keepends=True) if not
                      (line.startswith("ATOM  ") and line[21] == "A" and line[22:26].strip() in {"6", "7", "8"}))
    source_path = tmp_path / "source.pdb"
    source_path.write_text(source)
    request = build_cli_diffusion_request(source_path, _target_sequences(source_path, None), seeds=(7,),
                                         profile="protpardelle-1c-mps-local-v1")
    raw = "".join(line[:30] + f"{float(line[30:38]) + 1.5:8.3f}" + line[38:]
                   if line.startswith("ATOM  ") and line[21] == "A" and line[22:26].strip() == "8"
                   and line[12:16].strip() == "C" else line for line in protein.splitlines(keepends=True))
    good = _structure(protein, "native.pdb")
    calls = []
    def native_minimize(text, **kwargs):
        calls.append((text, kwargs))
        assert not kwargs["observed_anchor_restraints"]
        return SimpleNamespace(coordinates_angstrom={atom: good.atoms[atom].coordinates_angstrom
                                                     for atom in request.generated_atoms})
    monkeypatch.setattr(cohort, "refine_generated_region", native_minimize)
    loaded = request, source, raw, cohort.verify_experimental_origin(source, original)
    report = cohort.run_case(("synthetic-control", 7), loaded, tmp_path / "run", AnchorMovementLimits(0.5, 1))
    assert len(calls) == report["unique_native_minimizations"] == 1
    assert len(report["policies"]) == 3
    assert all(row["selected_stage"] == "generated-only" and row["research_pass"] for row in report["policies"])
    assert [row["minimizer_calls"][0]["reused_physical_minimization"] for row in report["policies"]] == [False, True, True]
