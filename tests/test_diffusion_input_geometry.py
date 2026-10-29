"""Bad observed geometry must never become immutable neural conditioning."""

import hashlib
import json
from dataclasses import replace
from pathlib import Path

import pytest
from Bio import SeqIO

from dvbfixer.diagnose.report import Finding, Severity
from dvbfixer.model.diffusion.input_geometry import (
    InputGeometryError,
    assess_input_geometry,
    require_input_geometry,
)
from dvbfixer.model.diffusion_cli import build_cli_diffusion_request


def _input(tmp_path):
    source = Path(__file__).parent / "fixtures/protein_context/1UBQ.pdb"
    sequence = str(next(SeqIO.parse(source, "pdb-seqres")).seq)
    path = tmp_path / "gap.pdb"
    path.write_text("".join(line for line in source.read_text().splitlines(keepends=True)
                            if not line.startswith("HETATM") and not
                            (line.startswith("ATOM  ") and 10 <= int(line[22:26]) <= 13)))
    request = build_cli_diffusion_request(path, {"A": sequence}, seeds=(7,), profile="protenix-v1-cuda")
    return path, request


def test_real_gap_is_not_a_stretched_observed_peptide(tmp_path):
    path, request = _input(tmp_path)
    before = path.read_bytes()
    findings = assess_input_geometry(request, before)
    assert not [f for f in findings if f.severity is Severity.ERROR]
    assert path.read_bytes() == before


def test_damaged_sidechain_blocks_sampling_and_retains_report_without_mutation(tmp_path, capsys):
    path, request = _input(tmp_path)
    lines = path.read_text().splitlines(keepends=True)
    for i, line in enumerate(lines):
        if line.startswith("ATOM  ") and line[22:26].strip() == "20" and line[12:16].strip() == "CB":
            lines[i] = line[:30] + f"{float(line[30:38]) + 4:8.3f}" + line[38:]
    damaged = "".join(lines).encode()
    request = replace(request, normalized_pdb=replace(request.normalized_pdb, sha256=hashlib.sha256(damaged).hexdigest()))
    with pytest.raises(InputGeometryError, match="Sampling was not started"):
        require_input_geometry(request, damaged, tmp_path)
    report = json.loads(next(tmp_path.glob("dvbfixer-input-geometry-*/report.json")).read_text())
    assert report["passed"] is False
    assert report["coordinates_modified"] is False
    assert any(f["category"] == "bond_length" and f["severity"] == "ERROR" for f in report["findings"])
    assert path.read_bytes() != damaged  # Original file was never touched.
    output = capsys.readouterr().out
    assert "ERROR: Diffusion input geometry:" in output
    assert "Input coordinates were not modified" in output
    assert "[bond_length]" in output
    assert "Full input geometry report:" in output


def test_input_audit_requires_exact_source_digest(tmp_path):
    path, request = _input(tmp_path)
    with pytest.raises(InputGeometryError, match="digest mismatch"):
        assess_input_geometry(request, path.read_bytes() + b"\n")


def test_warning_only_input_notifies_user_and_retains_complete_report(tmp_path, monkeypatch, capsys):
    from dvbfixer.model.diffusion import input_geometry

    path, request = _input(tmp_path)
    warning = Finding(Severity.WARNING, "bond_angle", "A", "20", "SER", "suspect observed angle")
    monkeypatch.setattr(input_geometry, "assess_input_geometry", lambda *_args: (warning,))
    require_input_geometry(request, path.read_bytes(), tmp_path)
    output = capsys.readouterr().out
    assert "WARNING: Diffusion input geometry: 0 error(s), 1 warning(s)" in output
    assert "A:20 SER [bond_angle]: suspect observed angle" in output
    report = json.loads(next(tmp_path.glob("dvbfixer-input-geometry-*/report.json")).read_text())
    assert report["passed"] is True
    assert len(report["findings"]) == 1


def test_clean_input_has_no_warning_banner_or_report_directory(tmp_path, monkeypatch, capsys):
    from dvbfixer.model.diffusion import input_geometry

    path, request = _input(tmp_path)
    monkeypatch.setattr(input_geometry, "assess_input_geometry", lambda *_args: ())
    require_input_geometry(request, path.read_bytes(), tmp_path)
    assert capsys.readouterr().out == ""
    assert not list(tmp_path.glob("dvbfixer-input-geometry-*"))


def test_pipeline_never_launches_runner_for_bad_observed_geometry(tmp_path, monkeypatch):
    from dvbfixer.model.diffusion import pipeline

    path, request = _input(tmp_path)
    original = path.read_bytes()
    lines = path.read_text().splitlines(keepends=True)
    for i, line in enumerate(lines):
        if line.startswith("ATOM  ") and line[22:26].strip() == "20" and line[12:16].strip() == "OG":
            lines[i] = line[:30] + f"{float(line[30:38]) + 4:8.3f}" + line[38:]
    damaged = "".join(lines).encode()
    bad_path = tmp_path / "damaged.pdb"
    bad_path.write_bytes(damaged)
    request = replace(request, normalized_pdb=replace(request.normalized_pdb, path=bad_path.name,
                                                     sha256=hashlib.sha256(damaged).hexdigest()))
    monkeypatch.setattr(pipeline, "run_diffusion_runner", lambda *_args, **_kwargs: pytest.fail("bad geometry reached runner"))
    with pytest.raises(pipeline.DiffusionPipelineError, match="Sampling was not started"):
        pipeline.run_diffusion_pipeline(request, ("runner",), source_root=tmp_path, work_parent=tmp_path,
                                        destination_bundle=tmp_path / "bundle")
    assert not (tmp_path / "bundle").exists()
    assert not list(tmp_path.glob(".dvbfixer-diffusion-run.*"))
    assert path.read_bytes() == original
    assert next(tmp_path.glob("dvbfixer-input-geometry-*/source.pdb")).read_bytes() == damaged
