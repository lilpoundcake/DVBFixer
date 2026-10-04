"""Tests for the narrow ``model --backend diffusion`` CLI slice."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from dvbfixer.model.cli import parse_args
from dvbfixer.model.diffusion.preflight import (
    AdapterPreflightCode,
    AdapterPreflightIssue,
    RunnerPreflightFacts,
    RunnerPreflightReport,
)
from dvbfixer.model.diffusion.scope import assess_diffusion_scope
from dvbfixer.model.diffusion_cli import (
    DiffusionCliError,
    build_cli_diffusion_request,
    run_diffusion_model,
)

FIXTURE = Path(__file__).parent / "fixtures" / "8cz8" / "8cz8_a_u.pdb"


def _gap_input(tmp_path: Path) -> Path:
    generated = {str(number) for number in range(65, 70)}
    selected = []
    for line in FIXTURE.read_text(encoding="utf-8").splitlines(keepends=True):
        if not line.startswith("ATOM  ") or line[21] != "C":
            continue
        number = line[22:26].strip()
        if number not in {str(value) for value in range(61, 74)} or number in generated:
            continue
        selected.append(line)
    path = tmp_path / "gap.pdb"
    path.write_text("".join(selected) + "TER\nEND\n", encoding="utf-8")
    return path


def test_diffusion_parser_requires_explicit_runtime_inputs(tmp_path: Path) -> None:
    checkpoint = tmp_path / "model.pt"
    checkpoint.write_bytes(b"checkpoint")
    args = parse_args(
        [
            "input.pdb",
            "--backend",
            "diffusion",
            "--diffusion-profile",
            "protpardelle-1c-mps",
            "--diffusion-runner",
            "runner",
            "--diffusion-checkpoint",
            str(checkpoint),
            "--diffusion-seed",
            "7",
            "--diffusion-seed",
            "11",
        ]
    )

    assert args.backend == "diffusion"
    assert args.diffusion_seeds == [7, 11]


def test_diffusion_parser_rejects_modeller_only_options() -> None:
    with pytest.raises(SystemExit):
        parse_args(
            [
                "input.pdb",
                "--backend",
                "diffusion",
                "--diffusion-profile",
                "protpardelle-1c-mps",
                "--diffusion-runner",
                "runner",
                "--diffusion-checkpoint",
                "model.pt",
                "--num-loops",
                "3",
            ]
        )


@pytest.mark.parametrize("timeout", ["nan", "inf", "0", "-1"])
def test_diffusion_parser_rejects_invalid_timeout(timeout: str) -> None:
    with pytest.raises(SystemExit):
        parse_args(
            [
                "input.pdb",
                "--backend",
                "diffusion",
                "--diffusion-profile",
                "protenix-v1-cuda",
                "--diffusion-runner",
                "runner",
                "--diffusion-checkpoint",
                "model.pt",
                "--diffusion-timeout",
                timeout,
            ]
        )


def test_request_builder_produces_admitted_one_gap_request(tmp_path: Path) -> None:
    input_path = _gap_input(tmp_path)
    request = build_cli_diffusion_request(
        input_path,
        {"C": "SNRFSGSKSGNTA"},
        seeds=(7,),
        profile="protpardelle-1c-mps",
    )

    assert request.gaps[0].target_interval.start == 4
    assert request.gaps[0].target_interval.stop == 9
    assert [residue.residue_number for residue in request.gaps[0].generated_residues] == [
        "65", "66", "67", "68", "69",
    ]
    assert request.candidate_count == 1
    assert assess_diffusion_scope(request, input_path.read_bytes()).supported


def test_failed_profile_preflight_blocks_inference(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_path = _gap_input(tmp_path)
    fasta = tmp_path / "target.fasta"
    fasta.write_text(">C\nSNRFSGSKSGNTA\n")
    checkpoint = tmp_path / "checkpoint.pt"
    checkpoint.write_bytes(b"checkpoint")
    report = RunnerPreflightReport(
        profile="protenix-v1-cuda",
        runner_protocol_version=4,
        contract_schema_version=4,
        facts=RunnerPreflightFacts(),
        issues=(
            AdapterPreflightIssue(
                AdapterPreflightCode.MISSING_CUDA,
                "required CUDA device is unavailable",
            ),
        ),
    )
    monkeypatch.setattr(
        "dvbfixer.model.diffusion_cli.invoke_runner_preflight",
        lambda **_kwargs: report,
    )
    monkeypatch.setattr(
        "dvbfixer.model.diffusion_cli.run_diffusion_pipeline",
        lambda *_args, **_kwargs: pytest.fail("inference must not run"),
    )
    args = SimpleNamespace(
        input=str(input_path),
        output=str(tmp_path / "output"),
        diffusion_checkpoint=str(checkpoint),
        diffusion_checkpoint_sha256=None,
        diffusion_profile="protenix-v1-cuda",
        diffusion_runner="runner",
        diffusion_timeout=1.0,
        diffusion_seeds=[7],
        diffusion_work_parent=None,
        fasta=str(fasta),
    )

    with pytest.raises(DiffusionCliError, match="missing-cuda"):
        run_diffusion_model(args)


def test_scope_rejection_precedes_runner_preflight(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_path = _gap_input(tmp_path)
    input_path.write_text(input_path.read_text().replace("ATOM  ", "HETATM", 1))
    fasta = tmp_path / "target.fasta"
    fasta.write_text(">C\nSNRFSGSKSGNTA\n")
    checkpoint = tmp_path / "checkpoint.pt"
    checkpoint.write_bytes(b"checkpoint")
    monkeypatch.setattr(
        "dvbfixer.model.diffusion_cli.invoke_runner_preflight",
        lambda **_kwargs: pytest.fail("preflight must follow scope admission"),
    )
    args = SimpleNamespace(
        input=str(input_path),
        output=str(tmp_path / "output"),
        diffusion_checkpoint=str(checkpoint),
        diffusion_checkpoint_sha256=None,
        diffusion_profile="protenix-v1-cuda",
        diffusion_runner="runner",
        diffusion_timeout=1.0,
        diffusion_seeds=[7],
        diffusion_work_parent=None,
        fasta=str(fasta),
    )

    with pytest.raises(DiffusionCliError):
        run_diffusion_model(args)
