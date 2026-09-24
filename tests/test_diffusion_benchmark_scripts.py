"""Tests for the internal diffusion benchmark reproduction scripts."""

from __future__ import annotations

import hashlib
import importlib.util
from pathlib import Path
from types import ModuleType

from dvbfixer.model.diffusion.contract import (
    DIFFUSION_SCHEMA_VERSION,
    ArtifactReference,
    BackendProvenance,
    DiffusionRequest,
    DiffusionStatus,
    RunnerCandidate,
    RunnerDiagnostics,
    RunnerResult,
)
from dvbfixer.model.diffusion.rfdiffusion_v1 import (
    build_synthetic_input,
    materialize_candidate,
)
from dvbfixer.model.diffusion.runner import DIFFUSION_RUNNER_PROTOCOL_VERSION
from dvbfixer.model.diffusion.scope import assess_diffusion_scope

ROOT = Path(__file__).resolve().parents[1]


def _load_script(name: str) -> ModuleType:
    path = ROOT / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_native_result(workspace: Path) -> None:
    request = DiffusionRequest.from_json((workspace / "request.json").read_text())
    candidate_path = workspace / "reference.pdb"
    candidate = RunnerCandidate(
        candidate_id="native-reference",
        seed=request.seeds[0],
        coordinate_artifact=ArtifactReference(
            "reference.pdb",
            hashlib.sha256(candidate_path.read_bytes()).hexdigest(),
        ),
        generated_atoms=request.generated_atoms,
        generated_residues=request.gaps[0].generated_residues,
        raw_backend_score=None,
        score_provenance="test-native-reference",
    )
    result = RunnerResult(
        schema_version=DIFFUSION_SCHEMA_VERSION,
        status=DiffusionStatus.SUCCESS,
        candidates=(candidate,),
        runner_diagnostics=RunnerDiagnostics(exit_code=0, timed_out=False),
        backend_provenance=BackendProvenance(
            backend="test",
            runner_protocol_version=DIFFUSION_RUNNER_PROTOCOL_VERSION,
            engine_repository="https://example.invalid/test",
            engine_revision="test",
        ),
    )
    (workspace / "result.json").write_text(result.to_json(), encoding="utf-8")


def _ordinal_backbone(reference: Path, chain: str) -> str:
    """Render reference backbone atoms in RFdiffusion's ordinal identity space."""
    residue_ordinals: dict[tuple[str, str], int] = {}
    output: list[str] = []
    serial = 1
    for line in reference.read_text().splitlines():
        if not line.startswith("ATOM  ") or line[21] != chain:
            continue
        identity = (line[22:26].strip(), line[26].strip())
        ordinal = residue_ordinals.setdefault(identity, len(residue_ordinals) + 1)
        if line[12:16].strip() not in {"N", "CA", "C", "O"}:
            continue
        padded = line.ljust(80)
        output.append(
            padded[:6]
            + f"{serial:5d}"
            + padded[11:21]
            + "A"
            + f"{ordinal:4d} "
            + padded[27:80]
            + "\n"
        )
        serial += 1
    output.append("END\n")
    return "".join(output)


def _ordinal_residues(reference: Path, chain: str) -> str:
    residue_ordinals: dict[tuple[str, str], int] = {}
    output: list[str] = []
    for line in reference.read_text().splitlines(keepends=True):
        if not line.startswith("ATOM  ") or line[21] != chain:
            output.append(line)
            continue
        identity = (line[22:26].strip(), line[26].strip())
        ordinal = residue_ordinals.setdefault(identity, len(residue_ordinals) + 1)
        output.append(line[:22] + f"{ordinal:4d} " + line[27:])
    return "".join(output)


def test_builder_and_analyzers_reproduce_native_reference(tmp_path: Path) -> None:
    builder = _load_script("build_diffusion_benchmark_request")
    pair_analyzer = _load_script("analyze_diffusion_benchmark_pair")
    modeller_analyzer = _load_script("analyze_modeller_benchmark")
    first = tmp_path / "first"
    second = tmp_path / "second"

    builder.build_workspace("8b01-chain-c-withheld-5", first, seed=7)
    builder.build_workspace("8b01-chain-c-withheld-5", second, seed=7)
    _write_native_result(first)
    _write_native_result(second)

    request = DiffusionRequest.from_json((first / "request.json").read_text())
    assert request.target_sequences[0].sequence[58:63] == "DNAKN"
    assert request.gaps[0].target_interval.start == 58
    assert request.gaps[0].target_interval.stop == 63

    pair = pair_analyzer.analyze(first, second)
    assert pair["validation_passed"] is True
    assert pair["quality"]["gap_backbone_rmsd_angstrom"] == 0.0
    assert pair["repeatability"]["classification"] == "deterministic"

    modeller = modeller_analyzer.analyze(
        first,
        (first / "reference.pdb", first / "reference.pdb"),
    )
    assert modeller["validation_pass_count"] == 2
    assert modeller["median_gap_backbone_rmsd_angstrom"] == 0.0
    assert modeller["median_fixed_heavy_rmsd_angstrom"] == 0.0


def test_builder_extracts_declared_target_chain_scope(tmp_path: Path) -> None:
    builder = _load_script("build_diffusion_benchmark_request")
    workspace = tmp_path / "glypro"

    builder.build_workspace("7x35-chain-a-glypro-7", workspace, seed=7)

    request = DiffusionRequest.from_json((workspace / "request.json").read_text())
    source = (workspace / "input/normalized.pdb").read_bytes()
    reference_lines = (workspace / "reference.pdb").read_text().splitlines()
    assert assess_diffusion_scope(request, source).supported is True
    assert request.target_sequences[0].sequence[133:140] == "PPGGPVP"
    assert all(
        not line.startswith(("ATOM  ", "HETATM")) or line[21] == "A"
        for line in reference_lines
    )
    assert not any(line.startswith("HETATM") for line in reference_lines)


def test_builder_extracts_target_and_partner_protein_scope(tmp_path: Path) -> None:
    builder = _load_script("build_diffusion_benchmark_request")
    workspace = tmp_path / "interface"

    builder.build_workspace("7x35-chain-a-interface-5", workspace, seed=7)

    request = DiffusionRequest.from_json((workspace / "request.json").read_text())
    source = (workspace / "input/normalized.pdb").read_bytes()
    reference_lines = (workspace / "reference.pdb").read_text().splitlines()
    fasta = (workspace / "target.fasta").read_text()
    assert assess_diffusion_scope(request, source).supported is True
    assert request.target_sequences[0].sequence[232:237] == "AWVPR"
    assert {
        line[21]
        for line in reference_lines
        if line.startswith(("ATOM  ", "HETATM"))
    } == {"A", "B"}
    assert not any(line.startswith("HETATM") for line in reference_lines)
    assert ">chain_A\n" in fasta
    assert ">chain_B\n" in fasta


def test_insertion_code_benchmark_restores_generated_identities(tmp_path: Path) -> None:
    builder = _load_script("build_diffusion_benchmark_request")
    modeller_analyzer = _load_script("analyze_modeller_benchmark")
    workspace = tmp_path / "insertion-codes"

    builder.build_workspace("7k8s-chain-h-insertion-codes-3", workspace, seed=7)

    request = DiffusionRequest.from_json((workspace / "request.json").read_text())
    source = (workspace / "input/normalized.pdb").read_text()
    generated = request.gaps[0].generated_residues
    assert assess_diffusion_scope(request, source.encode()).supported is True
    assert request.target_sequences[0].sequence[82:85] == "NSL"
    assert [(item.residue_number, item.insertion_code) for item in generated] == [
        ("82", "A"),
        ("82", "B"),
        ("82", "C"),
    ]

    synthetic = build_synthetic_input(request, source)
    assert synthetic.contig == "A1-82/3-3/A86-125"
    assert all(
        line[26] == " "
        for line in synthetic.pdb_text.splitlines()
        if line.startswith("ATOM  ")
    )

    candidate = materialize_candidate(
        request,
        source_text=source,
        synthetic_text=synthetic.pdb_text,
        raw_text=_ordinal_backbone(workspace / "reference.pdb", "H"),
    )
    represented = {
        (line[22:26].strip(), line[26].strip())
        for line in candidate.splitlines()
        if line.startswith("ATOM  ") and line[21] == "H"
    }
    assert {("82", ""), ("82", "A"), ("82", "B"), ("82", "C"), ("83", "")} <= represented

    modeller_path = workspace / "ordinal-modeller.pdb"
    modeller_path.write_text(_ordinal_residues(workspace / "reference.pdb", "H"))
    modeller = modeller_analyzer.analyze(workspace, (modeller_path,))
    assert modeller["validation_pass_count"] == 1
    assert modeller["median_gap_backbone_rmsd_angstrom"] == 0.0
    assert modeller["median_fixed_heavy_rmsd_angstrom"] == 0.0
    assert modeller["candidates"][0]["identity_normalized"] is True
