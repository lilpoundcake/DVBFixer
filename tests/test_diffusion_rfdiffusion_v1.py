"""CPU tests for the pinned RFdiffusion v1 benchmark adapter."""

from __future__ import annotations

import hashlib
import math
from dataclasses import replace
from io import StringIO
from pathlib import Path

import pytest
from openmm.app import PDBFile
from openmm.unit import angstrom

from dvbfixer.diagnose.report import Severity
from dvbfixer.diagnose.steric import clashes_python
from dvbfixer.model.diffusion.boundary_refinement import refine_generated_region
from dvbfixer.model.diffusion.contract import (
    DIFFUSION_SCHEMA_VERSION,
    ArtifactReference,
    AtomIdentity,
    DiffusionRequest,
    GapRegion,
    ResidueIdentity,
    SequencePlacement,
    TargetInterval,
    TargetSequence,
)
from dvbfixer.model.diffusion.rfdiffusion_v1 import (
    RFdiffusionAdapterError,
    build_synthetic_input,
    materialize_candidate,
    unsupported_request_reason,
)

FIXTURE = Path(__file__).parent / "fixtures" / "8cz8" / "8cz8_a_u.pdb"
GENERATED_NUMBERS = tuple(str(number) for number in range(65, 70))


def _residue(line: str) -> ResidueIdentity:
    return ResidueIdentity(line[21], line[22:26].strip(), line[26].strip())


def _atom(line: str) -> AtomIdentity:
    residue = _residue(line)
    return AtomIdentity(
        residue.chain,
        residue.residue_number,
        residue.insertion_code,
        line[12:16].strip(),
    )


def _heavy(line: str) -> bool:
    return (line[76:78].strip() or line[12:16].strip()[:1]).upper() != "H"


def _case() -> tuple[DiffusionRequest, str, str]:
    full_lines = FIXTURE.read_text(encoding="utf-8").splitlines(keepends=True)
    generated = {ResidueIdentity("C", number) for number in GENERATED_NUMBERS}
    source = "".join(
        line
        for line in full_lines
        if not (line.startswith(("ATOM  ", "HETATM")) and _residue(line) in generated)
    )
    generated_atoms = tuple(
        _atom(line)
        for line in full_lines
        if line.startswith(("ATOM  ", "HETATM"))
        and _residue(line) in generated
        and _heavy(line)
    )
    fixed_atoms = tuple(
        _atom(line)
        for line in source.splitlines()
        if line.startswith(("ATOM  ", "HETATM")) and _heavy(line)
    )
    source_bytes = source.encode("utf-8")
    left = ResidueIdentity("C", "64")
    right = ResidueIdentity("C", "70")
    generated_order = tuple(ResidueIdentity("C", number) for number in GENERATED_NUMBERS)
    request = DiffusionRequest(
        schema_version=DIFFUSION_SCHEMA_VERSION,
        normalized_pdb=ArtifactReference(
            "input/normalized.pdb", hashlib.sha256(source_bytes).hexdigest()
        ),
        target_sequences=(TargetSequence("C", "FSGSKSG"),),
        sequence_placements=(
            SequencePlacement(
                chain="C",
                target_length=7,
                observed_target_indices=(0, 6),
                observed_residues=(left, right),
            ),
        ),
        gaps=(
            GapRegion(
                chain="C",
                target_interval=TargetInterval(1, 6),
                left_anchor=left,
                right_anchor=right,
                generated_residues=generated_order,
                movable_junction_residues=(left, *generated_order, right),
            ),
        ),
        fixed_atoms=fixed_atoms,
        generated_atoms=generated_atoms,
        retained_explicit_links=(),
        candidate_count=1,
        seeds=(7,),
    )

    raw_lines: list[str] = []
    serial = 1
    target_number = 0
    for line in full_lines:
        if not line.startswith("ATOM  ") or line[21] != "C":
            continue
        residue_number = int(line[22:26])
        if not 64 <= residue_number <= 70 or line[12:16].strip() not in {"N", "CA", "C", "O"}:
            continue
        target_number = residue_number - 63
        padded = line.rstrip("\n").ljust(80)
        raw_lines.append(
            padded[:6]
            + f"{serial:5d}"
            + padded[11:17]
            + ("GLY" if 65 <= residue_number <= 69 else padded[17:20])
            + padded[20:21]
            + "A"
            + f"{target_number:4d} "
            + padded[27:]
            + "\n"
        )
        serial += 1
    raw_lines.append("END\n")
    return request, source, "".join(raw_lines)


def _case_with_partner() -> tuple[DiffusionRequest, str, str]:
    request, source, raw = _case()
    partner_lines: list[str] = []
    for line in FIXTURE.read_text(encoding="utf-8").splitlines(keepends=True):
        if (
            line.startswith("ATOM  ")
            and line[21] == "C"
            and line[22:26].strip() in {"71", "72"}
        ):
            partner_lines.append(line[:21] + "d" + line[22:])
    source_lines = source.splitlines(keepends=True)
    end_index = next(
        index for index, line in enumerate(source_lines) if line.startswith("END")
    )
    source_lines[end_index:end_index] = [*partner_lines, "TER\n"]
    source = "".join(source_lines)
    partner_atoms = tuple(_atom(line) for line in partner_lines if _heavy(line))
    request = replace(request, fixed_atoms=(*request.fixed_atoms, *partner_atoms))

    raw_lines = raw.splitlines(keepends=True)[:-1]
    serial = len(raw_lines) + 1
    for line in partner_lines:
        if line[12:16].strip() not in {"N", "CA", "C", "O"}:
            continue
        ordinal = int(line[22:26]) - 70
        padded = line.rstrip("\n").ljust(80)
        raw_lines.append(
            padded[:6]
            + f"{serial:5d}"
            + padded[11:21]
            + "B"
            + f"{request.sequence_placements[0].target_length + ordinal:4d} "
            + padded[27:]
            + "\n"
        )
        serial += 1
    raw_lines.append("END\n")
    return request, source, "".join(raw_lines)


def test_build_synthetic_input_uses_stable_target_ordinals() -> None:
    request, source, _raw = _case()

    synthetic = build_synthetic_input(request, source)

    assert synthetic.contig == "A1-1/5-5/A7-7"
    atom_lines = [line for line in synthetic.pdb_text.splitlines() if line.startswith("ATOM  ")]
    assert {line[21] for line in atom_lines} == {"A"}
    assert {line[22:26].strip() for line in atom_lines} == {"1", "7"}
    assert all(line[26] == " " for line in atom_lines)


def test_build_synthetic_input_includes_fixed_partner_as_receptor_context() -> None:
    request, source, _raw = _case_with_partner()

    synthetic = build_synthetic_input(request, source)

    assert synthetic.contig == "A1-1/5-5/A7-7/0 B1-2"
    assert len(synthetic.partner_chains) == 1
    partner = synthetic.partner_chains[0]
    assert (
        partner.source_chain,
        partner.synthetic_chain,
        partner.residue_count,
        partner.output_start,
    ) == ("d", "B", 2, 8)
    partner_lines = [
        line
        for line in synthetic.pdb_text.splitlines()
        if line.startswith("ATOM  ") and line[21] == "B"
    ]
    assert {line[22:26].strip() for line in partner_lines} == {"1", "2"}


def test_adapter_rejects_multiple_gaps_before_launch() -> None:
    request, _source, _raw = _case()
    request = DiffusionRequest(
        schema_version=request.schema_version,
        normalized_pdb=request.normalized_pdb,
        target_sequences=request.target_sequences,
        sequence_placements=request.sequence_placements,
        gaps=(request.gaps[0], request.gaps[0]),
        fixed_atoms=request.fixed_atoms,
        generated_atoms=request.generated_atoms,
        retained_explicit_links=(),
        candidate_count=1,
        seeds=(7,),
    )

    assert unsupported_request_reason(request) == (
        "RFdiffusion v1 benchmark adapter supports exactly one internal gap"
    )


def test_materialization_restores_identity_and_all_expected_heavy_atoms() -> None:
    request, source, raw = _case()
    synthetic = build_synthetic_input(request, source)

    candidate = materialize_candidate(
        request,
        source_text=source,
        synthetic_text=synthetic.pdb_text,
        raw_text=raw,
    )

    candidate_atoms = {
        _atom(line)
        for line in candidate.splitlines()
        if line.startswith(("ATOM  ", "HETATM")) and _heavy(line)
    }
    assert set(request.generated_atoms) <= candidate_atoms
    assert set(request.fixed_atoms) <= candidate_atoms
    assert not any(
        line.startswith("ATOM  ") and line[21] == "A"
        for line in candidate.splitlines()
    )


def test_materialization_uses_and_restores_partner_context() -> None:
    request, source, raw = _case_with_partner()
    synthetic = build_synthetic_input(request, source)

    candidate = materialize_candidate(
        request,
        source_text=source,
        synthetic_text=synthetic.pdb_text,
        raw_text=raw,
        partner_chains=synthetic.partner_chains,
    )

    source_partner_lines = {
        line
        for line in source.splitlines()
        if line.startswith("ATOM  ") and line[21] == "d"
    }
    assert source_partner_lines <= set(candidate.splitlines())


def test_materialization_rejects_missing_partner_context_output() -> None:
    request, source, raw = _case_with_partner()
    synthetic = build_synthetic_input(request, source)
    raw_without_partner = "".join(
        line
        for line in raw.splitlines(keepends=True)
        if not (line.startswith("ATOM  ") and line[21] == "B")
    )

    with pytest.raises(
        RFdiffusionAdapterError,
        match="missing partner-context backbone atom",
    ):
        materialize_candidate(
            request,
            source_text=source,
            synthetic_text=synthetic.pdb_text,
            raw_text=raw_without_partner,
            partner_chains=synthetic.partner_chains,
        )


def test_materialization_rejects_partner_context_drift() -> None:
    request, source, raw = _case_with_partner()
    synthetic = build_synthetic_input(request, source)
    shifted_lines: list[str] = []
    for line in raw.splitlines(keepends=True):
        if line.startswith("ATOM  ") and line[21] == "B":
            x = float(line[30:38]) + 4.0
            line = line[:30] + f"{x:8.3f}" + line[38:]
        shifted_lines.append(line)

    with pytest.raises(
        RFdiffusionAdapterError,
        match="partner context drifted after sampling",
    ):
        materialize_candidate(
            request,
            source_text=source,
            synthetic_text=synthetic.pdb_text,
            raw_text="".join(shifted_lines),
            partner_chains=synthetic.partner_chains,
        )


def test_boundary_refinement_preserves_fixed_records_and_generated_identity() -> None:
    request, source, raw = _case()
    synthetic = build_synthetic_input(request, source)
    unrefined = materialize_candidate(
        request,
        source_text=source,
        synthetic_text=synthetic.pdb_text,
        raw_text=raw,
    )

    result = refine_generated_region(
        unrefined,
        generated_residues=request.gaps[0].generated_residues,
        generated_atoms=request.generated_atoms,
        max_iterations=50,
        restart_count=0,
    )

    assert set(result.coordinates_angstrom) == set(request.generated_atoms)
    assert result.final_energy_kj_mol < result.initial_energy_kj_mol
    assert len(result.pre_coordinate_sha256) == 64
    assert len(result.post_coordinate_sha256) == 64


def test_boundary_refinement_same_seed_is_exactly_repeatable() -> None:
    request, source, raw = _case()
    synthetic = build_synthetic_input(request, source)
    unrefined = materialize_candidate(
        request,
        source_text=source,
        synthetic_text=synthetic.pdb_text,
        raw_text=raw,
    )

    first = refine_generated_region(
        unrefined,
        generated_residues=request.gaps[0].generated_residues,
        generated_atoms=request.generated_atoms,
        max_iterations=10,
        restart_count=1,
        random_seed=7,
    )
    second = refine_generated_region(
        unrefined,
        generated_residues=request.gaps[0].generated_residues,
        generated_atoms=request.generated_atoms,
        max_iterations=10,
        restart_count=1,
        random_seed=7,
    )

    assert first.pre_coordinate_sha256 == second.pre_coordinate_sha256
    assert first.post_coordinate_sha256 == second.post_coordinate_sha256


def test_materialization_with_refinement_keeps_source_atom_lines_exact() -> None:
    request, source, raw = _case()
    synthetic = build_synthetic_input(request, source)

    candidate = materialize_candidate(
        request,
        source_text=source,
        synthetic_text=synthetic.pdb_text,
        raw_text=raw,
        boundary_refinement=True,
        refinement_restart_count=0,
    )

    source_atom_lines = {
        line
        for line in source.splitlines()
        if line.startswith(("ATOM  ", "HETATM"))
    }
    candidate_lines = set(candidate.splitlines())
    assert source_atom_lines <= candidate_lines
    candidate_atoms = {
        _atom(line)
        for line in candidate.splitlines()
        if line.startswith(("ATOM  ", "HETATM")) and _heavy(line)
    }
    assert set(request.generated_atoms) <= candidate_atoms

    pdb = PDBFile(StringIO(candidate))
    lysine = next(
        residue
        for residue in pdb.topology.residues()
        if residue.chain.id == "C" and residue.id == "68"
    )
    atoms = {atom.name: atom for atom in lysine.atoms()}

    def distance(first: str, second: str) -> float:
        first_xyz = pdb.positions[atoms[first].index].value_in_unit(angstrom)
        second_xyz = pdb.positions[atoms[second].index].value_in_unit(angstrom)
        return math.dist(first_xyz, second_xyz)

    assert 1.3 < distance("CD", "CE") < 1.7
    assert 1.3 < distance("CE", "NZ") < 1.7

    lysine_clashes = []
    for finding in clashes_python(
        pdb.topology,
        pdb.positions,
        clash_warn_a=0.9,
        clash_error_a=0.9,
    ):
        partner = finding.extra.get("clash_partner", {})
        if finding.severity is Severity.ERROR and (
            (finding.chain, finding.resid) == ("C", "68")
            or (partner.get("chain"), partner.get("resid")) == ("C", "68")
        ):
            lysine_clashes.append(finding)
    assert not lysine_clashes
