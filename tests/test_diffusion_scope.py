"""Tests for fail-closed admission to the initial diffusion slice."""

from __future__ import annotations

import hashlib
from dataclasses import replace
from pathlib import Path

import pytest

from dvbfixer.model.diffusion.contract import (
    DIFFUSION_SCHEMA_VERSION,
    ArtifactReference,
    AtomIdentity,
    DiffusionRequest,
    ExplicitLink,
    GapRegion,
    ResidueIdentity,
    SequencePlacement,
    TargetInterval,
    TargetSequence,
)
from dvbfixer.model.diffusion.scope import (
    DiffusionScopeError,
    assess_diffusion_scope,
)

FIXTURE = Path(__file__).parent / "fixtures" / "8cz8" / "8cz8_a_u.pdb"
GENERATED_NUMBERS = tuple(str(number) for number in range(65, 70))


def _digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _atom_identity(line: str) -> AtomIdentity:
    return AtomIdentity(
        line[21],
        line[22:26].strip(),
        line[26].strip(),
        line[12:16].strip(),
    )


def _residue_identity(line: str) -> ResidueIdentity:
    return ResidueIdentity(line[21], line[22:26].strip(), line[26].strip())


def _source_and_request() -> tuple[bytes, DiffusionRequest]:
    source_lines = FIXTURE.read_text(encoding="utf-8").splitlines(keepends=True)
    generated_residues = {
        ResidueIdentity("C", number)
        for number in GENERATED_NUMBERS
    }
    selected_numbers = {str(number) for number in range(61, 74)}
    selected_lines = [
        line
        for line in source_lines
        if line.startswith(("ATOM  ", "HETATM"))
        and line[21] == "C"
        and line[22:26].strip() in selected_numbers
    ]
    source = (
        "".join(
            line
            for line in selected_lines
            if _residue_identity(line) not in generated_residues
        )
        + "TER\nEND\n"
    ).encode()
    source_residues: list[ResidueIdentity] = []
    seen_residues: set[ResidueIdentity] = set()
    fixed_atoms = tuple(
        _atom_identity(line)
        for line in source.decode().splitlines()
        if line.startswith(("ATOM  ", "HETATM"))
        and line[76:78].strip().upper() != "H"
    )
    for line in source.decode().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        residue = _residue_identity(line)
        if residue not in seen_residues:
            source_residues.append(residue)
            seen_residues.add(residue)
    generated_atoms = tuple(
        _atom_identity(line)
        for line in selected_lines
        if _residue_identity(line) in generated_residues
        and line[76:78].strip().upper() != "H"
    )
    left = ResidueIdentity("C", "64")
    generated = tuple(
        ResidueIdentity("C", number)
        for number in GENERATED_NUMBERS
    )
    right = ResidueIdentity("C", "70")
    target_sequence = "SNRFSGSKSGNTA"
    observed_indices = (0, 1, 2, 3, 9, 10, 11, 12)
    request = DiffusionRequest(
        schema_version=DIFFUSION_SCHEMA_VERSION,
        normalized_pdb=ArtifactReference("input/normalized.pdb", _digest(source)),
        target_sequences=(TargetSequence("C", target_sequence),),
        sequence_placements=(
            SequencePlacement(
                chain="C",
                target_length=len(target_sequence),
                observed_target_indices=observed_indices,
                observed_residues=tuple(source_residues),
            ),
        ),
        gaps=(
            GapRegion(
                chain="C",
                target_interval=TargetInterval(4, 9),
                left_anchor=left,
                right_anchor=right,
                generated_residues=generated,
                movable_junction_residues=(left, *generated, right),
            ),
        ),
        fixed_atoms=fixed_atoms,
        generated_atoms=generated_atoms,
        retained_explicit_links=(),
        candidate_count=2,
        seeds=(7, 11),
    )
    return source, request


def _replace_source(
    request: DiffusionRequest,
    source: bytes,
) -> DiffusionRequest:
    return replace(
        request,
        normalized_pdb=ArtifactReference(
            request.normalized_pdb.path,
            _digest(source),
        ),
    )


def test_scope_accepts_canonical_two_anchor_internal_gap() -> None:
    source, request = _source_and_request()

    admission = assess_diffusion_scope(request, source)

    assert admission.supported is True
    assert admission.reasons == ()


def test_scope_rejects_multiple_models_and_target_ambiguity() -> None:
    source, request = _source_and_request()
    first_end = source.rfind(b"END")
    coordinates = source[:first_end]
    multi_model = (
        b"MODEL        1\n"
        + coordinates
        + b"ENDMDL\nMODEL        2\n"
        + coordinates
        + b"ENDMDL\nEND\n"
    )
    multi_request = _replace_source(request, multi_model)

    admission = assess_diffusion_scope(multi_request, multi_model)

    assert admission.supported is False
    assert "multiple-models" in admission.reasons

    sequence = request.target_sequences[0].sequence
    noncanonical = replace(
        request,
        target_sequences=(
            TargetSequence("C", sequence[:6] + "X" + sequence[7:]),
        ),
    )
    admission = assess_diffusion_scope(noncanonical, source)
    assert admission.supported is False
    assert "non-canonical-target-sequence" in admission.reasons


def test_scope_rejects_adjacent_altloc_and_retained_heterogen() -> None:
    source, request = _source_and_request()
    text = source.decode()
    anchor_line = next(
        line
        for line in text.splitlines(keepends=True)
        if line.startswith("ATOM  ")
        and line[21] == "C"
        and line[22:26].strip() == "64"
        and line[12:16].strip() == "CA"
    )
    altloc_line = anchor_line[:16] + "A" + anchor_line[17:]
    altloc = text.replace(anchor_line, altloc_line).encode()
    altloc_request = _replace_source(request, altloc)

    admission = assess_diffusion_scope(altloc_request, altloc)

    assert admission.supported is False
    assert "adjacent-alternate-location-ambiguity" in admission.reasons

    second_altloc_line = anchor_line[:16] + "B" + anchor_line[17:]
    duplicate_altloc = text.replace(
        anchor_line,
        altloc_line + second_altloc_line,
    ).encode()
    duplicate_altloc_request = _replace_source(request, duplicate_altloc)
    admission = assess_diffusion_scope(
        duplicate_altloc_request,
        duplicate_altloc,
    )
    assert admission.supported is False
    assert "adjacent-alternate-location-ambiguity" in admission.reasons

    hetatm = text.replace("ATOM  ", "HETATM", 1).encode()
    hetatm_request = _replace_source(request, hetatm)
    admission = assess_diffusion_scope(hetatm_request, hetatm)
    assert admission.supported is False
    assert "retained-heterogen-or-noncanonical-residue" in admission.reasons


def test_scope_uses_heavy_atoms_for_fixed_mask() -> None:
    source, request = _source_and_request()
    lines = source.decode().splitlines(keepends=True)
    first_atom = next(
        line
        for line in lines
        if line.startswith("ATOM  ")
    )
    hydrogen = (
        first_atom[:6]
        + "99999"
        + first_atom[11:12]
        + " H1 "
        + first_atom[16:76]
        + " H"
        + first_atom[78:]
    )
    hydrogen_source = (
        "".join(
            line
            for line in lines
            if not line.startswith("TER")
            and not line.startswith("END")
        )
        + hydrogen
        + "TER\nEND\n"
    ).encode()
    hydrogen_request = _replace_source(request, hydrogen_source)

    admission = assess_diffusion_scope(
        hydrogen_request,
        hydrogen_source,
    )

    assert admission.supported is True


def test_scope_rejects_mask_and_sequence_mismatches() -> None:
    source, request = _source_and_request()
    bad_fixed = replace(request, fixed_atoms=request.fixed_atoms[1:])

    admission = assess_diffusion_scope(bad_fixed, source)

    assert admission.supported is False
    assert "fixed-atom-mask-mismatch" in admission.reasons

    bad_generated = replace(request, generated_atoms=request.generated_atoms[1:4])
    admission = assess_diffusion_scope(bad_generated, source)
    assert admission.supported is False
    assert "generated-atom-mask-mismatch" in admission.reasons

    sequence = request.target_sequences[0].sequence
    wrong_sequence = replace(
        request,
        target_sequences=(
            TargetSequence("C", "A" + sequence[1:]),
        ),
    )
    admission = assess_diffusion_scope(wrong_sequence, source)
    assert admission.supported is False
    assert "observed-sequence-placement-mismatch" in admission.reasons

    placement = request.sequence_placements[0]
    incomplete_coverage = replace(
        request,
        sequence_placements=(
            replace(
                placement,
                observed_target_indices=placement.observed_target_indices[:-1],
                observed_residues=placement.observed_residues[:-1],
            ),
        ),
    )
    admission = assess_diffusion_scope(incomplete_coverage, source)
    assert admission.supported is False
    assert "sequence-placement-gap-coverage-mismatch" in admission.reasons


def test_scope_rejects_external_covalent_link_near_gap() -> None:
    source, request = _source_and_request()
    link = ExplicitLink(
        AtomIdentity("C", "64", "", "CA"),
        AtomIdentity("C", "70", "", "CA"),
        "LINK",
    )
    link_line = (
        "LINK         CA  SER C  64                 CA  SER C  70\n"
    )
    linked_source = source.replace(b"TER\n", link_line.encode() + b"TER\n")
    linked = replace(
        _replace_source(request, linked_source),
        retained_explicit_links=(link,),
    )

    admission = assess_diffusion_scope(linked, linked_source)

    assert admission.supported is False
    assert "external-covalent-link-near-generated-region" in admission.reasons

    missing_atom_link = ExplicitLink(
        AtomIdentity("C", "64", "", "CA"),
        AtomIdentity("C", "65", "", "CA"),
        "LINK",
    )
    linked = replace(request, retained_explicit_links=(missing_atom_link,))
    admission = assess_diffusion_scope(linked, source)
    assert admission.supported is False
    assert "retained-explicit-link-atom-missing" in admission.reasons
    assert "retained-explicit-link-mismatch" in admission.reasons


def test_scope_rejects_unreported_source_explicit_link() -> None:
    source, request = _source_and_request()
    atoms = {
        atom.atom_name: atom
        for atom in request.fixed_atoms
        if atom.residue_number == "61"
    }
    assert "N" in atoms and "CA" in atoms
    source_lines = source.decode().splitlines(keepends=True)
    serials = {
        _atom_identity(line): int(line[6:11])
        for line in source_lines
        if line.startswith("ATOM  ")
    }
    conect = (
        f"CONECT{serials[atoms['N']]:5d}{serials[atoms['CA']]:5d}\n"
    ).encode()
    linked_source = source.replace(b"TER\n", conect + b"TER\n")
    linked_request = _replace_source(request, linked_source)

    admission = assess_diffusion_scope(linked_request, linked_source)

    assert admission.supported is False
    assert "retained-explicit-link-mismatch" in admission.reasons


def test_scope_preserves_case_sensitive_chain_identity() -> None:
    source, request = _source_and_request()
    lower_source = b"".join(
        line[:21] + b"d" + line[22:]
        if line.startswith((b"ATOM  ", b"HETATM"))
        else line
        for line in source.splitlines(keepends=True)
    )
    lower_fixed = tuple(
        AtomIdentity("d", atom.residue_number, atom.insertion_code, atom.atom_name)
        for atom in request.fixed_atoms
    )
    lower_generated = tuple(
        AtomIdentity("d", atom.residue_number, atom.insertion_code, atom.atom_name)
        for atom in request.generated_atoms
    )
    lower_observed_residues = tuple(
        ResidueIdentity(
            "d",
            residue.residue_number,
            residue.insertion_code,
        )
        for residue in request.sequence_placements[0].observed_residues
    )
    lower_left = ResidueIdentity("d", "64")
    lower_generated_residues = tuple(
        ResidueIdentity("d", number)
        for number in GENERATED_NUMBERS
    )
    lower_right = ResidueIdentity("d", "70")
    lower_request = DiffusionRequest(
        schema_version=DIFFUSION_SCHEMA_VERSION,
        normalized_pdb=ArtifactReference(
            request.normalized_pdb.path,
            _digest(lower_source),
        ),
        target_sequences=(
            TargetSequence("d", request.target_sequences[0].sequence),
        ),
        sequence_placements=(
            SequencePlacement(
                chain="d",
                target_length=request.sequence_placements[0].target_length,
                observed_target_indices=(
                    request.sequence_placements[0].observed_target_indices
                ),
                observed_residues=lower_observed_residues,
            ),
        ),
        gaps=(
            GapRegion(
                chain="d",
                target_interval=TargetInterval(4, 9),
                left_anchor=lower_left,
                right_anchor=lower_right,
                generated_residues=lower_generated_residues,
                movable_junction_residues=(
                    lower_left,
                    *lower_generated_residues,
                    lower_right,
                ),
            ),
        ),
        fixed_atoms=lower_fixed,
        generated_atoms=lower_generated,
        retained_explicit_links=(),
        candidate_count=request.candidate_count,
        seeds=request.seeds,
    )

    admission = assess_diffusion_scope(lower_request, lower_source)

    assert admission.supported is True


def test_scope_fails_closed_on_digest_or_malformed_coordinates() -> None:
    source, request = _source_and_request()

    with pytest.raises(DiffusionScopeError, match="SHA-256"):
        assess_diffusion_scope(request, source + b"REMARK changed\n")

    lines = source.splitlines(keepends=True)
    atom_index = next(
        index
        for index, line in enumerate(lines)
        if line.startswith(b"ATOM  ")
    )
    lines[atom_index] = (
        lines[atom_index][:30]
        + b"notfloat"
        + lines[atom_index][38:]
    )
    malformed = b"".join(lines)
    malformed_request = _replace_source(request, malformed)
    with pytest.raises(DiffusionScopeError, match="malformed PDB"):
        assess_diffusion_scope(malformed_request, malformed)
