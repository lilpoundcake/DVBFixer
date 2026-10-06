"""Tests for record-preserving diffusion candidate materialization."""

from __future__ import annotations

import pytest

from dvbfixer.model.diffusion.contract import (
    DIFFUSION_SCHEMA_VERSION,
    ArtifactReference,
    AtomIdentity,
    DiffusionContractError,
    DiffusionRequest,
    GapKind,
    GapRegion,
    ResidueIdentity,
    SequencePlacement,
    TargetInterval,
    TargetSequence,
)
from dvbfixer.model.diffusion.pdb_materialize import materialize_candidate_pdb


def _atom_line(
    serial: int,
    residue: int,
    x: float,
    *,
    occupancy: float,
    bfactor: float,
) -> str:
    return (
        f"ATOM  {serial:5d}  CA  ALA A{residue:4d}    "
        f"{x:8.3f}{2.0:8.3f}{3.0:8.3f}{occupancy:6.2f}{bfactor:6.2f}"
        "           C  \n"
    )


def _request() -> DiffusionRequest:
    left = ResidueIdentity("A", "10")
    generated = tuple(ResidueIdentity("A", str(number)) for number in range(11, 14))
    right = ResidueIdentity("A", "14")
    return DiffusionRequest(
        schema_version=DIFFUSION_SCHEMA_VERSION,
        normalized_pdb=ArtifactReference("input.pdb", "a" * 64),
        target_sequences=(TargetSequence("A", "AAAAA"),),
        sequence_placements=(
            SequencePlacement("A", 5, (0, 4), (left, right)),
        ),
        gaps=(
            GapRegion(
                "A",
                TargetInterval(1, 4),
                left,
                right,
                generated,
                (left, *generated, right),
            ),
        ),
        fixed_atoms=(
            AtomIdentity("A", "10", "", "CA"),
            AtomIdentity("A", "14", "", "CA"),
        ),
        generated_atoms=tuple(
            AtomIdentity("A", residue.residue_number, "", "CA")
            for residue in generated
        ),
        retained_explicit_links=(),
        candidate_count=1,
        seeds=(7,),
    )


def test_materializer_preserves_headers_metadata_anisou_and_links() -> None:
    request = _request()
    left = _atom_line(10, 10, 1.0, occupancy=0.55, bfactor=42.25)
    right = _atom_line(40, 14, 9.0, occupancy=0.65, bfactor=31.50)
    anisou = "ANISOU   10  CA  ALA A  10    10000  11000  12000    100    200    300       C  \n"
    source = (
        "HEADER    RECORD-PRESERVATION TEST\n"
        "SEQRES   1 A    5  ALA ALA ALA ALA ALA\n"
        + left
        + anisou
        + "REMARK 999 KEEP THIS RECORD\n"
        + right
        + "TER      41      ALA A  14\n"
        + "CONECT   10   40\nEND\n"
    )
    axis = (*request.fixed_atoms, *request.generated_atoms)
    coordinates = (
        (1.125, 2.25, 3.375),
        (9.125, 2.25, 3.375),
        (4.0, 5.0, 6.0),
        (5.0, 6.0, 7.0),
        (6.0, 7.0, 8.0),
    )

    rendered = materialize_candidate_pdb(
        source,
        axis,
        coordinates,
        ("ALA",) * len(axis),
        request,
        elements=("C",) * len(axis),
    )

    assert "HEADER    RECORD-PRESERVATION TEST\n" in rendered
    assert "SEQRES   1 A    5  ALA ALA ALA ALA ALA\n" in rendered
    assert anisou in rendered
    assert "REMARK 999 KEEP THIS RECORD\n" in rendered
    assert "CONECT   10   40\nEND\n" in rendered
    rendered_lines = rendered.splitlines(keepends=True)
    left_rendered = next(line for line in rendered_lines if line.startswith("ATOM     10"))
    right_rendered = next(line for line in rendered_lines if line.startswith("ATOM     40"))
    assert left_rendered[:30] == left[:30]
    assert left_rendered[54:] == left[54:]
    assert right_rendered[:30] == right[:30]
    assert right_rendered[54:] == right[54:]
    generated = [
        line
        for line in rendered_lines
        if line.startswith("ATOM  ") and line[22:26].strip() in {"11", "12", "13"}
    ]
    assert len(generated) == 3
    assert max(rendered_lines.index(line) for line in generated) < rendered_lines.index(
        right_rendered
    )
    assert {int(line[6:11]) for line in generated}.isdisjoint({10, 40, 41})


def test_materializer_preserves_unrepresented_fixed_atoms_from_source() -> None:
    request = _request()
    left = _atom_line(10, 10, 1.0, occupancy=1.0, bfactor=0.0)
    right = _atom_line(40, 14, 9.0, occupancy=1.0, bfactor=0.0)
    source = left + right + "END\n"

    rendered = materialize_candidate_pdb(
        source,
        request.generated_atoms,
        ((4.0, 5.0, 6.0), (5.0, 6.0, 7.0), (6.0, 7.0, 8.0)),
        ("ALA",) * 3,
        request,
        elements=("C",) * 3,
    )

    assert left in rendered
    assert right in rendered


@pytest.mark.parametrize(
    ("gap_kind", "generated_numbers", "expected_before_anchor"),
    [
        (GapKind.N_TERMINAL, (7, 8, 9), True),
        (GapKind.C_TERMINAL, (11, 12, 13), False),
    ],
)
def test_materializer_inserts_one_anchor_terminal_atoms(
    gap_kind: GapKind,
    generated_numbers: tuple[int, ...],
    expected_before_anchor: bool,
) -> None:
    anchor = ResidueIdentity("A", "10")
    generated = tuple(
        ResidueIdentity("A", str(number)) for number in generated_numbers
    )
    left_anchor = anchor if gap_kind is GapKind.C_TERMINAL else None
    right_anchor = anchor if gap_kind is GapKind.N_TERMINAL else None
    request = DiffusionRequest(
        schema_version=DIFFUSION_SCHEMA_VERSION,
        normalized_pdb=ArtifactReference("input.pdb", "a" * 64),
        target_sequences=(TargetSequence("A", "AAAA"),),
        sequence_placements=(
            SequencePlacement(
                "A",
                4,
                (3,) if gap_kind is GapKind.N_TERMINAL else (0,),
                (anchor,),
            ),
        ),
        gaps=(GapRegion(
            "A",
            TargetInterval(0, 3) if gap_kind is GapKind.N_TERMINAL else TargetInterval(1, 4),
            left_anchor,
            right_anchor,
            generated,
            (*generated, anchor),
            gap_kind,
        ),),
        fixed_atoms=(AtomIdentity("A", "10", "", "CA"),),
        generated_atoms=tuple(
            AtomIdentity("A", residue.residue_number, "", "CA")
            for residue in generated
        ),
        retained_explicit_links=(),
        candidate_count=1,
        seeds=(7,),
    )
    source = _atom_line(10, 10, 1.0, occupancy=1.0, bfactor=0.0) + "TER\nEND\n"
    rendered = materialize_candidate_pdb(
        source,
        request.generated_atoms,
        tuple((float(index), 2.0, 3.0) for index in range(3)),
        ("ALA",) * 3,
        request,
        elements=("C",) * 3,
    )
    atom_residues = [
        int(line[22:26])
        for line in rendered.splitlines()
        if line.startswith("ATOM  ")
    ]

    assert atom_residues == (
        [*generated_numbers, 10]
        if expected_before_anchor
        else [10, *generated_numbers]
    )
    assert rendered.index(f"A{generated_numbers[-1]:4d}") < rendered.index("TER")


def test_materializer_rejects_ter_between_gap_anchors() -> None:
    request = _request()
    source = (
        _atom_line(10, 10, 1.0, occupancy=1.0, bfactor=0.0)
        + "TER      11      ALA A  10\n"
        + _atom_line(40, 14, 9.0, occupancy=1.0, bfactor=0.0)
        + "END\n"
    )

    with pytest.raises(DiffusionContractError, match="TER between"):
        materialize_candidate_pdb(
            source,
            (*request.fixed_atoms, *request.generated_atoms),
            ((1.0, 2.0, 3.0),) * 5,
            ("ALA",) * 5,
            request,
            elements=("C",) * 5,
        )


def test_materializer_drops_stale_conect_and_master_records() -> None:
    request = _request()
    source = (
        _atom_line(10, 10, 1.0, occupancy=1.0, bfactor=0.0)
        + _atom_line(40, 14, 9.0, occupancy=1.0, bfactor=0.0)
        + "CONECT   10   40   99\n"
        + "MASTER        0    0    0    0    0    0    0    0    2    0    1    0\n"
        + "END\n"
    )

    rendered = materialize_candidate_pdb(
        source,
        (*request.fixed_atoms, *request.generated_atoms),
        ((1.0, 2.0, 3.0),) * 5,
        ("ALA",) * 5,
        request,
        elements=("C",) * 5,
    )

    assert "CONECT" not in rendered
    assert "MASTER" not in rendered
