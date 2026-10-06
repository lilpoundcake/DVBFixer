from __future__ import annotations

import hashlib
import json
from dataclasses import replace
from pathlib import Path

import pytest

from dvbfixer.ccd import (
    CcdError,
    CcdProvenance,
    ccd_to_heavy_atom_graph,
    load_ccd_component,
    parse_ccd_component,
)
from dvbfixer.domain.nonprotein_reconstruction import ComponentClass, ReconstructionStatus
from dvbfixer.pdb_component_reconstruction import (
    PdbComponentError,
    publish_pdb_bundle,
    reconstruct_pdb_component,
)


def _ccd(component_id: str = "LIG", *, ambiguous: str = "N") -> bytes:
    return f"""data_{component_id}
_chem_comp.id {component_id}
_chem_comp.name 'TEST LIGAND'
_chem_comp.type NON-POLYMER
_chem_comp.formula 'C2 N O'
_chem_comp.pdbx_formal_charge -1
_chem_comp.formula_weight 70.0
_chem_comp.pdbx_initial_date 2026-01-01
_chem_comp.pdbx_modified_date 2026-01-02
_chem_comp.pdbx_ambiguous_flag {ambiguous}
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.charge
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
_chem_comp_atom.pdbx_model_Cartn_x_ideal
_chem_comp_atom.pdbx_model_Cartn_y_ideal
_chem_comp_atom.pdbx_model_Cartn_z_ideal
C1 C1 C 0 N N R 0.0 0.0 0.0
C2 C2 C 0 N N N 1.0 0.0 0.0
N1 N1 N 0 N N N 0.0 1.0 0.0
O1 O1 O -1 N N N 0.0 0.0 1.0
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
C1 C2 SING N N
C1 N1 SING N N
C1 O1 SING N N
""".encode()


def _component(component_id: str = "LIG", *, ambiguous: str = "N"):
    content = _ccd(component_id, ambiguous=ambiguous)
    return parse_ccd_component(
        content,
        expected_component_id=component_id,
        provenance=CcdProvenance(
            "local-cif", f"/{component_id}.cif", hashlib.sha256(content).hexdigest()
        ),
    )


def _atom(
    serial: int,
    name: str,
    residue: str,
    chain: str,
    number: int,
    position: tuple[float, float, float],
    element: str,
    *,
    record: str = "HETATM",
) -> str:
    atom_field = f" {name:<3}" if len(name) < 4 else name
    x, y, z = position
    return (
        f"{record:<6}{serial:5d} {atom_field} {residue:>3} {chain}{number:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00          {element:>2}\n"
    )


def _pdb(component_id: str = "LIG", *, external_link: bool = False) -> bytes:
    lines = [
        "HEADER    COMPONENT TEST\n",
        _atom(1, "C1", component_id, "A", 10, (10.0, 20.0, 30.0), "C"),
        _atom(2, "C2", component_id, "A", 10, (11.0, 20.0, 30.0), "C"),
        _atom(3, "N1", component_id, "A", 10, (10.0, 21.0, 30.0), "N"),
    ]
    if external_link:
        lines.extend(
            [
                _atom(9, "SG", "CYS", "B", 20, (12.0, 20.0, 30.0), "S", record="ATOM"),
                "CONECT    1    9\n",
            ]
        )
    lines.append("END\n")
    return "".join(lines).encode("ascii")


def test_ccd_parser_exposes_authoritative_chemistry_and_heavy_graph() -> None:
    component = _component()
    graph = ccd_to_heavy_atom_graph(component)

    assert component.component_id == "LIG"
    assert component.name == "TEST LIGAND"
    assert component.provisional_class == "A"
    assert len(component.heavy_atoms) == 4
    assert graph.net_charge == -1
    assert graph.state_id.startswith("ccd-")
    assert [(bond.atom1, bond.atom2, bond.order) for bond in graph.bonds] == [
        ("C1", "C2", 1.0),
        ("C1", "N1", 1.0),
        ("C1", "O1", 1.0),
    ]


def test_ccd_information_preserves_unknown_cofactor_charges_without_guessing() -> None:
    content = _ccd("HEM").replace(
        b"_chem_comp.pdbx_formal_charge -1", b"_chem_comp.pdbx_formal_charge ?"
    ).replace(b"O1 O1 O -1", b"O1 O1 O ?")
    component = parse_ccd_component(
        content,
        expected_component_id="HEM",
        provenance=CcdProvenance(
            "local-cif", "/HEM.cif", hashlib.sha256(content).hexdigest()
        ),
    )

    assert component.formal_charge is None
    assert component.atoms[-1].formal_charge is None
    assert component.provisional_class == "C"
    with pytest.raises(CcdError, match="formal charges must be explicit"):
        ccd_to_heavy_atom_graph(component)


def test_ambiguous_ccd_record_is_not_a_reconstruction_authority() -> None:
    component = _component(ambiguous="Y")
    assert component.provisional_class == "D"
    with pytest.raises(CcdError, match="chemically ambiguous"):
        ccd_to_heavy_atom_graph(component)


def test_local_reconstruction_authority_requires_exact_digest_pin(tmp_path: Path) -> None:
    path = tmp_path / "LIG.cif"
    path.write_bytes(_ccd())

    with pytest.raises(CcdError, match="does not match"):
        load_ccd_component(path, "LIG", expected_sha256="0" * 64)

    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    assert load_ccd_component(path, "LIG", expected_sha256=digest).provenance.sha256 == digest


def test_online_discovery_record_cannot_be_used_for_reconstruction() -> None:
    component = _component()
    online = replace(
        component,
        provenance=replace(component.provenance, mode="online-ccd"),
    )

    with pytest.raises(PdbComponentError, match="online lookup is read-only"):
        reconstruct_pdb_component(
            _pdb(), online, model=1, chain_id="A", sequence_number="10"
        )


def test_class_a_reconstruction_preserves_source_atoms_and_adds_missing_heavy_atom() -> None:
    source = _pdb()
    bundle, output = reconstruct_pdb_component(
        source,
        _component(),
        model=1,
        chain_id="A",
        sequence_number="10",
    )

    assert output is not None
    assert bundle.result.status is ReconstructionStatus.SUCCEEDED
    assert bundle.result.component_class is ComponentClass.A
    assert bundle.added_atom_names == ("O1",)
    source_atoms = [line for line in source.splitlines() if line.startswith((b"ATOM  ", b"HETATM"))]
    output_atoms = [line for line in output.splitlines() if line.startswith((b"ATOM  ", b"HETATM"))]
    assert output_atoms[:3] == source_atoms
    assert output_atoms[3][12:16].strip() == b"O1"
    assert tuple(float(output_atoms[3][start : start + 8]) for start in (30, 38, 46)) == (
        10.0,
        20.0,
        31.0,
    )
    assert b"CONECT    4    1\n" in output


def test_explicit_external_conect_classifies_component_as_unsupported_class_b() -> None:
    bundle, output = reconstruct_pdb_component(
        _pdb(external_link=True),
        _component(),
        model=1,
        chain_id="A",
        sequence_number="10",
    )

    assert output is None
    assert bundle.result.status is ReconstructionStatus.UNSUPPORTED
    assert bundle.result.component_class is ComponentClass.B
    assert bundle.result.findings[0].code == "class-b-unsupported"


def test_protected_cofactor_is_refused_as_class_c_before_materialization() -> None:
    bundle, output = reconstruct_pdb_component(
        _pdb("FAD"),
        _component("FAD"),
        model=1,
        chain_id="A",
        sequence_number="10",
    )

    assert output is None
    assert bundle.result.status is ReconstructionStatus.UNSUPPORTED
    assert bundle.result.component_class is ComponentClass.C
    assert bundle.result.findings[0].code == "class-c-unsupported"


def test_environment_clash_prevents_materialization() -> None:
    source = _pdb().replace(
        b"END\n",
        _atom(9, "CA", "ALA", "B", 20, (10.0, 20.0, 31.0), "C", record="ATOM").encode()
        + b"END\n",
    )
    with pytest.raises(PdbComponentError, match="clashes"):
        reconstruct_pdb_component(
            source,
            _component(),
            model=1,
            chain_id="A",
            sequence_number="10",
        )


def test_successful_structure_and_provenance_publish_as_atomic_bundle(tmp_path: Path) -> None:
    bundle, output = reconstruct_pdb_component(
        _pdb(),
        _component(),
        model=1,
        chain_id="A",
        sequence_number="10",
    )
    assert output is not None

    destination = publish_pdb_bundle(bundle, output, tmp_path, "ligand-result")

    assert (destination / "structure.pdb").read_bytes() == output
    report = json.loads((destination / "provenance.json").read_text())
    assert report["output_sha256"] == hashlib.sha256(output).hexdigest()
    assert report["result"]["geometry_approved"] is False
    assert report["result"]["md_ready"] is False
    with pytest.raises(FileExistsError):
        publish_pdb_bundle(bundle, output, tmp_path, "ligand-result")
