from pathlib import Path

import pytest

from dvbfixer.structure_input import (
    StructureInputError,
    normalize_structure,
    normalized_command_inputs,
    run_with_normalized_inputs,
)


def _mmcif(*, chains: tuple[str, ...] = ("AA",)) -> str:
    rows = []
    atom_id = 1
    for chain_index, chain in enumerate(chains, 1):
        label = f"X{chain_index}"
        for atom_name, element, x in (
            ("N", "N", 0.0), ("CA", "C", 1.4), ("C", "C", 2.4), ("O", "O", 3.4),
        ):
            rows.append(
                f"ATOM {atom_id} {element} {atom_name} . ALA {label} {chain_index} 1 ? "
                f"{x:.1f} {chain_index:.1f} 0 1 10 ? 1 ALA {chain} {atom_name} 1"
            )
            atom_id += 1
    return """data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
""" + "\n".join(rows) + "\n#\n"


def test_mmcif_conversion_preserves_single_ids_and_maps_only_long_ids(tmp_path: Path):
    source = tmp_path / "input.cif"
    output = tmp_path / "internal.pdb"
    source.write_text(_mmcif(chains=("A", "LONG", "b")))

    converted = normalize_structure(source, output)

    assert converted.dialect == "pdbx/mmCIF"
    assert converted.chain_map == {"A": "A", "LONG": "B", "b": "b"}
    text = output.read_text()
    assert "CIF_CHAIN_MAP B LONG" in text
    assert {line[21] for line in text.splitlines() if line.startswith("ATOM")} == {"A", "B", "b"}


def test_small_cif_converts_fractional_asymmetric_unit_and_bonds(tmp_path: Path):
    source = tmp_path / "ligand.cif"
    output = tmp_path / "internal.pdb"
    source.write_text("""data_small
_cell_length_a 10
_cell_length_b 10
_cell_length_c 10
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'P 1'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
C1 C 0.1 0.2 0.3
O1 O 0.2 0.2 0.3
loop_
_geom_bond_atom_site_label_1
_geom_bond_atom_site_label_2
_geom_bond_distance
C1 O1 1.0
""")

    converted = normalize_structure(source, output)

    assert converted.dialect == "small-molecule CIF"
    text = output.read_text()
    assert "CRYST1   10.000" in text
    assert "       1.000   2.000   3.000" in text
    first_atom = next(line for line in text.splitlines() if line.startswith("HETATM"))
    assert first_atom[12:16].strip() == "C1"
    assert first_atom[17:20] == "LIG"
    assert first_atom[21] == "A"
    assert "CONECT" in text
    assert sum(line.startswith("HETATM") for line in text.splitlines()) == 2


def test_command_wrapper_keeps_original_default_output_and_translates_chain(tmp_path: Path):
    source = tmp_path / "complex.cif"
    source.write_text(_mmcif(chains=("AA",)))

    with normalized_command_inputs("prepare", [str(source), "--cap-chain", "AA"]) as argv:
        assert Path(argv[0]).suffix == ".pdb"
        assert Path(argv[0]).parent.parent.name.startswith(".dvbfixer_cif_")
        assert argv[argv.index("--cap-chain") + 1] == "A"
        assert argv[argv.index("-o") + 1] == str(tmp_path / "complex_prepared.pdb")


def test_command_wrapper_translates_fasta_headers(tmp_path: Path):
    source = tmp_path / "complex.cif"
    source.write_text(_mmcif(chains=("AA",)))
    fasta = tmp_path / "chains.fasta"
    fasta.write_text(">chain_AA description\nAAAA\n")

    with normalized_command_inputs(
        "model", [str(source), "--fasta", str(fasta)],
    ) as argv:
        mapped_fasta = Path(argv[argv.index("--fasta") + 1])
        assert mapped_fasta.read_text() == ">chain_A description\nAAAA\n"


def test_command_wrapper_translates_multi_chain_selector_grammars(tmp_path: Path):
    source = tmp_path / "complex.cif"
    source.write_text(_mmcif(chains=("LEFT", "RIGHT")))

    with normalized_command_inputs(
        "top", [str(source), "--ss", "LEFT:10:RIGHT:20"],
    ) as argv:
        assert argv[argv.index("--ss") + 1] == "A:10:B:20"

    with normalized_command_inputs(
        "transplant", [str(source), "--donor", str(source), "--align", "LEFT:RIGHT"],
    ) as argv:
        assert argv[argv.index("--align") + 1] == "A:B"


def test_homology_template_plan_cif_paths_and_chains_are_rewritten(tmp_path: Path):
    source = tmp_path / "template.cif"
    source.write_text(_mmcif(chains=("AA",)))
    plan = tmp_path / "plan.json"
    plan.write_text(
        '{"templates": [{"file": "template.cif", "chain": "AA", "targetChain": "T"}]}\n'
    )

    with normalized_command_inputs("homology", ["target.fasta", "--template-plan", str(plan)]) as argv:
        import json

        mapped_plan = Path(argv[argv.index("--template-plan") + 1])
        payload = json.loads(mapped_plan.read_text())
        template = payload["templates"][0]
        assert Path(template["file"]).suffix == ".pdb"
        assert template["chain"] == "A"
        assert template["targetChain"] == "T"


def test_rejects_unrecognized_cif_dialect(tmp_path: Path):
    source = tmp_path / "invalid.cif"
    source.write_text("data_invalid\n_chemical_formula_sum 'C H4'\n")
    with pytest.raises(StructureInputError, match="unrecognized CIF dialect"):
        normalize_structure(source, tmp_path / "out.pdb")


def test_existing_cif_output_path_is_not_mistaken_for_an_input(tmp_path: Path):
    source = tmp_path / "input.cif"
    source.write_text(_mmcif())
    output = tmp_path / "existing.cif"
    output.write_text("old output")

    with pytest.raises(StructureInputError, match="CIF output is not supported"):
        with normalized_command_inputs("prepare", [str(source), "-o", str(output)]):
            pass


def test_mmcif_conversion_retains_sequences_and_biological_assemblies(
    tmp_path: Path, biological_assembly_pdb: Path,
):
    import gemmi

    source_structure = gemmi.read_structure(str(biological_assembly_pdb))
    cif = tmp_path / "assembly.mmcif"
    source_structure.make_mmcif_document().write_file(str(cif))
    output = tmp_path / "assembly.pdb"

    normalize_structure(cif, output)

    lines = output.read_text().splitlines(keepends=True)
    assert sum(line.startswith("SEQRES") for line in lines) > 0
    from dvbfixer.biological_assembly import parse_biological_assemblies
    assemblies = parse_biological_assemblies(lines)
    assert set(assemblies) == {"1", "2", "3", "4"}


def test_existing_pdb_command_runs_unchanged_after_cif_normalization(tmp_path: Path):
    from dvbfixer.rename import main as rename_main

    source = tmp_path / "complex.cif"
    source.write_text(_mmcif(chains=("LONG",)))

    run_with_normalized_inputs("rename", rename_main, [str(source)])

    output = tmp_path / "complex_renamed.pdb"
    assert output.is_file()
    text = output.read_text()
    assert "REMARK 999 DVBFIXER CIF_CHAIN_MAP A LONG" in text
    assert all(line[21] == "A" for line in text.splitlines() if line.startswith("ATOM"))
