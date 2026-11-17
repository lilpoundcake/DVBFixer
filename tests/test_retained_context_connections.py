import json

import gemmi
import pytest
from scripts.audit_retained_context_connections import audit, deposited_facts, projection_link_facts
from scripts.mac_reconstruction_comparison import digest, save


def source_block(members="a,Z", endpoints=("a", "Z"), kind="metalc"):
    block = gemmi.cif.Block("TEST")
    block.set_mmcif_category("_entity_poly.", {"entity_id": ["1"], "type": ["polypeptide(L)"]})
    block.set_mmcif_category("_struct_asym.", {"id": ["a", "A", "Z"], "entity_id": ["1", "1", "2"]})
    block.set_mmcif_category("_pdbx_struct_assembly_gen.", {
        "assembly_id": ["1"], "asym_id_list": [members], "oper_expression": ["1"]})
    block.set_mmcif_category("_struct_conn.", {"id": ["connection1"], "conn_type_id": [kind],
        "ptnr1_label_asym_id": [endpoints[0]], "ptnr2_label_asym_id": [endpoints[1]],
        "ptnr1_symmetry": ["1_555"], "ptnr2_symmetry": ["2_555"]})
    return block


def atom(serial, chain, number, *, icode=" ", x=0., record="ATOM  ", name="CA", residue="ALA"):
    return (f"{record}{serial:5d} {name:^4s} {residue:>3s} {chain}{number:4d}{icode}   "
            f"{x:8.3f}{0.:8.3f}{0.:8.3f}{1.:6.2f}{20.:6.2f}           C\n")


def test_original_polymer_membership_not_selected_sequence_candidates():
    # A second protein outside assembly1 is not a nonprotein ligand.
    facts = deposited_facts(source_block(members="a", endpoints=("a", "A"), kind="disulf"))
    assert facts["connections"][0]["protein_endpoint_flags"] == [True, True]
    assert facts["source_cross_boundary_types"] == {}


def test_source_connection_outside_assembly_is_not_projection_loss():
    facts = deposited_facts(source_block(members="a", endpoints=("A", "Z")))
    assert facts["source_cross_boundary_types"] == {"metalc": 1}
    assert facts["assembly_member_label_cross_boundary_types"] == {}
    assert facts["outside_assembly_member_label_cross_boundary_types"] == {"metalc": 1}
    assert facts["connections"][0]["assembly_member_label_flags"] == [False, False]


def test_both_member_labels_do_not_prove_operator_pair_transport():
    facts = deposited_facts(source_block())
    assert facts["assembly_member_label_cross_boundary_types"] == {"metalc": 1}
    assert not facts["connections"][0]["actual_operator_pair_or_coordinate_mapping_established"]
    assert facts["connections"][0]["raw_annotation"]["ptnr2_symmetry"] == "2_555"


@pytest.mark.parametrize("change", ["member", "endpoint", "duplicate", "assembly"])
def test_invalid_source_identity_fails_closed(change):
    block = source_block()
    if change == "member":
        block.set_mmcif_category("_pdbx_struct_assembly_gen.", {
            "assembly_id": ["1"], "asym_id_list": ["UNKNOWN"], "oper_expression": ["1"]})
    elif change == "endpoint":
        block = source_block(endpoints=("missing", "Z"))
    elif change == "duplicate":
        block.set_mmcif_category("_struct_asym.", {"id": ["a", "a"], "entity_id": ["1", "1"]})
    else:
        block.set_mmcif_category("_pdbx_struct_assembly_gen.", {
            "assembly_id": ["2"], "asym_id_list": ["a"], "oper_expression": ["1"]})
    with pytest.raises(ValueError):
        deposited_facts(block)


def test_actual_projection_crossing_respects_case_and_insertion_code():
    protein = atom(1, "a", 82, icode="A") + atom(2, "A", 82, icode="A", x=5.)
    ligand = atom(3, "a", 90, record="HETATM", residue="LIG", x=2.)
    full = protein + ligand + "CONECT    1    3\nCONECT    3    1\n"
    facts = projection_link_facts(full, protein)
    pair, = facts["resolved_projection_boundary_pairs"]
    assert {tuple(pair[key][k] for k in ("chain", "residue_number", "insertion_code"))
            for key in ("first", "second")} == {("a", "82", "A"), ("a", "90", "")}
    assert facts["resolved_full_explicit_pairs"] == 1
    assert facts["raw_full_connection_record_counts"] == {"CONECT": 2}
    assert not facts["scientific_admission"]


def test_dangling_record_preserved_not_partially_rewritten_or_qualified():
    protein = atom(1, "a", 1)
    ligand = atom(2, "a", 2, record="HETATM", residue="LIG")
    record = "CONECT    1    2 9999"
    facts = projection_link_facts(protein + ligand + record + "\n", protein)
    assert facts["resolved_projection_boundary_pairs"] == []
    assert facts["raw_full_connection_records"] == [record]
    assert facts["unresolved_records_and_missing_annotations_not_qualified"]


def test_same_context_pair_is_not_a_protein_boundary_link():
    protein = atom(1, "A", 1)
    context = atom(2, "A", 2, record="HETATM", residue="LIG") + atom(3, "A", 3, record="HETATM", residue="LIG")
    facts = projection_link_facts(protein + context + "CONECT    2    3\n", protein)
    assert facts["resolved_full_explicit_pairs"] == 1
    assert facts["resolved_projection_boundary_pairs"] == []


def test_projection_changes_fail_closed():
    with pytest.raises(ValueError, match="changed original"):
        projection_link_facts(atom(1, "a", 1), atom(1, "a", 1, x=2.))
    with pytest.raises(ValueError, match="unknown atom"):
        projection_link_facts(atom(1, "a", 1), atom(1, "A", 1))


def test_retained_graph_loss_is_not_hidden_by_zero_projection_boundary_links():
    protein = atom(1, "a", 1) + atom(2, "a", 2)
    facts = projection_link_facts(protein + "CONECT    1    2\n", protein)
    assert facts["resolved_projection_boundary_pairs"] == []
    assert len(facts["lost_resolved_pairs_between_retained_atoms"]) == 1
    assert not facts["resolved_retained_pair_custody_verified"]


def test_extra_projection_link_is_reported_and_serial_changes_do_not_change_identity():
    protein = atom(1, "a", 1) + atom(2, "a", 2)
    facts = projection_link_facts(protein, protein + "CONECT    1    2\n")
    assert len(facts["new_resolved_projection_pairs"]) == 1
    assert not facts["resolved_retained_pair_custody_verified"]
    reserialized = atom(5, "a", 1) + atom(6, "a", 2) + "CONECT    5    6\n"
    facts = projection_link_facts(protein + "CONECT    1    2\n", reserialized)
    assert facts["resolved_retained_pair_custody_verified"]
    assert facts["resolved_projection_pairs"] == 1


def inputs(root):
    full, protein, source = [root / n for n in ("full.pdb", "protein.pdb", "source.cif")]
    protein.write_text(atom(1, "a", 1))
    full.write_text("REMARK 999 DVBfixer BIOLOGICAL ASSEMBLY 1\n" + protein.read_text())
    doc = gemmi.cif.Document()
    doc.add_copied_block(source_block())
    doc.write_file(str(source))
    ctx = {"entry_id": "TEST", "scope": "pending", "omitted_atoms": 0,
           "original_full_context_path": str(full), "original_projection_path": str(protein)}
    cons = {"entry_id": "TEST", "scope": "pending", "original_source_path": str(source),
            "original_source_sha256": digest(source)}
    paths = []
    for name, detail, bindings in (("context", ctx, {str(p): digest(p) for p in (full, protein)}),
                                   ("construct", cons, {str(source): digest(source)})):
        directory = root / name
        directory.mkdir()
        save(directory / "entry.json", detail)
        path = directory / "audit.json"
        save(path, {"rows": [{"entry_id": "TEST", "report": "entry.json"}], "bindings": bindings,
                    "reports_sha256": {"entry.json": digest(directory / "entry.json")}})
        paths.append(path)
    return paths


def test_complete_bound_audit_and_no_overwrite(tmp_path):
    context, construct = inputs(tmp_path)
    output = tmp_path / "output"
    audit(context, construct, output)
    data = json.loads((output / "connection-audit.json").read_text())
    assert data["summary"]["entries"] == 1
    # Label co-membership is kept separate from actual available coordinate pairs.
    assert data["summary"]["assembly_member_label_cross_boundary_annotations"] == 1
    assert data["summary"]["actual_resolved_projection_boundary_pairs"] == 0
    assert data["TSV_sha256"] == digest(output / "connection-review.tsv")
    with pytest.raises(FileExistsError):
        audit(context, construct, output)


@pytest.mark.parametrize("name", ["full.pdb", "source.cif", "context/entry.json"])
def test_changed_evidence_refused_before_output(tmp_path, name):
    context, construct = inputs(tmp_path)
    path = tmp_path / name
    path.write_text(path.read_text() + "\n")
    with pytest.raises(ValueError):
        audit(context, construct, tmp_path / "output")
    assert not (tmp_path / "output").exists()
