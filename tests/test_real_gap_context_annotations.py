import gemmi
from scripts.classify_real_gap_context import declared_connections


def test_deposited_coordination_covalent_and_disulfide_classes_stay_distinct():
    document = gemmi.cif.read_string("""data_test
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
one metalc a SER 17 A MG ?
two covale a CYS 20 B LIG ?
three disulf a CYS 30 a CYS 40
""")
    connections = declared_connections(document.sole_block(), {"a"})
    assert [c["depositor_connection_type"] for c in connections] == ["metalc", "covale", "disulf"]
    assert [c["crosses_protein_nonprotein_boundary"] for c in connections] == [True, True, False]
    assert connections[0]["first"]["label_chain"] == "a"
    assert connections[0]["second"]["label_chain"] == "A"
    assert all("not established" in c["interpretation"] for c in connections)


def test_absent_connection_table_does_not_invent_bonds():
    document = gemmi.cif.read_string("data_test\n_entry.id TEST\n")
    assert declared_connections(document.sole_block(), {"A"}) == []
