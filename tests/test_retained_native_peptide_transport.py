from types import SimpleNamespace

from scripts.preflight_retained_native_peptide_transport import junction_transport

from dvbfixer.model.diffusion.contract import AtomIdentity, ResidueIdentity


def test_native_graph_readback_distinguishes_source_mapping_from_covalent_transport():
    left, gap, right = (ResidueIdentity("a", str(n), "A") for n in (82, 83, 84))
    request = SimpleNamespace(gaps=(SimpleNamespace(left_anchor=left, right_anchor=right, generated_residues=(gap,)),))
    axis = (AtomIdentity("a", "82", "A", "C"), AtomIdentity("a", "83", "A", "N"),
            AtomIdentity("a", "83", "A", "C"), AtomIdentity("a", "84", "A", "N"))
    disconnected = junction_transport(request, axis, [(1, 2, 1)], ["A", "B", "B", "C"])
    assert disconnected["intended_internal_junctions"] == 2
    assert disconnected["native_bonds_present"] == 0 and disconnected["cross_native_chain_junctions"] == 2
    connected = junction_transport(request, axis, [(0, 1, 1), (1, 2, 1), (2, 3, 1)], ["A"]*4)
    assert connected["native_bonds_present"] == 2
    assert not connected["causality_of_raw_geometry_failure_established"]
