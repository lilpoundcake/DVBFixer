from scripts.audit_retained_1he8_segment_frames import backbone_edges, frame_screen

from dvbfixer.model.diffusion.contract import AtomIdentity


def test_explicit_graph_does_not_bridge_case_distinct_chains_or_numbering_gaps():
    atoms = [AtomIdentity("A", "82", "A", "C"), AtomIdentity("A", "83", "", "N"),
             AtomIdentity("a", "83", "", "N"), AtomIdentity("A", "84", "", "CA")]
    bonds = [(atoms[1], atoms[3])]
    assert backbone_edges(set(atoms), bonds) == tuple(bonds)
    assert backbone_edges(set(atoms), []) == ()


def test_segment_frame_reports_intersection_without_claiming_a_knot_or_gate():
    a, b = AtomIdentity("A", "1", "", "N"), AtomIdentity("A", "1", "", "CA")
    c, d = AtomIdentity("a", "1", "A", "N"), AtomIdentity("a", "1", "A", "CA")
    points = {a: (-1., 0., 0.), b: (1., 0., 0.), c: (0., -1., 0.), d: (0., 1., 0.)}
    edges = ((a, b), (c, d))
    result = frame_screen(points, edges, edges, frozenset({a, b}))
    assert len(result["alarms"]) == 1
    assert result["alarms"][0]["ownership"] == "generated-observed"
    assert not result["alarms"][0]["topological_defect_established"]
    assert not result["publication_gate"]
    # A truly absent bond is never invented to keep the intersection alarm.
    assert not frame_screen(points, edges[:1], edges[:1], frozenset())["alarms"]
