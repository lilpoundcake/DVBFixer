from scripts.reaudit_anchor_diagnostics import classify, difference, python_report


def contact(version="2", reverse=False):
    axis = [{"chain": "B", "resid": "161", "icode": "A", "atom": "OG"},
            {"chain": "b", "resid": "12", "icode": "", "atom": "O"}]
    return {"category": "clash", "severity": "ERROR", "message": "close contact",
            "details": {"atoms": list(reversed(axis)) if reverse else axis,
                        "engine": {"version": version}}}


def test_engine_metadata_does_not_create_new_contact():
    assert difference([contact("1")], [contact(reverse=True)]) == []
    assert len(difference([], [contact()])) == 1
    assert len(difference([contact()], [contact(), contact()])) == 1


def test_classification_keeps_chain_case_and_insertion_codes():
    assert classify(contact(), set(), set())["only_immutable_atoms"]
    assert not classify(contact(), {("B", "161", "", "OG")}, set())["touches_generated"]
    assert classify(contact(), {("B", "161", "A", "OG")}, set())["touches_generated"]
    assert not classify(contact(), {("B", "12", "", "O")}, set())["touches_generated"]
    assert classify(contact(), set(), {("b", "12", "", "O")})["touches_observed_anchor"]


def test_python_report_records_fixed_ser_water_clash(tmp_path):
    from openmm import Vec3
    from openmm.app import Topology, element
    from openmm.unit import nanometer

    topology = Topology()
    chain = topology.addChain("B")
    ser = topology.addResidue("SER", chain, "161")
    water = topology.addResidue("HOH", chain, "12")
    topology.addAtom("OG", element.oxygen, ser)
    topology.addAtom("O", element.oxygen, water)
    report = python_report(topology, [Vec3(0, 0, 0), Vec3(.1658382, 0, 0)]*nanometer,
                           tmp_path / "report.json")
    assert report["steric"]["engine_version"] == "2"
    assert report["steric"]["status"] == "complete"
    clashes = [f for f in report["findings"] if f["category"] == "clash"]
    assert len(clashes) == 1 and clashes[0]["severity"] == "ERROR"
