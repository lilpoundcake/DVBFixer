"""Residue-coherent altloc policies, metadata, and non-destructive CLI."""

import hashlib
import json
from pathlib import Path

import pytest

from dvbfixer.altloc_selection import AltlocError, select_altlocs
from dvbfixer.altlocs import main


def atom(serial, name="CB", label="A", occ="0.60", chain="C", number=82, icode="A", resname="SER"):
    return (f"ATOM  {serial:5d} {name:^4}{label or ' '}{resname:>3} {chain}{number:4d}{icode or ' '}   "
            f"{serial:8.3f}{2.:8.3f}{3.:8.3f}{occ:>6} 20.00           C  \n")


def source():
    return atom(1, "N", "") + atom(2, label="B", occ="0.40") + atom(3)


def select(text=None, **kwargs):
    text = text or source()
    return select_altlocs(text, source_sha256=hashlib.sha256(text.encode()).hexdigest(), **kwargs)


def test_occupancy_is_whole_state_and_preserves_coordinates():
    text = source() + atom(4, "OG", "A") + atom(5, "OG", "B", "0.90")
    output, report = select(text)
    assert report["residues"][0]["selected"] == "B"  # mean .65 versus .60
    assert [int(line[6:11]) for line in output.splitlines()] == [1, 2, 5]
    assert all(line[16] == " " for line in output.splitlines())
    assert output.splitlines()[1][30:] == text.splitlines()[1][30:]
    assert "0.40" in output  # do not normalize occupancy to 1


def test_first_preferred_manual_override_and_missing_preferred():
    assert select(policy="first")[1]["residues"][0]["selected"] == "B"
    assert select(policy="prefer", prefer="A")[1]["residues"][0]["selected"] == "A"
    assert select(policy="prefer", prefer="Z")[0] is None
    assert select(policy="prefer", prefer="Z", fallback="first")[1]["residues"][0]["selected"] == "B"
    _, report = select(policy="manual")
    key = report["residues"][0]["key"]
    assert report["unresolved"] == [key]
    assert select(policy="occupancy", overrides={key: "B"})[1]["residues"][0]["selected"] == "B"
    with pytest.raises(AltlocError, match="absent"):
        select(overrides={key: "Z"})


def test_ties_missing_occupancy_and_partial_inventory_are_visible():
    text = atom(1, label="B") + atom(2)
    _, report = select(text)
    assert report["residues"][0]["selected"] == "B"
    assert any("Equal occupancy" in warning for warning in report["warnings"])
    output, report = select(atom(1, occ="") + atom(2, label="B", name="OG"))
    assert output is None
    assert any("inventories differ" in warning for warning in report["warnings"])


def test_chain_case_insertion_model_occurrences_stay_distinct():
    text = "MODEL        1\n" + source() + atom(4, chain="c") + atom(5, number=82, icode="B") + "TER\n" + atom(6) + "ENDMDL\nMODEL        2\n" + source() + "ENDMDL\n"
    _, report = select(text)
    assert len(report["residues"]) == 5
    assert len({row["key"] for row in report["residues"]}) == 5
    assert {row["model"] for row in report["residues"]} == {1, 2}
    assert {row["occurrence"] for row in report["residues"]} == {1, 2}


def test_shared_collision_rejected_and_microheterogeneity_selected_whole():
    with pytest.raises(AltlocError, match="collision"):
        select(atom(1, label="") + atom(2))
    output, report = select(atom(1, resname="SER") + atom(2, label="B", resname="CYS", occ="0.30"))
    assert "SER" in output and "CYS" not in output
    assert len(report["residues"][0]["states"]) == 2


def test_anisou_conect_and_link_are_reconciled_without_serial_reassignment():
    text = source()
    anisou = "ANISOU" + text.splitlines()[2][6:27] + " " + "  10000" * 6 + "       C  \n"
    link = list(" " * 80)
    link[:6] = "LINK  "
    link[12:27] = text.splitlines()[1][12:27]  # removed B
    link[42:57] = text.splitlines()[0][12:27]
    text += anisou + "CONECT    1    2    3\n" + "".join(link) + "\nEND\n"
    output, report = select(text)
    assert "CONECT    1    3" in output
    assert "LINK" not in output
    assert "ANISOU" in output and next(l for l in output.splitlines() if l.startswith("ANISOU"))[16] == " "
    assert report["warnings"]


def test_cli_digest_plan_and_source_safety(tmp_path):
    path = tmp_path / "input.pdb"
    path.write_text(source())
    report_path = tmp_path / "analysis.json"
    main([str(path), "--policy", "manual", "--report-json", str(report_path)])
    report = json.loads(report_path.read_text())
    key = report["residues"][0]["key"]
    plan = tmp_path / "plan.json"
    plan.write_text(json.dumps({"schemaVersion": 1, "sourceSha256": report["sourceSha256"], "selections": {key: "B"}}))
    output = tmp_path / "output.pdb"
    main([str(path), "--selection-plan", str(plan), "-o", str(output)])
    assert "   2 " in output.read_text()
    assert path.read_text() == source()
    with pytest.raises(SystemExit, match="already exists"):
        main([str(path), "-o", str(output)])
    with pytest.raises(SystemExit, match="different"):
        main([str(path), "-o", str(path)])
    path.write_text(source() + "END\n")
    with pytest.raises(SystemExit, match="digest mismatch"):
        main([str(path), "--selection-plan", str(plan)])


def test_cli_cif_normalization(tmp_path):
    gemmi = pytest.importorskip("gemmi")
    pdb = gemmi.read_pdb_string(source())
    cif = tmp_path / "input.cif"
    pdb.make_mmcif_document().write_file(str(cif))
    output = tmp_path / "selected.pdb"
    main([str(cif), "-o", str(output)])
    assert output.exists()
    assert all(line[16] == " " for line in output.read_text().splitlines() if line.startswith("ATOM"))


def test_retained_link_labels_are_cleared_and_coordinates_not_rewritten():
    text = source()
    link = list(" " * 80)
    link[:6] = "LINK  "
    link[12:27] = text.splitlines()[2][12:27]  # retained A
    link[42:57] = text.splitlines()[0][12:27]
    output, _ = select(text + "".join(link) + "\n")
    record = next(line for line in output.splitlines() if line.startswith("LINK"))
    assert record[16] == record[46] == " "
    assert record[12:16].strip() == "CB"


def test_dangling_conect_does_not_retarget_serial_and_global_models_require_agreement():
    output, report = select(source() + "CONECT    1  999\n")
    assert "CONECT" not in output
    assert any("dangling" in warning for warning in report["warnings"])
    text = "MODEL        1\n" + source() + "ENDMDL\nMODEL        2\n" + source() + "ENDMDL\nCONECT    1    3\n"
    _, inventory = select(text, policy="manual")
    keys = [row["key"] for row in inventory["residues"]]
    with pytest.raises(AltlocError, match="incompatible selections across MODELs"):
        select(text, overrides={keys[0]: "A", keys[1]: "B"})


def test_cli_dry_run_never_writes_structure_and_duplicate_plan_keys_fail(tmp_path):
    path = tmp_path / "input.pdb"
    path.write_text(source())
    output = tmp_path / "selected.pdb"
    report = tmp_path / "report.json"
    main([str(path), "--dry-run", "-o", str(output), "--report-json", str(report)])
    assert not output.exists()
    assert json.loads(report.read_text())["outputWritten"] is False
    plan = tmp_path / "plan.json"
    plan.write_text('{"schemaVersion":1,"schemaVersion":1}')
    with pytest.raises(SystemExit, match="duplicate JSON key"):
        main([str(path), "--selection-plan", str(plan), "-o", str(output)])
    assert not output.exists()


def test_report_publication_failure_rolls_back_only_new_output(tmp_path, monkeypatch):
    import dvbfixer.altlocs as cli

    path = tmp_path / "input.pdb"
    path.write_text(source())
    output = tmp_path / "selected.pdb"
    report = tmp_path / "report.json"
    original = cli._write_new

    def fail_report(target, content):
        if target == report:
            raise OSError("injected report failure")
        return original(target, content)

    monkeypatch.setattr(cli, "_write_new", fail_report)
    with pytest.raises(SystemExit, match="injected report failure"):
        main([str(path), "-o", str(output), "--report-json", str(report)])
    assert not output.exists()
    assert path.read_text() == source()


@pytest.mark.parametrize("name", ["protein_context/1A3N.pdb", "regressions/2VLN.pdb", "regressions/2VLQ.pdb", "regressions/1DQJ.pdb"])
def test_tracked_structures_select_states_preserve_coordinates_and_remove_altlocs(name):
    path = Path(__file__).parent / "fixtures" / name
    text = path.read_text()
    output, report = select(text)
    assert not report["unresolved"]
    assert report["residues"]
    before = {int(line[6:11]): line for line in text.splitlines() if line.startswith(("ATOM  ", "HETATM"))}
    seen = set()
    for line in output.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        assert line[16] == " "
        identity = (line[21:27], line[12:16])
        assert identity not in seen
        seen.add(identity)
        assert line[30:] == before[int(line[6:11])][30:]
