from types import SimpleNamespace

import numpy as np
import pytest
from scripts import audit_anchor_calibration_grid as audit

from dvbfixer.model.diffusion.contract import AtomIdentity, GapKind, ResidueIdentity
from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates


def line(serial, name, number, point, record="ATOM", resname="GLY"):
    return (f"{record:6s}{serial:5d} {name:^4s} {resname:3s} a{number:4d}    "
            + "".join(f"{v:8.3f}" for v in point)
            + f"{1.:6.2f}{10.:6.2f}          {name[0]:>2s}  \n")


def test_evaluation_view_removes_only_context_not_actual_anchor_motion(tmp_path):
    rows = []
    identities = []
    for number in (1, 2, 3):
        for name in ("N", "CA", "C", "O"):
            rows.append(line(len(rows)+1, name, number, (number*3., len(rows), 0.)))
            identities.append(AtomIdentity("a", str(number), "", name))
    source = "".join(row for row in rows if row[22:26].strip() != "2")+"END\n"
    residues = [ResidueIdentity("a", str(i)) for i in (1, 2, 3)]
    request = SimpleNamespace(target_sequences=(SimpleNamespace(chain="a", sequence="GGG"),),
        sequence_placements=(SimpleNamespace(chain="a", observed_target_indices=(0, 2),
                                            observed_residues=(residues[0], residues[2])),),
        gaps=(SimpleNamespace(chain="a", gap_kind=GapKind.INTERNAL, target_interval=SimpleNamespace(start=1, stop=2),
                              generated_residues=(residues[1],), left_anchor=residues[0], right_anchor=residues[2]),),
        fixed_atoms=tuple(a for a in identities if a.residue_number != "2"),
        generated_atoms=tuple(a for a in identities if a.residue_number == "2"))
    water = line(13, "O", 9, (12., 9., 3.), record="HETATM", resname="HOH")
    before = "".join(rows)+water+"END\n"
    anchor = AtomIdentity("a", "1", "", "CA")
    full = rewrite_generated_coordinates(before, {anchor: np.array([3.25, 1., 0.])})
    final = tmp_path / "full.pdb"
    final.write_text(full)
    view = audit.evaluation_view(request, source, final, {anchor})
    view_path = tmp_path / "view.pdb"
    view_path.write_text(view)
    complete, projection = audit.atoms(final), audit.atoms(view_path)
    assert len(complete) == 13 and len(projection) == 12
    assert all(row["xyz"] == complete[atom]["xyz"] for atom, row in projection.items())
    assert projection[audit.key(anchor)]["xyz"] == (3.25, 1., 0.)
    assert "HETATM" not in view
    audit.record_audit(before, full, {audit.key(anchor)})
    with pytest.raises(ValueError, match="record count"):
        audit.record_audit(before, full.replace(water, ""), {audit.key(anchor)})


def test_cli_baseline_difference_preserves_multiplicity_and_severity():
    warning = {"severity": "WARNING", "category": "clash", "chain": "a", "resid": "82A"}
    error = {**warning, "severity": "ERROR"}
    assert audit.finding_difference([warning], [warning]) == []
    assert audit.finding_difference([warning], [warning, warning]) == [warning]
    assert audit.finding_difference([warning], [error]) == [error]


def test_native_excursion_is_not_hidden_by_serialized_rounding():
    identities = {("a", "1", "", name) for name in ("N", "CA", "C", "O")}
    source = {a: {"xyz": (0., 0., 0.)} for a in identities}
    native = {a: (.100001, 0., 0.) for a in identities}
    serialized = {a: (.100, 0., 0.) for a in identities}
    assert not audit.movement(identities, source, native, (.1, .1, .5))[0]["within_limits"]
    assert audit.movement(identities, source, serialized, (.1, .1, .5))[0]["within_limits"]


def test_summary_distinguishes_missing_measurement_from_zero():
    from scripts.summarize_anchor_calibration import make_row, maximum

    assert maximum(None, "heavy_max") is None
    assert maximum([], "heavy_max") == 0.
    assert maximum([{"heavy_max": .5}], "heavy_max") == .5
    audited = {"entry": "case", "cell": {"id": "fixed", "window": 0, "envelope": None, "strength": None},
               "final_candidate_available": False, "rejection": ["global_chirality_failure"]}
    row = make_row({"shared_preparation_seconds": 2.}, {}, audited, {"elapsed_seconds": 3., "peak_rss_bytes": 100})
    assert row["charged_seconds"] == 5.
    assert row["max_anchor_heavy_displacement"] is None
    assert row["internal_generated_issue_count"] is None
    assert row["current_cli_errors"] is None
