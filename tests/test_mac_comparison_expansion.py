import hashlib
from pathlib import Path

import pytest
from scripts.mac_reconstruction_comparison import modeller_alignment, target_axis
from scripts.prepare_mac_comparison_expansion import CASES, prepare

from dvbfixer.model.diffusion.contract import DiffusionRequest, GapKind


@pytest.mark.parametrize("index,expected_gaps", [(0, 1), (1, 2)])
def test_new_masks_preserve_observations_and_exact_ownership(tmp_path, index, expected_gaps):
    spec = CASES[index]
    fixture = Path(spec["fixture"])
    before = hashlib.sha256(fixture.read_bytes()).hexdigest()
    output = tmp_path / "case"
    manifest = prepare(spec, output)
    request = DiffusionRequest.from_json((output / "evaluation-mask.json").read_text())
    assert len(request.gaps) == expected_gaps
    assert set(request.fixed_atoms).isdisjoint(request.generated_atoms)
    assert manifest["source_admissible"] and manifest["reference_baseline_accepted"]
    assert manifest["sampling_executed"] is False and manifest["public_capability_enabled"] is False
    expected_hidden = {(chain, str(number)) for chain, first, last in spec["hidden"] for number in range(first, last+1)}
    source_lines = (output / "source.pdb").read_text().splitlines(keepends=True)
    assert all((line[21], line[22:26].strip()) not in expected_hidden for line in source_lines if line.startswith("ATOM  "))
    original = set(fixture.read_text().splitlines(keepends=True))
    assert all(line in original for line in source_lines if line.startswith("ATOM  "))
    assert {(a.chain, a.residue_number) for a in request.generated_atoms} == expected_hidden
    assert len(target_axis(request)["A"]) == 72
    if index == 0:
        assert request.gaps[0].gap_kind is GapKind.N_TERMINAL
        assert request.gaps[0].left_anchor is None
        assert modeller_alignment(request).splitlines()[2].startswith("----")
    assert hashlib.sha256(fixture.read_bytes()).hexdigest() == before
    with pytest.raises(FileExistsError):
        prepare(spec, output)


def test_partner_partial_occupancy_is_not_silently_repaired(tmp_path):
    with pytest.raises(ValueError, match="partial occupancy"):
        prepare(CASES[2], tmp_path / "refused")
    assert not (tmp_path / "refused/evaluation-mask.json").exists()


def test_save_failure_does_not_leave_partial_json(tmp_path):
    from scripts.mac_reconstruction_comparison import save
    path = tmp_path / "opaque.json"
    with pytest.raises(TypeError):
        save(path, {"opaque": object()})
    assert not path.exists()
