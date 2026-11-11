import json

import pytest
from scripts.prepare_mac_comparison_expansion import CASES, prepare
from scripts.prepare_mac_expansion_inputs import prepare as encode


@pytest.mark.parametrize("index", [0, 1])
def test_private_expansion_templates_contain_only_observations(tmp_path, index):
    case, output = tmp_path / "case", tmp_path / "input"
    prepare(CASES[index], case)
    receipt = encode(case, output)
    hidden = {i for _, first, last in CASES[index]["hidden"] for i in range(first, last+1)}
    for name in ("source.pdb", "boltz-2/template.pdb"):
        rows = (output / name).read_text().splitlines()
        assert all(int(line[22:26]) not in hidden for line in rows if line.startswith("ATOM  "))
    template = json.loads((output / "protenix-v1/template.json").read_text())[0]
    assert not hidden & {i+1 for i in template["queryIndices"]}
    assert len(template["queryIndices"]) == 72-len(hidden)
    for line in (output / "rfd3/conditioning.pdb").read_text().splitlines():
        if line.startswith("ATOM  ") and int(line[22:26]) in hidden:
            assert [float(line[i:i+8]) for i in (30, 38, 46)] == [0., 0., 0.]
    assert receipt["atomic_group"] == list(range(index+1))
    assert receipt["public_capability_enabled"] is False
    assert receipt["checkpoint_loaded"] is False
    with pytest.raises(FileExistsError):
        encode(case, output)
