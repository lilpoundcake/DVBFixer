import json
import shutil

from dvbfixer.homology_plan import materialize_template_plan


def _atom(serial, residue, number, x):
    return (f"ATOM  {serial:5d}  CA  {residue} A{number:4d}    "
            f"{x:8.3f}   0.000   0.000  1.00 20.00           C  \n")


def test_template_plan_builds_one_fitted_mosaic(tmp_path, monkeypatch):
    left = tmp_path / "left.pdb"
    right = tmp_path / "right.pdb"
    left.write_text(_atom(1, "ALA", 1, 1) + _atom(2, "GLY", 2, 2))
    right.write_text(_atom(1, "ALA", 1, 31) + _atom(2, "GLY", 2, 32))

    def fake_fit(specs, output, fit_dir, **kwargs):
        fit_dir.mkdir(parents=True)
        fitted = []
        for index, spec in enumerate(specs, 1):
            source_text, chain = spec.rsplit(":", 1)
            target = fit_dir / f"template_{index}_{chain}_fit.pdb"
            shutil.copy2(source_text, target)
            fitted.append(target)
        output.write_text("test\n")
        return fitted

    monkeypatch.setattr("dvbfixer.homology_plan.run_biopython_superposition", fake_fit)
    plan = tmp_path / "plan.json"
    plan.write_text(json.dumps({
        "templates": [
            {"id": "left", "path": str(left), "chain": "A", "targetChain": "A"},
            {"id": "right", "path": str(right), "chain": "A", "targetChain": "A"},
        ],
        "alignmentGroups": [{
            "chainId": "A",
            "rows": [
                {"id": "A", "kind": "target", "sequence": "AG"},
                {"id": "left", "kind": "template", "templateId": "left", "sequence": "AG"},
                {"id": "right", "kind": "template", "templateId": "right", "sequence": "AG"},
            ],
            "masks": {"left": [{"start": 0, "end": 1}], "right": [{"start": 1, "end": 2}]},
            "maskModes": {"left": "ranges", "right": "ranges"},
        }],
    }))
    templates, alignment = materialize_template_plan(plan, tmp_path / "work")
    assert len(templates) == 1
    mosaic = (tmp_path / "work" / "selected_template_mosaic.pdb").read_text()
    assert "   1.000" in mosaic
    assert "  32.000" in mosaic
    assert ">P1;selected_template_mosaic" in open(alignment).read()


def test_vh_vl_groups_get_distinct_pdb_chain_ids(tmp_path):
    template = tmp_path / "complex.pdb"
    template.write_text(
        _atom(1, "ALA", 1, 1).replace(" A   1", " H   1") +
        _atom(2, "GLY", 1, 2).replace(" A   1", " L   1")
    )
    plan = tmp_path / "plan.json"
    plan.write_text(json.dumps({
        "templates": [
            {"id": "heavy", "path": str(template), "chain": "H", "targetChain": "VH"},
            {"id": "light", "path": str(template), "chain": "L", "targetChain": "VL"},
        ],
        "alignmentGroups": [
            {"chainId": "VH", "rows": [
                {"id": "VH", "kind": "target", "sequence": "A"},
                {"id": "heavy", "kind": "template", "templateId": "heavy", "sequence": "A"},
            ], "masks": {}, "maskModes": {"heavy": "all"}},
            {"chainId": "VL", "rows": [
                {"id": "VL", "kind": "target", "sequence": "G"},
                {"id": "light", "kind": "template", "templateId": "light", "sequence": "G"},
            ], "masks": {}, "maskModes": {"light": "all"}},
        ],
    }))
    templates, _ = materialize_template_plan(plan, tmp_path / "work")
    atom_lines = [line for line in open(templates[0]) if line.startswith("ATOM")]
    assert [line[21] for line in atom_lines] == ["H", "L"]
