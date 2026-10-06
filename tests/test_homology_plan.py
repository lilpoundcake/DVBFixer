import json
import shutil

import pytest

from dvbfixer.homology_plan import (
    materialize_template_plan,
    prepare_mosaic_diffusion_inputs,
)


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
    coverage = json.loads(
        (tmp_path / "work" / "selected_template_mosaic.coverage.json").read_text()
    )
    assert coverage["maskConvention"] == "zero-based-half-open"
    assert coverage["overlapPolicy"] == "earlier-template-wins"
    assert [item["templateId"] for item in coverage["groups"][0]["coverage"]] == [
        "left",
        "right",
    ]
    assert [item["pdbResidue"]["residueNumber"] for item in coverage["groups"][0]["coverage"]] == [
        "1",
        "2",
    ]


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


def test_mosaic_diffusion_inputs_preserve_ownership_and_chain_mapping(tmp_path):
    coverage = tmp_path / "selected_template_mosaic.coverage.json"
    coverage.write_text(json.dumps({
        "schemaVersion": 1,
        "groups": [{
            "targetChain": "VH",
            "pdbChain": "H",
            "targetSequence": "AGST",
            "coverage": [
                {
                    "targetIndex": index,
                    "targetResidue": residue,
                    "covered": index in {0, 3},
                    "templateId": "left" if index == 0 else "right" if index == 3 else None,
                    "templateResidue": residue if index in {0, 3} else None,
                    "pdbResidue": {
                        "chain": "H",
                        "residueNumber": str(index + 1),
                        "insertionCode": "",
                    },
                }
                for index, residue in enumerate("AGST")
            ],
        }],
    }))

    fasta, ownership = prepare_mosaic_diffusion_inputs(
        coverage,
        [("VH", "AGST")],
        tmp_path,
    )

    assert fasta.read_text() == ">H\nAGST\n"
    assert [(item.chain, item.target_index, item.template_id) for item in ownership] == [
        ("H", 0, "left"),
        ("H", 3, "right"),
    ]


def test_mosaic_diffusion_rejects_covered_substitution(tmp_path):
    coverage = tmp_path / "selected_template_mosaic.coverage.json"
    coverage.write_text(json.dumps({
        "schemaVersion": 1,
        "groups": [{
            "targetChain": "A",
            "pdbChain": "A",
            "targetSequence": "AGGG",
            "coverage": [
                {
                    "targetIndex": index,
                    "targetResidue": "A" if index == 0 else "G",
                    "covered": index != 2,
                    "templateId": "template" if index != 2 else None,
                    "templateResidue": "S" if index == 1 else "A" if index == 0 else "G" if index == 3 else None,
                    "pdbResidue": {
                        "chain": "A",
                        "residueNumber": str(index + 1),
                        "insertionCode": "",
                    },
                }
                for index in range(4)
            ],
        }],
    }))

    with pytest.raises(ValueError, match="covered substitutions"):
        prepare_mosaic_diffusion_inputs(coverage, [("A", "AGGG")], tmp_path)
