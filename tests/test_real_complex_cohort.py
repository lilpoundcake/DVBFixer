"""Protocol and independent deposited-gap evidence, without checkpoint sampling."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


def script(name):
    path = Path(__file__).resolve().parents[1] / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def atom(chain, number, residue="ALA"):
    return f"ATOM      1  CA  {residue} {chain}{number:4d}       0.000   0.000   0.000  1.00 20.00           C  \n"


def missing(chain, number, residue="GLY"):
    return "REMARK 465" + " " * 5 + f"{residue} {chain} {number:5d} \n"


def request(chain="A"):
    return {"target_sequences": [{"chain": chain, "sequence": "AGA"}], "gaps": [{
        "chain": chain, "target_interval": {"start": 1, "stop": 2}, "gap_kind": "internal",
        "generated_residues": [{"chain": chain, "residue_number": "2", "insertion_code": ""}],
    }]}


def test_approved_cohort_keeps_ten_structures_and_default_seed():
    module = script("prepare_real_complex_cohort")
    assert len(module.SOURCES) == 10
    assert module.SEEDS == (7, 101, 102, 103, 104)
    assert module.PROFILE == "protpardelle-1c-cc94-mps-v1"
    assert all(len(digest) == 64 for digest in module.SOURCES.values())


def test_real_missing_annotation_not_an_artificial_deletion():
    module = script("audit_real_complex_gaps")
    source = atom("A", 1) + atom("A", 3) + missing("A", 2)
    assert module.audit(source, source, request())["passed"]
    with pytest.raises(ValueError, match="not unambiguously deposited"):
        module.audit(source.replace(missing("A", 2), ""), source, request())


def test_assembly_copy_maps_back_to_deposited_chain():
    module = script("audit_real_complex_gaps")
    source = atom("A", 1) + atom("A", 3) + missing("A", 2)
    assembled = atom("b", 1) + atom("b", 3)
    report = module.audit(source, assembled, request("b"))
    assert report["assembly_chain_sources"] == {"b": ["A"]}


def test_missing_annotation_sequence_must_match_target():
    module = script("audit_real_complex_gaps")
    source = atom("A", 1) + atom("A", 3) + missing("A", 2, "VAL")
    with pytest.raises(ValueError, match="sequence differ"):
        module.audit(source, source, request())


def test_chain_identifiers_remain_case_sensitive():
    module = script("audit_real_complex_gaps")
    source = atom("A", 1) + atom("A", 3) + missing("a", 2)
    with pytest.raises(ValueError, match="not unambiguously deposited"):
        module.audit(source, source, request())
