from dataclasses import dataclass
from types import SimpleNamespace

import numpy as np
import pytest
from scripts.audit_retained_6z38_restoration import motion
from scripts.run_retained_6z38_restoration import compare_findings, observation_audit

from dvbfixer.model.diffusion.contract import AtomIdentity


def observations():
    fixed = AtomIdentity("A", "4", "", "CA")
    generated = AtomIdentity("A", "1", "", "CA")
    request = SimpleNamespace(fixed_atoms=(fixed,), generated_atoms=(generated,))
    source = {fixed: ("THR", "C", np.array([1., 2., 3.]))}
    native = {**source, generated: ("THR", "C", np.array([2., 3., 4.]))}
    return request, source, native, fixed


def test_native_exact_observations_pass():
    request, source, native, _ = observations()
    assert observation_audit(request, source, native)["native_observations_exact_at_PDB_precision"]


def test_even_rigid_native_motion_is_not_hidden_by_fitting():
    request, source, native, fixed = observations()
    native[fixed] = ("THR", "C", source[fixed][2] + [0.001, 0., 0.])
    report = observation_audit(request, source, native)
    assert not report["native_observations_exact_at_PDB_precision"]
    assert report["max_motion_angstrom"] == pytest.approx(.001)


def test_native_missing_atom_fails_closed():
    request, source, native, fixed = observations()
    del native[fixed]
    with pytest.raises(ValueError, match="ownership"):
        observation_audit(request, source, native)


def test_native_chemistry_change_is_not_coordinate_preservation():
    request, source, native, fixed = observations()
    native[fixed] = ("ALA", "C", source[fixed][2])
    with pytest.raises(ValueError, match="chemistry"):
        observation_audit(request, source, native)


@dataclass
class Finding:
    category: str
    measurement: float


def test_diagnostics_preserve_duplicates_and_changed_measurements():
    baseline = [Finding("source", 1.)]
    result = compare_findings(baseline, [*baseline, *baseline, Finding("source", 2.)])
    assert len(result["unchanged_source_findings"]) == 1
    assert len(result["new_or_changed_findings"]) == 2


def test_independent_readback_reports_permutations_without_excusing_custody():
    first, second = ("A", "5", "", "OD1"), ("A", "5", "", "OD2")
    source = {first: {"xyz": (0., 0., 0.)}, second: {"xyz": (1., 0., 0.)}}
    native = {first: source[second], second: source[first]}
    result = motion(source, native)
    assert result["changed_named_observations"] == 2
    assert result["residues_with_exact_coordinate_permutations"] == [["A", "5", ""]]
    assert result["coordinate_permutation_is_not_named_atom_custody"]
