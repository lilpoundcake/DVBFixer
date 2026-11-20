import json
from dataclasses import asdict
from types import SimpleNamespace

import numpy as np
import pytest
from openmm.app import PDBFile
from openmm.unit import nanometer
from scripts import refine_1non_generated_only as refine
from scripts.audit_1non_generated_refinement import audit_prepared

from dvbfixer.model.diffusion.contract import AtomIdentity


def inputs(monkeypatch):
    fixed = tuple(AtomIdentity("A", str(i), "", "CA") for i in range(5040))
    generated = tuple(AtomIdentity("d", str(i), "", "CA") for i in range(556))
    request = SimpleNamespace(fixed_atoms=fixed, generated_atoms=generated, target_sequences=())
    def key(a):
        return a.chain, a.residue_number, a.insertion_code, a.atom_name
    source = {key(a): {"line": "ATOM  fixed", "element": "C"} for a in fixed}
    proposal = {**source, **{key(a): {"line": "ATOM  generated", "element": "C"} for a in generated}}
    record = SimpleNamespace(added_atoms=[{"chain": a.chain, "resid": a.residue_number,
                                          "icode": a.insertion_code, "atom": a.atom_name} for a in generated])
    monkeypatch.setattr(refine, "target_axis", lambda request: {})
    return request, source, proposal, record


def test_exact_axis_and_dat_required(monkeypatch):
    request, source, proposal, record = inputs(monkeypatch)
    refine.validate_inputs(request, source, proposal, record)
    record.added_atoms[0]["chain"] = "D"
    with pytest.raises(ValueError, match="ownership"):
        refine.validate_inputs(request, source, proposal, record)


def test_observed_change_and_context_refused(monkeypatch):
    request, source, proposal, record = inputs(monkeypatch)
    a = next(iter(source))
    proposal[a] = {"line": "ATOM  moved", "element": "C"}
    with pytest.raises(ValueError, match="observed"):
        refine.validate_inputs(request, source, proposal, record)
    proposal[a] = source[a]
    proposal[("W", "1", "", "O")] = {"line": "HETATM water", "element": "O"}
    with pytest.raises(ValueError, match="complete retained"):
        refine.validate_inputs(request, source, proposal, record)


def test_phase_budget_is_bounded():
    assert [p[-1] for p in refine.PHASES] == [100, 75, 75, 1000]
    assert refine.PHASES[-1][1] == 1 | (1 << 31)


def test_independent_prepared_heavy_custody_detects_motion(tmp_path):
    case, prepared = tmp_path / "case", tmp_path / "prepared"
    case.mkdir()
    prepared.mkdir()
    text = "ATOM      1  CA  GLY d  82       1.000   2.000   3.000  1.00 20.00           C  \nEND\n"
    path = case / "proposal.pdb"
    path.write_text(text)
    raw = PDBFile(str(path))
    a = next(raw.topology.atoms())
    row = {"index": 0, "identity": asdict(refine.identity(a)), "element": "C", "temporary_addition": False}
    (prepared / "temporary-axis.json").write_text(json.dumps({"atoms": [row]}))
    points = np.asarray(raw.positions.value_in_unit(nanometer))
    np.savez_compressed(prepared / "initial-state.npz", positions_nm=points)
    assert audit_prepared(case, prepared)["raw_to_prepared_heavy_bits_exact"]
    points[0, 0] += 1.e-10
    np.savez_compressed(prepared / "initial-state.npz", positions_nm=points)
    with pytest.raises(ValueError, match="custody"):
        audit_prepared(case, prepared)
