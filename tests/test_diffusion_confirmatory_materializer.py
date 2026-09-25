"""Offline tests for confirmatory coordinate screening and materialization."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from types import ModuleType

import pytest

from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.scope import assess_diffusion_scope

ROOT = Path(__file__).resolve().parents[1]


def _load_script() -> ModuleType:
    name = "materialize_diffusion_confirmatory_cohort"
    path = ROOT / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _mmcif(
    *,
    residues: int = 16,
    models: tuple[int, ...] = (1,),
    residue_name: str = "ALA",
    auth_start: int = 1,
    heterogens_near_every_residue: bool = False,
    altloc_index: int | None = None,
    missing_indices: frozenset[int] = frozenset(),
    include_hydrogens: bool = False,
) -> bytes:
    rows: list[str] = []
    atom_id = 1
    one_letter = {"ALA": "A", "MSE": "M"}[residue_name]
    del one_letter
    for model in models:
        for index in range(residues):
            if index in missing_indices:
                continue
            base = index * 2.0
            atoms = [
                ("N", "N", base, 0.0, 0.0),
                ("CA", "C", base + 1.0, 0.0, 0.0),
                ("C", "C", base + 1.0, 1.0, 0.0),
                ("O", "O", base + 1.0, 2.0, 0.0),
                ("CB", "C", base + 1.0, 0.0, -1.0),
            ]
            if include_hydrogens:
                atoms.append(("H", "H", base - 0.5, 0.0, 0.0))
            for atom_name, element, x, y, z in atoms:
                altloc = "A" if altloc_index == index and atom_name == "CA" else "."
                rows.append(
                    f"ATOM {atom_id} {element} {atom_name} {altloc} {residue_name} L 1 "
                    f"{index + 1} ? {x:.3f} {y:.3f} {z:.3f} 1 10 ? "
                    f"{auth_start + index} {residue_name} a {atom_name} {model}"
                )
                atom_id += 1
            if heterogens_near_every_residue:
                rows.append(
                    f"HETATM {atom_id} ZN ZN . ZN H 2 . ? {base + 1.0:.3f} 0.0 1.0 "
                    f"1 10 2 {index + 1} ZN Z ZN {model}"
                )
                atom_id += 1
    header = """data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
"""
    return (header + "\n".join(rows) + "\n#\n").encode()


def _case(
    *, url: str = "https://files.rcsb.org/download/TEST.cif", length: int = 5
) -> dict[str, object]:
    return {
        "pdb_id": "TEST",
        "entity_id": "1",
        "polymer_entity_id": "TEST_1",
        "independence_group": "rcsb-30pct-representative:TEST_1",
        "label_asym_id": "L",
        "auth_asym_id": "a",
        "sequence": "A" * 16,
        "sequence_length": 16,
        "planned_gap_length": length,
        "request_sequence_basis": "observed-coordinate-sequence",
        "mmcif_url": url,
    }


def _cohort(path: Path, case: dict[str, object]) -> None:
    case = {**case, "screening_index": 0}
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "eligible_backends": ["protenix-v1"],
                "screening_pool_count": 1,
                "final_minimum_independence_groups": 1,
                "final_maximum_independence_groups": 1,
                "request_sequence_basis": "observed-coordinate-sequence",
                "screening_order_locked": True,
                "cases": [case],
            }
        )
        + "\n"
    )


def test_materializes_verified_workspace_and_resumes_without_network(tmp_path: Path) -> None:
    script = _load_script()
    cohort = tmp_path / "cohort.json"
    output = tmp_path / "out"
    _cohort(cohort, _case())
    calls: list[str] = []

    def fetch(url: str) -> bytes:
        calls.append(url)
        return _mmcif()

    report = script.materialize_cohort(cohort, output, fetch=fetch)

    assert report["summary"] == {
        "total": 1,
        "accepted": 1,
        "rejected": 0,
        "operational_errors": 0,
        "pending": 0,
        "selected_for_inference": 1,
        "reserve_accepted": 0,
        "final_minimum": 1,
        "final_maximum": 1,
        "accepted_count_sufficient": True,
        "screening_complete": True,
        "confirmatory_sample_valid": True,
    }
    record = report["cases"][0]
    workspace = output / "cases/test/workspace"
    request = DiffusionRequest.from_json((workspace / "request.json").read_text())
    source = (workspace / "input/normalized.pdb").read_bytes()
    assert assess_diffusion_scope(request, source).supported is True
    assert record["source"]["sha256"] == script._sha256(_mmcif())
    assert record["independence_group"] == "rcsb-30pct-representative:TEST_1"
    assert record["leakage_resolved"] is True
    assert record["selected_for_inference"] is True
    assert record["selection_status"] == "selected-for-inference"
    assert record["mask"]["length"] == 5
    assert record["mask"]["sequence"] == "AAAAA"
    assert (output / "cases/test/.complete.json").is_file()
    assert set(path.name for path in workspace.iterdir()) == {
        "input", "reference.pdb", "target.fasta", "request.json", "source.json"
    }
    assert calls == ["https://files.rcsb.org/download/TEST.cif"]

    def no_fetch(_url: str) -> bytes:
        raise AssertionError("completed case attempted a network request")

    resumed = script.materialize_cohort(cohort, output, fetch=no_fetch)
    assert resumed == report


def test_natural_entity_gaps_are_omitted_from_observed_coordinate_request(
    tmp_path: Path,
) -> None:
    script = _load_script()
    cohort = tmp_path / "cohort.json"
    output = tmp_path / "out"
    case = _case()
    case["sequence"] = "A" * 20
    case["sequence_length"] = 20
    _cohort(cohort, case)

    report = script.materialize_cohort(
        cohort,
        output,
        fetch=lambda _url: _mmcif(
            residues=20,
            missing_indices=frozenset({3, 4}),
        ),
    )

    assert report["summary"]["confirmatory_sample_valid"] is True
    workspace = output / "cases/test/workspace"
    request = DiffusionRequest.from_json((workspace / "request.json").read_text())
    source = (workspace / "input/normalized.pdb").read_bytes()
    provenance = json.loads((workspace / "source.json").read_text())
    assert request.target_sequences[0].sequence == "A" * 18
    assert request.sequence_placements[0].target_length == 18
    assert assess_diffusion_scope(request, source).supported is True
    assert provenance["request_sequence_basis"] == "observed-coordinate-sequence"
    assert provenance["observed_entity_indices_zero_based"] == [
        0, 1, 2, *range(5, 20)
    ]
    assert provenance["observed_coordinate_residue_count"] == 18
    assert all(
        residue.residue_number not in {"4", "5"}
        for residue in request.sequence_placements[0].observed_residues
    )


def test_explicit_hydrogens_are_excluded_from_workspace_atom_contract(
    tmp_path: Path,
) -> None:
    script = _load_script()
    cohort = tmp_path / "cohort.json"
    output = tmp_path / "out"
    _cohort(cohort, _case())

    report = script.materialize_cohort(
        cohort,
        output,
        fetch=lambda _url: _mmcif(include_hydrogens=True),
    )

    assert report["summary"]["confirmatory_sample_valid"] is True
    workspace = output / "cases/test/workspace"
    request = DiffusionRequest.from_json((workspace / "request.json").read_text())
    normalized = (workspace / "input/normalized.pdb").read_bytes()
    reference = (workspace / "reference.pdb").read_text().splitlines()

    def atom_identities(lines: list[str]) -> set[tuple[str, str, str, str]]:
        return {
            (line[21], line[22:26].strip(), line[26].strip(), line[12:16].strip())
            for line in lines
            if line.startswith(("ATOM  ", "HETATM"))
        }

    normalized_lines = normalized.decode().splitlines()
    assert all(script._is_heavy(line) for line in normalized_lines if line.startswith("ATOM  "))
    assert all(script._is_heavy(line) for line in reference if line.startswith("ATOM  "))
    assert atom_identities(normalized_lines) == {
        (atom.chain, atom.residue_number, atom.insertion_code, atom.atom_name)
        for atom in request.fixed_atoms
    }
    assert atom_identities(reference) == {
        (atom.chain, atom.residue_number, atom.insertion_code, atom.atom_name)
        for atom in (*request.fixed_atoms, *request.generated_atoms)
    }
    assert assess_diffusion_scope(request, normalized).supported is True


def test_report_marks_accepted_count_sufficiency_only_after_full_screen() -> None:
    script = _load_script()
    records = [
        {"status": "accepted", "selected_for_inference": True},
        {"status": "accepted", "selected_for_inference": True},
        {"status": "rejected"},
    ]

    sufficient = script._report(
        Path("cohort.json"),
        "0" * 64,
        records,
        expected_count=3,
        final_minimum=2,
        final_maximum=3,
    )
    insufficient = script._report(
        Path("cohort.json"),
        "0" * 64,
        records,
        expected_count=3,
        final_minimum=3,
        final_maximum=3,
    )

    assert sufficient["summary"]["accepted_count_sufficient"] is True
    assert sufficient["summary"]["confirmatory_sample_valid"] is True
    assert insufficient["summary"]["accepted_count_sufficient"] is False
    assert insufficient["summary"]["confirmatory_sample_valid"] is False


def test_final_selection_uses_locked_order_and_marks_reserve(tmp_path: Path) -> None:
    script = _load_script()
    records = [
        {
            "case_id": f"case{index}",
            "screening_index": index,
            "status": "accepted",
            "selected_for_inference": False,
            "selection_status": "awaiting-full-screen",
        }
        for index in range(4)
    ]
    for record in records:
        (tmp_path / "cases" / record["case_id"]).mkdir(parents=True)

    finalized = script._finalize_inference_selection(
        records,
        tmp_path,
        "0" * 64,
        expected_count=4,
        final_maximum=3,
    )

    assert [record["selected_for_inference"] for record in finalized] == [
        True, True, True, False
    ]
    assert finalized[-1]["selection_status"] == "reserve-not-selected"
    for record in finalized:
        stored = json.loads(
            (tmp_path / "cases" / record["case_id"] / "accepted.json").read_text()
        )
        assert stored["selected_for_inference"] == record["selected_for_inference"]
    report = script._report(
        Path("cohort.json"),
        "0" * 64,
        finalized,
        expected_count=4,
        final_minimum=2,
        final_maximum=3,
    )
    assert report["summary"]["accepted"] == 4
    assert report["summary"]["selected_for_inference"] == 3
    assert report["summary"]["reserve_accepted"] == 1
    assert report["summary"]["confirmatory_sample_valid"] is True


def test_rejects_multiple_models_with_atomic_marker(tmp_path: Path) -> None:
    script = _load_script()
    cohort = tmp_path / "cohort.json"
    output = tmp_path / "out"
    _cohort(cohort, _case())

    report = script.materialize_cohort(
        cohort,
        output,
        fetch=lambda _url: _mmcif(models=(1, 2)),
    )

    record = report["cases"][0]
    assert record["status"] == "rejected"
    assert record["reason_codes"] == ["multiple-models"]
    assert not (output / "cases/test/workspace").exists()
    assert (output / "cases/test/.complete.json").is_file()


def test_rejects_noncanonical_target_and_non_https_source(tmp_path: Path) -> None:
    script = _load_script()
    cohort = tmp_path / "cohort.json"

    _cohort(cohort, _case())
    noncanonical = script.materialize_cohort(
        cohort,
        tmp_path / "noncanonical",
        fetch=lambda _url: _mmcif(residue_name="MSE"),
    )
    assert noncanonical["cases"][0]["reason_codes"] == ["noncanonical-residue"]

    _cohort(cohort, _case(url="http://example.test/TEST.cif"))
    insecure = script.materialize_cohort(
        cohort,
        tmp_path / "insecure",
        fetch=lambda _url: (_ for _ in ()).throw(AssertionError("must not fetch")),
    )
    assert insecure["cases"][0]["reason_codes"] == ["non-https-source"]


def test_rejects_altloc_and_heterogen_context(tmp_path: Path) -> None:
    script = _load_script()
    cohort = tmp_path / "cohort.json"
    _cohort(cohort, _case())

    altloc = script.materialize_cohort(
        cohort,
        tmp_path / "altloc",
        fetch=lambda _url: _mmcif(altloc_index=8),
    )
    assert altloc["cases"][0]["reason_codes"] == [
        "alternate-location-ambiguity-in-target-chain"
    ]

    heterogen = script.materialize_cohort(
        cohort,
        tmp_path / "heterogen",
        fetch=lambda _url: _mmcif(heterogens_near_every_residue=True),
    )
    record = heterogen["cases"][0]
    assert record["reason_codes"] == ["no-clean-internal-mask"]
    assert record["evidence"]["local_rejection_counts"]["heteroatom-context-near-mask"] > 0


def test_rejects_covalent_link_context_for_every_mask(tmp_path: Path) -> None:
    import gemmi

    script = _load_script()
    source = tmp_path / "links.cif"
    source.write_bytes(_mmcif())
    case = _case()
    structure = script._load_structure(source)
    _matches, observed = script._select_instance(structure, case)
    for residue_number in range(1, 17):
        connection = gemmi.Connection()
        connection.partner1.chain_name = "a"
        connection.partner1.res_id.seqid = gemmi.SeqId(residue_number, " ")
        connection.partner1.atom_name = "CB"
        connection.partner2.chain_name = "Z"
        connection.partner2.res_id.seqid = gemmi.SeqId(1, " ")
        connection.partner2.atom_name = "ZN"
        structure.connections.append(connection)

    with pytest.raises(script.ScreeningRejection) as error:
        script._select_mask(structure, observed, case)

    assert error.value.code == "no-clean-internal-mask"
    assert error.value.evidence["local_rejection_counts"][
        "covalent-link-context-near-mask"
    ] > 0


def test_fixed_column_overflow_is_scientific_rejection(tmp_path: Path) -> None:
    script = _load_script()
    cohort = tmp_path / "cohort.json"
    _cohort(cohort, _case())

    report = script.materialize_cohort(
        cohort,
        tmp_path / "overflow",
        fetch=lambda _url: _mmcif(auth_start=10000),
    )

    assert report["cases"][0]["reason_codes"] == ["pdb-fixed-column-overflow"]
