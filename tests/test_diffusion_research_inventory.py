"""Validate the predeclared diffusion research inventory and benchmark policy."""

from __future__ import annotations

import hashlib
import tomllib
from pathlib import Path
from typing import Any

from dvbfixer.model.diffusion.rfdiffusion_v1 import RFDIFFUSION_ENVIRONMENT_SHA256

REPO_ROOT = Path(__file__).resolve().parent.parent
INVENTORY = REPO_ROOT / "docs/research/diffusion-gap-reconstruction-inventory.toml"
FIXTURE_MANIFEST = REPO_ROOT / "tests/fixtures/MANIFEST.sha256"
RFDIFFUSION_ENVIRONMENT = REPO_ROOT / "deploy/rfdiffusion-v1/environment.yml"

CANONICAL_RESIDUES = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


def _load_inventory() -> dict[str, Any]:
    with INVENTORY.open("rb") as handle:
        return tomllib.load(handle)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fixture_manifest() -> dict[str, str]:
    entries: dict[str, str] = {}
    for line in FIXTURE_MANIFEST.read_text().splitlines():
        digest, relative = line.split("  ", 1)
        entries[relative] = digest
    return entries


def _fasta_sequences(path: Path) -> dict[str, str]:
    sequences: dict[str, list[str]] = {}
    current: str | None = None
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            current = line[1:].split(maxsplit=1)[0]
            assert current not in sequences, f"duplicate FASTA record: {current}"
            sequences[current] = []
            continue
        assert current is not None, f"sequence before FASTA header in {path}"
        sequences[current].append(line)
    return {name: "".join(parts) for name, parts in sequences.items()}


def _pdb_residues(path: Path, chain: str) -> list[tuple[str, str, str, set[str]]]:
    residues: list[tuple[str, str, str, set[str]]] = []
    keys: dict[tuple[str, str, str], int] = {}
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM  ") or len(line) < 27 or line[21] != chain:
            continue
        key = (line[22:26].strip(), line[26], line[17:20].strip())
        index = keys.get(key)
        if index is None:
            keys[key] = len(residues)
            residues.append((*key, {line[16]}))
        else:
            residues[index][3].add(line[16])
    return residues


def _identity(number: str, icode: str, residue_name: str, chain: str) -> str:
    encoded_icode = icode if icode.strip() else "_"
    return f"{chain}:{number}:{encoded_icode}:{residue_name}"


def test_diffusion_inventory_declares_versioned_scope_and_thresholds() -> None:
    inventory = _load_inventory()

    assert inventory["schema_version"] == 1
    assert inventory["status"] == "research"
    assert inventory["scope"]["public_cli_enabled"] is False
    assert inventory["scope"]["default_backend"] == "modeller"
    assert inventory["scope"]["automatic_fallback"] is False
    assert inventory["scope"]["initial_gap_lengths"] == [3, 12]

    thresholds = inventory["thresholds"]
    assert thresholds["version"] == 1
    assert thresholds["fixed_heavy_atom_rmsd_angstrom_max"] == 0.01
    assert thresholds["fixed_heavy_atom_displacement_angstrom_max"] == 0.03
    assert thresholds["detectable_d_ca_max"] == 0
    assert thresholds["junction_cn_angstrom_min"] == 1.20
    assert thresholds["junction_cn_angstrom_max"] == 1.45
    assert thresholds["generated_backbone_break_angstrom_max"] == 1.80
    assert thresholds["same_seed_repeat_rmsd_angstrom_max"] == 0.01


def test_diffusion_benchmark_cases_match_reviewed_fixtures() -> None:
    inventory = _load_inventory()
    fixture_manifest = _fixture_manifest()
    case_ids: set[str] = set()
    strata: set[str] = set()

    for case in inventory["benchmark_cases"]:
        case_id = case["id"]
        assert case_id not in case_ids
        case_ids.add(case_id)
        strata.add(case["stratum"])

        fixture_relative = case["fixture"]
        fixture = REPO_ROOT / fixture_relative
        assert fixture.is_file(), f"missing benchmark fixture for {case_id}"
        assert fixture_manifest[fixture_relative] == case["fixture_sha256"]
        assert _sha256(fixture) == case["fixture_sha256"]

        sequence_relative = case["sequence_fixture"]
        if sequence_relative:
            sequence_fixture = REPO_ROOT / sequence_relative
            assert fixture_manifest[sequence_relative] == case["sequence_fixture_sha256"]
            assert _sha256(sequence_fixture) == case["sequence_fixture_sha256"]

        assert case["leakage_status"] in {"resolved", "unresolved", "not_applicable"}
        assert case["leakage_reason"]

        if case["classification"] != "supported":
            assert case["classification"] == "unsupported"
            assert case["unsupported_reason"]
            continue

        start, stop = case["target_interval_zero_based_half_open"]
        assert 3 <= stop - start <= 12
        assert case["left_anchor"] and case["right_anchor"]
        assert len(case["generated_residues"]) == stop - start

        sequences = _fasta_sequences(REPO_ROOT / sequence_relative)
        fasta_key = f"8CZ8_{case['target_chain']}"
        sequence = sequences[fasta_key]
        assert sequence[start:stop] == case["expected_sequence"]

        residues = _pdb_residues(fixture, case["target_chain"])
        selected = residues[start:stop]
        assert len(selected) == stop - start
        assert all(altlocs == {" "} for _, _, _, altlocs in residues[start - 1:stop + 1])
        assert "".join(CANONICAL_RESIDUES[name] for _, _, name, _ in selected) == case["expected_sequence"]
        assert [
            _identity(number, icode, name, case["target_chain"])
            for number, icode, name, _ in selected
        ] == case["generated_residues"]

        left = residues[start - 1]
        right = residues[stop]
        assert _identity(*left[:3], case["target_chain"]) == case["left_anchor"]
        assert _identity(*right[:3], case["target_chain"]) == case["right_anchor"]

        pdb_text = fixture.read_text()
        assert pdb_text.count("\nMODEL ") + int(pdb_text.startswith("MODEL ")) <= 1
        assert not any(line.startswith("HETATM") for line in pdb_text.splitlines())
        assert not any(line.startswith(("LINK  ", "CONECT", "SSBOND")) for line in pdb_text.splitlines())

    assert {"withheld-internal-3-5", "withheld-internal-6-12"} <= strata
    assert {"terminal-one-anchor", "multiple-models", "retained-heterogens"} <= strata


def test_diffusion_engine_inventory_is_fail_closed() -> None:
    inventory = _load_inventory()
    engine_ids: set[str] = set()

    for engine in inventory["engines"]:
        engine_id = engine["id"]
        assert engine_id not in engine_ids
        engine_ids.add(engine_id)
        assert engine["repository"].startswith("https://github.com/")
        assert len(engine["revision"]) == 40
        int(engine["revision"], 16)
        assert engine["source_license"]
        assert engine["source_license_url"].startswith("https://github.com/")
        assert engine["training_cutoff"]
        assert engine["hardware"]
        assert engine["distribution_status"] in {"blocked", "comparator-only"}
        assert engine["blocking_reason"]

        if engine["distribution_status"] == "blocked":
            unresolved_artifacts = [
                engine["checkpoint_sha256"],
                engine["container_digest"],
                engine["package_lock_hash"],
            ]
            assert not all(unresolved_artifacts), (
                f"{engine_id} cannot become distributable merely by filling metadata; "
                "a reviewed status change is required"
            )

    assert engine_ids == {"rfdiffusion-v1", "protenix-v1", "boltz-2", "patchr"}
    assert {entry["id"] for entry in inventory["excluded_engines"]} == {
        "framedipt",
        "chroma",
        "rfdiffusion2",
    }


def test_rfdiffusion_smoke_evidence_and_environment_are_pinned() -> None:
    inventory = _load_inventory()
    engine = next(item for item in inventory["engines"] if item["id"] == "rfdiffusion-v1")
    smoke = inventory["rfdiffusion_v1_smoke"]

    assert engine["checkpoint_sha256"] == (
        "0fcf7d7c32b4848030aca3a051e6768de194616f96ba6c38186351a33bfc6eca"
    )
    assert engine["package_lock_hash"] == _sha256(RFDIFFUSION_ENVIRONMENT)
    assert RFDIFFUSION_ENVIRONMENT_SHA256 == engine["package_lock_hash"]
    assert smoke["same_seed_repeat_coordinate_rmsd_angstrom"] == 0.0
    assert smoke["validation_pass_count"] == 0
    assert smoke["failed_gate"] == "junction-peptide-connectivity"
    assert smoke["observed_peak_vram_bytes"] > 0
    assert smoke["withheld_10_detectable_d_ca"] == 2
    assert "d-ca-chirality" in smoke["withheld_10_failed_gates"]

    refined = inventory["rfdiffusion_v1_boundary_refinement"]
    modeller = inventory["modeller_comparator"]
    assert refined["per_step_reinjection"] is False
    assert refined["withheld_5_validation_passed"] is True
    assert refined["withheld_10_validation_passed"] is True
    assert refined["withheld_5_fixed_heavy_rmsd_angstrom"] == 0.0
    assert refined["withheld_10_fixed_heavy_rmsd_angstrom"] == 0.0
    assert refined["withheld_5_detectable_d_ca"] == 0
    assert refined["withheld_10_detectable_d_ca"] == 0
    assert refined["withheld_5_severe_steric_overlaps"] == 0
    assert refined["withheld_10_severe_steric_overlaps"] == 0
    assert refined["repeatability_status"] == "open"
    assert (
        refined["withheld_5_gap_backbone_rmsd_angstrom"]
        < modeller["withheld_5_median_gap_backbone_rmsd_angstrom"]
    )
    assert (
        refined["withheld_10_gap_backbone_rmsd_angstrom"]
        < modeller["withheld_10_median_gap_backbone_rmsd_angstrom"]
    )
