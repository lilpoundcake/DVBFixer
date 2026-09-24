"""Validate the predeclared diffusion research inventory and benchmark policy."""

from __future__ import annotations

import hashlib
import math
import os
import subprocess
import tomllib
from pathlib import Path
from typing import Any

from dvbfixer.model.diffusion.rfdiffusion_v1 import (
    PARTNER_CONTEXT_BACKBONE_DISPLACEMENT_MAX_ANGSTROM,
    PARTNER_CONTEXT_BACKBONE_RMSD_MAX_ANGSTROM,
    RFDIFFUSION_ENVIRONMENT_SHA256,
)
from dvbfixer.model.diffusion.sampler import (
    SamplerCapabilities,
    SamplingAblationMode,
    assess_sampler_conformance,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
INVENTORY = REPO_ROOT / "docs/research/diffusion-gap-reconstruction-inventory.toml"
FIXTURE_MANIFEST = REPO_ROOT / "tests/fixtures/MANIFEST.sha256"
RFDIFFUSION_ENVIRONMENT = REPO_ROOT / "deploy/rfdiffusion-v1/environment.yml"
RFDIFFUSION_ADAPTER_ENVIRONMENT = (
    REPO_ROOT / "deploy/rfdiffusion-v1/adapter-environment.yml"
)
RFDIFFUSION_DOCKERFILE = REPO_ROOT / "deploy/rfdiffusion-v1/Dockerfile"
RFDIFFUSION_ENTRYPOINT = REPO_ROOT / "deploy/rfdiffusion-v1/container-entrypoint.sh"
PROTENIX_HOOK_PATCH = REPO_ROOT / "deploy/protenix-v1/per-step-callback.patch"
PROTENIX_HOOK_SMOKE = REPO_ROOT / "deploy/protenix-v1/hook-smoke.py"
PROTENIX_REINJECTION = REPO_ROOT / "deploy/protenix-v1/reinjection.py"
PROTENIX_CHECKPOINT_SMOKE = REPO_ROOT / "deploy/protenix-v1/checkpoint_gap_smoke.py"
PROTENIX_REFINEMENT_SMOKE = REPO_ROOT / "deploy/protenix-v1/refine_candidate.py"
BOLTZ_HOOK_PATCH = REPO_ROOT / "deploy/boltz-2/per-step-callback.patch"
BOLTZ_HOOK_SMOKE = REPO_ROOT / "deploy/boltz-2/hook-smoke.py"

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
    assert thresholds["version"] == 2
    assert thresholds["fixed_heavy_atom_rmsd_angstrom_max"] == 0.01
    assert thresholds["fixed_heavy_atom_displacement_angstrom_max"] == 0.03
    assert thresholds["detectable_d_ca_max"] == 0
    assert thresholds["junction_cn_angstrom_min"] == 1.20
    assert thresholds["junction_cn_angstrom_max"] == 1.45
    assert thresholds["generated_backbone_break_angstrom_max"] == 1.80
    assert thresholds["same_seed_repeat_rmsd_angstrom_max"] == 0.01
    assert (
        thresholds["partner_context_backbone_rmsd_angstrom_max"]
        == PARTNER_CONTEXT_BACKBONE_RMSD_MAX_ANGSTROM
    )
    assert (
        thresholds["partner_context_backbone_displacement_angstrom_max"]
        == PARTNER_CONTEXT_BACKBONE_DISPLACEMENT_MAX_ANGSTROM
    )
    partner_history = inventory["threshold_history"]["partner_context_v2"]
    assert partner_history["previous_backbone_displacement_angstrom_max"] == 1.5
    assert partner_history["revised_backbone_displacement_angstrom_max"] == 2.0
    assert "1.708 Angstrom" in partner_history["triggering_observation"]


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
        fixture_start, fixture_stop = case["fixture_interval_zero_based_half_open"]
        assert 3 <= stop - start <= 12
        assert fixture_stop - fixture_start == stop - start
        assert case["request_sequence_basis"] == "observed-coordinate-sequence"
        assert case["left_anchor"] and case["right_anchor"]
        assert len(case["generated_residues"]) == stop - start

        sequences = _fasta_sequences(REPO_ROOT / sequence_relative)
        fasta_key = case["sequence_record"]
        sequence = sequences[fasta_key]
        assert sequence[start:stop] == case["expected_sequence"]

        residues = _pdb_residues(fixture, case["target_chain"])
        selected = residues[fixture_start:fixture_stop]
        assert len(selected) == stop - start
        assert all(
            altlocs == {" "}
            for _, _, _, altlocs in residues[fixture_start - 1:fixture_stop + 1]
        )
        assert "".join(CANONICAL_RESIDUES[name] for _, _, name, _ in selected) == case["expected_sequence"]
        assert [
            _identity(number, icode, name, case["target_chain"])
            for number, icode, name, _ in selected
        ] == case["generated_residues"]

        left = residues[fixture_start - 1]
        right = residues[fixture_stop]
        assert _identity(*left[:3], case["target_chain"]) == case["left_anchor"]
        assert _identity(*right[:3], case["target_chain"]) == case["right_anchor"]

        pdb_text = fixture.read_text()
        assert pdb_text.count("\nMODEL ") + int(pdb_text.startswith("MODEL ")) <= 1
        structure_scope = case.get("request_structure_scope", "full-fixture")
        if structure_scope == "full-fixture":
            assert not any(line.startswith("HETATM") for line in pdb_text.splitlines())
            assert not any(
                line.startswith(("LINK  ", "CONECT", "SSBOND"))
                for line in pdb_text.splitlines()
            )
        elif structure_scope == "target-protein-chain-only":
            local_residues = {
                (number, icode.strip())
                for number, icode, _name, _altlocs in residues[
                    fixture_start - 1 : fixture_stop + 1
                ]
            }
            local_coordinates: list[tuple[float, float, float]] = []
            non_target_coordinates: list[tuple[float, float, float]] = []
            heterogen_coordinates: list[tuple[float, float, float]] = []
            for line in pdb_text.splitlines():
                record = line[:6].strip()
                if record not in {"ATOM", "HETATM"}:
                    continue
                coordinate = (
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                )
                identity = (line[22:26].strip(), line[26].strip())
                if (
                    record == "ATOM"
                    and line[21] == case["target_chain"]
                    and identity in local_residues
                ):
                    local_coordinates.append(coordinate)
                elif record == "ATOM" and line[21] != case["target_chain"]:
                    non_target_coordinates.append(coordinate)
                elif record == "HETATM":
                    heterogen_coordinates.append(coordinate)
            minimum_non_target = min(
                math.dist(local, other)
                for local in local_coordinates
                for other in non_target_coordinates
            )
            minimum_heterogen = min(
                math.dist(local, other)
                for local in local_coordinates
                for other in heterogen_coordinates
            )
            assert abs(
                minimum_non_target
                - case["minimum_candidate_anchor_non_target_atom_distance_angstrom"]
            ) < 0.01
            assert abs(
                minimum_heterogen
                - case["minimum_candidate_anchor_heterogen_distance_angstrom"]
            ) < 0.01
            assert minimum_non_target > 5.0
            assert minimum_heterogen > 5.0
        else:
            assert structure_scope == "target-and-context-protein-chains"
            context_chains = set(case["context_chains"])
            assert context_chains
            assert case["target_chain"] not in context_chains
            generated_identities = {
                (identity.split(":")[1], identity.split(":")[2].replace("_", ""))
                for identity in case["generated_residues"]
            }
            generated_coordinates: list[tuple[float, float, float]] = []
            partner_coordinates: list[tuple[float, float, float]] = []
            heterogen_coordinates = []
            for line in pdb_text.splitlines():
                record = line[:6].strip()
                if record not in {"ATOM", "HETATM"}:
                    continue
                coordinate = (
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                )
                residue_identity = (line[22:26].strip(), line[26].strip())
                if (
                    record == "ATOM"
                    and line[21] == case["target_chain"]
                    and residue_identity in generated_identities
                ):
                    generated_coordinates.append(coordinate)
                elif record == "ATOM" and line[21] in context_chains:
                    partner_coordinates.append(coordinate)
                elif record == "HETATM":
                    heterogen_coordinates.append(coordinate)
            minimum_partner = min(
                math.dist(generated, partner)
                for generated in generated_coordinates
                for partner in partner_coordinates
            )
            minimum_heterogen = min(
                math.dist(generated, heterogen)
                for generated in generated_coordinates
                for heterogen in heterogen_coordinates
            )
            assert abs(
                minimum_partner
                - case["minimum_generated_partner_atom_distance_angstrom"]
            ) < 0.01
            assert abs(
                minimum_heterogen
                - case["minimum_generated_heterogen_distance_angstrom"]
            ) < 0.01
            assert minimum_partner < 3.5
            assert minimum_heterogen > 5.0

    assert {"withheld-internal-3-5", "withheld-internal-6-12"} <= strata
    assert {
        "withheld-independent-regular-loop-3-5",
        "withheld-independent-regular-loop-6-12",
    } <= strata
    assert {
        "withheld-difficult-glypro-3-5",
        "withheld-difficult-glypro-6-12",
    } <= strata
    assert "withheld-interface-adjacent-3-5" in strata
    assert "withheld-antibody-insertion-codes-3-5" in strata
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


def test_all_atom_hook_spikes_record_checkpoint_backed_control() -> None:
    inventory = _load_inventory()
    spike = inventory["all_atom_hook_spike"]

    assert spike["selected_sampler"] == "protenix-v1-maintained-patch"
    assert (
        spike["decision"]
        == "protenix-initial-three-way-ablation-passed-broader-evidence-pending"
    )
    assert (REPO_ROOT / spike["evidence"]).is_file()

    patch = inventory["protenix_v1_hook_patch"]
    assert patch["source_revision"] == (
        "85767b811c40ed46e73a9b39519cf6bfca8701ba"
    )
    assert patch["patch_sha256"] == _sha256(PROTENIX_HOOK_PATCH)
    assert REPO_ROOT / patch["patch"] == PROTENIX_HOOK_PATCH
    assert REPO_ROOT / patch["smoke"] == PROTENIX_HOOK_SMOKE
    assert REPO_ROOT / patch["reinjection"] == PROTENIX_REINJECTION
    assert REPO_ROOT / patch["checkpoint_smoke"] == PROTENIX_CHECKPOINT_SMOKE
    assert patch["checkpoint_smoke_sha256"] == _sha256(PROTENIX_CHECKPOINT_SMOKE)
    assert REPO_ROOT / patch["refinement_smoke"] == PROTENIX_REFINEMENT_SMOKE
    assert patch["refinement_smoke_sha256"] == _sha256(PROTENIX_REFINEMENT_SMOKE)
    assert patch["callback_stage"] == "post-update-before-next-step"
    assert patch["stable_atom_axis"] is True
    assert patch["torch_native_weighted_kabsch"] is True
    assert patch["exact_fixed_overwrite_cpu_smoke"] is True
    assert patch["callback_count"] == patch["expected_callback_count"] == 200
    assert patch["real_checkpoint_smoke"] is True
    assert patch["status"] == "initial-three-way-ablation-complete"

    smoke = inventory["protenix_v1_checkpoint_smoke"]
    assert smoke["status"] == "passed-base-model-only"
    assert smoke["model_parameters"] == 368_484_735
    assert smoke["residue_count"] == 20
    assert smoke["atom_count"] == 168
    assert smoke["diffusion_steps"] == 2
    assert smoke["per_step_callback_exercised"] is False
    assert len(smoke["candidate_sha256"]) == 64
    assert len(smoke["summary_sha256"]) == 64

    gap_smoke = inventory["protenix_v1_gap_hook_smoke"]
    assert gap_smoke["status"] == "hook-passed-raw-candidate-failed-junction-gate"
    assert gap_smoke["model_atom_count"] == 2118
    assert gap_smoke["written_atom_count"] == 2117
    assert gap_smoke["callback_count"] == gap_smoke["expected_callback_count"] == 200
    assert gap_smoke["maximum_post_projection_error_angstrom"] == 0.0
    assert gap_smoke["fixed_heavy_rmsd_angstrom"] == 0.0
    assert gap_smoke["fixed_heavy_max_displacement_angstrom"] == 0.0
    assert gap_smoke["validation_passed"] is False
    assert gap_smoke["failed_gates"] == ["junction-peptide-connectivity"]
    assert len(gap_smoke["candidate_sha256"]) == 64
    assert len(gap_smoke["summary_sha256"]) == 64
    assert gap_smoke["candidate_sha256"] == gap_smoke[
        "same_seed_repeat_candidate_sha256"
    ]
    assert gap_smoke["same_seed_repeat_coordinate_rmsd_angstrom"] == 0.0
    assert gap_smoke["repeatability_status"] == "passed"

    refinement = inventory["protenix_v1_boundary_refinement"]
    assert refinement["status"] == "passed-initial-ablation-case"
    assert refinement["validation_passed"] is True
    assert refinement["fixed_heavy_rmsd_angstrom"] == 0.0
    assert refinement["fixed_heavy_max_displacement_angstrom"] == 0.0
    assert refinement["detectable_d_ca"] == 0
    assert refinement["severe_steric_overlaps"] == 0
    assert refinement["repeatability_status"] == "passed"
    assert refinement["candidate_sha256"] == refinement["repeat_candidate_sha256"]
    assert refinement["repeat_coordinate_rmsd_angstrom"] == 0.0
    assert refinement["repeat_max_displacement_angstrom"] == 0.0
    assert (
        refinement["gap_backbone_rmsd_angstrom"]
        < refinement["modeller_best_gap_backbone_rmsd_angstrom"]
    )

    ablation = inventory["protenix_v1_three_way_ablation"]
    assert ablation["status"] == "passed-initial-case"
    assert ablation["template_only_validation_passed"] is False
    assert ablation["template_only_callback_count"] == 0
    assert ablation["template_only_fixed_heavy_rmsd_angstrom"] > 0.01
    assert ablation["reinjection_validation_passed"] is False
    assert ablation["reinjection_callback_count"] == 200
    assert ablation["reinjection_fixed_heavy_rmsd_angstrom"] == 0.0
    assert ablation["reinjection_repeat_coordinate_rmsd_angstrom"] == 0.0
    assert ablation["refined_validation_passed"] is True
    assert ablation["refined_fixed_heavy_rmsd_angstrom"] == 0.0
    assert ablation["refined_repeat_coordinate_rmsd_angstrom"] == 0.0

    protenix = next(
        engine for engine in inventory["engines"] if engine["id"] == "protenix-v1"
    )
    assert protenix["checkpoint_url"] == (
        "https://protenix.tos-cn-beijing.volces.com/checkpoint/"
        "protenix_base_default_v1.0.0.pt"
    )
    assert protenix["checkpoint_sha256"] == (
        "2b7d5a8b30494514fc47fd2271a16260528cdba170ba09cc112fdecd8f85ec04"
    )

    patch_text = PROTENIX_HOOK_PATCH.read_text()
    update_offset = patch_text.index("x_l = x_noisy + step_scale_eta")
    callback_offset = patch_text.index("updated = step_callback")
    assert callback_offset > update_offset

    engines = {engine["id"]: engine for engine in inventory["engines"]}
    protenix = engines["protenix-v1"]
    protenix_capabilities = SamplerCapabilities(
        mutable_state_each_step=protenix["mutable_state_each_step"],
        identity_mapping_each_step=protenix["identity_mapping_each_step"],
        fixed_coordinate_overwrite_each_step=protenix[
            "fixed_coordinate_overwrite_each_step"
        ],
        localized_boundary_refinement=protenix["localized_boundary_refinement"],
    )
    assert assess_sampler_conformance(
        protenix_capabilities, SamplingAblationMode.REINJECTION
    ).supported
    assert assess_sampler_conformance(
        protenix_capabilities,
        SamplingAblationMode.REINJECTION_BOUNDARY_REFINEMENT,
    ).supported

    boltz = engines["boltz-2"]
    assert boltz["checkpoint_sha256"] == (
        "090e82ac8c92f5e943fa1b39e7410a44027bea7243c0bbb3caa67a77fc1428e1"
    )
    boltz_capabilities = SamplerCapabilities(
        mutable_state_each_step=boltz["mutable_state_each_step"],
        identity_mapping_each_step=boltz["identity_mapping_each_step"],
        fixed_coordinate_overwrite_each_step=boltz[
            "fixed_coordinate_overwrite_each_step"
        ],
        localized_boundary_refinement=boltz["localized_boundary_refinement"],
    )
    decision = assess_sampler_conformance(
        boltz_capabilities, SamplingAblationMode.REINJECTION
    )
    assert not decision.supported
    assert decision.missing_capabilities == (
        "identity-mapping-each-step",
    )

    boltz_patch = inventory["boltz_v2_hook_patch"]
    assert boltz_patch["patch_sha256"] == _sha256(BOLTZ_HOOK_PATCH)
    assert boltz_patch["smoke_sha256"] == _sha256(BOLTZ_HOOK_SMOKE)
    assert boltz_patch["cpu_callback_count"] == 2
    assert boltz_patch["gpu_callback_count"] == 2
    assert boltz_patch["maximum_post_projection_error_angstrom"] == 0.0
    assert boltz_patch["final_fixed_coordinates_exact"] is True

    boltz_smoke = inventory["boltz_v2_checkpoint_smoke"]
    assert boltz_smoke["strict_checkpoint_load"] is True
    assert boltz_smoke["model_parameters"] == 506_724_992
    assert boltz_smoke["failed_examples"] == 0
    assert boltz_smoke["candidate_sha256"] == boltz_smoke[
        "repeat_candidate_sha256"
    ]
    assert boltz_smoke["per_step_callback_exercised"] is False

    broader = inventory["protenix_v1_broader_gap_evidence"]
    assert broader["status"] == "passed-three-additional-single-chain-cases"
    assert broader["case_count"] == 3
    assert broader["all_fixed_heavy_rmsd_angstrom"] == 0.0
    assert broader["all_raw_repeat_coordinate_rmsd_angstrom"] == 0.0
    assert broader["all_refined_repeat_coordinate_rmsd_angstrom"] == 0.0
    assert broader["case_8cz8_10_refined_validation_passed"] is True
    assert broader["case_7x35_7_raw_validation_passed"] is False
    assert broader["case_7x35_7_refined_validation_passed"] is True
    assert broader["case_7k8s_3_refined_validation_passed"] is True


def test_rfdiffusion_smoke_evidence_and_environment_are_pinned() -> None:
    inventory = _load_inventory()
    engine = next(item for item in inventory["engines"] if item["id"] == "rfdiffusion-v1")
    smoke = inventory["rfdiffusion_v1_smoke"]

    assert engine["checkpoint_sha256"] == (
        "0fcf7d7c32b4848030aca3a051e6768de194616f96ba6c38186351a33bfc6eca"
    )
    assert engine["package_lock_hash"] == _sha256(RFDIFFUSION_ENVIRONMENT)
    assert RFDIFFUSION_ENVIRONMENT_SHA256 == engine["package_lock_hash"]
    assert engine["container_base_digest"] == (
        "sha256:008e06cd8432eb558faa4738a092f30b38dd8db3137a5dd3fca57374a790825b"
    )
    assert engine["container_recipe_sha256"] == _sha256(RFDIFFUSION_DOCKERFILE)
    assert engine["container_adapter_environment_sha256"] == _sha256(
        RFDIFFUSION_ADAPTER_ENVIRONMENT
    )
    assert engine["container_entrypoint_sha256"] == _sha256(RFDIFFUSION_ENTRYPOINT)
    assert engine["container_embeds_checkpoint"] is False
    assert engine["container_digest"] == ""
    dockerfile = RFDIFFUSION_DOCKERFILE.read_text()
    adapter_environment = RFDIFFUSION_ADAPTER_ENVIRONMENT.read_text()
    entrypoint = RFDIFFUSION_ENTRYPOINT.read_text()
    assert engine["container_base_digest"] in dockerfile
    assert engine["revision"] in dockerfile
    assert engine["checkpoint_sha256"] in entrypoint
    assert "COPY Base_ckpt.pt" not in dockerfile
    assert "python=3.11.16" in adapter_environment
    assert "python=3.13" not in adapter_environment
    assert "python=3.14" not in adapter_environment
    assert "numpy==2.4.6" in adapter_environment
    assert "scipy==1.17.1" in adapter_environment
    assert "propka" not in adapter_environment.lower()
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
    assert refined["revision"] == "dvbfixer-openmm-boundary-refinement-v3"
    assert refined["repeatability_status"] == "passed-initial-subset"
    assert refined["withheld_5_repeat_coordinate_rmsd_angstrom"] == 0.0
    assert refined["withheld_10_repeat_coordinate_rmsd_angstrom"] == 0.0
    assert (
        refined["withheld_5_gap_backbone_rmsd_angstrom"]
        < modeller["withheld_5_median_gap_backbone_rmsd_angstrom"]
    )
    assert (
        refined["withheld_10_gap_backbone_rmsd_angstrom"]
        < modeller["withheld_10_median_gap_backbone_rmsd_angstrom"]
    )

    independent = inventory["independent_regular_loop_evidence"]
    assert independent["fixture"] == "tests/fixtures/numbering/8b01_a_b.pdb"
    assert independent["withheld_5_validation_passed"] is True
    assert independent["withheld_10_validation_passed"] is True
    assert independent["withheld_5_repeat_coordinate_rmsd_angstrom"] == 0.0
    assert independent["withheld_10_repeat_coordinate_rmsd_angstrom"] == 0.0
    assert independent["withheld_5_modeller_validation_pass_count"] == 0
    assert independent["withheld_10_modeller_validation_pass_count"] == 0
    assert independent["modeller_comparison_passed"] is True
    assert independent["repeatability_status"] == "passed"
    assert (
        independent["withheld_5_gap_backbone_rmsd_angstrom"]
        < independent["withheld_5_modeller_median_gap_backbone_rmsd_angstrom"]
    )
    assert (
        independent["withheld_10_gap_backbone_rmsd_angstrom"]
        < independent["withheld_10_modeller_median_gap_backbone_rmsd_angstrom"]
    )

    glypro = inventory["difficult_glypro_loop_evidence"]
    assert glypro["request_structure_scope"] == "target-protein-chain-only"
    assert glypro["glypro_5_validation_passed"] is True
    assert glypro["glypro_7_validation_passed"] is True
    assert glypro["glypro_5_repeat_coordinate_rmsd_angstrom"] == 0.0
    assert glypro["glypro_7_repeat_coordinate_rmsd_angstrom"] == 0.0
    assert glypro["glypro_5_modeller_validation_pass_count"] == 0
    assert glypro["glypro_7_modeller_validation_pass_count"] == 0
    assert glypro["modeller_median_comparison_passed"] is True
    assert (
        glypro["glypro_5_gap_backbone_rmsd_angstrom"]
        < glypro["glypro_5_modeller_median_gap_backbone_rmsd_angstrom"]
    )
    assert (
        glypro["glypro_7_gap_backbone_rmsd_angstrom"]
        < glypro["glypro_7_modeller_median_gap_backbone_rmsd_angstrom"]
    )
    assert (
        glypro["glypro_7_gap_backbone_rmsd_angstrom"]
        > glypro["glypro_7_modeller_best_gap_backbone_rmsd_angstrom"]
    )

    interface = inventory["interface_adjacent_evidence"]
    assert interface["case"] == "7x35-chain-a-interface-5"
    assert interface["context_chains"] == ["B"]
    assert interface["validation_passed"] is True
    assert interface["repeatability_status"] == "passed"
    assert interface["repeat_coordinate_rmsd_angstrom"] == 0.0
    assert interface["fixed_heavy_rmsd_angstrom"] == 0.0
    assert interface["detectable_d_ca"] == 0
    assert interface["severe_steric_overlaps"] == 0
    assert interface["modeller_validation_pass_count"] == 0
    assert interface["modeller_median_comparison_passed"] is True
    assert (
        interface["gap_backbone_rmsd_angstrom"]
        < interface["modeller_median_gap_backbone_rmsd_angstrom"]
    )
    assert interface["partner_backbone_rmsd_angstrom"] < (
        PARTNER_CONTEXT_BACKBONE_RMSD_MAX_ANGSTROM
    )
    assert interface["partner_backbone_max_displacement_angstrom"] < (
        PARTNER_CONTEXT_BACKBONE_DISPLACEMENT_MAX_ANGSTROM
    )

    insertion = inventory["insertion_code_evidence"]
    assert insertion["case"] == "7k8s-chain-h-insertion-codes-3"
    assert insertion["generated_residues"] == [
        "H:82:A:ASN",
        "H:82:B:SER",
        "H:82:C:LEU",
    ]
    assert insertion["identity_restoration_passed"] is True
    assert insertion["validation_passed"] is True
    assert insertion["repeatability_status"] == "passed"
    assert insertion["repeat_coordinate_rmsd_angstrom"] == 0.0
    assert insertion["fixed_heavy_rmsd_angstrom"] == 0.0
    assert insertion["detectable_d_ca"] == 0
    assert insertion["severe_steric_overlaps"] == 0
    assert insertion["modeller_validation_pass_count"] == 0
    assert insertion["modeller_median_comparison_passed"] is True
    assert (
        insertion["gap_backbone_rmsd_angstrom"]
        < insertion["modeller_median_gap_backbone_rmsd_angstrom"]
    )


def test_rfdiffusion_container_entrypoint_fails_closed_on_checkpoint(
    tmp_path: Path,
) -> None:
    missing = tmp_path / "missing.pt"
    environment = {**os.environ, "DVBFIXER_RF_CHECKPOINT": str(missing)}
    result = subprocess.run(
        ("sh", str(RFDIFFUSION_ENTRYPOINT)),
        capture_output=True,
        text=True,
        check=False,
        env=environment,
    )
    assert result.returncode == 2
    assert "missing mounted RFdiffusion checkpoint" in result.stderr

    invalid = tmp_path / "invalid.pt"
    invalid.write_bytes(b"not the pinned checkpoint")
    environment["DVBFIXER_RF_CHECKPOINT"] = str(invalid)
    result = subprocess.run(
        ("sh", str(RFDIFFUSION_ENTRYPOINT)),
        capture_output=True,
        text=True,
        check=False,
        env=environment,
    )
    assert result.returncode == 2
    assert "checkpoint digest mismatch" in result.stderr
