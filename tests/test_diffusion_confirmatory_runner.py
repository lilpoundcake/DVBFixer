"""Tests for resumable confirmatory-cohort orchestration and aggregation."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _load_runner() -> ModuleType:
    path = ROOT / "scripts/run_diffusion_confirmatory_cohort.py"
    spec = importlib.util.spec_from_file_location("diffusion_confirmatory_runner", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_stage_resume_rejects_corrupt_output_and_reruns(tmp_path: Path) -> None:
    runner = _load_runner()
    source = tmp_path / "input.txt"
    output = tmp_path / "output.txt"
    counter = tmp_path / "counter.txt"
    source.write_text("input", encoding="utf-8")
    code = (
        "from pathlib import Path; "
        f"p=Path({str(counter)!r}); n=int(p.read_text())+1 if p.exists() else 1; "
        "p.write_text(str(n)); "
        f"Path({str(output)!r}).write_text('result-'+str(n))"
    )
    stage = runner.StageSpec(
        "fake",
        (sys.executable, "-c", code),
        (source,),
        (output,),
        (output,),
    )

    resumed, first = runner.run_stage(tmp_path, stage)
    assert not resumed
    assert first["exit_status"] == 0
    assert first["inputs"]["input.txt"] == _digest(source)
    assert first["outputs"]["output.txt"] == _digest(output)
    assert first["resources"]["wall_time_seconds"] >= 0
    assert first["resources"]["peak_ram_bytes"] >= 0

    resumed, _ = runner.run_stage(tmp_path, stage)
    assert resumed
    assert counter.read_text() == "1"

    output.write_text("corrupt", encoding="utf-8")
    resumed, third = runner.run_stage(tmp_path, stage)
    assert not resumed
    assert third["exit_status"] == 0
    assert counter.read_text() == "2"
    assert output.read_text() == "result-2"


def test_stage_failure_is_atomic_and_never_marked_complete(tmp_path: Path) -> None:
    runner = _load_runner()
    source = tmp_path / "input.txt"
    source.write_text("input", encoding="utf-8")
    stage = runner.StageSpec(
        "broken",
        (sys.executable, "-c", "raise SystemExit(9)"),
        (source,),
        (tmp_path / "missing.txt",),
        (tmp_path / "missing.txt",),
    )

    resumed, record = runner.run_stage(tmp_path, stage)

    assert not resumed
    assert record["exit_status"] == 9
    assert not (tmp_path / ".confirmatory/stages/broken.json").exists()
    failure = json.loads(
        (tmp_path / ".confirmatory/failures/broken.json").read_text(encoding="utf-8")
    )
    assert failure["command"] == list(stage.command)
    assert failure["outputs"] == {}


def test_config_preserves_symlinked_virtualenv_interpreter(tmp_path: Path) -> None:
    runner = _load_runner()
    base_python = tmp_path / "base-python"
    base_python.write_text("#!/bin/sh\n", encoding="utf-8")
    virtualenv = tmp_path / "protenix-env/bin"
    virtualenv.mkdir(parents=True)
    linked_python = virtualenv / "python"
    linked_python.symlink_to(base_python)
    config_path = tmp_path / "config.json"
    config_path.write_text(
        json.dumps(
            {
                "protenix_checkout": str(tmp_path / "protenix"),
                "protenix_python": str(linked_python),
                "kalign": str(tmp_path / "kalign"),
                "checkpoint": str(tmp_path / "protenix.pt"),
                "dvbfixer_python": str(tmp_path / "dvbfixer-python"),
                "modeller_python": str(tmp_path / "modeller-python"),
                "boltz_checkout": str(tmp_path / "boltz"),
                "boltz_python": str(tmp_path / "boltz-env/bin/python"),
                "boltz_cache": str(tmp_path / "boltz-cache"),
                "boltz_checkpoint": str(tmp_path / "boltz2_conf.ckpt"),
            }
        ),
        encoding="utf-8",
    )

    config = runner._load_config(config_path)
    specs = runner._stage_specs(tmp_path / "workspace", config, tmp_path / "complete.json")

    assert config.protenix_python == linked_python.absolute()
    assert config.protenix_python != linked_python.resolve()
    assert {spec.name: spec for spec in specs}["protenix-inference"].command[0] == str(
        linked_python.absolute()
    )


def test_stage_recreates_declared_output_directory_after_cleanup(tmp_path: Path) -> None:
    runner = _load_runner()
    source = tmp_path / "input.txt"
    source.write_text("input", encoding="utf-8")
    output_dir = tmp_path / "owned-output"
    output_dir.mkdir()
    (output_dir / "stale.txt").write_text("stale", encoding="utf-8")
    output = output_dir / "result.txt"
    stage = runner.StageSpec(
        "needs-parent",
        (
            sys.executable,
            "-c",
            f"from pathlib import Path; Path({str(output)!r}).write_text('result')",
        ),
        (source,),
        (output,),
        (output_dir,),
        directories_to_create=(output_dir,),
    )

    resumed, record = runner.run_stage(tmp_path, stage)

    assert not resumed
    assert record["exit_status"] == 0
    assert output.read_text(encoding="utf-8") == "result"
    assert not (output_dir / "stale.txt").exists()


def test_directory_tree_digest_is_deterministic_and_rejects_symlinks(tmp_path: Path) -> None:
    runner = _load_runner()
    first = tmp_path / "first"
    second = tmp_path / "second"
    (first / "nested").mkdir(parents=True)
    (second / "nested").mkdir(parents=True)
    (first / "z.txt").write_text("z", encoding="utf-8")
    (first / "nested/a.txt").write_text("a", encoding="utf-8")
    (second / "nested/a.txt").write_text("a", encoding="utf-8")
    (second / "z.txt").write_text("z", encoding="utf-8")

    assert runner._directory_tree_sha256(first) == runner._directory_tree_sha256(second)

    (first / "linked.txt").symlink_to(first / "z.txt")
    with pytest.raises(runner.CohortRunnerError, match="contains a symlink"):
        runner._directory_tree_sha256(first)


def test_unchanged_directory_digest_does_not_reread_file_contents(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    runner = _load_runner()
    cache = tmp_path / "cache"
    cache.mkdir()
    artifact = cache / "artifact.bin"
    artifact.write_bytes(b"version-one")
    original_sha256 = runner._sha256
    content_reads: list[Path] = []

    def tracking_sha256(path: Path) -> str:
        if path == artifact:
            content_reads.append(path)
        return original_sha256(path)

    monkeypatch.setattr(runner, "_sha256", tracking_sha256)
    runner._directory_content_sha256.cache_clear()

    first = runner._directory_tree_sha256(cache)
    second = runner._directory_tree_sha256(cache)

    assert first == second
    assert content_reads == [artifact]

    original_metadata = artifact.stat()
    artifact.write_bytes(b"version-two")
    os.utime(
        artifact,
        ns=(original_metadata.st_atime_ns, original_metadata.st_mtime_ns + 1),
    )
    third = runner._directory_tree_sha256(cache)

    assert third != first
    assert content_reads == [artifact, artifact]


def test_directory_input_mutation_invalidates_stage_resume(tmp_path: Path) -> None:
    runner = _load_runner()
    cache = tmp_path / "cache"
    cache.mkdir()
    cache_file = cache / "artifact.bin"
    cache_file.write_bytes(b"version-one")
    output = tmp_path / "output.txt"
    counter = tmp_path / "counter.txt"
    code = (
        "from pathlib import Path; "
        f"p=Path({str(counter)!r}); n=int(p.read_text())+1 if p.exists() else 1; "
        "p.write_text(str(n)); "
        f"Path({str(output)!r}).write_text('result-'+str(n))"
    )
    stage = runner.StageSpec(
        "directory-input",
        (sys.executable, "-c", code),
        (cache,),
        (output,),
        (output,),
    )

    resumed, _ = runner.run_stage(tmp_path, stage)
    assert not resumed
    resumed, _ = runner.run_stage(tmp_path, stage)
    assert resumed
    assert counter.read_text(encoding="utf-8") == "1"

    original_metadata = cache_file.stat()
    cache_file.write_bytes(b"version-two")
    os.utime(
        cache_file,
        ns=(original_metadata.st_atime_ns, original_metadata.st_mtime_ns + 1),
    )
    resumed, record = runner.run_stage(tmp_path, stage)

    assert not resumed
    assert record["exit_status"] == 0
    assert counter.read_text(encoding="utf-8") == "2"
    assert output.read_text(encoding="utf-8") == "result-2"


def _materialized_case(root: Path) -> Path:
    workspace = root / "case-1"
    (workspace / "input").mkdir(parents=True)
    for relative, content in (
        ("request.json", "{}\n"),
        ("input/normalized.pdb", "END\n"),
        ("reference.pdb", "END\n"),
        ("target.fasta", ">chain_A\nA\n"),
    ):
        (workspace / relative).write_text(content, encoding="utf-8")
    marker = {
        "status": "complete",
        "case_id": "case-1",
        "screening_index": 0,
        "selected_for_inference": True,
        "independence_group": "structure-1",
        "leakage_resolved": True,
        "artifacts": {
            relative: {"sha256": _digest(workspace / relative)}
            for relative in runner_files()
        },
    }
    (workspace / "materialization.json").write_text(json.dumps(marker), encoding="utf-8")
    (workspace / "confirmatory/protenix-refined").mkdir(parents=True)
    (workspace / "confirmatory/protenix-refined/summary.json").write_text(
        json.dumps(
            {
                "validation": {"passed": True},
                "quality": {
                    "fixed_heavy_rmsd_angstrom": 0.0,
                    "gap_backbone_rmsd_angstrom": 0.5,
                },
            }
        ),
        encoding="utf-8",
    )
    (workspace / "confirmatory/boltz-refined").mkdir(parents=True)
    (workspace / "confirmatory/boltz-refined/summary.json").write_text(
        json.dumps(
            {
                "validation": {"passed": True},
                "quality": {
                    "fixed_heavy_rmsd_angstrom": 0.0,
                    "gap_backbone_rmsd_angstrom": 0.75,
                },
            }
        ),
        encoding="utf-8",
    )
    (workspace / "confirmatory/modeller-analysis.json").write_text(
        json.dumps(
            {
                "candidates": [
                    {
                        "validation_passed": True,
                        "fixed_heavy_rmsd_angstrom": 0.25,
                        "gap_backbone_rmsd_angstrom": 1.0,
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    stages = workspace / ".confirmatory/stages"
    stages.mkdir(parents=True)
    for name, vram in (
        ("protenix-input", None),
        ("protenix-inference", 4096),
        ("protenix-reference-refinement", None),
        ("boltz-input", None),
        ("boltz-inference", 8192),
        ("boltz-reference-refinement", None),
        ("modeller-comparator", None),
        ("analysis", None),
    ):
        (stages / f"{name}.json").write_text(
            json.dumps(
                {
                    "resources": {
                        "wall_time_seconds": 1.0,
                        "peak_ram_bytes": 1024,
                        "peak_vram_bytes": vram,
                    }
                }
            ),
            encoding="utf-8",
        )
    return workspace


def runner_files() -> tuple[str, ...]:
    return "request.json", "input/normalized.pdb", "reference.pdb", "target.fasta"


def test_runner_excludes_reserve_materialization(tmp_path: Path) -> None:
    runner = _load_runner()
    workspace = _materialized_case(tmp_path)
    assert runner._workspaces(tmp_path, "materialization.json") == [workspace]
    marker_path = workspace / "materialization.json"
    marker = json.loads(marker_path.read_text(encoding="utf-8"))
    marker["selected_for_inference"] = False
    marker_path.write_text(json.dumps(marker), encoding="utf-8")

    with pytest.raises(runner.CohortRunnerError, match="not selected"):
        runner.load_materialization(workspace, "materialization.json")
    assert runner._workspaces(tmp_path, "materialization.json") == []


def test_aggregator_emits_intervals_failures_and_non_applicable_modeller_vram(
    tmp_path: Path,
) -> None:
    runner = _load_runner()
    _materialized_case(tmp_path)
    output = tmp_path / "aggregate.json"

    result = runner.aggregate(
        tmp_path,
        "materialization.json",
        output,
        bootstrap_samples=20,
    )

    assert output.is_file()
    completion = output.with_suffix(".json.complete.json")
    assert completion.is_file()
    assert json.loads(completion.read_text(encoding="utf-8"))["stage"] == "aggregation"
    assert result["preregistered_thresholds"] == {
        "minimum_independence_groups": 180,
        "maximum_independence_groups": 300,
        "minimum_both_valid_independence_groups": 30,
    }
    comparison = result["comparison"]
    assert comparison["candidate"]["backend_id"] == "protenix-v1"
    assert comparison["baseline"]["backend_id"] == "modeller-10.8"
    assert comparison["candidate"]["validation_pass_rate_interval"]["confidence_level"] == 0.95
    assert comparison["pass_rate_difference_interval"]["confidence_level"] == 0.95
    assert comparison["candidate"]["resource_reporting_complete"] is True
    assert comparison["baseline"]["resource_reporting_complete"] is True
    assert comparison["baseline"]["gpu_memory_applicable"] is False
    assert comparison["baseline"]["peak_vram_bytes"] is None
    assert result["failure_rates"] == {
        "protenix-v1": 0.0,
        "modeller-10.8": 0.0,
        "boltz-2-proxy-only": 0.0,
    }
    assert result["resources_by_case_and_stage"]["case-1"]["protenix-inference"][
        "peak_vram_bytes"
    ] == 4096
    proxy = result["proxy_only"]
    assert proxy["exclusion_reason"] == "unresolved-training-leakage"
    assert proxy["backend_summary"]["eligible"] is False
    assert proxy["backend_summary"]["ineligibility_reasons"] == [
        "unresolved-training-leakage"
    ]
    assert proxy["backend_summary"]["resource_reporting_complete"] is True
    assert "decision" not in proxy
    assert proxy["descriptive_comparisons"]["protenix-v1"] == {
        "comparison_backend_id": "protenix-v1",
        "paired_case_count": 1,
        "both_passed_count": 1,
        "proxy_only_passed_count": 0,
        "comparison_only_passed_count": 0,
        "neither_passed_count": 0,
        "proxy_failure_count": 0,
        "comparison_failure_count": 0,
        "both_valid_rmsd_count": 1,
        "proxy_median_valid_gap_backbone_rmsd_angstrom": 0.75,
        "comparison_median_valid_gap_backbone_rmsd_angstrom": 0.5,
        "median_paired_rmsd_difference_angstrom": 0.25,
    }
    assert result["case_results"]["boltz-2-proxy-only"][0][
        "leakage_resolved"
    ] is False
    assert result["resources_by_case_and_stage"]["case-1"]["boltz-inference"][
        "peak_vram_bytes"
    ] == 8192

    first_completion = completion.read_bytes()
    rendered = json.loads(output.read_text(encoding="utf-8"))
    assert runner.aggregate(
        tmp_path,
        "materialization.json",
        output,
        bootstrap_samples=20,
    ) == rendered
    assert completion.read_bytes() == first_completion


def test_stage_specs_use_distinct_refinements_and_pinned_boltz_command(tmp_path: Path) -> None:
    runner = _load_runner()
    config = runner.RunnerConfig(
        protenix_checkout=tmp_path / "protenix",
        protenix_python=tmp_path / "protenix-python",
        kalign=tmp_path / "kalign",
        checkpoint=tmp_path / "protenix.pt",
        dvbfixer_python=tmp_path / "dvbfixer-python",
        modeller_python=tmp_path / "modeller-python",
        boltz_checkout=tmp_path / "boltz",
        boltz_python=tmp_path / "boltz-python",
        boltz_cache=tmp_path / "boltz-cache",
        boltz_checkpoint=tmp_path / "boltz2_conf.ckpt",
    )
    specs = runner._stage_specs(tmp_path / "workspace", config, tmp_path / "complete.json")
    by_name = {spec.name: spec for spec in specs}

    assert "reference-refinement" not in by_name
    assert "protenix-reference-refinement" in by_name
    assert "boltz-reference-refinement" in by_name
    assert by_name["modeller-comparator"].directories_to_create == (
        tmp_path / "workspace/confirmatory/modeller",
    )
    assert by_name["protenix-input"].directories_to_create == ()
    assert by_name["protenix-inference"].directories_to_create == ()
    assert by_name["boltz-input"].directories_to_create == ()
    assert by_name["boltz-inference"].directories_to_create == ()
    command = by_name["boltz-inference"].command
    assert command[0] == str(config.boltz_python)
    assert command[command.index("--cache") + 1] == str(config.boltz_cache)
    assert command[command.index("--checkpoint") + 1] == str(config.boltz_checkpoint)
    assert command[command.index("--recycling-steps") + 1] == "1"
    assert command[command.index("--sampling-steps") + 1] == "200"
    assert not any(path.name == "reinjection.py" for path in by_name["boltz-inference"].inputs)
    assert config.expected_boltz_revision == "b1ebfc46ecf57f5414e0d1a6f9027bbb122c53bc"
    assert config.expected_boltz_checkpoint_sha256 == (
        "090e82ac8c92f5e943fa1b39e7410a44027bea7243c0bbb3caa67a77fc1428e1"
    )


def test_verify_config_checks_boltz_revision_patch_checkpoint_and_cache(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    runner = _load_runner()
    protenix_checkout = tmp_path / "protenix"
    boltz_checkout = tmp_path / "boltz"
    protenix_checkout.mkdir()
    boltz_checkout.mkdir()
    executables: list[Path] = []
    for name in ("protenix-python", "kalign", "dvbfixer-python", "modeller-python", "boltz-python"):
        path = tmp_path / name
        path.write_text("#!/bin/sh\n", encoding="utf-8")
        path.chmod(0o755)
        executables.append(path)
    protenix_checkpoint = tmp_path / "protenix.pt"
    boltz_checkpoint = tmp_path / "boltz2_conf.ckpt"
    protenix_checkpoint.write_bytes(b"protenix")
    boltz_checkpoint.write_bytes(b"boltz")
    cache = tmp_path / "boltz-cache"
    (cache / "mols").mkdir(parents=True)
    (cache / "ccd.pkl").write_bytes(b"ccd")
    (cache / "mols/component.pkl").write_bytes(b"mol")
    config = runner.RunnerConfig(
        protenix_checkout=protenix_checkout,
        protenix_python=executables[0],
        kalign=executables[1],
        checkpoint=protenix_checkpoint,
        dvbfixer_python=executables[2],
        modeller_python=executables[3],
        boltz_checkout=boltz_checkout,
        boltz_python=executables[4],
        boltz_cache=cache,
        boltz_checkpoint=boltz_checkpoint,
        expected_checkpoint_sha256=_digest(protenix_checkpoint),
        expected_boltz_checkpoint_sha256=_digest(boltz_checkpoint),
    )
    commands: list[tuple[str, ...]] = []

    def fake_run(command: tuple[str, ...], **_kwargs: object) -> SimpleNamespace:
        commands.append(command)
        stdout = ""
        if command[-2:] == ("rev-parse", "HEAD"):
            stdout = (
                runner.BOLTZ_REVISION
                if str(boltz_checkout) in command
                else runner.PROTENIX_REVISION
            )
        return SimpleNamespace(returncode=0, stdout=stdout)

    monkeypatch.setattr(runner.subprocess, "run", fake_run)

    verified = runner.verify_config(config)

    assert verified["boltz_revision"] == runner.BOLTZ_REVISION
    assert verified["boltz_checkpoint_sha256"] == _digest(boltz_checkpoint)
    assert any(
        command[-1].endswith("deploy/boltz-2/per-step-callback.patch")
        and "--reverse" in command
        for command in commands
    )


def test_boltz_failure_does_not_block_confirmatory_analysis(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    runner = _load_runner()
    workspace = _materialized_case(tmp_path)
    config = runner.RunnerConfig(
        protenix_checkout=tmp_path,
        protenix_python=Path(sys.executable),
        kalign=Path(sys.executable),
        checkpoint=tmp_path / "protenix.pt",
        dvbfixer_python=Path(sys.executable),
        modeller_python=Path(sys.executable),
        boltz_checkout=tmp_path,
        boltz_python=Path(sys.executable),
        boltz_cache=tmp_path,
        boltz_checkpoint=tmp_path / "boltz.pt",
        materialization_marker="materialization.json",
    )
    names = (
        "protenix-input",
        "protenix-inference",
        "protenix-reference-refinement",
        "boltz-input",
        "boltz-inference",
        "boltz-reference-refinement",
        "modeller-comparator",
        "analysis",
    )
    specs = [runner.StageSpec(name, (name,), (), (), ()) for name in names]
    calls: list[str] = []

    def fake_stage(_workspace: Path, spec: object) -> tuple[bool, dict[str, object]]:
        calls.append(spec.name)
        status = 1 if spec.name == "boltz-inference" else 0
        return False, {"exit_status": status}

    monkeypatch.setattr(runner, "_stage_specs", lambda *_args: specs)
    monkeypatch.setattr(runner, "_resolve_refinement_digest", lambda spec: spec)
    monkeypatch.setattr(runner, "run_stage", fake_stage)

    assert runner.run_case(workspace, config) is False
    assert "boltz-reference-refinement" not in calls
    assert "analysis" in calls
    case = json.loads((workspace / ".confirmatory/case.json").read_text(encoding="utf-8"))
    assert case["confirmatory_complete"] is True
    assert case["proxy_complete"] is False


def test_zero_exit_with_output_error_blocks_analysis(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    runner = _load_runner()
    workspace = _materialized_case(tmp_path)
    config = runner.RunnerConfig(
        protenix_checkout=tmp_path,
        protenix_python=Path(sys.executable),
        kalign=Path(sys.executable),
        checkpoint=tmp_path / "protenix.pt",
        dvbfixer_python=Path(sys.executable),
        modeller_python=Path(sys.executable),
        boltz_checkout=tmp_path,
        boltz_python=Path(sys.executable),
        boltz_cache=tmp_path,
        boltz_checkpoint=tmp_path / "boltz.pt",
        materialization_marker="materialization.json",
    )
    names = (
        "protenix-input",
        "protenix-inference",
        "protenix-reference-refinement",
        "boltz-input",
        "boltz-inference",
        "boltz-reference-refinement",
        "modeller-comparator",
        "analysis",
    )
    specs = [runner.StageSpec(name, (name,), (), (), ()) for name in names]
    calls: list[str] = []

    def fake_stage(_workspace: Path, spec: object) -> tuple[bool, dict[str, object]]:
        calls.append(spec.name)
        if spec.name == "modeller-comparator":
            return False, {"exit_status": 0, "error": "missing output"}
        return False, {"exit_status": 0}

    monkeypatch.setattr(runner, "_stage_specs", lambda *_args: specs)
    monkeypatch.setattr(runner, "_resolve_refinement_digest", lambda spec: spec)
    monkeypatch.setattr(runner, "run_stage", fake_stage)

    assert runner.run_case(workspace, config) is False
    assert "analysis" not in calls


def test_failed_boltz_counts_only_as_proxy_failure(tmp_path: Path) -> None:
    runner = _load_runner()
    workspace = _materialized_case(tmp_path)
    (workspace / "confirmatory/boltz-refined/summary.json").unlink()
    refinement_marker = workspace / ".confirmatory/stages/boltz-reference-refinement.json"
    resources = json.loads(refinement_marker.read_text(encoding="utf-8"))["resources"]
    refinement_marker.unlink()
    failures = workspace / ".confirmatory/failures"
    failures.mkdir()
    (failures / "boltz-reference-refinement.json").write_text(
        json.dumps({"exit_status": 1, "resources": resources}),
        encoding="utf-8",
    )

    result = runner.aggregate(
        tmp_path,
        "materialization.json",
        tmp_path / "aggregate.json",
        bootstrap_samples=20,
    )

    assert result["comparison"]["both_passed_count"] == 1
    assert result["failure_rates"]["protenix-v1"] == 0.0
    assert result["failure_rates"]["modeller-10.8"] == 0.0
    assert result["failure_rates"]["boltz-2-proxy-only"] == 1.0
    assert result["proxy_only"]["backend_summary"]["resource_reporting_complete"] is True
    assert result["proxy_only"]["descriptive_comparisons"]["protenix-v1"][
        "comparison_only_passed_count"
    ] == 1


def test_backend_result_rejects_fabricated_non_applicable_vram() -> None:
    from dvbfixer.model.diffusion.benchmark import BackendCaseResult

    with pytest.raises(ValueError, match="not applicable"):
        BackendCaseResult(
            "modeller",
            "case",
            "group",
            True,
            True,
            True,
            True,
            1.0,
            peak_vram_bytes=1,
            peak_ram_bytes=2,
            gpu_memory_applicable=False,
        )
