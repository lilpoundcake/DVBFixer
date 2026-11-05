"""Final installer outcomes remain truthful and actionable."""

import subprocess
from types import SimpleNamespace

import pytest

from dvbfixer import install, install_progress, install_result


@pytest.mark.parametrize(("code", "preview", "cancelled", "expected"), [
    (0, False, False, "SUCCESS"),
    (0, True, False, "PREVIEW COMPLETED"),
    (0, False, True, "CANCELLED"),
    (130, False, False, "CANCELLED"),
    (2, False, False, "ACTION REQUIRED"),
    (1, False, False, "ERROR"),
])
def test_final_status(code, preview, cancelled, expected, capsys):
    install_result.finish(install_result.InstallResult(preview=preview, cancelled=cancelled), code)
    output = capsys.readouterr().out
    assert f"{expected} —" in output
    if expected != "SUCCESS":
        assert "SUCCESS —" not in output


def test_main_always_finishes_and_restores_context(monkeypatch, capsys):
    monkeypatch.setattr(install, "_main", lambda argv: 0)
    assert install.main([]) == 0
    assert capsys.readouterr().out.count("=== Installation result ===") == 1
    assert install_result.current.get() is None


def test_help_is_not_an_installation(monkeypatch, capsys):
    with pytest.raises(SystemExit) as exc:
        install.main(["--help"])
    assert exc.value.code == 0
    assert "Installation result" not in capsys.readouterr().out
    assert install_result.current.get() is None


def test_interrupt_is_not_success(monkeypatch, capsys):
    def interrupted(argv):
        raise KeyboardInterrupt()

    monkeypatch.setattr(install, "_main", interrupted)
    assert install.main([]) == 130
    output = capsys.readouterr().out
    assert "CANCELLED" in output and "SUCCESS" not in output


def test_usage_error_is_not_missing_license(capsys):
    with pytest.raises(SystemExit) as exc:
        install.main(["--not-an-installer-option"])
    assert exc.value.code == 2
    output = capsys.readouterr().out
    assert "ERROR —" in output
    assert "Invalid installer arguments" in output
    assert "ACTION REQUIRED" not in output


def test_nested_failure_keeps_actual_stage_and_log(tmp_path, monkeypatch, capsys):
    log = tmp_path / "solver.log"
    log.write_text("Solving\nERROR: No space left on device\n")

    def failed(argv):
        try:
            with install_progress.stage("Provision backend", interval=60):
                with install_progress.stage("Install dependencies", log, follow=False):
                    install_result.started()
                    raise subprocess.CalledProcessError(1, ["manager", "create"])
        except subprocess.SubprocessError as exc:
            install_result.failure(exc)
            with install_progress.stage("Recovery", follow=False):
                pass
            return 1

    monkeypatch.setattr(install, "_main", failed)
    assert install.main([]) == 1
    output = capsys.readouterr().out.split("=== Installation result ===")[1]
    assert "Uncompleted step: Install dependencies" in output
    assert "No space left" in output and "Free disk space" in output
    assert f"Step log: {log}" in output
    assert "no automatic rollback" in output


@pytest.mark.parametrize("code, finding, remedy", [
    (2, "action-required: modeller: Configure MODELLER license", "protected key descriptor"),
    (1, "error: pdbfixer-resources: missing soft.xml", "do not copy individual"),
])
def test_doctor_findings_survive_subsequent_checks(tmp_path, monkeypatch, capsys, code, finding, remedy):
    log = tmp_path / "doctor.log"
    log.write_text(f"MISSING optional feature\n{finding}\n")

    def checked(argv):
        install_result.readiness(code, log)
        with install_progress.stage("Check Desktop launcher", follow=False):
            pass
        return code

    monkeypatch.setattr(install, "_main", checked)
    assert install.main([]) == code
    output = capsys.readouterr().out.split("=== Installation result ===")[1]
    assert finding in output and remedy in output
    assert "SUCCESS" not in output
    assert "MISSING optional feature" not in output
    assert f"Step log: {log}" in output


def test_runtime_missing_checkpoint_is_action_required(monkeypatch, capsys):
    def provisioned(argv):
        install_result.required(["Supply the checkpoint before resuming"])
        return 2

    monkeypatch.setattr(install, "_main", provisioned)
    assert install.main([]) == 2
    output = capsys.readouterr().out
    assert "ACTION REQUIRED" in output and "managed download" in output


def test_missing_weights_does_not_mislabel_successful_versions_or_infer_network_failure(tmp_path, monkeypatch, capsys):
    log = tmp_path / "versions-refinement.log"
    log.write_text('{"packages": ["connection-pool", "certificate-bundle", "node"]}')

    def provisioned(argv):
        install_result.stage_started("versions-refinement", log)
        install_result.required(["Provide cc94_epoch3100.pth with --checkpoint"])
        return 2

    monkeypatch.setattr(install, "_main", provisioned)
    assert install.main([]) == 2
    output = capsys.readouterr().out
    assert "Pending: required setup actions" in output
    assert "Uncompleted step" not in output and "Step log:" not in output
    assert "network/proxy" not in output and "GUI build log" not in output


@pytest.mark.parametrize("scope, expected", [("package_only", "package only"), ("gui_only", "GUI only")])
def test_narrow_scope_success_does_not_claim_full_profile(scope, expected, monkeypatch, capsys):
    def applied(argv):
        args = SimpleNamespace(prefix=None, dry_run=False, runtime=None, gui_only=False,
                               package_only=False, extend=False, profile="desktop")
        setattr(args, scope, True)
        install_result.configure(args)
        return 0

    monkeypatch.setattr(install, "_main", applied)
    assert install.main([]) == 0
    output = capsys.readouterr().out
    assert expected in output
    assert "scientific qualification" not in output
