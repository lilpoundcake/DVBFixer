"""Interactive waiting is not work, and only the innermost stage reports activity."""

import threading

import pytest

from dvbfixer import install_progress, install_ui
from dvbfixer.install_progress import stage, waiting_for_input


def test_only_innermost_stage_reports_activity(capsys):
    with stage("outer provisioning", interval=0.01):
        with stage("inner publication", interval=0.01):
            threading.Event().wait(0.08)
            output = capsys.readouterr().out
            assert "Still running: inner publication" in output
            assert "Still running: outer provisioning" not in output
    assert install_progress._current.get() is None


@pytest.mark.parametrize("answer", ["1", "0", None])
def test_confirmation_wait_silences_all_enclosing_stages(monkeypatch, capsys, answer):
    def respond(prompt):
        assert prompt == "Choice: "
        capsys.readouterr()
        threading.Event().wait(0.08)
        assert capsys.readouterr().out == ""
        if answer is None:
            raise KeyboardInterrupt()
        return answer

    monkeypatch.setattr("builtins.input", respond)
    with stage("Provision additional backend", interval=0.01):
        with stage("Publish validated backend launcher", interval=0.01):
            if answer == "0":
                with pytest.raises(EOFError):
                    install_ui.confirm("Replace the existing launcher?")
            elif answer is None:
                with pytest.raises(KeyboardInterrupt):
                    install_ui.confirm("Replace the existing launcher?")
            else:
                assert install_ui.confirm("Replace the existing launcher?")
            threading.Event().wait(0.08)
            resumed = capsys.readouterr().out
            assert "Still running: Publish validated backend launcher" in resumed
            assert "Still running: Provision additional backend" not in resumed
    assert install_progress._current.get() is None


def test_wait_time_is_excluded_even_for_nested_prompts(monkeypatch, capsys):
    clock = [0.0]
    monkeypatch.setattr(install_progress.time, "monotonic", lambda: clock[0])
    with stage("publication", interval=60):
        with waiting_for_input():
            clock[0] += 1000
            with waiting_for_input():
                clock[0] += 1000
            clock[0] += 131
    assert "Completed: publication (0.0s active)" in capsys.readouterr().out
    assert install_progress._current.get() is None


def test_prompt_and_log_following_do_not_interleave(tmp_path, monkeypatch, capsys):
    log = tmp_path / "publication.log"
    log.write_text("")

    def respond(prompt):
        capsys.readouterr()
        log.write_text("retained diagnostic\n")
        threading.Event().wait(0.08)
        assert capsys.readouterr().out == ""
        return "2"

    monkeypatch.setattr("builtins.input", respond)
    with stage("publication", log, interval=0.01):
        assert not install_ui.confirm("Replace?")
        threading.Event().wait(0.08)
    assert "retained diagnostic" in capsys.readouterr().out
    assert log.read_text() == "retained diagnostic\n"


def test_independent_progress_contexts_do_not_silence_each_other(capsys):
    with stage("waiting parent", interval=0.01), waiting_for_input():
        worker = threading.Thread(target=lambda: _independent_work())
        worker.start()
        worker.join(timeout=2)
        assert not worker.is_alive()
        output = capsys.readouterr().out
        assert "Still running: independent work" in output
        assert "Still running: waiting parent" not in output


def _independent_work():
    with stage("independent work", interval=0.01):
        threading.Event().wait(0.08)
