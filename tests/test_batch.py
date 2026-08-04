from argparse import Namespace
from pathlib import Path

import pytest

from dvbfixer.batch import extract_batch_options, run_directory


def test_extract_batch_options_preserves_command_arguments():
    options, remaining = extract_batch_options(
        ["--input-dir", "structures", "--no-solvent", "--recursive"]
    )
    assert options.input_dir == "structures"
    assert options.recursive is True
    assert options.fail_fast is False
    assert remaining == ["--no-solvent"]


def test_directory_runs_each_pdb_and_preserves_subdirectories(tmp_path: Path):
    source = tmp_path / "structures"
    nested = source / "nested"
    nested.mkdir(parents=True)
    (source / "a.pdb").write_text("END\n")
    (nested / "b.ent").write_text("END\n")
    (source / "ignore.txt").write_text("not a structure\n")
    output = tmp_path / "results"
    calls = []

    run_directory(
        "zbs",
        calls.append,
        Namespace(
            input_dir=str(source),
            output_dir=str(output),
            recursive=True,
            fail_fast=False,
        ),
        ["--no-solvent"],
    )

    assert calls == [
        [str(source / "a.pdb"), "--no-solvent", "-o", str(output / "a_zbs.pdb")],
        [
            str(nested / "b.ent"),
            "--no-solvent",
            "-o",
            str(output / "nested" / "b_zbs.pdb"),
        ],
    ]


def test_directory_rejects_ambiguous_output_option(tmp_path: Path):
    (tmp_path / "a.pdb").write_text("END\n")
    with pytest.raises(SystemExit, match="Use --output-dir"):
        run_directory(
            "prepare",
            lambda argv: None,
            Namespace(
                input_dir=str(tmp_path),
                output_dir=None,
                recursive=False,
                fail_fast=False,
            ),
            ["-o", "one.pdb"],
        )


def test_directory_continues_by_default_and_prints_clear_summary(tmp_path: Path, capsys):
    (tmp_path / "a.pdb").write_text("END\n")
    (tmp_path / "b.pdb").write_text("END\n")
    calls = []

    def fail(argv):
        calls.append(Path(argv[0]).name)
        raise SystemExit(1)

    with pytest.raises(SystemExit) as caught:
        run_directory(
            "zbs", fail,
            Namespace(input_dir=str(tmp_path), output_dir=None,
                      recursive=False, fail_fast=False),
            ["--no-solvent"],
        )
    assert caught.value.code == 1
    assert calls == ["a.pdb", "b.pdb"]
    output = capsys.readouterr().out
    assert "Batch mode completed: 0 succeeded, 2 failed, 0 not processed" in output
    assert "Batch failed for" not in output


def test_fail_fast_stops_after_first_failure(tmp_path: Path):
    (tmp_path / "a.pdb").write_text("END\n")
    (tmp_path / "b.pdb").write_text("END\n")
    calls = []

    def fail(argv):
        calls.append(Path(argv[0]).name)
        raise SystemExit(2)

    with pytest.raises(SystemExit):
        run_directory(
            "prepare", fail,
            Namespace(input_dir=str(tmp_path), output_dir=None,
                      recursive=False, fail_fast=True),
            [],
        )
    assert calls == ["a.pdb"]
