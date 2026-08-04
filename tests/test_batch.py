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
            continue_on_error=False,
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
                continue_on_error=False,
            ),
            ["-o", "one.pdb"],
        )
