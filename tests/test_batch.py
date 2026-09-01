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


@pytest.mark.parametrize(
    "argv",
    [
        ["input.pdb", "--output", "output.pdb"],
        ["input.pdb", "--log", "command-specific-value"],
        ["--input", "command-specific-value"],
        ["input.pdb", "--fail", "command-specific-value"],
    ],
)
def test_extract_batch_options_does_not_abbreviate_global_flags(argv):
    options, remaining = extract_batch_options(argv)
    assert options.input_dir is None
    assert options.output_dir is None
    assert options.log_file is None
    assert options.fail_fast is False
    assert remaining == argv


def test_directory_runs_each_pdb_and_preserves_subdirectories(tmp_path: Path):
    source = tmp_path / "structures"
    nested = source / "nested"
    nested.mkdir(parents=True)
    (source / "a.pdb").write_text("END\n")
    (nested / "b.ent").write_text("END\n")
    (source / "c.cif").write_text("data_c\n")
    (nested / "d.mmcif").write_text("data_d\n")
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
        [str(source / "c.cif"), "--no-solvent", "-o", str(output / "c_zbs.pdb")],
        [
            str(nested / "b.ent"),
            "--no-solvent",
            "-o",
            str(output / "nested" / "b_zbs.pdb"),
        ],
        [
            str(nested / "d.mmcif"),
            "--no-solvent",
            "-o",
            str(output / "nested" / "d_zbs.pdb"),
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
    assert "command exited with status 1" in output
    assert "How to fix:" in output


def test_directory_summary_retains_command_error_and_explains_fasta_chain_mismatch(
    tmp_path: Path, capsys,
):
    (tmp_path / "8cxi_t_b.pdb").write_text("END\n")

    def fail(_argv):
        print(
            "Error: FASTA missing sequences for chain(s): A. "
            "FASTA has: B, D, E. PDB has: A.",
            file=__import__("sys").stderr,
        )
        raise SystemExit(1)

    with pytest.raises(SystemExit, match="1"):
        run_directory(
            "zbs", fail,
            Namespace(input_dir=str(tmp_path), output_dir=None,
                      recursive=False, fail_fast=False),
            ["--fasta", "chains.fasta"],
        )

    output = capsys.readouterr()
    combined = output.out + output.err
    assert "Cause: FASTA missing sequences for chain(s): A" in combined
    assert "Rename/add the FASTA header for every listed PDB chain" in combined
    assert "Chain IDs are case-sensitive" in combined


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


def test_diagnose_exit_one_is_findings_not_execution_failure(tmp_path: Path, capsys):
    (tmp_path / "a.pdb").write_text("END\n")
    with pytest.raises(SystemExit) as caught:
        run_directory(
            "diagnose", lambda argv: (_ for _ in ()).throw(SystemExit(1)),
            Namespace(input_dir=str(tmp_path), output_dir=None,
                      recursive=False, fail_fast=False),
            [],
        )
    assert caught.value.code == 1  # retains diagnose's CI-friendly semantics
    output = capsys.readouterr().out
    assert "1 with ERROR findings, 0 execution failures" in output
    assert "FAILED:" not in output
