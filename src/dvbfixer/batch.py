"""Shared directory-input support for single-structure commands."""

from __future__ import annotations

import argparse
from collections.abc import Callable, Sequence
from pathlib import Path

OUTPUT_SUFFIXES = {
    "split": "_split.pdb",
    "renumber": "_renum.pdb",
    "model": "_model.pdb",
    "pull": "_pulled.pdb",
    "prepare": "_prepared.pdb",
    "minimize": "_minimized.pdb",
    "protonate": "_prot.pdb",
    "rename": "_canon.pdb",
    "convert": "_converted.pdb",
    "conect": "_conect.pdb",
    "puppet": "_puppet.pdb",
    "diagnose": "_diagnose.txt",
    "zbs": "_zbs.pdb",
}


def extract_batch_options(argv: Sequence[str]) -> tuple[argparse.Namespace, list[str]]:
    """Remove global batch options while leaving command options untouched."""
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--input-dir")
    parser.add_argument("--output-dir")
    parser.add_argument("--recursive", action="store_true")
    parser.add_argument("--fail-fast", action="store_true")
    return parser.parse_known_args(list(argv))


def run_directory(
    command: str,
    command_main: Callable[[list[str]], object],
    options: argparse.Namespace,
    command_argv: list[str],
) -> None:
    """Run a normal single-input command once for every structure in a folder."""
    if command not in OUTPUT_SUFFIXES:
        supported = ", ".join(sorted(OUTPUT_SUFFIXES))
        raise SystemExit(
            f"dvbfixer {command} does not support --input-dir. Supported commands: {supported}"
        )
    if "-o" in command_argv or "--output" in command_argv:
        raise SystemExit("Use --output-dir, not -o/--output, with --input-dir")

    input_dir = Path(options.input_dir).expanduser()
    if not input_dir.is_dir():
        raise SystemExit(f"Input directory does not exist or is not a directory: {input_dir}")
    output_dir = Path(options.output_dir or f"{input_dir.name}_dvbfixer").expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)

    globber = input_dir.rglob if options.recursive else input_dir.glob
    extensions = {".pdb", ".ent"}
    if command == "split":
        extensions.add(".gro")
    inputs = sorted(p for p in globber("*") if p.is_file() and p.suffix.lower() in extensions)
    if not inputs:
        raise SystemExit(f"No supported structure files found in {input_dir}")

    failures: list[tuple[Path, str]] = []
    diagnostic_findings: list[tuple[Path, Path]] = []
    successes = 0
    policy = "stop at first failure" if options.fail_fast else "continue after failures"
    print(f"Batch mode: run '{command}' independently on {len(inputs)} structure(s)")
    print(f"Output directory: {output_dir}")
    print(f"Failure policy: {policy}")
    for index, input_path in enumerate(inputs, 1):
        relative = input_path.relative_to(input_dir)
        destination_dir = output_dir / relative.parent
        destination_dir.mkdir(parents=True, exist_ok=True)
        suffix = OUTPUT_SUFFIXES[command]
        if command == "diagnose" and "--format" in command_argv:
            pos = command_argv.index("--format")
            if pos + 1 < len(command_argv) and command_argv[pos + 1] == "json":
                suffix = "_diagnose.json"
        output_path = destination_dir / f"{input_path.stem}{suffix}"
        print(f"[{index}/{len(inputs)}] {relative}")
        try:
            command_main([str(input_path), *command_argv, "-o", str(output_path)])
        except SystemExit as exc:
            code = exc.code if isinstance(exc.code, int) else (0 if exc.code is None else 1)
            if code == 0:
                successes += 1
                continue
            if command == "diagnose" and code == 1:
                diagnostic_findings.append((relative, output_path))
                print(f"  FINDINGS: {relative} has ERROR-severity findings "
                      f"(report: {output_path})")
                if options.fail_fast:
                    break
                continue
            reason = (f"command exit status {code}" if isinstance(exc.code, int)
                      else str(exc.code))
            failures.append((relative, reason))
            print(f"  FAILED: {relative} ({reason}; see error above)")
        except Exception as exc:  # isolate each structure from failures in the others
            reason = f"{type(exc).__name__}: {exc}"
            failures.append((relative, reason))
            print(f"  FAILED: {relative} ({reason})")
        else:
            successes += 1
        if failures and options.fail_fast:
            break

    if failures:
        processed = successes + len(diagnostic_findings) + len(failures)
        print(f"Batch mode completed: {successes} succeeded, "
              f"{len(failures)} failed, {len(inputs) - processed} not processed.")
        print("Failed structures:")
        for path, reason in failures:
            print(f"  - {path}: {reason}")
        raise SystemExit(1)
    if diagnostic_findings:
        processed = successes + len(diagnostic_findings)
        print(f"Batch diagnose completed: {successes} clean, "
              f"{len(diagnostic_findings)} with ERROR findings, "
              f"0 execution failures, {len(inputs) - processed} not processed.")
        print("Structures with ERROR findings:")
        for path, report_path in diagnostic_findings:
            print(f"  - {path} (report: {report_path})")
        raise SystemExit(1)
    print(f"Batch mode completed: {successes} succeeded, 0 failed.")
