"""Structure-based multiple alignment using Modeller SALIGN."""

from __future__ import annotations

import argparse
import os
import shutil
import tempfile
from pathlib import Path


def parse_template_spec(spec: str) -> tuple[Path, str | None]:
    """Parse ``PDB[:CHAIN]`` while still accepting paths containing colons."""
    direct = Path(spec).expanduser()
    if direct.exists():
        return direct.resolve(), None
    path_text, separator, chain = spec.rpartition(":")
    candidate = Path(path_text).expanduser()
    if separator and chain and candidate.exists():
        return candidate.resolve(), chain
    raise ValueError(f"template not found: {spec}")


def extract_chain(input_path: Path, chain: str, output_path: Path) -> None:
    """Write one PDB chain plus harmless global records to ``output_path``."""
    atom_count = 0
    lines: list[str] = []
    for line in input_path.read_text().splitlines(keepends=True):
        if line.startswith(("ATOM  ", "HETATM", "ANISOU", "TER   ")):
            if len(line) > 21 and line[21].strip() == chain:
                lines.append(line)
                if line.startswith(("ATOM  ", "HETATM")):
                    atom_count += 1
        elif line.startswith(("HEADER", "TITLE ", "REMARK", "CRYST1")):
            lines.append(line)
    if not atom_count:
        raise ValueError(f"chain {chain!r} has no atoms in {input_path}")
    lines.append("END\n")
    output_path.write_text("".join(lines))


def run_salign(specs: list[str], output: Path, fit_dir: Path | None = None,
               fit_atoms: str = "CA", rms_cutoff: float = 3.5,
               verbose: bool = False) -> list[Path]:
    if len(specs) < 2:
        raise ValueError("structural alignment requires at least two templates")
    parsed = [parse_template_spec(spec) for spec in specs]
    output.parent.mkdir(parents=True, exist_ok=True)
    if fit_dir:
        fit_dir.mkdir(parents=True, exist_ok=True)

    try:
        from modeller import Alignment, Environ, Model, log
    except ImportError as exc:
        raise RuntimeError("Modeller is required for structural alignment") from exc
    if verbose:
        log.verbose()
    else:
        log.none()

    fitted: list[Path] = []
    with tempfile.TemporaryDirectory(prefix="dvbfixer_salign_") as workdir_text:
        workdir = Path(workdir_text)
        codes: list[str] = []
        staged: list[Path] = []
        for index, (source, chain) in enumerate(parsed, start=1):
            code = f"template_{index}_{chain or 'all'}"
            destination = workdir / f"{code}.pdb"
            if chain:
                extract_chain(source, chain, destination)
            else:
                shutil.copy2(source, destination)
            codes.append(code)
            staged.append(destination)

        previous = Path.cwd()
        try:
            os.chdir(workdir)
            env = Environ()
            env.io.atom_files_directory = [str(workdir)]
            env.io.hetatm = True
            alignment = Alignment(env)
            for code, path in zip(codes, staged):
                model = Model(env, file=str(path))
                alignment.append_model(model, align_codes=code, atom_files=str(path))
            alignment.salign(
                alignment_type="tree", auto_overhang=True,
                gap_penalties_1d=(-450, -50), output="",
            )
            alignment.salign(
                alignment_type="progressive", align_block=1,
                feature_weights=(0, 0, 0, 0, 1, 0), fit=True,
                fit_atoms=fit_atoms, rms_cutoff=rms_cutoff,
                improve_alignment=True, write_fit=fit_dir is not None,
                output="QUALITY" if verbose else "",
            )
            local_alignment = workdir / "structural_alignment.pir"
            alignment.write(file=str(local_alignment), alignment_format="PIR")
        finally:
            os.chdir(previous)

        if not local_alignment.exists():
            raise RuntimeError("Modeller SALIGN did not produce an alignment")
        shutil.copy2(local_alignment, output)
        if fit_dir:
            for code in codes:
                candidate = workdir / f"{code}_fit.pdb"
                if candidate.exists():
                    destination = fit_dir / candidate.name
                    shutil.copy2(candidate, destination)
                    fitted.append(destination)
    return fitted


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="dvbfixer salign",
        description="Create a structure-based multiple alignment with Modeller SALIGN.",
    )
    io = parser.add_argument_group("Input / output")
    io.add_argument("template", nargs="+", help="Template PDBs as PATH or PATH:CHAIN (at least two)")
    io.add_argument("-o", "--output", default="structural_alignment.pir",
                    help="Output PIR alignment (default: structural_alignment.pir)")
    io.add_argument("--fit-dir", help="Optional directory for fitted/superposed PDB files")
    alignment = parser.add_argument_group("Alignment")
    alignment.add_argument("--fit-atoms", default="CA", help="Atoms used for fitting (default: CA)")
    alignment.add_argument("--rms-cutoff", type=float, default=3.5,
                           help="RMS cutoff in angstroms (default: 3.5)")
    diagnostics = parser.add_argument_group("Diagnostics")
    diagnostics.add_argument("-v", "--verbose", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    try:
        fitted = run_salign(
            args.template, Path(args.output).resolve(),
            Path(args.fit_dir).resolve() if args.fit_dir else None,
            args.fit_atoms, args.rms_cutoff, args.verbose,
        )
    except (ValueError, RuntimeError) as exc:
        raise SystemExit(f"ERROR: {exc}") from exc
    except Exception as exc:
        raise SystemExit(f"ERROR: Modeller SALIGN failed: {exc}") from exc
    print(f"Wrote {Path(args.output).resolve()}")
    if fitted:
        print(f"Wrote {len(fitted)} fitted structure(s) to {Path(args.fit_dir).resolve()}")


if __name__ == "__main__":
    main()
