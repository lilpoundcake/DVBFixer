"""Sequence-guided structural superposition with Biopython or Modeller SALIGN."""

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


def run_biopython_superposition(specs: list[str], output: Path, fit_dir: Path | None = None,
                                fit_atoms: str = "CA", msa_engine: str = "auto",
                                verbose: bool = False) -> list[Path]:
    """Sequence-guided multiple structural superposition without Modeller.

    Protein sequences are aligned by the configured DVBfixer MSA engine, then
    Biopython's SVD Superimposer fits every template onto the first template
    using corresponding atoms from the aligned residues.
    """
    if len(specs) < 2:
        raise ValueError("structural alignment requires at least two templates")
    if fit_atoms != "CA":
        raise ValueError("the biopython engine currently supports --fit-atoms CA")
    from Bio.PDB import PDBIO, PDBParser, Superimposer
    from Bio.PDB.Polypeptide import is_aa
    from Bio.SeqUtils import seq1

    from dvbfixer.msa import run_alignment

    parsed = [parse_template_spec(spec) for spec in specs]
    output.parent.mkdir(parents=True, exist_ok=True)
    if fit_dir:
        fit_dir.mkdir(parents=True, exist_ok=True)
    fitted: list[Path] = []
    with tempfile.TemporaryDirectory(prefix="dvbfixer_biosalign_") as workdir_text:
        workdir = Path(workdir_text)
        staged: list[Path] = []
        codes: list[str] = []
        chains: list[str] = []
        for index, (source, chain) in enumerate(parsed, start=1):
            if not chain:
                raise ValueError("the biopython engine requires PATH:CHAIN for every template")
            code = f"template_{index}_{chain}"
            destination = workdir / f"{code}.pdb"
            if source.suffix.lower() == ".pdb":
                extract_chain(source, chain, destination)
            else:
                from Bio.PDB import MMCIFParser
                structure = MMCIFParser(QUIET=True).get_structure(code, str(source))
                model = next(structure.get_models())
                if chain not in model:
                    raise ValueError(f"chain {chain!r} not found in {source}")
                io = PDBIO()
                io.set_structure(model[chain])
                io.save(str(destination))
            staged.append(destination)
            codes.append(code)
            chains.append(chain)

        parser = PDBParser(QUIET=True)
        structures = [parser.get_structure(code, str(file)) for code, file in zip(codes, staged)]
        residues = [[residue for residue in structure.get_residues()
                     if is_aa(residue, standard=False) and fit_atoms in residue]
                    for structure in structures]
        sequences = ["".join(seq1(residue.resname, custom_map={"MSE": "M"}, undef_code="X")
                             for residue in rows) for rows in residues]
        fasta = workdir / "structures.fasta"
        fasta.write_text("".join(f">{code}\n{sequence}\n" for code, sequence in zip(codes, sequences)))
        aligned_fasta = workdir / "structures_aligned.fasta"
        selected_engine = run_alignment(fasta, aligned_fasta, msa_engine, "fasta", verbose)
        aligned: dict[str, str] = {}
        current = ""
        for raw in aligned_fasta.read_text().splitlines():
            if raw.startswith(">"):
                current = raw[1:].split()[0]
                aligned[current] = ""
            elif current:
                aligned[current] += raw.strip()

        reference_row = aligned[codes[0]]
        reference_by_column: list[object | None] = []
        residue_index = 0
        for character in reference_row:
            reference_by_column.append(None if character == "-" else residues[0][residue_index])
            if character != "-":
                residue_index += 1

        for index, (code, structure, row_residues) in enumerate(zip(codes, structures, residues)):
            if index:
                mobile_by_column: list[object | None] = []
                residue_index = 0
                for character in aligned[code]:
                    mobile_by_column.append(None if character == "-" else row_residues[residue_index])
                    if character != "-":
                        residue_index += 1
                pairs = [(ref, mobile) for ref, mobile in zip(reference_by_column, mobile_by_column)
                         if ref is not None and mobile is not None]
                if len(pairs) < 3:
                    raise RuntimeError(f"fewer than three aligned CA atoms for {code}")
                superimposer = Superimposer()
                superimposer.set_atoms([pair[0][fit_atoms] for pair in pairs],
                                       [pair[1][fit_atoms] for pair in pairs])
                superimposer.apply(list(structure.get_atoms()))
                if verbose:
                    print(f"{code}: {len(pairs)} CA pairs, RMSD {superimposer.rms:.3f} Å ({selected_engine})")
            if fit_dir:
                destination = fit_dir / f"{code}_fit.pdb"
                io = PDBIO()
                io.set_structure(structure)
                io.save(str(destination))
                fitted.append(destination)

        lines: list[str] = []
        for code, chain in zip(codes, chains):
            lines.extend((f">P1;{code}", f"structureX:{code}::{chain}::{chain}::::",
                          aligned[code] + "*"))
        output.write_text("\n".join(lines) + "\n")
    return fitted


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="dvbfixer salign",
        description="Align and superpose multiple structures with Biopython or Modeller SALIGN.",
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
    alignment.add_argument("--engine", choices=["biopython", "modeller"], default="biopython",
                           help="Structural fitting engine (default: biopython)")
    alignment.add_argument("--msa-engine", choices=["auto", "mafft", "muscle", "clustalo"], default="auto",
                           help="Sequence correspondence engine used by biopython (default: auto)")
    diagnostics = parser.add_argument_group("Diagnostics")
    diagnostics.add_argument("-v", "--verbose", action="store_true")
    from dvbfixer.batch import add_runtime_help
    add_runtime_help(parser)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    try:
        runner = run_biopython_superposition if args.engine == "biopython" else run_salign
        kwargs = {"fit_atoms": args.fit_atoms, "verbose": args.verbose}
        if args.engine == "biopython":
            kwargs["msa_engine"] = args.msa_engine
        else:
            kwargs["rms_cutoff"] = args.rms_cutoff
        fitted = runner(args.template, Path(args.output).resolve(),
                        Path(args.fit_dir).resolve() if args.fit_dir else None, **kwargs)
    except (ValueError, RuntimeError) as exc:
        raise SystemExit(f"ERROR: {exc}") from exc
    except Exception as exc:
        raise SystemExit(f"ERROR: structural alignment failed: {exc}") from exc
    print(f"Wrote {Path(args.output).resolve()}")
    if fitted:
        print(f"Wrote {len(fitted)} fitted structure(s) to {Path(args.fit_dir).resolve()}")


if __name__ == "__main__":
    main()
