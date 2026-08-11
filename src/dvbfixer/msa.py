"""Multiple protein-sequence alignment through common external engines."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import tempfile
from pathlib import Path

ENGINES = {"mafft": "mafft", "muscle": "muscle", "clustalo": "clustalo"}


def available_engines() -> dict[str, str | None]:
    """Return executable paths for the supported MSA engines."""
    return {name: shutil.which(executable) for name, executable in ENGINES.items()}


def _read_fasta(path: Path, min_records: int = 2) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    identifier: str | None = None
    chunks: list[str] = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            if identifier is not None:
                records.append((identifier, "".join(chunks)))
            identifier = line[1:].split()[0]
            if not identifier:
                raise ValueError("FASTA headers must contain an identifier")
            chunks = []
        elif identifier is None:
            raise ValueError("FASTA sequence found before the first header")
        else:
            chunks.append("".join(line.split()).upper())
    if identifier is not None:
        records.append((identifier, "".join(chunks)))
    if len(records) < min_records:
        raise ValueError(f"multiple alignment requires at least {min_records} FASTA records")
    if len({identifier for identifier, _ in records}) != len(records):
        raise ValueError("FASTA identifiers must be unique")
    allowed = set("ABCDEFGHIKLMNPQRSTVWXYZJUO-*?")
    for identifier, sequence in records:
        if not sequence:
            raise ValueError(f"FASTA record {identifier!r} is empty")
        invalid = sorted(set(sequence) - allowed)
        if invalid:
            raise ValueError(f"FASTA record {identifier!r} contains invalid symbols: {''.join(invalid)}")
    return records


def _validate_alignment(path: Path, expected: list[tuple[str, str]]) -> list[tuple[str, str]]:
    aligned = _read_fasta(path)
    if {name for name, _ in aligned} != {name for name, _ in expected}:
        raise RuntimeError("alignment engine changed the FASTA record identifiers")
    lengths = {len(sequence) for _, sequence in aligned}
    if len(lengths) != 1:
        raise RuntimeError("alignment engine produced rows of different lengths")
    expected_map = {name: sequence.replace("-", "").replace("*", "") for name, sequence in expected}
    for name, sequence in aligned:
        if sequence.replace("-", "").replace("*", "") != expected_map[name]:
            raise RuntimeError(f"alignment engine changed the ungapped sequence for {name!r}")
    return aligned


def _write_fasta(records: list[tuple[str, str]], output: Path) -> None:
    lines: list[str] = []
    for identifier, sequence in records:
        lines.append(f">{identifier}")
        lines.extend(sequence[index:index + 80] for index in range(0, len(sequence), 80))
    output.write_text("\n".join(lines) + "\n")


def _write_pir(records: list[tuple[str, str]], output: Path) -> None:
    lines: list[str] = []
    for identifier, sequence in records:
        lines.extend((f">P1;{identifier}", f"sequence:{identifier}::::::::"))
        terminated = sequence.rstrip("*") + "*"
        lines.extend(terminated[index:index + 75] for index in range(0, len(terminated), 75))
    output.write_text("\n".join(lines) + "\n")


def run_alignment(input_path: Path, output_path: Path, engine: str = "auto",
                  output_format: str = "fasta", verbose: bool = False,
                  templates: list[str] | None = None) -> str:
    expected = _read_fasta(input_path, min_records=1 if templates else 2)
    if templates:
        from dvbfixer.homology import get_template_chains
        from dvbfixer.salign import parse_template_spec

        for index, spec in enumerate(templates, start=1):
            pdb_path, selected_chain = parse_template_spec(spec)
            chains = get_template_chains(pdb_path)
            if selected_chain:
                if selected_chain not in chains:
                    raise ValueError(f"chain {selected_chain!r} not found in {pdb_path}")
                chains = {selected_chain: chains[selected_chain]}
            for chain, sequence in chains.items():
                expected.append((f"template_{index}_{pdb_path.stem}_{chain}", sequence))
    if len(expected) < 2:
        raise ValueError("multiple alignment requires at least two total sequences")
    paths = available_engines()
    if engine == "auto":
        engine = next((name for name in ENGINES if paths[name]), "")
        if not engine:
            raise RuntimeError("no MSA engine found; install MAFFT, MUSCLE 5, or Clustal Omega")
    executable = paths.get(engine)
    if not executable:
        raise RuntimeError(f"{engine} executable is not available")

    with tempfile.TemporaryDirectory(prefix="dvbfixer_msa_") as workdir:
        engine_input = Path(workdir) / "input.fasta"
        _write_fasta(expected, engine_input)
        raw_output = Path(workdir) / "aligned.fasta"
        if engine == "mafft":
            command = [executable, "--auto", str(engine_input)]
            result = subprocess.run(command, capture_output=True, text=True, check=False)
            if result.returncode == 0:
                raw_output.write_text(result.stdout)
        elif engine == "muscle":
            command = [executable, "-align", str(engine_input), "-output", str(raw_output)]
            result = subprocess.run(command, capture_output=True, text=True, check=False)
        else:
            command = [executable, "-i", str(engine_input), "-o", str(raw_output),
                       "--force", "--outfmt", "fa"]
            result = subprocess.run(command, capture_output=True, text=True, check=False)
        if verbose:
            print("Running:", " ".join(command))
            if result.stderr.strip():
                print(result.stderr.rstrip())
        if result.returncode != 0:
            detail = (result.stderr or result.stdout).strip()
            raise RuntimeError(f"{engine} failed with exit code {result.returncode}: {detail}")
        if not raw_output.exists() or raw_output.stat().st_size == 0:
            raise RuntimeError(f"{engine} did not produce an alignment")
        records = _validate_alignment(raw_output, expected)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_format == "pir":
        _write_pir(records, output_path)
    else:
        _write_fasta(records, output_path)
    return engine


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="dvbfixer msa",
        description="Align protein FASTA records with MAFFT, MUSCLE 5, or Clustal Omega.",
    )
    io = parser.add_argument_group("Input / output")
    io.add_argument("input", nargs="?", help="Input FASTA containing at least two records")
    io.add_argument("-o", "--output", help="Output alignment (default: <input>_aligned.fasta)")
    alignment = parser.add_argument_group("Alignment")
    alignment.add_argument("--engine", choices=["auto", *ENGINES], default="auto",
                           help="Alignment engine (default: auto; prefers MAFFT)")
    alignment.add_argument("--format", choices=["fasta", "pir"], default="fasta",
                           help="Output alignment format (default: fasta)")
    alignment.add_argument("--template", action="append", default=[],
                           help="Append sequences extracted from PDB[:CHAIN] (repeatable)")
    diagnostics = parser.add_argument_group("Diagnostics")
    diagnostics.add_argument("--list-engines", action="store_true",
                             help="Print detected alignment engines and exit")
    diagnostics.add_argument("-v", "--verbose", action="store_true")
    from dvbfixer.batch import add_runtime_help
    add_runtime_help(parser)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    if args.list_engines:
        for name, executable in available_engines().items():
            print(f"{name}: {executable or 'MISSING'}")
        return
    if not args.input:
        raise SystemExit("ERROR: input FASTA is required unless --list-engines is used")
    input_path = Path(args.input).resolve()
    if not input_path.exists():
        raise SystemExit(f"input FASTA not found: {input_path}")
    suffix = ".pir" if args.format == "pir" else ".fasta"
    output_path = Path(args.output).resolve() if args.output else input_path.with_name(
        f"{input_path.stem}_aligned{suffix}")
    try:
        selected = run_alignment(
            input_path, output_path, args.engine, args.format, args.verbose, args.template,
        )
    except (ValueError, RuntimeError) as exc:
        raise SystemExit(f"ERROR: {exc}") from exc
    print(f"Aligned with {selected}: {output_path}")


if __name__ == "__main__":
    main()
