"""CLI for authoritative fixed-anchor Class A component reconstruction."""

from __future__ import annotations

import argparse
from pathlib import Path

from dvbfixer.ccd import CcdError, load_ccd_component
from dvbfixer.pdb_component_reconstruction import (
    PdbComponentError,
    publish_pdb_bundle,
    reconstruct_pdb_component,
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="dvbfixer reconstruct-component",
        description=(
            "Reconstruct missing heavy atoms in one isolated Class A HETATM "
            "component from an explicitly selected local wwPDB CCD CIF."
        ),
    )
    parser.add_argument("input", help="Input PDB (source is never modified)")
    parser.add_argument("component", help="Exact CCD/PDB component ID")
    parser.add_argument("--ccd-cif", required=True, help="Pinned local CCD component CIF")
    parser.add_argument("--ccd-sha256", required=True,
                        help="Required lowercase SHA-256 pin for --ccd-cif")
    parser.add_argument("--chain", required=True, help="Exact case-sensitive one-character chain ID")
    parser.add_argument("--residue", required=True, help="Exact PDB residue sequence number")
    parser.add_argument("--icode", default="", help="Exact insertion code (default: blank)")
    parser.add_argument("--model", type=int, default=1, help="Exact MODEL serial (default: 1)")
    parser.add_argument("--occurrence", type=int, default=1,
                        help="Occurrence when an exact residue identity repeats (default: 1)")
    parser.add_argument("--ph", type=float, default=7.0,
                        help="Recorded pH; CCD state is not changed (default: 7.0)")
    parser.add_argument("--output-root", default=".", help="Existing publication directory")
    parser.add_argument("--bundle-name", required=True,
                        help="Name of the new atomic output bundle directory")
    from dvbfixer.batch import add_runtime_help

    add_runtime_help(parser)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    try:
        if len(args.chain) != 1:
            raise PdbComponentError("--chain must be exactly one case-sensitive character")
        if len(args.icode) > 1:
            raise PdbComponentError("--icode must be blank or one character")
        if args.model < 1 or args.occurrence < 1:
            raise PdbComponentError("--model and --occurrence must be positive")
        input_path = Path(args.input).expanduser().resolve(strict=True)
        component = load_ccd_component(
            Path(args.ccd_cif),
            args.component,
            expected_sha256=args.ccd_sha256,
        )
        bundle, output = reconstruct_pdb_component(
            input_path.read_bytes(),
            component,
            model=args.model,
            chain_id=args.chain,
            sequence_number=args.residue,
            insertion_code=args.icode,
            occurrence=args.occurrence,
            ph=args.ph,
        )
        if output is None:
            finding = bundle.result.findings[0] if bundle.result.findings else None
            detail = f" [{finding.code}] {finding.message}" if finding else ""
            raise PdbComponentError(
                f"reconstruction {bundle.result.status.value} as Class "
                f"{bundle.result.component_class.value}:{detail}"
            )
        destination = publish_pdb_bundle(
            bundle,
            output,
            Path(args.output_root),
            args.bundle_name,
        )
    except (CcdError, PdbComponentError, FileExistsError, OSError) as exc:
        raise SystemExit(f"ERROR: {exc}") from exc
    print(f"Published {destination}")
    print(f"Added heavy atoms: {len(bundle.added_atom_names)}")
    if bundle.added_atom_names:
        print("Atoms: " + ", ".join(bundle.added_atom_names))
