#!/usr/bin/env python3
"""Run the internal MODELLER comparator with a locked diffusion placement."""

from __future__ import annotations

import argparse
from pathlib import Path

from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.pipeline import main as model_main
from dvbfixer.model.pipeline import parse_fasta


def _locked_placements(
    request: DiffusionRequest,
    fasta_path: Path,
) -> dict[str, tuple[int, ...]]:
    fasta_sequences = parse_fasta(fasta_path)
    target_sequences = {
        target.chain: target.sequence for target in request.target_sequences
    }
    if fasta_sequences != target_sequences:
        raise ValueError(
            "diffusion request target sequences do not exactly match the MODELLER FASTA"
        )
    return {
        placement.chain: placement.observed_target_indices
        for placement in request.sequence_placements
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("request", type=Path)
    parser.add_argument("model_args", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    if not args.model_args:
        parser.error("MODELLER command arguments are required")
    return args


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    request = DiffusionRequest.from_json(args.request.read_text(encoding="utf-8"))
    try:
        fasta_index = args.model_args.index("--fasta")
        fasta_path = Path(args.model_args[fasta_index + 1])
    except (ValueError, IndexError) as exc:
        raise ValueError("MODELLER comparator requires --fasta PATH") from exc
    model_main(
        args.model_args,
        observed_target_indices_by_chain=_locked_placements(request, fasta_path),
    )


if __name__ == "__main__":
    main()
