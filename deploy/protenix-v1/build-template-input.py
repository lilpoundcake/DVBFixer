"""Build a single-chain Protenix template input from a DiffusionRequest."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import gemmi

from dvbfixer.model.diffusion.contract import DiffusionRequest


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_input(request_path: Path, output_dir: Path) -> Path:
    request_path = request_path.resolve()
    request = DiffusionRequest.from_json(request_path.read_text(encoding="utf-8"))
    if len(request.target_sequences) != 1 or len(request.gaps) != 1:
        raise ValueError("Protenix hook smoke requires one target chain and one gap")
    target = request.target_sequences[0]
    placement = request.sequence_placements[0]
    if target.chain != placement.chain or request.gaps[0].chain != target.chain:
        raise ValueError("target, placement, and gap chains must match")
    if any(atom.chain != target.chain for atom in (*request.fixed_atoms, *request.generated_atoms)):
        raise ValueError("Protenix hook smoke accepts target-chain-only structures")

    pdb_path = (request_path.parent / request.normalized_pdb.path).resolve()
    if _sha256(pdb_path) != request.normalized_pdb.sha256:
        raise ValueError("normalized PDB digest mismatch")
    structure = gemmi.read_structure(str(pdb_path))
    structure.setup_entities()
    mmcif = structure.make_mmcif_document().as_string()

    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=False)
    template_path = output_dir / "template.json"
    template_path.write_text(
        json.dumps(
            [
                {
                    "mmcif": mmcif,
                    "queryIndices": list(placement.observed_target_indices),
                    "templateIndices": list(range(len(placement.observed_target_indices))),
                }
            ],
            separators=(",", ":"),
        ),
        encoding="utf-8",
    )
    input_path = output_dir / "input.json"
    input_path.write_text(
        json.dumps(
            [
                {
                    "name": "dvbfixer_protenix_gap_smoke",
                    "sequences": [
                        {
                            "proteinChain": {
                                "sequence": target.sequence,
                                "count": 1,
                                "templatesPath": str(template_path),
                            }
                        }
                    ],
                    "modelSeeds": list(request.seeds),
                }
            ],
            separators=(",", ":"),
        ),
        encoding="utf-8",
    )
    return input_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("request", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    print(build_input(args.request, args.output_dir))


if __name__ == "__main__":
    main()
