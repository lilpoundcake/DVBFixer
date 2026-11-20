"""Build a shared-frame Boltz-2 protein template from a DiffusionRequest."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from atom_mapping import ONE_TO_THREE, target_residue_map

from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _seqres_lines(chain: str, sequence: str) -> list[str]:
    names = [ONE_TO_THREE[symbol] for symbol in sequence]
    return [
        f"SEQRES {row:3d} {chain} {len(names):4d}  {' '.join(names[start:start + 13])}\n"
        for row, start in enumerate(range(0, len(names), 13), start=1)
    ]


def build_input(request_path: Path, output_dir: Path) -> Path:
    request_path = request_path.resolve()
    request = DiffusionRequest.from_json(request_path.read_text(encoding="utf-8"))
    targets = request.target_sequences
    for target in targets:
        target_residue_map(request,target.chain)
    if request.candidate_count != 1 or len(request.seeds) != 1:
        raise ValueError("Boltz checkpoint smoke requires one candidate and one seed")

    source_path = (request_path.parent / request.normalized_pdb.path).resolve()
    if _sha256(source_path) != request.normalized_pdb.sha256:
        raise ValueError("normalized PDB digest mismatch")
    source_lines = source_path.read_text(encoding="ascii").splitlines(keepends=True)
    observed_residues = {r for p in request.sequence_placements for r in p.observed_residues}
    source_atoms: set[AtomIdentity] = set()
    source_residues = set()
    for line in source_lines:
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if line[16:17] not in {"", " ", "A"}:
            continue
        identity = AtomIdentity(
            line[21:22], line[22:26].strip(), line[26:27].strip(), line[12:16].strip()
        )
        if identity in source_atoms:
            raise ValueError("normalized PDB contains duplicate atom identities")
        source_atoms.add(identity)
        source_residues.add((identity.chain, identity.residue_number, identity.insertion_code))
    expected_observed = {
        (residue.chain, residue.residue_number, residue.insertion_code)
        for residue in observed_residues
    }
    if source_residues != expected_observed:
        raise ValueError("normalized PDB residues do not match the request placement")
    if set(request.fixed_atoms) != source_atoms:
        raise ValueError("normalized PDB atom identities do not match fixed_atoms")

    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=False)
    template_path = output_dir / "template.pdb"
    body = [line for line in source_lines if not line.startswith("SEQRES")]
    template_path.write_text("".join([*(line for t in targets for line in _seqres_lines(t.chain,t.sequence)), *body]), encoding="ascii")

    input_path = output_dir / "input.yaml"
    input_path.write_text(
        json.dumps(
            {
                "sequences": [
                    {
                        "protein": {
                            "id": target.chain,
                            "sequence": target.sequence,
                            "msa": "empty",
                        }
                    } for target in targets
                ],
                "templates": [
                    {
                        "pdb": str(template_path),
                        "chain_id": [t.chain for t in targets] if len(targets)>1 else targets[0].chain,
                        "template_id": [f"{t.chain}1" for t in targets] if len(targets)>1 else f"{targets[0].chain}1",
                    }
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
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
