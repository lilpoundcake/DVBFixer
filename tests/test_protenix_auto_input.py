"""Production Protenix template construction at the full/motif token boundary.

These are input-transport tests, not checkpoint-backed sampling acceptance.
"""

import importlib.util
from dataclasses import replace
from pathlib import Path

import pytest
from Bio import SeqIO

from dvbfixer.model.diffusion.protein_context import (
    ProteinContextError,
    prepare_protenix_request,
    sampling_entities,
)
from dvbfixer.model.diffusion_cli import build_cli_diffusion_request

ROOT = Path(__file__).resolve().parents[1]


def make_input(directory: Path, total: int):
    fixture = ROOT / "tests/fixtures/protein_context/1UBQ.pdb"
    sequence = str(next(SeqIO.parse(fixture, "pdb-seqres")).seq)
    template = [line for line in fixture.read_text().splitlines(keepends=True) if line.startswith("ATOM  ")]
    sequences = {}
    lines = []
    serial = 1
    for index, start in enumerate(range(0, total, len(sequence))):
        chain = chr(ord("A") + index)
        length = min(total - start, len(sequence))
        sequences[chain] = sequence[:length]
        for line in template:
            number = int(line[22:26])
            if number > length or (chain == "A" and 10 <= number <= 13):
                continue
            # Distant translated copies provide a deterministic size control,
            # not a claim that these chains form a biological assembly.
            lines.append(f"{line[:6]}{serial:5d}{line[11:21]}{chain}{line[22:30]}"
                         f"{float(line[30:38]) + index * 80:8.3f}{line[38:]}")
            serial += 1
        lines.append("TER\n")
    path = directory / "input.pdb"
    path.write_text("".join(lines) + "END\n")
    request = build_cli_diffusion_request(path, sequences, seeds=(7,), profile="protenix-v1-cuda")
    return path, request


@pytest.mark.parametrize("total", [152, 512, 513])
def test_auto_full_or_motif_builds_real_protenix_templates(tmp_path, total):
    import json

    path, request = make_input(tmp_path, total)
    original = path.read_bytes()
    selected = prepare_protenix_request(request, original, "auto")
    axis = sum(len(entity.target.sequence) for entity in sampling_entities(selected))
    if total <= 512:
        assert selected.protein_context.policy.partner_mode == "full"
        assert axis == total
        assert selected == prepare_protenix_request(request, original, "full")
    else:
        assert selected.protein_context.policy.target_flanks == 8
        assert axis < 512
        assert selected == prepare_protenix_request(request, original, "local")
        with pytest.raises(ProteinContextError, match="supported budget is 512"):
            prepare_protenix_request(request, original, "full")
    request_path = tmp_path / "request.json"
    request_path.write_text(selected.to_json())
    spec = importlib.util.spec_from_file_location("protenix_template_acceptance", ROOT / "deploy/protenix-v1/build-template-input.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    model_input = module.build_input(request_path, tmp_path / "engine-input")
    payload, = json.loads(model_input.read_text())
    assert sum(len(row["proteinChain"]["sequence"]) for row in payload["sequences"]) == axis
    assert payload["modelSeeds"] == [7]
    for row in payload["sequences"]:
        template, = json.loads(Path(row["proteinChain"]["templatesPath"]).read_text())
        assert len(template["queryIndices"]) == len(template["templateIndices"])
        assert template["mmcif"]
    assert path.read_bytes() == original


def test_invalid_full_input_does_not_silently_fall_back_to_motif(tmp_path):
    import hashlib

    path, request = make_input(tmp_path, 152)
    source = b"".join(line for line in path.read_bytes().splitlines(keepends=True) if not (
        line.startswith(b"ATOM") and chr(line[21]) == "B" and line[22:26].strip() == b"20" and line[12:16].strip() == b"CB"))
    request = replace(request, normalized_pdb=replace(request.normalized_pdb, sha256=hashlib.sha256(source).hexdigest()),
                      fixed_atoms=tuple(atom for atom in request.fixed_atoms if not (
                          atom.chain == "B" and atom.residue_number == "20" and atom.atom_name == "CB")))
    with pytest.raises(ProteinContextError, match="heavy-atom-incomplete"):
        prepare_protenix_request(request, source, "auto")
