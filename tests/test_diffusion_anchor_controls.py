"""A complete experimental conformation is a no-minimization control on 1UBQ."""

from pathlib import Path

import pytest

pytest.importorskip("openmm")
pytest.importorskip("pdbfixer")

from dvbfixer.model.diffusion.anchor_acceptance import _assessment, _structure
from dvbfixer.model.diffusion.pdb_metadata import reconcile_pdb_metadata
from dvbfixer.model.diffusion.refinement_stages import run_refinement_stages
from dvbfixer.model.diffusion_cli import _target_sequences, build_cli_diffusion_request


@pytest.mark.parametrize("experimental_header", [True, False])
def test_1ubq_native_gap_control_is_unchanged_without_minimization(tmp_path, experimental_header):
    original = (Path(__file__).parent / "fixtures/protein_context/1UBQ.pdb").read_text()
    # A declared synthetic gap (A6-A8), not a diffusion accuracy benchmark. The
    # raw control restores the experimentally observed atom coordinates exactly.
    protein = "".join(line for line in original.splitlines(keepends=True)
                      if not line.startswith(("HETATM", "CONECT"))
                      and (experimental_header or not line.startswith("EXPDTA")))
    source_text = "".join(line for line in protein.splitlines(keepends=True)
                          if not (line.startswith("ATOM  ") and line[21] == "A" and
                                  line[22:26].strip() in {"6", "7", "8"}))
    path = tmp_path / "source.pdb"
    path.write_text(source_text)
    request = build_cli_diffusion_request(path, _target_sequences(path, None), seeds=(7,),
                                         profile="protpardelle-1c-mps-local-v1")
    raw = reconcile_pdb_metadata(protein, request).text
    source = _structure(source_text, "source.pdb")
    def assess(text):
        assessment, strict = _assessment(request, source, _structure(text, "candidate.pdb"))
        assert strict.passed, strict.hard_gate_failures
        return assessment
    def forbidden(*_args):
        pytest.fail("clean experimental control must not be minimized or warned")
    result = run_refinement_stages(raw, source_text, request, assess=assess,
                                   minimize=forbidden, notify=forbidden)
    assert result.selected_stage == "raw"
    assert result.coordinate_text == raw
    assert len(result.attempts) == 1
