"""Versioned union of the existing anchor shells and a finite CA–CA capsule."""

import hashlib
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from dvbfixer.model.diffusion.contract import ArtifactReference, BackendOption
from dvbfixer.model.diffusion.protein_context import (
    ANCHOR_SEGMENT_CONTEXT,
    ANCHOR_SEGMENT_OPTION,
    CUDA_CC89_CONTEXT,
    CUDA_CC89_PROFILE,
    MPS_CC89_PROFILE,
    TERMINAL_SHELL_CONTEXT,
    TERMINAL_SHELL_OPTION,
    ProteinContextError,
    _point_segment_distances,
    prepare_mps_cc89_request,
    select_protein_context,
    validate_context_source,
    validate_mps_cc89_request,
)
from tests.test_protein_context import context_request


def test_segment_distances_clamp_to_endpoints_and_handle_zero_length():
    points = np.array([[5., 10., 0.], [-5., 0., 0.], [15., 0., 0.]])
    np.testing.assert_allclose(_point_segment_distances(points, np.zeros(3), np.array([10., 0., 0.])), [10., 5., 5.])
    np.testing.assert_allclose(_point_segment_distances(points, np.zeros(3), np.zeros(3)), np.linalg.norm(points, axis=1))


def test_segment_adds_midpoint_partner_without_removing_anchor_shell_atoms():
    request, pdb = context_request(partner_y={1: 4, 2: 100, 3: 100})
    lines = []
    for line in pdb.decode().splitlines(keepends=True):
        if line.startswith('ATOM'):
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            if line[21] == 'D' and int(line[22:26]) == 5:
                x += 48  # CA endpoints at 0 and 60 Å.
            if line[21] == 'd' and int(line[22:26]) == 2:
                x += 27
                y -= 90.1  # 9.9 Å from the segment, far from both anchors.
            line = line[:30] + f'{x:8.3f}{y:8.3f}{z:8.3f}' + line[54:]
        lines.append(line)
    pdb = ''.join(lines).encode()
    request = replace(request, normalized_pdb=ArtifactReference('input.pdb', hashlib.sha256(pdb).hexdigest()))
    old = select_protein_context(request, pdb, replace(CUDA_CC89_CONTEXT, sequence_padding=0))
    new = select_protein_context(request, pdb, replace(ANCHOR_SEGMENT_CONTEXT, sequence_padding=0))
    old_atoms = {item.atom for item in old.protein_context.fixed_conditioning}
    new_atoms = {item.atom for item in new.protein_context.fixed_conditioning}
    assert old_atoms < new_atoms
    assert not any(atom.chain == 'd' and atom.residue_number == '2' for atom in old_atoms)
    assert any(atom.chain == 'd' and atom.residue_number == '2' and atom.insertion_code == 'A' for atom in new_atoms)
    assert not any(atom.chain == 'd' and atom.residue_number == '3' for atom in new_atoms)
    validate_context_source(new, pdb)


def test_corrected_profile_requires_explicit_segment_option_and_preserves_default():
    request, pdb = context_request()
    request = replace(request, backend_options=(BackendOption('profile', MPS_CC89_PROFILE),))
    baseline = prepare_mps_cc89_request(request, pdb, 'auto')
    assert baseline.protein_context.policy == CUDA_CC89_CONTEXT
    expanded = prepare_mps_cc89_request(replace(request, backend_options=(*request.backend_options, ANCHOR_SEGMENT_OPTION)), pdb, 'auto')
    assert expanded.protein_context.policy == ANCHOR_SEGMENT_CONTEXT
    validate_mps_cc89_request(expanded)
    validate_context_source(expanded, pdb)
    with pytest.raises(ProteinContextError, match='frozen source-only'):
        validate_mps_cc89_request(replace(expanded, backend_options=baseline.backend_options))


@pytest.mark.parametrize('profile', [CUDA_CC89_PROFILE, MPS_CC89_PROFILE, 'protpardelle-1c-mps', 'protenix-v1-cuda'])
def test_new_cli_requests_declare_union_only_for_corrected_cc89(tmp_path, profile):
    from dvbfixer.model.diffusion_cli import build_cli_diffusion_request

    _, pdb = context_request()
    source = tmp_path / 'input.pdb'
    source.write_bytes(pdb)
    request = build_cli_diffusion_request(source, {'D': 'AAAAA', 'd': 'AAA'}, seeds=(7,), profile=profile)
    assert (ANCHOR_SEGMENT_OPTION in request.backend_options) == (profile in {CUDA_CC89_PROFILE, MPS_CC89_PROFILE})
    assert (TERMINAL_SHELL_OPTION in request.backend_options) == (profile in {CUDA_CC89_PROFILE, MPS_CC89_PROFILE})


def test_terminal_views_keep_anchor_shells_and_frozen_requests_do_not_change():
    from tests.test_diffusion_independent_gaps import prepared

    frozen, pdb = prepared()
    frozen_bytes = frozen.to_json()
    assert prepare_mps_cc89_request(frozen, pdb, 'auto').to_json() == frozen_bytes
    request = replace(frozen, gap_contexts=None, backend_options=(*frozen.backend_options, ANCHOR_SEGMENT_OPTION))
    union = prepare_mps_cc89_request(request, pdb, 'auto')
    assert union.gap_contexts[-1].policy == ANCHOR_SEGMENT_CONTEXT
    assert frozen.gap_contexts[-1].spans == union.gap_contexts[-1].spans
    assert frozen.gap_contexts[-1].fixed_conditioning == union.gap_contexts[-1].fixed_conditioning
    for old, new in zip(frozen.gap_contexts, union.gap_contexts, strict=True):
        assert {item.atom for item in old.fixed_conditioning} <= {item.atom for item in new.fixed_conditioning}
    validate_context_source(union, pdb)


def test_union_budget_failure_never_discards_contacts():
    request, pdb = context_request()
    with pytest.raises(ProteinContextError, match='no neighbors were silently dropped'):
        select_protein_context(request, pdb, replace(ANCHOR_SEGMENT_CONTEXT, residue_budget=6))


def test_new_5j13_internal_request_includes_the_previously_omitted_surface(tmp_path):
    from dvbfixer.model.diffusion_cli import _target_sequences, build_cli_diffusion_request

    source = Path(__file__).parent / 'fixtures/protein_context/5j13_selected.pdb'
    normalized = tmp_path / 'protein.pdb'
    normalized.write_text(''.join(line for line in source.read_text().splitlines(keepends=True)
        if not line.startswith(('HETATM', 'CONECT', 'LINK  '))))
    request = build_cli_diffusion_request(normalized, _target_sequences(normalized, None),
        seeds=(7,), profile=MPS_CC89_PROFILE, no_terminal=True)
    request = prepare_mps_cc89_request(request, normalized.read_bytes(), 'auto')
    assert len(request.gaps) == 4
    assert sum(len(g.generated_residues) for g in request.gaps) == 27
    for gap, plan in zip(request.gaps, request.gap_contexts, strict=True):
        assert sum(span.stop - span.start for span in plan.spans) <= 512
        if gap.chain == 'A' and len(gap.generated_residues) == 13:
            # This fixture keeps original author numbers; renumbered A77/79/85
            # correspond to the same selected source identities below.
            legacy = select_protein_context(replace(request, gaps=(gap,), gap_contexts=None,
                generated_atoms=tuple(a for a in request.generated_atoms
                    if any(a.chain == r.chain and a.residue_number == r.residue_number
                        and a.insertion_code == r.insertion_code for r in gap.generated_residues))),
                normalized.read_bytes(), CUDA_CC89_CONTEXT)
            old_atoms = {item.atom for item in legacy.protein_context.fixed_conditioning}
            new_atoms = {item.atom for item in plan.fixed_conditioning}
            assert old_atoms < new_atoms
            assert len(new_atoms - old_atoms) == 84
    validate_context_source(request, normalized.read_bytes())


def test_terminal20_requires_recorded_options_and_does_not_change_internal_selection():
    from tests.test_diffusion_independent_gaps import prepared

    old, pdb = prepared()
    baseline = prepare_mps_cc89_request(replace(old, gap_contexts=None,
        backend_options=(*old.backend_options, ANCHOR_SEGMENT_OPTION)), pdb, 'auto')
    expanded = prepare_mps_cc89_request(replace(baseline, gap_contexts=None,
        backend_options=(*baseline.backend_options, TERMINAL_SHELL_OPTION)), pdb, 'auto')
    assert expanded.gap_contexts[0].fixed_conditioning == baseline.gap_contexts[0].fixed_conditioning
    assert expanded.gap_contexts[0].spans == baseline.gap_contexts[0].spans
    assert all(plan.policy == TERMINAL_SHELL_CONTEXT for plan in expanded.gap_contexts)
    for prior, plan in zip(baseline.gap_contexts, expanded.gap_contexts, strict=True):
        assert {item.atom for item in prior.fixed_conditioning} <= {item.atom for item in plan.fixed_conditioning}
    validate_context_source(expanded, pdb)
    invalid = replace(old, gap_contexts=None, backend_options=(*old.backend_options, TERMINAL_SHELL_OPTION))
    with pytest.raises(ProteinContextError, match='matching source geometry'):
        prepare_mps_cc89_request(invalid, pdb, 'auto')


def test_5j13_terminal20_context_fits_budget_and_matches_preflight_counts(tmp_path):
    from dvbfixer.model.diffusion_cli import _target_sequences, build_cli_diffusion_request

    source = Path(__file__).parent / 'fixtures/protein_context/5j13_selected.pdb'
    normalized = tmp_path / 'protein.pdb'
    normalized.write_text(''.join(line for line in source.read_text().splitlines(keepends=True)
        if not line.startswith(('HETATM', 'CONECT', 'LINK  '))))
    original = build_cli_diffusion_request(normalized, _target_sequences(normalized, None), seeds=(7,), profile=MPS_CC89_PROFILE)
    assert TERMINAL_SHELL_OPTION in original.backend_options
    expanded = prepare_mps_cc89_request(original, normalized.read_bytes(), 'auto')
    counts = [sum(span.stop - span.start for span in plan.spans) for plan in expanded.gap_contexts]
    assert counts == [167, 77, 85, 267, 310, 184, 309, 67, 51, 180]
    assert sum(len(gap.generated_residues) for gap in expanded.gaps) == 124
    assert expanded.fixed_atoms == original.fixed_atoms
    validate_context_source(expanded, normalized.read_bytes())


def test_new_fasta_request_uses_default_terminal_union_without_extra_options(tmp_path):
    from dvbfixer.model.diffusion_cli import _target_sequences, build_cli_diffusion_request

    fixture = Path(__file__).parent / 'fixtures/protein_context/5j13_selected.pdb'
    fasta = fixture.with_suffix('.fasta')
    before_pdb, before_fasta = fixture.read_bytes(), fasta.read_bytes()
    source = tmp_path / 'input.pdb'
    source.write_text(''.join(line for line in before_pdb.decode().splitlines(keepends=True)
        if not line.startswith(('HETATM', 'CONECT', 'LINK  '))))
    request = build_cli_diffusion_request(source, _target_sequences(source, str(fasta)), seeds=(7,), profile=MPS_CC89_PROFILE)
    assert TERMINAL_SHELL_OPTION in request.backend_options
    prepared = prepare_mps_cc89_request(request, source.read_bytes(), 'auto')
    assert {target.chain: len(target.sequence) for target in prepared.target_sequences} == {'A': 134, 'B': 213, 'C': 226}
    assert len(prepared.gaps) == 7
    assert sum(len(gap.generated_residues) for gap in prepared.gaps) == 42
    assert [sum(span.stop - span.start for span in plan.spans) for plan in prepared.gap_contexts] == [154, 77, 85, 267, 67, 51, 174]
    assert all(plan.policy == TERMINAL_SHELL_CONTEXT for plan in prepared.gap_contexts)
    assert prepare_mps_cc89_request(prepared, source.read_bytes(), 'auto').to_json() == prepared.to_json()
    assert fixture.read_bytes() == before_pdb and fasta.read_bytes() == before_fasta
    validate_context_source(prepared, source.read_bytes())
