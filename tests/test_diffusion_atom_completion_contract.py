"""Independent atom-completion protocol must not widen schema-4 runner support."""

import hashlib
import json
import os
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import pytest

from dvbfixer.model.diffusion.atom_completion import plan_atom_completion
from dvbfixer.model.diffusion.atom_completion_contract import (
    AtomCompletionCapabilities,
    AtomCompletionRequest,
    audit_completion_request_source,
    require_completion_capability,
)
from dvbfixer.model.diffusion.atom_completion_materialize import materialize_atom_completion
from dvbfixer.model.diffusion.contract import (
    ArtifactReference,
    AtomIdentity,
    DiffusionContractError,
    DiffusionRequest,
    GapKind,
    GapRegion,
    ResidueIdentity,
    SequencePlacement,
    TargetInterval,
    TargetSequence,
    TemplateOwnership,
)
from dvbfixer.model.diffusion.protein_atoms import canonical_heavy_atom_identities


@pytest.fixture
def completion_case():
    residue = ResidueIdentity("d", "82", "A")
    fixed = tuple(AtomIdentity("d", "82", "A", name) for name in ("N", "CA", "CB", "C", "O"))
    source = ("".join(
        f"ATOM  {serial:5d} {atom.atom_name:>4s} LYS d  82A   "
        f"{serial * 1.1:8.3f}{0:8.3f}{0:8.3f}  1.00 20.00          {atom.atom_name[0]:>2s}\n"
        for serial, atom in enumerate(fixed, 1)
    ) + "END\n").encode()
    digest = hashlib.sha256(source).hexdigest()
    targets = (TargetSequence("d", "K"),)
    placements = (SequencePlacement("d", 1, (0,), (residue,)),)
    masks = plan_atom_completion(source, source_sha256=digest, target_sequences=targets,
                                 sequence_placements=placements, fixed_atoms=fixed)
    request = AtomCompletionRequest(1, ArtifactReference("input.pdb", digest), targets, placements,
                                    (), masks, (7,), "experimental-test-only")
    return source, request


def _capability(request, **kwargs):
    return AtomCompletionCapabilities(1, request.profile, sidechain=True, max_target_residues=10,
                                      max_generated_atoms=100, max_candidates=2, **kwargs)


def test_gapless_roundtrip_and_source_audit(completion_case):
    source, request = completion_case
    assert AtomCompletionRequest.from_json(request.to_json()) == request
    assert request.to_dict()["schema_version"] == 1
    assert len(request.generated_atoms) == 4
    audit_completion_request_source(request, source)
    require_completion_capability(request, _capability(request))
    with pytest.raises(DiffusionContractError):
        DiffusionRequest.from_json(request.to_json())
    with pytest.raises(ValueError, match="digest mismatch"):
        audit_completion_request_source(request, source + b"\n")


def test_mixed_request_requires_complete_gap_ownership(completion_case):
    source, request = completion_case
    anchor = request.sequence_placements[0].observed_residues[0]
    generated = ResidueIdentity("d", "83")
    gap = GapRegion("d", TargetInterval(1, 2), anchor, None, (generated,), (anchor, generated), GapKind.C_TERMINAL)
    atoms = canonical_heavy_atom_identities(generated, "GLY")
    mixed = replace(request, target_sequences=(TargetSequence("d", "KG"),),
                    sequence_placements=(replace(request.sequence_placements[0], target_length=2),),
                    gaps=(gap,), masks=replace(request.masks, gap_generated_atoms=atoms))
    audit_completion_request_source(mixed, source)
    with pytest.raises(DiffusionContractError, match="mixed"):
        require_completion_capability(mixed, _capability(mixed, whole_gaps=True))
    require_completion_capability(mixed, _capability(mixed, whole_gaps=True, mixed_gaps=True))
    with pytest.raises(DiffusionContractError, match="exactly every canonical"):
        replace(mixed, masks=replace(mixed.masks, gap_generated_atoms=atoms[:-1]))
    with pytest.raises(DiffusionContractError, match="adjacent observed"):
        replace(mixed, gaps=(replace(gap, left_anchor=ResidueIdentity("d", "81"),
                                    movable_junction_residues=(ResidueIdentity("d", "81"), generated)),))


@pytest.mark.parametrize("change,match", [
    ({"schema_version": 4}, "schema_version"),
    ({"profile": " "}, "explicit runner profile"),
    ({"seeds": ()}, "unique seeds"),
    ({"seeds": (7, 7)}, "unique seeds"),
    ({"seeds": (-1,)}, "32-bit"),
    ({"seeds": (2**32,)}, "32-bit"),
    ({"seeds": (True,)}, "32-bit"),
])
def test_invalid_request_refused(completion_case, change, match):
    with pytest.raises(DiffusionContractError, match=match):
        replace(completion_case[1], **change)


def test_masks_cannot_forge_target_identity_or_source_binding(completion_case):
    _, request = completion_case
    partial, = request.masks.partial_residues
    with pytest.raises(DiffusionContractError, match="observed target placement"):
        replace(request, masks=replace(request.masks, partial_residues=(replace(partial, target_index=1),)))
    with pytest.raises(DiffusionContractError, match="source digest"):
        replace(request, masks=replace(request.masks, source_sha256="a" * 64))
    with pytest.raises(DiffusionContractError, match="differs from target sequence"):
        replace(request, target_sequences=(TargetSequence("d", "R"),))
    with pytest.raises(DiffusionContractError, match="every target position"):
        replace(request, target_sequences=(TargetSequence("d", "KK"),),
                sequence_placements=(replace(request.sequence_placements[0], target_length=2),))
    with pytest.raises(DiffusionContractError, match="template ownership"):
        replace(request, template_ownership=(TemplateOwnership("d", 1, partial.residue, "template"),))


def test_capability_defaults_deny_and_limits_are_enforced(completion_case):
    _, request = completion_case
    default = AtomCompletionCapabilities(1, request.profile)
    assert AtomCompletionCapabilities.from_json(default.to_json()) == default
    with pytest.raises(DiffusionContractError, match="does not support"):
        require_completion_capability(request, default)
    for change in ({"max_target_residues": 0}, {"max_generated_atoms": 3}, {"max_candidates": 0}):
        with pytest.raises(DiffusionContractError, match="resource bounds"):
            require_completion_capability(request, replace(_capability(request), **change))
    with pytest.raises(DiffusionContractError, match="profile mismatch"):
        require_completion_capability(request, replace(_capability(request), profile="another"))
    for change in ({"sidechain": 1}, {"max_target_residues": True}, {"max_candidates": -1}):
        with pytest.raises(DiffusionContractError):
            replace(default, **change)


@pytest.mark.parametrize("mutation", [
    lambda data: data.update(unknown=True),
    lambda data: data.update(schema_version=True),
    lambda data: data["masks"].update(unknown=True),
    lambda data: data.update(seeds=[True]),
])
def test_strict_wire_decoder_rejects_unknown_fields_and_wrong_types(completion_case, mutation):
    data = completion_case[1].to_dict()
    mutation(data)
    with pytest.raises(DiffusionContractError):
        AtomCompletionRequest.from_json(json.dumps(data))


def test_contract_import_is_stdlib_only():
    root = Path(__file__).resolve().parents[1]
    env = dict(os.environ, PYTHONPATH=str(root / "src"))
    subprocess.run([sys.executable, "-S", "-c", """
import sys
from dvbfixer.model.diffusion.atom_completion_contract import AtomCompletionRequest
assert not any(name in sys.modules for name in ('openmm', 'numpy', 'Bio', 'torch'))
"""], env=env, check=True, capture_output=True, text=True)


def test_materializer_inserts_only_missing_atoms_and_preserves_source_bytes(completion_case):
    source, request = completion_case
    axis = request.generated_atoms
    result = materialize_atom_completion(source, request, axis, [(1, 2, 3)] * len(axis))
    observations = source.splitlines(keepends=True)[:-1]
    assert result.splitlines(keepends=True)[:len(observations)] == observations
    generated = result.splitlines()[len(observations):-1]
    assert {line[12:16].strip().decode() for line in generated} == {"CG", "CD", "CE", "NZ"}
    assert all(line[17:27] == b"LYS d  82A" for line in generated)
    assert len({int(line[6:11]) for line in result.splitlines() if line.startswith(b"ATOM")}) == 9
    with pytest.raises(DiffusionContractError, match="exactly generated"):
        materialize_atom_completion(source, request, axis + request.masks.fixed_atoms[:1], [(1, 2, 3)] * 5)
    with pytest.raises(DiffusionContractError, match="exactly generated"):
        materialize_atom_completion(source, request, axis[:-1], [(1, 2, 3)] * 3)


@pytest.mark.parametrize("xyz", [(float("nan"), 0, 0), (float("inf"), 0, 0), (10000, 0, 0), (-1000, 0, 0), (1, 2)])
def test_materializer_rejects_bad_coordinates(completion_case, xyz):
    source, request = completion_case
    with pytest.raises(DiffusionContractError, match="coordinate"):
        materialize_atom_completion(source, request, request.generated_atoms, [xyz] * 4)


def _rebind(source, request):
    digest = hashlib.sha256(source).hexdigest()
    return replace(request, normalized_pdb=replace(request.normalized_pdb, sha256=digest),
                   masks=replace(request.masks, source_sha256=digest))


def test_crlf_anisou_and_valid_conect_are_preserved(completion_case):
    source, request = completion_case
    lines = source.splitlines(keepends=True)
    anisou = b"ANISOU" + lines[-2][6:27] + b"    100    100    100      0      0      0       O\n"
    source = (b"".join(lines[:-1]) + anisou + b"CONECT    1    2\nEND\n").replace(b"\n", b"\r\n")
    request = _rebind(source, request)
    result = materialize_atom_completion(source, request, request.generated_atoms, [(1, 2, 3)] * 4)
    assert anisou.replace(b"\n", b"\r\n") in result
    assert b"CONECT    1    2\r\n" in result
    assert result.count(b"\n") == result.count(b"\r\n")
    output = result.splitlines()
    anisou_index = next(i for i, line in enumerate(output) if line.startswith(b"ANISOU"))
    assert output[anisou_index - 1] == lines[-2].rstrip(b"\n")
    assert output[anisou_index + 1][12:16].strip() == b"CG"


@pytest.mark.parametrize("change,match", [
    (lambda source: b"MODEL        1\n" + source, "single-model"),
    (lambda source: source.replace(b"ATOM      2", b"ATOM      1"), "serials"),
    (lambda source: source.replace(b"END\n", b"CONECT    1 9999\nEND\n"), "CONECT"),
    (lambda source: source.replace(b"END\n", b"ANISOU    9  CG  LYS d  82A\nEND\n"), "ANISOU"),
    (lambda source: source.replace(b"ATOM      2", b"TER\nATOM      2"), "split or repeated"),
])
def test_ambiguous_source_layout_refused_before_insertion(completion_case, change, match):
    source, request = completion_case
    source = change(source)
    request = _rebind(source, request)
    with pytest.raises(DiffusionContractError, match=match):
        materialize_atom_completion(source, request, request.generated_atoms, [(1, 2, 3)] * 4)


def test_mixed_materialization_reuses_real_gap_without_fixed_updates(completion_case):
    source, request = completion_case
    anchor = request.sequence_placements[0].observed_residues[0]
    generated = ResidueIdentity("d", "83")
    gap = GapRegion("d", TargetInterval(1, 2), anchor, None, (generated,), (anchor, generated), GapKind.C_TERMINAL)
    request = replace(request, target_sequences=(TargetSequence("d", "KG"),),
                      sequence_placements=(replace(request.sequence_placements[0], target_length=2),),
                      gaps=(gap,), masks=replace(request.masks, gap_generated_atoms=canonical_heavy_atom_identities(generated, "GLY")))
    result = materialize_atom_completion(source, request, request.generated_atoms, [(1, 2, 3)] * 8)
    assert all(line in result for line in source.splitlines(keepends=True) if line.startswith(b"ATOM"))
    atoms = [line for line in result.splitlines() if line.startswith(b"ATOM")]
    assert len(atoms) == 13
    assert [line[22:27] for line in atoms] == [b"  82A"] * 9 + [b"  83 "] * 4
    assert len({int(line[6:11]) for line in atoms}) == 13


def test_stale_master_removed_and_post_end_records_refused(completion_case):
    source, request = completion_case
    source = source.replace(b"END\n", b"MASTER        0    0    0    0    0    0    0    0    5    0    0    0\nEND\n")
    request = _rebind(source, request)
    result = materialize_atom_completion(source, request, request.generated_atoms, [(1, 2, 3)] * 4)
    assert b"MASTER" not in result
    source += b"REMARK invalid record after END\n"
    with pytest.raises(DiffusionContractError, match="after END"):
        materialize_atom_completion(source, _rebind(source, request), request.generated_atoms, [(1, 2, 3)] * 4)


def test_backbone_support_is_not_inferred_from_sidechain_capability(completion_case):
    source, request = completion_case
    source = b"".join(source.splitlines(keepends=True)[1:])  # Remove observed N, never rebuild it here.
    digest = hashlib.sha256(source).hexdigest()
    masks = plan_atom_completion(source, source_sha256=digest, target_sequences=request.target_sequences,
                                 sequence_placements=request.sequence_placements, fixed_atoms=request.masks.fixed_atoms[1:])
    request = replace(request, normalized_pdb=replace(request.normalized_pdb, sha256=digest), masks=masks)
    with pytest.raises(DiffusionContractError, match="backbone"):
        require_completion_capability(request, _capability(request))
    require_completion_capability(request, _capability(request, backbone=True))
