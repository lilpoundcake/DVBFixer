"""Atom inventories must never turn a partial residue into a whole-residue gap."""

import hashlib
from dataclasses import replace

import pytest

from dvbfixer.model.diffusion.atom_completion import (
    AtomCompletionError,
    CompletionKind,
    PartialResidueMask,
    audit_atom_completion_masks,
    plan_atom_completion,
    plan_request_atom_completion,
)
from dvbfixer.model.diffusion.contract import (
    ArtifactReference,
    AtomIdentity,
    DiffusionRequest,
    GapKind,
    GapRegion,
    ResidueIdentity,
    SequencePlacement,
    TargetInterval,
    TargetSequence,
)
from dvbfixer.model.diffusion.scope import CANONICAL_HEAVY_ATOMS, canonical_heavy_atom_identities


def _source(residues):
    lines = []
    fixed = []
    for residue, name, atom_names in residues:
        for atom_name in atom_names:
            atom = AtomIdentity(residue.chain, residue.residue_number, residue.insertion_code, atom_name)
            serial = len(lines) + 1
            lines.append(f"ATOM  {serial:5d} {atom_name:>4s} {name:3s} {residue.chain}"
                         f"{int(residue.residue_number):4d}{residue.insertion_code or ' '}   "
                         f"{serial * 1.1:8.3f}{0:8.3f}{0:8.3f}  1.00 20.00          {atom_name[0]:>2s}\n")
            fixed.append(atom)
    return ("".join(lines) + "END\n").encode(), tuple(fixed)


def _request():
    first, second, third = ResidueIdentity("D", "82"), ResidueIdentity("D", "82", "A"), ResidueIdentity("d", "82")
    source, fixed = _source([(r, "LYS", ("N", "CA", "CB", "C", "O")) for r in (first, second, third)])
    gap = ResidueIdentity("D", "83")
    request = DiffusionRequest(schema_version=4,
        normalized_pdb=ArtifactReference("input.pdb", hashlib.sha256(source).hexdigest()),
        target_sequences=(TargetSequence("D", "KKG"), TargetSequence("d", "K")),
        sequence_placements=(SequencePlacement("D", 3, (0, 1), (first, second)), SequencePlacement("d", 1, (0,), (third,))),
        gaps=(GapRegion("D", TargetInterval(2, 3), second, None, (gap,), (second, gap), GapKind.C_TERMINAL),),
        fixed_atoms=fixed, generated_atoms=canonical_heavy_atom_identities(gap, "GLY"), retained_explicit_links=(),
        candidate_count=1, seeds=(7,))
    return source, request


def test_inventory_keeps_case_insertion_and_whole_gap_ownership_separate():
    source, request = _request()
    before = request.to_json()
    plan = plan_request_atom_completion(request, source)
    assert [item.residue for item in plan.partial_residues] == [ResidueIdentity("D", "82"), ResidueIdentity("D", "82", "A"), ResidueIdentity("d", "82")]
    assert [item.target_index for item in plan.partial_residues] == [0, 1, 0]
    assert all(item.kind is CompletionKind.SIDECHAIN for item in plan.partial_residues)
    assert all({a.atom_name for a in item.missing_atoms} == {"CG", "CD", "CE", "NZ"} for item in plan.partial_residues)
    assert plan.fixed_atoms == request.fixed_atoms
    assert plan.gap_generated_atoms == request.generated_atoms
    assert set(plan.missing_atoms).isdisjoint(request.fixed_atoms + request.generated_atoms)
    assert request.to_json() == before  # No frozen request/mask mutation.
    audit_atom_completion_masks(plan, request, source)


@pytest.mark.parametrize("name", sorted(CANONICAL_HEAVY_ATOMS))
def test_gapless_inventory_uses_the_shared_canonical_atom_catalog(name):
    from dvbfixer.model.diffusion.scope import _CANONICAL_ONE_TO_THREE

    residue = ResidueIdentity("A", "10", "B")
    canonical = CANONICAL_HEAVY_ATOMS[name]
    source, fixed = _source([(residue, name, canonical[:-1])])
    one = next(key for key, value in _CANONICAL_ONE_TO_THREE.items() if value == name)
    plan = plan_atom_completion(source, source_sha256=hashlib.sha256(source).hexdigest(),
        target_sequences=(TargetSequence("A", one),), sequence_placements=(SequencePlacement("A", 1, (0,), (residue,)),),
        fixed_atoms=fixed)
    partial, = plan.partial_residues
    assert partial.missing_atoms == (AtomIdentity("A", "10", "B", canonical[-1]),)
    assert set(partial.observed_atoms + partial.missing_atoms) == set(canonical_heavy_atom_identities(residue, name))


def test_backbone_damage_is_classified_not_granted_sampling_capability():
    residue = ResidueIdentity("A", "1")
    source, fixed = _source([(residue, "ALA", ("CA", "CB", "C", "O"))])
    plan = plan_atom_completion(source, source_sha256=hashlib.sha256(source).hexdigest(),
        target_sequences=(TargetSequence("A", "A"),), sequence_placements=(SequencePlacement("A", 1, (0,), (residue,)),),
        fixed_atoms=fixed)
    partial, = plan.partial_residues
    assert partial.kind is CompletionKind.BACKBONE
    assert not hasattr(plan, "supported")


def test_independent_audit_rejects_forged_target_ownership():
    source, request = _request()
    plan = plan_request_atom_completion(request, source)
    changed = replace(plan.partial_residues[0], target_index=5)
    forged = replace(plan, partial_residues=(changed, *plan.partial_residues[1:]))
    with pytest.raises(AtomCompletionError, match="independent source inventory"):
        audit_atom_completion_masks(forged, request, source)


def test_source_digest_and_exact_observed_mask_are_required():
    source, request = _request()
    with pytest.raises(AtomCompletionError, match="digest mismatch"):
        plan_request_atom_completion(request, source + b"\n")
    with pytest.raises(AtomCompletionError, match="every observed source heavy atom"):
        plan_request_atom_completion(replace(request, fixed_atoms=request.fixed_atoms[:-1]), source)
    upper = replace(request, normalized_pdb=replace(request.normalized_pdb, sha256=request.normalized_pdb.sha256.upper()))
    assert plan_request_atom_completion(upper, source) == plan_request_atom_completion(request, source)


def test_neither_substitutions_nor_altlocs_are_silently_resolved():
    source, request = _request()
    with pytest.raises(AtomCompletionError, match="substitution"):
        plan_request_atom_completion(replace(request, target_sequences=(TargetSequence("D", "RKG"), request.target_sequences[1])), source)
    altered = source[:16] + b"A" + source[17:]
    renamed = replace(request, normalized_pdb=replace(request.normalized_pdb, sha256=hashlib.sha256(altered).hexdigest()))
    with pytest.raises(AtomCompletionError, match="alternate locations"):
        plan_request_atom_completion(renamed, altered)


def test_plan_cannot_remove_observed_atoms_or_claim_a_whole_gap_residue():
    source, request = _request()
    plan = plan_request_atom_completion(request, source)
    with pytest.raises(AtomCompletionError, match="ownership must be disjoint"):
        replace(plan, gap_generated_atoms=(*plan.gap_generated_atoms, plan.missing_atoms[0]))
    with pytest.raises(AtomCompletionError, match="must not also be owned"):
        replace(plan, gap_generated_atoms=(*plan.gap_generated_atoms, AtomIdentity("D", "82", "", "OXT")))
    with pytest.raises(AtomCompletionError, match="must remain fixed"):
        replace(plan, fixed_atoms=plan.fixed_atoms[1:])


def test_partial_mask_cannot_include_duplicate_or_noncanonical_atoms():
    residue = ResidueIdentity("A", "1")
    canonical = canonical_heavy_atom_identities(residue, "ALA")
    with pytest.raises(AtomCompletionError, match="duplicate"):
        PartialResidueMask(residue, 0, "ALA", (canonical[0], canonical[0]), canonical[1:])
    with pytest.raises(AtomCompletionError, match="partition"):
        PartialResidueMask(residue, 0, "ALA", canonical[:-1], (AtomIdentity("d", "1", "", "O"),))


@pytest.mark.parametrize("atom_name", ["XX", "OXT"])
def test_unknown_atoms_and_nonterminal_oxt_are_not_silently_deleted(atom_name):
    source, request = _request()
    changed = source.replace(b"  CB ", f"{atom_name:>4s} ".encode(), 1)
    altered = replace(request, normalized_pdb=replace(request.normalized_pdb, sha256=hashlib.sha256(changed).hexdigest()),
                      fixed_atoms=(request.fixed_atoms[0], request.fixed_atoms[1],
                                   replace(request.fixed_atoms[2], atom_name=atom_name), *request.fixed_atoms[3:]))
    with pytest.raises(AtomCompletionError, match="unsupported heavy atoms or nonterminal OXT"):
        plan_request_atom_completion(altered, changed)


def test_terminal_oxt_is_retained_but_not_part_of_the_canonical_completion_mask():
    residue = ResidueIdentity("A", "1")
    source, fixed = _source([(residue, "LYS", ("N", "CA", "C", "O", "CB", "OXT"))])
    plan = plan_atom_completion(source, source_sha256=hashlib.sha256(source).hexdigest(),
        target_sequences=(TargetSequence("A", "K"),), sequence_placements=(SequencePlacement("A", 1, (0,), (residue,)),),
        fixed_atoms=fixed)
    assert AtomIdentity("A", "1", "", "OXT") in plan.fixed_atoms
    assert all(atom.atom_name != "OXT" for atom in plan.missing_atoms)


def test_complete_input_requires_no_completion_and_empty_target_is_not_an_admission():
    residue = ResidueIdentity("A", "1")
    source, fixed = _source([(residue, "ALA", CANONICAL_HEAVY_ATOMS["ALA"])])
    plan = plan_atom_completion(source, source_sha256=hashlib.sha256(source).hexdigest(),
        target_sequences=(TargetSequence("A", "A"),), sequence_placements=(SequencePlacement("A", 1, (0,), (residue,)),),
        fixed_atoms=fixed)
    assert not plan.partial_residues
    assert not plan.missing_atoms
    assert plan.fixed_atoms == fixed
    with pytest.raises(AtomCompletionError, match="requires observed protein targets"):
        plan_atom_completion(source, source_sha256=hashlib.sha256(source).hexdigest(),
                             target_sequences=(), sequence_placements=(), fixed_atoms=fixed)


def test_heterogen_observations_are_preserved_without_claiming_backend_support():
    source, request = _request()
    extra, extra_atoms = _source([(ResidueIdentity("Z", "99"), "LIG", ("C1",))])
    heterogen = extra.decode().replace("ATOM  ", "HETATM").replace("    1 ", " 9999 ", 1).encode()
    source = source.replace(b"END\n", heterogen)
    altered = replace(request, normalized_pdb=replace(request.normalized_pdb, sha256=hashlib.sha256(source).hexdigest()),
                      fixed_atoms=(*request.fixed_atoms, *extra_atoms))
    plan = plan_request_atom_completion(altered, source)
    assert extra_atoms[0] in plan.fixed_atoms
    assert len(plan.partial_residues) == 3


def test_malformed_source_and_duplicate_target_ownership_fail_closed():
    source, request = _request()
    malformed = source.replace(b"  1.100", b"    nan", 1)
    with pytest.raises(AtomCompletionError, match="invalid atom-completion source"):
        plan_request_atom_completion(replace(request, normalized_pdb=replace(request.normalized_pdb,
            sha256=hashlib.sha256(malformed).hexdigest())), malformed)
    plan = plan_request_atom_completion(request, source)
    duplicate = replace(plan.partial_residues[1], target_index=0)
    with pytest.raises(AtomCompletionError, match="target positions must be unique"):
        replace(plan, partial_residues=(plan.partial_residues[0], duplicate, plan.partial_residues[2]))


def test_original_5j13_inventory_matches_pdbfixer_without_rebuilding_atoms(tmp_path, fixtures_root):
    from pdbfixer import PDBFixer

    from dvbfixer.model.diffusion_cli import build_cli_diffusion_request
    from dvbfixer.model.pipeline import parse_fasta

    original = fixtures_root / "protein_context/5j13_selected.pdb"
    before = original.read_bytes()
    # Inventory-only protein projection: no prepared/template-built atoms and
    # no assertion that the original SO4 input is admitted for sampling.
    projection = tmp_path / "protein-observations.pdb"
    projection.write_bytes(b"".join(line for line in before.splitlines(keepends=True)
                                   if line.startswith(b"ATOM  ")) + b"END\n")
    targets = parse_fasta(str(original.with_suffix(".fasta")))
    request = build_cli_diffusion_request(projection, targets, seeds=(7,), profile="protpardelle-1c-cc94-mps-v1")
    plan = plan_request_atom_completion(request, projection.read_bytes())
    fixer = PDBFixer(filename=str(original))
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    independent = {AtomIdentity(residue.chain.id, residue.id, residue.insertionCode.strip(), atom.name)
                   for residue, atoms in fixer.missingAtoms.items() for atom in atoms}
    assert set(plan.missing_atoms) == independent
    assert len(plan.partial_residues) == 20
    assert len(plan.missing_atoms) == 62
    assert len(request.gaps) == 7
    assert sum(len(gap.generated_residues) for gap in request.gaps) == 42
    assert original.read_bytes() == before
