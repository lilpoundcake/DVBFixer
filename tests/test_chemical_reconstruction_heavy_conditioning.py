"""Heavy-conditioning transport is not source provenance or geometry acceptance."""
from __future__ import annotations

from dataclasses import replace

import pytest
from scripts.chemical_reconstruction_all20_controls import (
    HEAVY_POLICY,
    PAIRS,
    compile_H,
    sidechain_heavy_mask,
    synthetic_input,
)

from dvbfixer.chemical_reconstruction.admission import admit
from dvbfixer.chemical_reconstruction.conformer_search import search_conformers
from dvbfixer.chemical_reconstruction.contracts import (
    ChemistryKind,
    Component,
    ContractError,
    Observation,
    OpaqueContextPolicy,
    Operation,
    Origin,
    Ownership,
    PreparedContextRequest,
    Unit,
    digest_record,
)
from dvbfixer.chemical_reconstruction.heavy_conditioning import (
    HeavyConditionedRequest,
    condition_on_heavy,
)
from dvbfixer.chemical_reconstruction.observed_peptide import compile_peptide_request_plan
from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef


def parent_case(source='SER', *, fractional=False):
    full, authority, incoming = synthetic_input(source, dict(PAIRS)[source], {})
    generated = frozenset(o.identity for o in full.observations if o.identity.component == authority.identity
                          and o.identity.name in sidechain_heavy_mask(authority))
    original = replace(full, observations=tuple(o for o in full.observations if o.identity not in generated),
                       ownership=Ownership(generated))
    if fractional:
        atom = Observation(ExactAtomRef(ComponentInstanceRef(1, 'X', '999'), 'Q'), 'C', (50., 50., 50.), Unit.ANGSTROM, .5)
        observations = (*original.observations, atom)
        original = PreparedContextRequest(
            original.source_sha256, original.operation,
            (*original.components, Component(atom.identity.component, ChemistryKind.UNKNOWN)), observations,
            original.targets, original.ownership, original.connections, original.constraints, original.ph,
            Origin('prepared-observed-peptide-v1', 'a'*64),
            context_policy=OpaqueContextPolicy(original.source_sha256, 'a'*64, frozenset({atom.identity}), digest_record(observations)),
        )
    plan = compile_peptide_request_plan(original, HEAVY_POLICY)
    roots = tuple(sorted(o.identity for o in original.observations if o.identity.component in original.targets))
    result = search_conformers(original, plan, authority, roots)
    assert result.selected_index is not None
    checkpoint = result.attempts[result.selected_index].rounded
    return original, authority, incoming, checkpoint


def test_numeric_conditioning_is_distinct_from_source_and_preserves_both_roles():
    original, _, _, checkpoint = parent_case()
    request = condition_on_heavy(original, checkpoint)
    assert isinstance(request, HeavyConditionedRequest) and not isinstance(request, PreparedContextRequest)
    assert request.original is original and request.checkpoint is checkpoint
    assert request.observations == original.observations + checkpoint.additions
    assert request.source_sha256 != original.source_sha256 and request.digest != original.digest
    assert not request.ownership.generated and request.original.ownership.generated
    assert any(f.code == 'completion-requires-generated-mask' for f in admit(request).findings)


@pytest.mark.parametrize('source', ['SER', 'THR', 'PRO', 'HIS'])
def test_existing_complete_H_compilation_preserves_conditioning_through_replace(source):
    original, authority, incoming, checkpoint = parent_case(source)
    request = condition_on_heavy(original, checkpoint)
    view = compile_H(request, authority, incoming, {})
    assert isinstance(view.request, HeavyConditionedRequest)
    assert view.request.original is original and view.request.checkpoint is checkpoint
    assert view.request.operation == Operation.HYDROGENS
    assert view.request.observations is request.observations
    assert not original.ownership.generated & view.request.ownership.generated
    assert {a.name for a in view.request.ownership.generated} == {a.name for a in authority.state.graph.atoms if a.element == 'H'}
    assert not any(f.code == 'generated-element-conflicts-with-operation' for f in admit(view.request).findings)


@pytest.mark.parametrize('mutation', ['partial', 'extra', 'nonrounded', 'element', 'unit', 'occupancy', 'binding', 'raw'])
def test_checkpoint_mutations_fail_closed(mutation):
    original, _, _, checkpoint = parent_case()
    first = checkpoint.additions[0]
    changes = {
        'nonrounded': {'position': (first.position[0]+.00001, *first.position[1:])},
        'element': {'element': 'N' if first.element == 'C' else 'C'},
        'unit': {'unit': Unit.NANOMETER}, 'occupancy': {'occupancy': .5},
    }
    if mutation == 'partial':
        checkpoint = replace(checkpoint, additions=checkpoint.additions[:-1])
    elif mutation == 'extra':
        checkpoint = replace(checkpoint, additions=(*checkpoint.additions, original.observations[0]))
    elif mutation == 'binding':
        checkpoint = replace(checkpoint, request_digest='0'*64)
    elif mutation == 'raw':
        checkpoint = replace(checkpoint, additions=tuple(replace(o, position=tuple(x+.00001 for x in o.position)) for o in checkpoint.additions))
    else:
        checkpoint = replace(checkpoint, additions=(replace(first, **changes[mutation]), *checkpoint.additions[1:]))
    with pytest.raises(ContractError):
        condition_on_heavy(original, checkpoint)


@pytest.mark.parametrize('mutation', ['observed', 'origin', 'source', 'ph', 'constraints', 'nested', 'partial-H', 'heavy-owned', 'mobile'])
def test_request_role_or_source_mutation_is_rejected(mutation):
    original, authority, incoming, checkpoint = parent_case()
    request = condition_on_heavy(original, checkpoint)
    if mutation == 'observed':
        changes = {'observations': tuple(reversed(request.observations))}
    elif mutation == 'origin':
        changes = {'origin': original.origin}
    elif mutation == 'source':
        changes = {'source_sha256': original.source_sha256}
    elif mutation == 'ph':
        changes = {'ph': 5.}
    elif mutation == 'constraints':
        assert original.constraints
        changes = {'constraints': ()}
    elif mutation == 'nested':
        changes = {'original': request}
    elif mutation == 'partial-H':
        view = compile_H(request, authority, incoming, {})
        changes = {'operation': Operation.HYDROGENS, 'ownership': Ownership(frozenset(list(view.request.ownership.generated)[:-1]))}
    elif mutation == 'heavy-owned':
        changes = {'ownership': original.ownership}
    else:
        changes = {'ownership': Ownership(frozenset(), mobile=frozenset({request.observations[0].identity}))}
    with pytest.raises(ContractError):
        replace(request, **changes)


def test_fractional_scope_is_inherited_only_from_original_prepared_proof():
    original, authority, incoming, checkpoint = parent_case(fractional=True)
    request = condition_on_heavy(original, checkpoint)
    view = compile_H(request, authority, incoming, {})
    assert request.allowed_fractional_atoms == original.context_policy.fractional_atoms
    assert not any(f.code == 'partial-occupancy-or-altloc-requires-explicit-selection' for f in admit(view.request).findings)
    atom = next(iter(request.allowed_fractional_atoms))
    with pytest.raises(ContractError):
        replace(view.request, ownership=Ownership(view.request.ownership.generated | {atom}))
    with pytest.raises(ContractError):
        replace(request, observations=tuple(replace(o, occupancy=1.) if o.identity == atom else o for o in request.observations))


def test_conditioning_does_not_acquire_prepared_source_transport_provenance():
    from types import SimpleNamespace

    from dvbfixer.chemical_reconstruction.source_preparation import audit_source_transport

    original, _, _, checkpoint = parent_case()
    request = condition_on_heavy(original, checkpoint)
    # Boundary refuses on DTO type before touching source bytes or probe details.
    report = audit_source_transport(None, SimpleNamespace(request=request))
    assert not report.passed and report.findings == ('prepared-context-request-required',)


def test_even_numerically_equal_source_coordinate_serialization_is_bound():
    original, _, _, checkpoint = parent_case()
    first = original.observations[0]
    # 0 and 0.0 compare equal in Python, but are distinct source records/digests.
    original = replace(original, observations=(replace(first, position=(0., *first.position[1:])), *original.observations[1:]))
    checkpoint = replace(checkpoint, request_digest=original.digest)
    request = condition_on_heavy(original, checkpoint)
    altered = replace(request.observations[0], position=(0, *first.position[1:]))
    assert altered == request.observations[0]
    with pytest.raises(ContractError, match='observation-transport'):
        replace(request, observations=(altered, *request.observations[1:]))


def test_generated_heavy_or_partial_H_cannot_be_smuggled_into_generic_operation():
    original, authority, incoming, checkpoint = parent_case()
    view = compile_H(condition_on_heavy(original, checkpoint), authority, incoming, {})
    with pytest.raises(ContractError, match='complete-H-mask'):
        replace(view.request, ownership=Ownership(original.ownership.generated | view.request.ownership.generated))
    with pytest.raises(ContractError, match='complete-H-mask'):
        replace(view.request, operation=Operation.REPACK)


@pytest.mark.parametrize('settings', [(("foreign", 'a'*64),),
                                    (("peptide_state_declaration", 'a'*64),),
                                    (("complete_H_policy", 'unversioned'), ("complete_H_authority", 'a'*64))])
def test_only_existing_declarative_compiler_metadata_is_transportable(settings):
    original, _, _, checkpoint = parent_case()
    request = condition_on_heavy(original, checkpoint)
    with pytest.raises(ContractError):
        replace(request, origin=replace(request.origin, settings=settings))
