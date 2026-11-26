"""Synthetic measurement/transport controls, not numerical R1 cohort execution."""
from __future__ import annotations

import hashlib
import json
import math
from dataclasses import replace

import pytest
from scripts.chemical_reconstruction_all20_controls import load_parent, synthetic_input
from scripts.chemical_reconstruction_ccd_controls import _json

from dvbfixer.chemical_reconstruction import cycle_audit
from dvbfixer.chemical_reconstruction.contracts import (
    ChemistryKind,
    Component,
    Constraint,
    ConstraintKind,
    Observation,
    OpaqueContextPolicy,
    Operation,
    Origin,
    Ownership,
    PreparedContextRequest,
    ReconstructionRequest,
    Unit,
    digest_record,
)
from dvbfixer.chemical_reconstruction.cycle_audit import (
    CycleGeometryAudit,
    CyclePairAudit,
    audit_cycle_geometry,
    audit_cycle_pair,
    measure_cycle_geometry,
    measure_cycle_pair,
    round_cycle_additions,
)
from dvbfixer.chemical_reconstruction.cycle_contract import compile_cycle_problem
from dvbfixer.chemical_reconstruction.cycle_initialization import initialize_cycle_problem
from dvbfixer.chemical_reconstruction.cycle_source import (
    CycleSourceMetadata,
    audit_cycle_source,
    bind_cycle_source,
)
from dvbfixer.chemical_reconstruction.metadata_chirality import assert_metadata_all_l
from dvbfixer.chemical_reconstruction.observed_peptide import _parse_source
from dvbfixer.chemical_reconstruction.research_decode import decode_research_record
from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.pdb_component_reconstruction import _component_ref


def source_fixture(parent='PRO', mask=('CD',), extra=()):
    """Synthetic full observations -> actual PDB bytes -> exact masked request.

    Witness additions come from synthetic input, not a solver or held-out native
    scientific trial. Adapt the synthetic occurrence to an actual single-residue
    PDB identity; separate helper tests exercise repeated occurrences/models.
    """
    full, _, _ = synthetic_input(parent, 'HIP' if parent == 'HIS' else parent, {})
    old = next(iter(full.targets))
    target = replace(old, occurrence=1)
    def ref(a):
        return replace(a, component=target) if a.component == old else a
    observations = tuple(replace(o, identity=ref(o.identity)) for o in full.observations) + tuple(extra)
    components = tuple(replace(c, identity=target) if c.identity == old else c for c in full.components)
    components += tuple(Component(c, ChemistryKind.UNKNOWN) for c in dict.fromkeys(o.identity.component for o in extra))
    lines = []
    for serial, o in enumerate(observations, 1):
        c = o.identity.component
        name = parent if c == target else 'ALA' if c.chain_id == 'A' else 'UNK'
        x, y, z = o.angstrom
        lines.append(f'ATOM  {serial:5d} {o.identity.name:>4s} {name:>3s} {c.chain_id}{int(c.sequence_number):4d}{c.insertion_code:1s}   '
                     f'{x:8.3f}{y:8.3f}{z:8.3f}{o.occupancy:6.2f}{0.:6.2f}          {o.element:>2s}\n')
    prepared = (''.join(lines)+'END\n').encode('ascii')
    parsed_lines, atoms = _parse_source(prepared)
    measured = tuple(Observation(ExactAtomRef(_component_ref(a), a.name), a.element, a.position, Unit.ANGSTROM,
        float(parsed_lines[a.line_index][54:60])) for a in atoms)
    generated = frozenset(ExactAtomRef(target, n) for n in mask)
    request = replace(full, source_sha256=hashlib.sha256(prepared).hexdigest(), components=components,
        targets=frozenset({target}), observations=tuple(o for o in measured if o.identity not in generated),
        ownership=Ownership(generated), constraints=tuple(replace(q, atoms=tuple(ref(a) for a in q.atoms)) for q in full.constraints),
        connections=tuple(replace(c, link=replace(c.link, component_atom=ref(c.link.component_atom), other_atom=ref(c.link.other_atom)))
                          for c in full.connections), origin=Origin('synthetic-cycle-source-test-v1', hashlib.sha256(prepared).hexdigest()))
    fractional = frozenset(o.identity for o in request.observations if o.occupancy != 1)
    if fractional:
        request = PreparedContextRequest(request.source_sha256, request.operation, request.components, request.observations,
            request.targets, request.ownership, request.connections, request.constraints, request.ph,
            Origin('prepared-observed-peptide-v1', 'b'*64),
            context_policy=OpaqueContextPolicy(request.source_sha256, 'b'*64, fractional, digest_record(request.observations)))
    problem = compile_cycle_problem(request, tuple(ExactAtomRef(target, n) for n in ('CA', 'N', 'C')))
    metadata = bind_cycle_source(problem, prepared)
    witness = tuple(o for o in measured if o.identity in generated)
    return problem, prepared, metadata, witness


@pytest.mark.parametrize('parent,mask', [('PRO', ('CD',)), ('PRO', ('CB', 'CD')), ('PRO', ('CB', 'CG', 'CD')),
    ('PHE', ('CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ')),
    ('TYR', ('CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH')),
    ('TRP', ('CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'))])
def test_synthetic_witness_measurement_covers_all_terms_and_not_constructor_approval(parent, mask):
    problem, prepared, metadata, witness = source_fixture(parent, mask)
    report = measure_cycle_geometry(problem, prepared, metadata, witness, phase='external-witness')
    assert report.passed and report.global_l_status == 'passed'
    assert len(report.terms) == len(problem.coverage.metrics)+len(problem.source.constraints)+len(problem.extra_constraints)
    assert {t.role for t in report.terms} == {'fixed-only', 'generated-involving'}
    assert not report.constructor_approved and not report.md_ready and not report.source_coverage_verified
    audit_cycle_geometry(problem, prepared, metadata, witness, report)


def test_actual_names_and_masked_bytes_are_bound_without_retaining_withheld_records():
    problem, prepared, metadata, witness = source_fixture()
    audit_cycle_source(problem, prepared, metadata)
    _, atoms = _parse_source(metadata.masked_pdb.encode('ascii'))
    assert not {ExactAtomRef(_component_ref(a), a.name) for a in atoms} & problem.source.ownership.generated
    assert dict(metadata.residue_names)[next(iter(problem.source.targets))] == 'PRO'
    assert metadata.observations_digest == digest_record(problem.source.observations)
    assert not metadata.source_coverage_verified
    assert len(witness) == 1 and 'CD' not in {a.name for a in atoms if a.chain_id == 'a'}


def test_self_consistent_relabelled_source_cannot_authenticate_against_original_bytes():
    problem, prepared, metadata, witness = source_fixture()
    changed = metadata.masked_pdb.replace(' PRO a', ' GLY a')
    names = tuple((c, 'GLY' if c in problem.source.targets else n) for c, n in metadata.residue_names)
    forged = replace(metadata, masked_pdb=changed, masked_sha256=hashlib.sha256(changed.encode()).hexdigest(), residue_names=names)
    with pytest.raises(ValueError, match='readback-binding'):
        measure_cycle_geometry(problem, prepared, forged, witness, phase='raw')


def test_actual_wrong_protein_parent_cannot_exploit_the_global_L_skip_list():
    problem, prepared, _, _ = source_fixture()
    changed = prepared.replace(b' PRO a', b' GLY a')
    request = replace(problem.source, source_sha256=hashlib.sha256(changed).hexdigest())
    other = compile_cycle_problem(request, problem.anchors)
    with pytest.raises(ValueError, match='target-parent-name-mismatch'):
        bind_cycle_source(other, changed)


def test_variant_parent_transport_preserves_actual_label_and_explicit_state():
    problem, prepared, _, witness = source_fixture('HIS', ('CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'))
    changed = prepared.replace(b' HIS a', b' HIP a')
    request = replace(problem.source, source_sha256=hashlib.sha256(changed).hexdigest())
    other = compile_cycle_problem(request, problem.anchors)
    metadata = bind_cycle_source(other, changed)
    target = next(iter(request.targets))
    assert dict(metadata.residue_names)[target] == 'HIP'
    assert other.source.components is problem.source.components
    assert measure_cycle_geometry(other, changed, metadata, witness, phase='external-witness').passed


def test_unknown_metadata_element_is_unmeasurable_not_a_global_L_pass():
    extra = Observation(ExactAtomRef(ComponentInstanceRef(1, 'z', '999'), 'Q'), 'X', (50., 50., 50.), Unit.ANGSTROM)
    problem, prepared, metadata, witness = source_fixture(extra=(extra,))
    report = measure_cycle_geometry(problem, prepared, metadata, witness, phase='raw')
    assert not report.passed and report.global_l_status == 'unmeasurable'
    assert 'metadata-topology-unknown-element' in report.global_l_detail


def test_withheld_coordinates_are_not_used_by_the_geometry_measurement():
    problem, prepared, metadata, witness = source_fixture()
    original = measure_cycle_geometry(problem, prepared, metadata, witness, phase='external-witness')
    lines = prepared.decode().splitlines(keepends=True)
    for index, line in enumerate(lines):
        if line.startswith('ATOM') and line[12:16].strip() == 'CD' and line[21] == 'a':
            lines[index] = line[:30]+f'{70.:8.3f}{70.:8.3f}{70.:8.3f}'+line[54:]
    changed = ''.join(lines).encode()
    request = replace(problem.source, source_sha256=hashlib.sha256(changed).hexdigest())
    other = compile_cycle_problem(request, problem.anchors)
    newer = bind_cycle_source(other, changed)
    report = measure_cycle_geometry(other, changed, newer, witness, phase='external-witness')
    assert newer.masked_pdb == metadata.masked_pdb
    assert newer.removed_records_digest != metadata.removed_records_digest
    assert report.terms == original.terms and report.findings == original.findings
    assert report.severe_overlap_pairs == original.severe_overlap_pairs and report.global_l_status == original.global_l_status


def test_isolated_resolved_free_molecule_is_not_restricted_to_protein_inputs():
    authority = load_parent('PRO')
    identity = ComponentInstanceRef(1, 'x', '12', 'A')
    full = tuple(Observation(ExactAtomRef(identity, a.name), a.element,
        tuple(float(f'{v:.3f}') for v in a.ideal_position), Unit.ANGSTROM)
        for a in authority.state.graph.atoms if a.element not in {'H', 'D'})
    lines = [f'ATOM  {i:5d} {o.identity.name:>4s} PRO x  12A   '
             f'{o.position[0]:8.3f}{o.position[1]:8.3f}{o.position[2]:8.3f}{1.:6.2f}{0.:6.2f}          {o.element:>2s}\n'
             for i, o in enumerate(full, 1)]
    prepared = (''.join(lines)+'END\n').encode()
    generated = frozenset({ExactAtomRef(identity, 'CD')})
    request = ReconstructionRequest(hashlib.sha256(prepared).hexdigest(), Operation.HEAVY,
        (Component(identity, ChemistryKind.ISOLATED, (authority.state,)),),
        tuple(o for o in full if o.identity not in generated), frozenset({identity}), Ownership(generated),
        constraints=authority.constraints(identity, include_hydrogens=False, minimum_volume=.05, planar_tolerance_angstrom=.05))
    problem = compile_cycle_problem(request, tuple(ExactAtomRef(identity, n) for n in ('CA', 'N', 'C')))
    metadata = bind_cycle_source(problem, prepared)
    witness = tuple(o for o in full if o.identity in generated)
    report = measure_cycle_geometry(problem, prepared, metadata, witness, phase='external-witness')
    assert report.passed and report.global_l_status == 'passed'


def test_degenerate_generated_geometry_is_retained_as_findings_not_zero_penalty():
    problem, prepared, metadata, witness = source_fixture()
    cg = next(o.position for o in problem.source.observations if o.identity.name == 'CG')
    collapsed = (replace(witness[0], position=cg),)
    report = measure_cycle_geometry(problem, prepared, metadata, collapsed, phase='raw')
    assert not report.passed and 'cycle-degenerate-metric' in report.findings
    assert any(t.error for t in report.terms)


@pytest.mark.parametrize('mutation', ['source', 'mask', 'name', 'metadata', 'observations'])
def test_source_and_metadata_tampering_fail_before_physical_measurement(mutation, monkeypatch):
    problem, prepared, metadata, witness = source_fixture()
    monkeypatch.setattr(cycle_audit, 'assert_metadata_all_l', lambda *a: pytest.fail('bad source reached global-L'))
    with pytest.raises(ValueError):
        if mutation == 'source':
            measure_cycle_geometry(problem, prepared+b'REMARK changed\n', metadata, witness, phase='raw')
        elif mutation == 'mask':
            forged = replace(metadata, masked_pdb=prepared.decode(), masked_sha256=hashlib.sha256(prepared).hexdigest())
            measure_cycle_geometry(problem, prepared, forged, witness, phase='raw')
        elif mutation == 'name':
            replace(metadata, residue_names=tuple((c, 'GLY') for c, _ in metadata.residue_names))
        elif mutation == 'metadata':
            measure_cycle_geometry(problem, prepared, replace(metadata, removed_records_digest='0'*64), witness, phase='raw')
        else:
            measure_cycle_geometry(problem, prepared, replace(metadata, observations_digest='0'*64), witness, phase='raw')


@pytest.mark.parametrize('mutation', ['empty', 'duplicate', 'extra', 'wrong-name', 'element', 'unit', 'occupancy'])
def test_bad_additions_are_not_spliced_or_globally_approved(mutation):
    problem, prepared, metadata, witness = source_fixture()
    if mutation == 'empty':
        values = ()
    elif mutation == 'duplicate':
        values = witness+witness
    elif mutation == 'extra':
        values = witness+(problem.source.observations[0],)
    else:
        fields = {'identity': replace(witness[0].identity, name='FAKE')} if mutation == 'wrong-name' else {
            'element': 'O'} if mutation == 'element' else {'unit': Unit.NANOMETER} if mutation == 'unit' else {'occupancy': .5}
        values = (replace(witness[0], **fields),)
    report = measure_cycle_geometry(problem, prepared, metadata, values, phase='raw')
    assert not report.passed and not report.terms and report.global_l_status == 'not-measured'


def test_full_opaque_context_overlap_is_not_cropped_or_excluded_by_a_guessed_bond():
    _, _, _, witness = source_fixture()
    extra = Observation(ExactAtomRef(ComponentInstanceRef(1, 'z', '999'), 'Q'), 'C', witness[0].position, Unit.ANGSTROM)
    problem, prepared, metadata, witness = source_fixture(extra=(extra,))
    report = measure_cycle_geometry(problem, prepared, metadata, witness, phase='raw')
    assert not report.passed and 'cycle-generated-context-minimum-separation' in report.findings
    assert any(extra.identity in (a, b) and d == 0 for a, b, d in report.severe_overlap_pairs)
    assert report.minimum_eligible_separation_angstrom == 0


@pytest.mark.parametrize('clashing', [False, True])
def test_fractional_context_is_preserved_and_still_in_separation_measurements(clashing):
    _, _, _, witness = source_fixture()
    point = witness[0].position if clashing else (50., 50., 50.)
    extra = Observation(ExactAtomRef(ComponentInstanceRef(1, 'z', '999'), 'Q'), 'C', point, Unit.ANGSTROM, .5)
    problem, prepared, metadata, witness = source_fixture(extra=(extra,))
    policy = problem.source.context_policy
    report = measure_cycle_geometry(problem, prepared, metadata, witness, phase='raw')
    assert report.passed != clashing and problem.source.context_policy is policy
    assert next(o for o in problem.source.observations if o.identity == extra.identity).occupancy == .5


@pytest.mark.parametrize('kind', ['finite-exclusion-upper', 'nanometer-distance'])
def test_declared_constraint_semantics_and_units_are_not_dropped(kind):
    problem, prepared, _, witness = source_fixture()
    target = next(iter(problem.source.targets))
    atoms = ExactAtomRef(target, 'N'), ExactAtomRef(target, 'CD')
    extra = (Constraint(ConstraintKind.EXCLUSION, atoms, Unit.ANGSTROM, .1, .2) if kind == 'finite-exclusion-upper'
             else Constraint(ConstraintKind.DISTANCE, atoms, Unit.NANOMETER, .14, .16))
    other = compile_cycle_problem(problem.source, problem.anchors, extra_constraints=(extra,))
    metadata = bind_cycle_source(other, prepared)
    report = measure_cycle_geometry(other, prepared, metadata, witness, phase='raw')
    assert report.terms[-1].value == pytest.approx(math.dist(*[next(o.position for o in (*problem.source.observations, *witness)
        if o.identity == a) for a in atoms])/(10 if kind == 'nanometer-distance' else 1))
    assert report.terms[-1].passed == (kind == 'nanometer-distance')
    if kind == 'finite-exclusion-upper':
        assert not report.passed and 'cycle-declared-constraint' in report.findings


def test_known_one_two_and_one_three_pairs_are_excluded_without_losing_other_context():
    problem, prepared, metadata, witness = source_fixture()
    report = measure_cycle_geometry(problem, prepared, metadata, witness, phase='raw')
    target = next(iter(problem.source.targets))
    excluded = {'CD', 'N', 'CG', 'CA', 'CB'}
    expected = [o for o in problem.source.observations if o.identity.component != target or o.identity.name not in excluded]
    # Incoming C is a declared N partner and therefore 1-3 from CD.
    expected = [o for o in expected if o.identity.component.chain_id != 'A']
    assert report.eligible_pair_count == len(expected)


def test_preexisting_D_in_retained_context_blocks_a_complete_geometry_pass():
    problem, _, _, witness = source_fixture()
    target = next(iter(problem.source.targets))
    atoms = {o.identity.name: o for o in (*problem.source.observations, *witness) if o.identity.component == target}
    # Reflect all coordinates in a plane: preserve context bond lengths, reverse C-alpha parity.
    context = ComponentInstanceRef(1, 'z', '999')
    extra = tuple(Observation(ExactAtomRef(context, n), atoms[n].element,
        (50.-atoms[n].position[0], 50.+atoms[n].position[1], 50.+atoms[n].position[2]), Unit.ANGSTROM)
        for n in ('N', 'CA', 'C', 'CB'))
    problem, prepared, metadata, witness = source_fixture(extra=extra)
    report = measure_cycle_geometry(problem, prepared, metadata, witness, phase='raw')
    assert report.global_l_status == 'rejected' and not report.passed
    assert 'cycle-global-L-rejected' in report.findings


def test_generated_D_is_detected_independently_from_declared_stereo_measurements():
    problem, prepared, metadata, witness = source_fixture(mask=('CB', 'CG', 'CD'))
    positions = {o.identity.name: o.position for o in problem.source.observations if o.identity.component in problem.source.targets}
    from dvbfixer.chemical_reconstruction.measurement import _cross, _dot, _sub, _unit
    normal = _unit(_cross(_sub(positions['N'], positions['CA']), _sub(positions['C'], positions['CA'])))
    reflected = tuple(replace(o, position=tuple(o.position[k]-2*_dot(_sub(o.position, positions['CA']), normal)*normal[k]
                                                for k in range(3))) for o in witness)
    report = measure_cycle_geometry(problem, prepared, metadata, reflected, phase='raw')
    assert report.global_l_status == 'rejected' and not report.passed


def test_raw_rounded_pair_is_measurement_only_and_does_not_require_a_passing_start():
    problem, prepared, metadata, witness = source_fixture(mask=('CB', 'CG', 'CD'))
    initializer = initialize_cycle_problem(problem)
    initial = measure_cycle_geometry(problem, prepared, metadata, initializer.starts[0].additions, phase='initial')
    assert not initial.passed
    rounded = round_cycle_additions(witness)
    report = measure_cycle_pair(problem, prepared, metadata, initializer, 0, witness, rounded)
    assert report.passed and not report.constructor_approved and not report.md_ready
    assert not initializer.solver_ready
    audit_cycle_pair(problem, prepared, metadata, initializer, witness, rounded, report)


def test_rounding_transport_mismatch_is_not_accepted_as_a_second_raw_pose():
    problem, prepared, metadata, witness = source_fixture()
    initializer = initialize_cycle_problem(problem)
    changed = (replace(witness[0], position=(witness[0].position[0]+.001, *witness[0].position[1:])),)
    report = measure_cycle_pair(problem, prepared, metadata, initializer, 0, witness, changed)
    assert not report.passed and not report.rounding_matches and 'cycle-rounding-transport' in report.findings


def test_rounded_cap_crossing_is_rejected_without_clipping():
    problem, prepared, metadata, _ = source_fixture()
    initializer = initialize_cycle_problem(problem)
    start = initializer.starts[0].additions[0]
    for sign in (-1, 1):
        raw = (replace(start, position=(start.position[0]+sign*(8.-1e-7), *start.position[1:])),)
        rounded = round_cycle_additions(raw)
        if math.dist(rounded[0].position, start.position) > 8:
            break
    else:
        pytest.fail('adversarial synthetic rounding boundary not constructed')
    report = measure_cycle_pair(problem, prepared, metadata, initializer, 0, raw, rounded)
    assert report.raw_maximum_movement_angstrom < 8 < report.rounded_maximum_movement_angstrom
    assert 'cycle-per-atom-movement-cap' in report.findings and not report.passed


def test_measurement_and_readback_call_no_constructor_optimizer_or_linear_solve(monkeypatch):
    import numpy as np
    from scipy import optimize

    from dvbfixer.chemical_reconstruction import construction, cycle_initialization

    problem, prepared, metadata, witness = source_fixture()
    initializer = initialize_cycle_problem(problem)
    def forbidden(*args, **kwargs):
        pytest.fail('independent measurement constructed coordinates')
    for module, name in ((construction, 'propose'), (cycle_initialization, 'initialize_cycle_problem'),
                         (cycle_initialization, '_direct_solve'), (cycle_initialization, '_make_start'),
                         (np.linalg, 'solve'), (optimize, 'least_squares')):
        monkeypatch.setattr(module, name, forbidden)
    report = measure_cycle_pair(problem, prepared, metadata, initializer, 0, witness, round_cycle_additions(witness))
    audit_cycle_pair(problem, prepared, metadata, initializer, witness, round_cycle_additions(witness), report)


@pytest.mark.parametrize('key', ['geometry', 'pairs', 'global', 'source', 'movement', 'rounding'])
def test_saved_measurements_cannot_be_relabelled_or_numerically_tampered(key):
    problem, prepared, metadata, witness = source_fixture()
    initializer = initialize_cycle_problem(problem)
    report = measure_cycle_pair(problem, prepared, metadata, initializer, 0, witness, round_cycle_additions(witness))
    with pytest.raises(ValueError):
        if key == 'geometry':
            altered = replace(report, raw=replace(report.raw, terms=report.raw.terms[:-1]))
        elif key == 'pairs':
            altered = replace(report, raw=replace(report.raw, eligible_pair_count=report.raw.eligible_pair_count+1))
        elif key == 'global':
            altered = replace(report, raw=replace(report.raw, global_l_status='rejected'))
        elif key == 'source':
            altered = replace(report, raw=replace(report.raw, source_metadata_digest='0'*64))
        elif key == 'movement':
            altered = replace(report, raw_maximum_movement_angstrom=0.)
        else:
            altered = replace(report, rounding_matches=False)
        audit_cycle_pair(problem, prepared, metadata, initializer, witness, round_cycle_additions(witness), altered)


def test_strict_transport_has_standalone_roles_and_external_source_audit_is_required():
    problem, prepared, metadata, witness = source_fixture()
    report = measure_cycle_geometry(problem, prepared, metadata, witness, phase='external-witness')
    restored_metadata = decode_research_record(json.dumps(_json(metadata)).encode(), CycleSourceMetadata)
    restored = decode_research_record(json.dumps(_json(report)).encode(), CycleGeometryAudit)
    audit_cycle_geometry(problem, prepared, restored_metadata, witness, restored)
    assert digest_record(restored) == digest_record(report)
    initializer = initialize_cycle_problem(problem)
    pair = measure_cycle_pair(problem, prepared, metadata, initializer, 0, witness, round_cycle_additions(witness))
    restored_pair = decode_research_record(json.dumps(_json(pair)).encode(), CyclePairAudit)
    audit_cycle_pair(problem, prepared, metadata, initializer, witness, round_cycle_additions(witness), restored_pair)


def test_metadata_topology_keeps_case_models_insertions_occurrences_and_late_additions(monkeypatch):
    import dvbfixer.ffutils.geometry as geometry

    problem, _, metadata, witness = source_fixture(mask=('CB', 'CG', 'CD'))
    target = next(iter(problem.source.targets))
    atoms = tuple(o for o in (*problem.source.observations, *witness) if o.identity.component == target)
    components = [ComponentInstanceRef(1, 'D', '82'), ComponentInstanceRef(1, 'd', '82', 'A'),
                  ComponentInstanceRef(2, 'D', '82'), ComponentInstanceRef(1, 'D', '82', occurrence=2)]
    values = tuple(replace(o, identity=replace(o.identity, component=c)) for c in components for o in atoms)
    values = (*values[1:], values[0])  # A previously seen residue receives an atom last.
    captured = []
    def inspect(topology, positions):
        captured.extend((r.chain.id, r.id, r.insertionCode, len(list(r.atoms()))) for r in topology.residues())
        assert len(positions) == len(values)
    monkeypatch.setattr(geometry, 'assert_all_l', inspect)
    assert_metadata_all_l(values, dict.fromkeys(components, 'PRO'))
    assert len(captured) == 4 and {row[0] for row in captured} == {'D', 'd'}
    assert all(row[3] == len(atoms) for row in captured)


@pytest.mark.parametrize('kind', ['unavailable', 'programming', 'unrelated-module'])
def test_environmental_failure_and_unexpected_errors_are_not_fake_passes(kind, monkeypatch):
    problem, prepared, metadata, witness = source_fixture()
    def failure(*args):
        if kind == 'unavailable':
            raise ModuleNotFoundError('missing OpenMM', name='openmm')
        if kind == 'unrelated-module':
            raise ModuleNotFoundError('unexpected dependency', name='other')
        raise RuntimeError('unexpected bug')
    monkeypatch.setattr(cycle_audit, 'assert_metadata_all_l', failure)
    if kind == 'unavailable':
        report = measure_cycle_geometry(problem, prepared, metadata, witness, phase='raw')
        assert report.global_l_status == 'unavailable' and not report.passed
    else:
        with pytest.raises(ModuleNotFoundError if kind == 'unrelated-module' else RuntimeError):
            measure_cycle_geometry(problem, prepared, metadata, witness, phase='raw')
