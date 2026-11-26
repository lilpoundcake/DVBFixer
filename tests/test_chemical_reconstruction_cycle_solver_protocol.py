"""R1 frozen metadata only: these tests do not execute or qualify its constructor."""
from __future__ import annotations

import hashlib
import json
from dataclasses import asdict
from pathlib import Path

import pytest

from dvbfixer.chemical_reconstruction.observed_peptide import _parse_source, select_interior_target

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / 'docs/plans/chemical-reconstruction-cycle-solver-v1.md'
LEDGER = ROOT / 'docs/plans/chemical-reconstruction-cycle-solver-cohort-v1.json'
PROTOCOL_SHA256 = '88cb776bdf3bb35053d9e463307d8d36ca8b765bab8b18e8d1350bd87ff3bd7b'
LEDGER_SHA256 = 'f6c097f6c14d1394b2fc0b9c64984fd7a35eb488b7e207d73137b8407f695075'


def ledger():
    return json.loads(LEDGER.read_bytes())


def test_R1_protocol_and_ledger_are_frozen_before_numerical_work():
    assert hashlib.sha256(PROTOCOL.read_bytes()).hexdigest() == PROTOCOL_SHA256
    assert hashlib.sha256(LEDGER.read_bytes()).hexdigest() == LEDGER_SHA256
    data = ledger()
    assert data['baseline_commit'] == '4565ac8'
    assert data['stage'] == 'prospective-protocol-only-no-R1-coordinate-execution'
    assert not data['source_admission_or_geometry_evaluated']
    assert not data['implementation_ready'] and not data['coordinate_execution_completed']
    assert not data['production_or_general_qualification']
    assert not data['H_parent_or_public_capability_enabled']


@pytest.mark.parametrize('source_id', ['1TIM', '1CRN', '2CI2'])
@pytest.mark.parametrize('parent', ['PRO', 'PHE', 'TYR', 'TRP'])
def test_identity_only_enrollment_including_absence_does_not_generate_coordinates(source_id, parent, monkeypatch):
    from dvbfixer.chemical_reconstruction import construction, metric_projection

    def forbidden(*args, **kwargs):
        pytest.fail('protocol metadata validation called a coordinate constructor')
    monkeypatch.setattr(construction, 'propose', forbidden)
    monkeypatch.setattr(construction, '_internal', forbidden)
    monkeypatch.setattr(construction, '_template', forbidden)
    monkeypatch.setattr(metric_projection, '_solve', forbidden)
    data = ledger()
    row = next(r for r in data['sources'] if r['id'] == source_id)
    source = (ROOT / row['path']).read_bytes()
    assert hashlib.sha256(source).hexdigest() == row['sha256']
    _, atoms = _parse_source(source)
    chain = next(a.chain_id for a in atoms if a.record == 'ATOM  ')
    assert chain == row['chain']  # Chain IDs never normalized.
    expected_number = row['targets'][parent]
    if expected_number is None:
        assert (source_id, parent) == ('1CRN', 'TRP')
        with pytest.raises(ValueError, match=row['absent_target_reasons'][parent]):
            select_interior_target(source, parent, chain=chain)
    else:
        selected = select_interior_target(source, parent, chain=chain)
        assert asdict(selected) == dict(data['identity_defaults'], chain_id=chain, sequence_number=expected_number)
    # Identity presence/absence, not source admission, geometry or native quality.


def test_all_seven_primary_subsets_and_ring_profiles_keep_30_cell_denominator():
    data = ledger()
    assert len(data['sources']) == data['source_systems'] == 3
    profiles = {p['parent']: p for p in data['profiles']}
    pro = profiles['PRO']
    wanted = {'CB', 'CG', 'CD'}
    masks = {frozenset(m) for m in pro['masks']}
    assert len(masks) == 7 and all(m and m <= wanted for m in masks)
    assert frozenset({'CB', 'CD'}) in masks  # Separate generated islands, same problem.
    assert [p['parent'] for p in data['profiles']] == ['PRO', 'PHE', 'TYR', 'TRP']
    assert all(p['state'] == p['parent'] for p in data['profiles'])  # No implicit HIS state.
    cells = [(s['id'], p['parent'], tuple(m)) for s in data['sources'] for p in data['profiles'] for m in p['masks']]
    assert len(set(cells)) == len(cells) == data['declared_cells'] == 30
    assert sum(p == 'PRO' for _, p, _ in cells) == data['primary_PRO_cells'] == 21
    assert sum(p != 'PRO' for _, p, _ in cells) == data['secondary_aromatic_cells'] == 9
    absent = [(s['id'], p['parent']) for s in data['sources'] for p in data['profiles']
              for _ in p['masks'] if s['targets'][p['parent']] is None]
    assert absent == [('1CRN', 'TRP')] and data['metadata_absent_cells'] == 1
    assert 'OH' in profiles['TYR']['masks'][0] and 'CH2' in profiles['TRP']['masks'][0]


def test_old_nine_cells_are_retrospective_regressions_not_new_independent_sources():
    old = json.loads((ROOT / 'docs/plans/chemical-reconstruction-comparison-cohort-v1.json').read_bytes())
    data = ledger()
    pro_masks = next(p['masks'] for p in data['profiles'] if p['parent'] == 'PRO')
    assert all(m in pro_masks for m in old['masks'])
    for enrolled, previous in zip(data['sources'], old['sources'], strict=True):
        assert all(enrolled[k] == previous[k] for k in ('id', 'path', 'sha256', 'chain'))
        assert enrolled['targets']['PRO'] == previous['residue_number']
        assert data['identity_defaults']['insertion_code'] == previous['insertion_code']
    assert data['comparison']['old_nine_comparator_repeats_are_new_independent_data'] is False


def test_numerical_settings_are_declared_not_executed_or_native_selected():
    data = ledger()
    start, solve = data['initializer'], data['solver']
    assert start['branches'] == [[sign, k] for sign in (1, -1) for k in range(4)]
    assert not start['external_start_input'] and not start['native_PDBFixer_rotamer_or_fitted_ideal_pose_input']
    assert not start['require_initial_geometry_pass'] and not start['deduplicate_or_success_short_circuit']
    assert solve['max_nfev_per_branch'] == 600
    assert solve['maximum_optimizer_residual_calls_per_branch_including_finite_difference'] == 50000
    assert solve['raw_and_rounded_per_atom_movement_cap_angstrom'] == 8.0
    assert solve['observation_movement_allowed_angstrom'] == 0.0
    assert solve['maximum_numerical_residual_scalars_including_context'] == 65536
    assert not solve['regularizer_forcefield_reflection_restarts_or_fallback']
    assert not data['selection']['native_reference_used']
    assert data['comparison']['H_or_old_projection_calls'] == 0
    assert data['comparison']['required_pass_or_win_rate'] is None


def test_historical_document_and_R0_contract_pins_are_not_replaced():
    pins = ledger()['baseline_pins']
    for name, expected in pins.items():
        assert len(expected) == 64 and all(c in '0123456789abcdef' for c in expected)
        # Source hashes record the named baseline, not a ban on separately declared
        # future factoring/codec additions. Historical docs and R0 remain immutable.
        if name.startswith('docs/') or name.endswith('/cycle_contract.py'):
            assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == expected


def test_authorities_are_pinned_existing_assets_not_bundled_new_tables():
    for name, expected in ledger()['authority_pins'].items():
        path = ROOT / name
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected
        if path.name == 'manifest.json':
            for component, item in json.loads(path.read_bytes())['components'].items():
                assert hashlib.sha256((path.parent / f'{component}.cif').read_bytes()).hexdigest() == item['sha256']
