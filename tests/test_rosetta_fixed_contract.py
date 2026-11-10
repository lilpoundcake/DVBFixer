"""Atom-fixed score selection and topology diagnostics never relax old gates."""
from __future__ import annotations

import math
from types import SimpleNamespace

import pytest
from scripts.research_rosetta_fixed_scoring import physical_snapshot
from scripts.rosetta_fixed_contract import (
    choose_scored_own,
    contacts,
    identity,
    import_audit,
    indexed,
    short_paths,
)
from scripts.summarize_rosetta_fixed_scoring import source_records_preserved


def atom(name, xyz=(0., 0., 0.), *, chain='a', icode='A', resid=82):
    return dict(chain=chain, resid=resid, icode=icode, resname='VAL', atom=name, element='C', xyz=list(xyz))


def test_exact_import_audit_rejects_sub_tolerance_motion_identity_and_topology_changes():
    a, b = atom('CG1'), atom('C', (3., 0., 0.))
    edge = (identity(a), identity(b))
    assert import_audit([a, b], [a, b], [edge], [edge], [], [], [])['accepted']
    tiny_move = atom('CG1', (1e-10, 0., 0.))
    assert not import_audit([a, b], [tiny_move, b], [edge], [edge], [], [], [])['accepted']
    for changed in (atom('CG1', chain='A'), atom('CG1', icode=' ')):
        assert not import_audit([a, b], [changed, b], [], [], [], [], [])['accepted']
    assert not import_audit([a, b], [a, b], [], [edge], [], [], [])['accepted']
    assert not import_audit([a, b], [a, b], [edge], [], [], [], [])['accepted']
    assert not import_audit([a, b], [a, b], [], [], [edge], [], [])['accepted']


def test_only_predeclared_extra_atom_is_accepted_and_is_reported():
    a, extra = atom('C'), atom('OXT', (1., 0., 0.))
    assert not import_audit([a], [a, extra], [], [], [], [], [])['accepted']
    result = import_audit([a], [a, extra], [], [], [], [], [identity(extra)])
    assert result['accepted'] and result['scorer_owned_extras'] == [identity(extra)]


@pytest.mark.parametrize('bad', [math.nan, math.inf, True])
def test_invalid_or_duplicate_snapshot_refuses(bad):
    with pytest.raises(ValueError, match='coordinates'):
        indexed([atom('C', (bad, 0., 0.))])
    with pytest.raises(ValueError, match='duplicate'):
        indexed([atom('C'), atom('C')])


def test_native_like_1_4_clash_is_classified_not_suppressed_and_disconnected_is_unresolved():
    rows = [atom('CG1', (2.762756, 0., 0.)), atom('CB', (3., 1., 0.)),
            atom('CA', (1.5, 1., 0.)), atom('C'), atom('X', (2.8, .1, 0.), chain='A')]
    keys = {a['atom']: identity(a) for a in rows}
    edges = [(keys[a], keys[b]) for a, b in [('CG1', 'CB'), ('CB', 'CA'), ('CA', 'C')]]
    results = contacts(rows, [keys['CG1']], edges, {'C': 1.7})
    native = next(r for r in results if r['right'] == keys['C'])
    assert native['topology'] == 'known_1_4' and native['partial_graph_bond_distance'] == 3
    assert native['old_finding_retained'] and not native['both_generated']
    unknown = next(r for r in results if r['right'] == keys['X'])
    assert unknown['topology'] == 'no_short_path_in_partial_graph'
    assert unknown['partial_graph_bond_distance'] is None
    assert not any(r['right'] in {keys['CB'], keys['CA']} for r in results)


def test_shortest_path_wins_over_longer_walk_and_real_new_new_collision_is_retained_once():
    rows = [atom('CG1'), atom('CG2', (.1, 0., 0.))]
    a, b = [identity(r) for r in rows]
    assert len(contacts(rows, [a, b], [], {'C': 1.7})) == 1
    assert contacts(rows, [a, b], [], {'C': 1.7})[0]['both_generated']
    x, y = identity(atom('CB')), identity(atom('CA'))
    assert short_paths(a, [(a, x), (x, y), (y, b), (a, b)])[b] == 1
    with pytest.raises(ValueError, match='endpoint absent'):
        contacts(rows, [a], [(a, x), (x, b)], {'C': 1.7})


def candidate(index, score, *, common=True, source=True, status='scored'):
    return dict(index=index, total_score=score, old_source_accepted=source,
                common_accepted=common, score_status=status)


def test_low_score_cannot_admit_a_rejected_candidate_or_use_a_truth_field():
    rows = [candidate(0, -100., common=False), candidate(1, -50., source=False),
            candidate(2, 10.), candidate(3, 20.)]
    rows[3]['native_rmsd'] = 0.
    assert choose_scored_own(rows)['selected_index'] == 2
    assert choose_scored_own(rows[:2])['status'] == 'no_old_gate_feasible_candidate'


@pytest.mark.parametrize('score', [None, math.nan, math.inf, True])
def test_missing_feasible_score_blocks_selection_instead_of_cherry_picking(score):
    assert choose_scored_own([candidate(0, score), candidate(1, 1.)])['selected_index'] is None
    assert choose_scored_own([candidate(0, 0., status='import_contract_blocker'), candidate(1, 1.)])['status'] == 'selection_blocked_incomplete_scores'


def test_ties_are_stable_by_index_and_duplicate_indices_refuse():
    assert choose_scored_own([candidate(2, 1.), candidate(0, 1.)])['selected_index'] == 0
    with pytest.raises(ValueError, match='duplicate'):
        choose_scored_own([candidate(0, 1.), candidate(0, 2.)])


@pytest.mark.parametrize('is_virtual', [True, False])
def test_virtual_inventory_uses_atom_type_not_the_nv_name(is_virtual):
    residue = SimpleNamespace(natoms=lambda: 1, nheavyatoms=lambda: 1,
                              atom_type=lambda a: SimpleNamespace(is_virtual=lambda: is_virtual, name=lambda: 'VIRT' if is_virtual else 'Nbb'),
                              xyz=lambda a: SimpleNamespace(x=1., y=2., z=3.),
                              name3=lambda: 'PRO', atom_name=lambda a: 'NV')
    info = SimpleNamespace(chain=lambda i: 'a', number=lambda i: 82, icode=lambda i: 'A')
    pose = SimpleNamespace(total_residue=lambda: 1, residue=lambda i: residue, pdb_info=lambda: info)
    physical, virtual = physical_snapshot(pose)
    assert len(physical) == int(not is_virtual) and len(virtual) == int(is_virtual)
    if is_virtual:
        assert virtual[0]['virtual_by_rosetta_type'] and virtual[0]['rosetta_atom_type'] == 'VIRT'


def test_source_record_check_detects_noncoordinate_edits_and_extra_headers():
    line = b'ATOM      1  C   VAL A  39       1.000   2.000   3.000  1.00 20.00           C  \n'
    added = b'ATOM      2  CB  VAL A  39       1.500   2.000   3.000  1.00  0.00           C  \n'
    source = line+b'TER\nEND\n'
    candidate = line+added+b'TER\nEND\n'
    generated = [('A', 39, ' ', 'VAL', 'CB')]
    assert source_records_preserved(source, candidate, generated)
    assert not source_records_preserved(source, candidate.replace(b'20.00', b'21.00'), generated)
    assert not source_records_preserved(source, candidate+b'REMARK unapproved\n', generated)
