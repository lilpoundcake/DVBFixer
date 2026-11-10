"""Outcome-blind whole-mask inventory and seen/unseen-variant stratification."""
import copy

import pytest
from scripts import research_ncac_expansion as module


def rows():
    return [dict(id='first-partial', mask=['CG1', 'CG2'], native_rmsd=0., common_accepted=True),
            dict(id='seen-whole', mask=['CB', 'CG1', 'CG2'], native_rmsd=99., common_accepted=False),
            dict(id='new-singleton', mask=['OG1'], native_rmsd=0., common_accepted=True),
            dict(id='new-whole', mask=['CB', 'OG1', 'CG2'], native_rmsd=99., common_accepted=False)]


def test_all_whole_masks_retained_in_original_order_regardless_of_prior_failures():
    source = rows()
    before = copy.deepcopy(source)
    selected = module.whole_cases(source, ['seen-whole'])
    assert [c['id'] for c in selected] == ['seen-whole', 'new-whole']
    assert [c['ncac_previously_executed'] for c in selected] == [True, False]
    assert source == before


def test_accuracy_and_common_acceptance_do_not_control_enrollment():
    source = rows()
    first = module.whole_cases(source, ['seen-whole'])
    for row in source:
        row['native_rmsd'] = -999.
        row['common_accepted'] = not row['common_accepted']
    second = module.whole_cases(source, ['seen-whole'])
    assert [(c['id'], c['ncac_previously_executed']) for c in first] == [(c['id'], c['ncac_previously_executed']) for c in second]


def test_previous_case_omission_and_duplicate_whole_identity_fail_closed():
    with pytest.raises(ValueError, match='omits previously tested'):
        module.whole_cases(rows(), ['absent'])
    with pytest.raises(ValueError, match='duplicates'):
        module.whole_cases(rows()+[rows()[1]], ['seen-whole'])


def test_empty_matched_denominator_is_none_not_an_invented_zero():
    from scripts.audit_ncac_expansion import median

    assert median([]) is None
    assert median([1., 3.]) == 2.


def test_accuracy_reports_rounded_coordinates_and_generated_only_symmetry_without_selection():
    from scripts.audit_ncac_expansion import accuracy

    from dvbfixer.domain.structure_identity import ComponentInstanceRef
    from dvbfixer.pdb_component_reconstruction import _format_generated_atom
    from dvbfixer.protein_completion.branched_source import plan_branched_source
    from dvbfixer.protein_completion.service import complete_sidechain
    from dvbfixer.protein_completion.source import _insert_candidates

    target = ComponentInstanceRef(1, 'I', '22')
    points = {'N': (-1.46, 0., 0.), 'CA': (0., 0., 0.), 'C': (.55, 1.43, 0.), 'O': (1.76, 1.66, 0.)}
    source = ''.join('ATOM  '+_format_generated_atom(serial=i, name=n, element=n[0], position=xyz,
                                                   residue_name='VAL', chain_id='I', sequence_number='22', insertion_code='')[6:]
                     for i, (n, xyz) in enumerate(points.items(), 1)).encode()+b'TER\nEND\n'
    request = plan_branched_source(source, target).request
    result = complete_sidechain(request)
    candidate = result.selected
    pdb = _insert_candidates(source, (request,), (candidate,))
    truth = points | {a.atom.name: a.xyz for a in candidate.additions}
    truth['CG1'], truth['CG2'] = truth['CG2'], truth['CG1']
    before = copy.deepcopy(candidate)
    report = accuracy(request, candidate.additions, pdb, truth)
    assert report['raw']['exact_name_rmsd_angstrom'] > 0
    assert report['raw']['generated_only_val_symmetry_rmsd_angstrom'] == 0
    assert report['rounded']['generated_only_val_symmetry_rmsd_angstrom'] < .001
    assert report['native_reporting_only'] and candidate == before
