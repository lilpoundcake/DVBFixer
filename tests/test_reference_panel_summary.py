"""Panel denominator control; fake audit numbers are not model evidence."""
from scripts.summarize_reference_reconstruction_panel import summarize_rows


def test_panel_keeps_failures_and_does_not_score_only_accepted_seeds():
    rows = [{'seed': 7, 'returncode': 1},
        {'seed': 42, 'returncode': 0, 'audit': {'heldout_generated_rmsd_angstrom': 1.,
            'assessments': {'source-restored': {'geometry_passed': True, 'admission_failures': []}}}},
        {'seed': 101, 'returncode': 0, 'audit': {'heldout_generated_rmsd_angstrom': 9.,
            'assessments': {'source-restored': {'geometry_passed': False, 'admission_failures': []}}}},
        {'seed': 102, 'returncode': 0}]
    result = summarize_rows(rows)
    assert result['declared'] == 4
    assert result['executed_successfully'] == 3
    assert result['independently_audited'] == 2
    assert result['source_restored_geometry_passes'] == 1
    assert result['technical_or_audit_failures'] == [7, 102]
    assert result['all_audited_heldout_heavy_rmsd_median'] == 5.
    assert result['per_seed'][0]['heldout_heavy_rmsd'] is None


def test_panel_has_no_zero_error_imputation_when_every_output_failed():
    result = summarize_rows([{'seed': 7, 'returncode': 1}])
    assert result['all_audited_heldout_heavy_rmsd_median'] is None
    assert result['source_restored_geometry_passes'] == 0


def test_restored_boltz_geometry_does_not_count_as_native_source_preservation():
    result = summarize_rows([{'seed': 7, 'returncode': 0, 'audit': {
        'heldout_generated_rmsd_angstrom': 1., 'assessments': {
            'aligned-native': {'geometry_passed': True, 'admission_failures': ['fixed-heavy-atom-rmsd']},
            'source-restored': {'geometry_passed': True, 'admission_failures': []}}}}])
    assert result['source_restored_geometry_passes'] == 1
    assert result['native_source_and_geometry_passes'] == 0
    assert result['per_seed'][0]['native_source_geometry_pass'] is False
