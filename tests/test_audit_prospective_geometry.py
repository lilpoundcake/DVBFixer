"""Original-frame accuracy never mixes relaxed coordinates with original anchors."""
import pytest
from scripts.audit_prospective_geometry import target_errors


def test_all_generated_atom_errors_use_exact_output_coordinates_without_alignment():
    truth=dict(CB=(0.,0.,0.),OG1=(1.,0.,0.),CG2=(0.,1.,0.))
    moved={n:(p[0]+2,p[1],p[2]) for n,p in truth.items()}
    result=target_errors(moved,truth,tuple(truth),'THR')
    assert result['exact_name_rmsd_angstrom']==2.
    assert not result['aligned'] and not result['used_for_selection']


def test_unavailable_truth_is_an_explicit_refusal_not_a_partial_accuracy_score():
    with pytest.raises(ValueError,match='lacks'):
        target_errors(dict(CB=(0.,0.,0.)),{},('CB',),'THR')
