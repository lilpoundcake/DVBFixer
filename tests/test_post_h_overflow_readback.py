"""Fixed readback envelopes, no Context/native trajectory."""
from scripts.research_post_h_overflow_readback import compare


def test_original_envelopes_are_not_widened_to_hide_gradient_error():
    policy = dict(gradient_absolute_bound=1e-3, gradient_relative_bound=2e-5,
                  energy_absolute_bound=1e-3, energy_relative_bound=2e-6,
                  penalty_energy_absolute_bound=1e-3, penalty_energy_relative_bound=2e-5)
    stats = {'system energy': 12., 'restraint energy': 0., 'restraint strength': 0.}
    good = compare([[1., 2., 3.]], 12., [[0., 0., 0.]], [1., 2., 3.], stats, [], policy)
    bad = compare([[1., 2., 3.]], 12., [[0., 0., 0.]], [1., 2., 3.01], stats, [], policy)
    assert good['original_numerical_envelopes_passed']
    assert not bad['original_numerical_envelopes_passed']
    assert bad['gradient_maximum_error'] > 1e-3
