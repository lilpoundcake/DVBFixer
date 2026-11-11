from pathlib import Path

import pytest
from scripts.run_rfd3_context_development_panel import sampling_arguments


def test_context_input_prefix_is_not_rewritten_as_output_directory():
    inputs = Path('/project/context-input')
    work = Path('/project/panel/seed-42')
    result = sampling_arguments(['design', 'inputs=/project/context-input/rfd3-input.json',
        'out_dir=/project/context/predictions', 'seed=7', 'ckpt_path=/assets/model.ckpt'],
        inputs=inputs, work=work, seed=42)
    assert 'inputs=/project/context-input/rfd3-input.json' in result
    assert 'out_dir=/project/panel/seed-42/run/predictions' in result
    assert 'seed=42' in result
    assert 'ckpt_path=/assets/model.ckpt' in result
    assert all('run-input' not in argument for argument in result)


def test_context_arguments_refuse_ambiguous_original_inputs():
    with pytest.raises(ValueError, match='ambiguous'):
        sampling_arguments(['inputs=one', 'inputs=two', 'out_dir=result'],
                           inputs=Path('/input'), work=Path('/output'), seed=42)
