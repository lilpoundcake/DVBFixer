"""Offline accounting does not invent process success or rerun a worker."""
import gzip
import json

from scripts.recover_hydrogen_guard_summary import assemble_case


def test_missing_second_index_cannot_be_reported_as_success(tmp_path):
    for number in (1, 2, 3):
        folder = tmp_path/f'phase-{number}'
        folder.mkdir()
        (folder/'result.json').write_text(json.dumps(dict(returned=True, callbacks=0)))
        (folder/'path.json.gz').write_bytes(gzip.compress(json.dumps(dict(returned=True, callbacks=[])).encode()))
        (folder/'returned-state.xml').write_text('state verified separately in recovery')
    result = assemble_case(1, dict(label='unseeded', mode='Reference'), tmp_path, None)
    assert result['status'] == 'phases-returned-exit-unrecorded'
    assert result['returncode'] is None


def test_partial_phase_remains_partial(tmp_path):
    result = assemble_case(1, dict(label='unseeded', mode='Reference'), tmp_path, None)
    assert result['status'] == 'partial-exit-unrecorded'
    assert result['returncode'] is None
