"""Exercise multi-path exclusive-write accounting with no scientific worker."""
import json

from scripts.research_post_hydrogen_native import run


def test_two_workers_have_separate_durable_records_and_one_final_index(tmp_path, monkeypatch):
    import scripts.research_post_hydrogen_native as module

    plan = dict(cases=[dict(label='unseeded', mode=m, static_index=i) for i, m in enumerate(('captured', 'Reference'))],
                protocol=dict(total_seconds=7200, per_path_seconds=1200))
    path = tmp_path/'plan.json'
    path.write_text(json.dumps(plan))
    calls = []

    def fake_invoke(plan, plan_path, index, folder, timeout):
        calls.append(index)
        (folder/'run.log').write_text('fake worker, no Context')
        return dict(status='completed' if index == 0 else 'execution-error', returncode=0 if index == 0 else 7)

    monkeypatch.setattr(module, 'check', lambda p: None)
    monkeypatch.setattr(module, 'invoke', fake_invoke)
    root = tmp_path/'run'
    run(path, root)
    result = json.loads((root/'evidence.json').read_text())
    assert calls == [0, 1]
    assert [r['returncode'] for r in result['cases']] == [0, 7]
    for i, row in enumerate(result['cases']):
        assert json.loads((root/f'path-{i:02d}-result.json').read_text()) == row


def test_worker_launch_error_is_retained_and_does_not_drop_next_case(tmp_path, monkeypatch):
    import scripts.research_post_hydrogen_native as module

    plan = dict(cases=[dict(label='x', mode='Reference', static_index=i) for i in range(2)],
                protocol=dict(total_seconds=7200, per_path_seconds=1200))
    path = tmp_path/'plan.json'
    path.write_text(json.dumps(plan))

    def failing_invoke(*args):
        raise OSError('launch blocked')

    monkeypatch.setattr(module, 'check', lambda p: None)
    monkeypatch.setattr(module, 'invoke', failing_invoke)
    root = tmp_path/'run'
    run(path, root)
    rows = json.loads((root/'evidence.json').read_text())['cases']
    assert len(rows) == 2
    assert all(r['status'] == 'orchestration-error' and r['returncode'] is None for r in rows)
