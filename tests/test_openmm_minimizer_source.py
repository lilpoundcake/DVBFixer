"""Source provenance is read-only and never attests a compiled binary."""
import io
import json
from contextlib import nullcontext

import pytest
from scripts import capture_openmm_minimizer_source as collector


def test_source_notices_hashes_and_declared_revision_are_retained(tmp_path, monkeypatch):
    openmm = pytest.importorskip('openmm')
    monkeypatch.setattr(openmm.version, 'git_revision', '1'*40)
    urls = []
    def response(url, timeout):
        urls.append(url)
        assert timeout == 60
        return nullcontext(io.BytesIO(b'/* Copyright retained. Synthetic test source. */\n'))
    monkeypatch.setattr(collector.urllib.request, 'urlopen', response)
    output = tmp_path/'source'
    collector.capture(output)
    record = json.loads((output/'manifest.json').read_bytes())
    assert len(urls) == len(collector.FILES) == 5
    assert all('/'+'1'*40+'/' in url for url in urls)
    assert all((output/f['path']).exists() for f in record['files'])
    assert record['upstream_notices_retained']
    assert not record['runtime_binary_attestation'] and not record['installed_runtime_modified']


def test_unversioned_runtime_refuses_before_creating_a_source_tree(tmp_path, monkeypatch):
    openmm = pytest.importorskip('openmm')
    monkeypatch.setattr(openmm.version, 'git_revision', 'main')
    output = tmp_path/'source'
    with pytest.raises(ValueError, match='exact upstream revision'):
        collector.capture(output)
    assert not output.exists()


def test_missing_source_notice_does_not_publish_a_provenance_manifest(tmp_path, monkeypatch):
    openmm = pytest.importorskip('openmm')
    monkeypatch.setattr(openmm.version, 'git_revision', '1'*40)
    monkeypatch.setattr(collector.urllib.request, 'urlopen', lambda *a, **k:nullcontext(io.BytesIO(b'not source')))
    output = tmp_path/'source'
    with pytest.raises(ValueError, match='unexpected upstream'):
        collector.capture(output)
    assert output.exists() and not (output/'manifest.json').exists()
