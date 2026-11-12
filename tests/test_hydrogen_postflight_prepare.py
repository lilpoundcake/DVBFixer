"""A geometry postcondition failure is not permission to resample hydrogens."""
import pytest


def test_prepare_geometry_error_does_not_fallback_or_publish(tmp_path, monkeypatch):
    pytest.importorskip('openmm')
    pytest.importorskip('pdbfixer')
    from openmm.app import Modeller
    from scripts.research_ser_hg_capture import synthetic_input

    from dvbfixer.ffutils import geometry
    from dvbfixer.prepare import main

    calls = []

    def add_hydrogens(model, *args, **kwargs):
        calls.append('native-H')

    error = geometry.HydrogenGeometryError([])

    def failed_repair(*args, **kwargs):
        raise error

    monkeypatch.setattr(Modeller, 'addHydrogens', add_hydrogens)
    monkeypatch.setattr(geometry, 'repair_misplaced_hydrogens', failed_repair)
    source, output = tmp_path/'input.pdb', tmp_path/'output.pdb'
    raw = synthetic_input()
    source.write_bytes(raw)
    with pytest.raises(geometry.HydrogenGeometryError) as caught:
        main([str(source), '-o', str(output), '--no-propka', '--no-protassign'])
    assert caught.value is error
    assert calls == ['native-H']
    assert not output.exists()
    assert not output.with_suffix('.dat').exists()
    assert source.read_bytes() == raw
