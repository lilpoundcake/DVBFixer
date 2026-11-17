from types import SimpleNamespace

from scripts.export_retained_return_review import generated_record, native_provenance

from dvbfixer.ffutils.dat import DatRecord
from dvbfixer.model.diffusion.contract import AtomIdentity


def test_private_sidecar_contains_only_generated_atoms_with_case_and_icode(tmp_path):
    known, generated = AtomIdentity("A", "82", "A", "CA"), AtomIdentity("a", "82", "A", "CA")
    template = "ATOM  {serial:5d}  CA  ALA {chain}  82A      1.000   2.000   3.000  1.00 20.00           C\n"
    text = template.format(serial=1, chain="A")+template.format(serial=2, chain="a")
    request = SimpleNamespace(generated_atoms=(generated,), fixed_atoms=(known,))
    path = tmp_path / "generated.dat"
    generated_record(request, text).save(path)
    record = DatRecord.load(path)
    assert len(record.added_atoms) == 1
    assert record.added_atoms[0] == {"chain": "a", "resid": "82", "icode": "A", "resname": "ALA", "atom": "CA", "element": "C"}


def test_review_never_relabels_native_modeller_seed_as_optimizer_or_neural_seed():
    result = native_provenance({"native_backend": "MODELLER", "native_seed": 42}, SimpleNamespace(seeds=(42,)), 11)
    assert result["native_backend"] == "MODELLER"
    assert result["native_seed"] == 42 and result["optimizer_seed"] == 11
    assert "MODELLER proposal completed" in result["native_full_inference_status"]
