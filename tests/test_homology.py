from types import SimpleNamespace

import pytest

from dvbfixer.homology import _make_model_with_loop_fallback


class _BaseModel:
    outputs = [{"failure": None, "name": "model.pdb", "molpdf": 1.0}]

    def __init__(self, env, **kwargs):
        self.env = env
        self.kwargs = kwargs
        self.loop = SimpleNamespace()


def test_no_loop_alignment_falls_back_to_automodel():
    class NoLoopModel(_BaseModel):
        def make(self):
            raise RuntimeError("No loops detected for refinement: you must redefine select_loop_atoms")

    class AutoModel(_BaseModel):
        def make(self):
            self.made = True

    model, used_loops = _make_model_with_loop_fallback(
        object(), AutoModel, NoLoopModel, "alignment.pir", ("template",),
        "target", 2, False, "fast",
    )
    assert isinstance(model, AutoModel)
    assert model.made
    assert used_loops is False


def test_unrelated_modeller_error_is_not_hidden():
    class BrokenLoopModel(_BaseModel):
        def make(self):
            raise RuntimeError("Sequence difference between alignment and pdb")

    with pytest.raises(RuntimeError, match="Sequence difference"):
        _make_model_with_loop_fallback(
            object(), _BaseModel, BrokenLoopModel, "alignment.pir", "template",
            "target", 1, False, "fast",
        )
