from types import SimpleNamespace

from scripts.inspect_retained_loop_context_contracts import profile_view

from dvbfixer.model.diffusion.contract import BackendOption


def test_profile_view_is_immutable_and_preserves_preparation_binding():
    from dataclasses import dataclass

    @dataclass(frozen=True)
    class Request:
        backend_options: tuple
        fixed_atoms: object

    original = Request((BackendOption("profile", "old"), BackendOption("input_preparation.revision", "bound")),
                       SimpleNamespace(observations="unchanged"))
    view = profile_view(original, "new")
    assert original.backend_options[0].value == "old"
    assert view.fixed_atoms is original.fixed_atoms
    assert [o.value for o in view.backend_options if o.name == "profile"] == ["new"]
    assert original.backend_options[1] in view.backend_options
