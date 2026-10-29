"""Legacy spelling is accepted but not advertised to CLI/GUI users."""

import argparse

from dvbfixer.model.context_options import add_context_argument


def test_public_choices_hide_none_without_breaking_legacy_parser():
    parser = argparse.ArgumentParser()
    add_context_argument(parser)
    action = next(a for a in parser._actions if a.dest == "diffusion_context")
    assert list(action.choices) == ["auto", "local", "full"]
    assert "{auto,local,full}" in parser.format_help()
    assert "none" not in parser.format_help()
    assert parser.parse_args([]).diffusion_context == "auto"
    assert parser.parse_args(["--diffusion-context", "none"]).diffusion_context == "none"
