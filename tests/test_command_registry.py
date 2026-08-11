from __future__ import annotations

from dvbfixer.batch import OUTPUT_SUFFIXES
from dvbfixer.cli import COMMANDS
from dvbfixer.command_registry import (
    COMMAND_BY_NAME,
    COMMAND_REGISTRY,
    validate_command_registry,
)


def test_command_registry_is_complete_and_unique() -> None:
    validate_command_registry()
    names = [command.name for command in COMMAND_REGISTRY]
    assert len(names) == len(set(names))
    assert COMMANDS == {
        command.name: command.description for command in COMMAND_REGISTRY
    }
    assert COMMAND_BY_NAME == {
        command.name: command for command in COMMAND_REGISTRY
    }
    assert all(command.module.startswith("dvbfixer.") for command in COMMAND_REGISTRY)
    assert all(command.category for command in COMMAND_REGISTRY)
    assert all(command.success_codes for command in COMMAND_REGISTRY)
    assert "glycam" not in COMMAND_BY_NAME


def test_batch_suffixes_are_derived_from_command_registry() -> None:
    assert OUTPUT_SUFFIXES == {
        command.name: command.batch_output_suffix
        for command in COMMAND_REGISTRY
        if command.batch_output_suffix is not None
    }


def test_generators_use_the_same_registry_metadata() -> None:
    from scripts import gen_cli_reference, gen_gui_spec

    assert gen_cli_reference.COMMANDS == [
        (command.name, command.module) for command in COMMAND_REGISTRY
    ]
    assert gen_gui_spec.MODULES == {
        command.name: command.module for command in COMMAND_REGISTRY
    }
    assert gen_gui_spec.CATEGORIES == {
        command.name: command.category for command in COMMAND_REGISTRY
    }
    assert gen_gui_spec.BATCH == {
        command.name for command in COMMAND_REGISTRY if command.batch
    }
    assert gen_gui_spec.OUTPUT_EXTENSIONS == {
        command.name: command.output_extension for command in COMMAND_REGISTRY
    }
    assert gen_gui_spec.OUTPUT_MODES == {
        command.name: command.output_mode for command in COMMAND_REGISTRY
    }
