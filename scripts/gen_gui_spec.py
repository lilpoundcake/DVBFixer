#!/usr/bin/env python
"""Generate the DVBfixer GUI command schema directly from argparse parsers."""

from __future__ import annotations

import argparse
import importlib
import json
from pathlib import Path
from unittest.mock import patch

from dvbfixer.cli import COMMANDS

ROOT = Path(__file__).resolve().parent.parent
OUTPUT = ROOT / "gui" / "server" / "generated-dvbfixer-spec.ts"
MODULES = {
    "split": "dvbfixer.split_chains", "convert": "dvbfixer.glycam",
    "model": "dvbfixer.model", "prepare": "dvbfixer.prepare",
    "minimize": "dvbfixer.minimize", "top": "dvbfixer.top",
    "diagnose": "dvbfixer.diagnose",
}
CATEGORIES = {
    "split": "Structure preparation", "renumber": "Structure preparation",
    "model": "Structure preparation", "prepare": "Structure preparation",
    "pull": "Refinement", "minimize": "Refinement", "protonate": "Refinement",
    "homology": "Modeling & alignment", "msa": "Modeling & alignment",
    "salign": "Modeling & alignment", "convert": "Glycoprotein preparation",
    "transplant": "Glycoprotein preparation", "top": "Topology & chemistry",
    "parametrize": "Topology & chemistry", "cluster": "Analysis",
    "diagnose": "Analysis", "rename": "Utilities", "conect": "Utilities",
    "puppet": "Utilities", "doctor": "Utilities", "zbs": "Pipeline",
}
BATCH = {
    "split", "renumber", "model", "pull", "prepare", "minimize", "protonate",
    "rename", "convert", "conect", "puppet", "diagnose", "zbs",
}
ARTIFACT_DESTS = {
    "input", "topology", "trajectory", "acceptor", "donor", "graft", "template",
    "fasta", "alignment", "dat", "gaussian_log",
}
OUTPUT_EXTENSIONS = {
    "split": ".pdb", "renumber": ".pdb", "model": ".pdb", "pull": ".pdb",
    "prepare": ".pdb", "minimize": ".pdb", "protonate": ".pdb", "rename": ".pdb",
    "transplant": ".pdb", "puppet": ".pdb", "convert": ".pdb", "conect": ".pdb",
    "homology": ".pdb", "salign": ".pir", "msa": ".fasta", "diagnose": ".txt",
    "zbs": ".pdb", "top": ".top", "cluster": "", "parametrize": "",
}
OUTPUT_MODES = {
    "doctor": "stdout", "homology": "prefix", "cluster": "directory",
    "parametrize": "directory",
}


class _Captured(Exception):
    def __init__(self, parser: argparse.ArgumentParser):
        self.parser = parser


def _capture_parser(module_name: str) -> argparse.ArgumentParser:
    module = importlib.import_module(module_name)

    def capture(parser, args=None, namespace=None):  # noqa: ANN001, ARG001
        raise _Captured(parser)

    try:
        with patch.object(argparse.ArgumentParser, "parse_args", capture):
            module.parse_args([])
    except _Captured as captured:
        return captured.parser
    raise RuntimeError(f"could not capture parser for {module_name}")


def _plain_default(value):  # noqa: ANN001, ANN201
    if value is None or value is argparse.SUPPRESS:
        return None
    if isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, (list, tuple)) and all(isinstance(item, (str, int, float)) for item in value):
        return list(value)
    return None


def _field(action: argparse.Action, group: str) -> dict:
    option_strings = list(action.option_strings)
    positive = next((item for item in option_strings if item.startswith("--") and not item.startswith("--no-")), None)
    flag = positive or next((item for item in option_strings if item.startswith("--")), None)
    flag = flag or (option_strings[-1] if option_strings else action.dest)
    choices = list(action.choices) if action.choices is not None else None
    is_bool = isinstance(action, (argparse._StoreTrueAction, argparse._StoreFalseAction,
                                  argparse.BooleanOptionalAction))
    if is_bool:
        field_type = "bool"
    elif choices:
        field_type = "select"
    elif action.type in (int, float):
        field_type = "number"
    elif action.dest in ARTIFACT_DESTS:
        field_type = "artifact"
    else:
        field_type = "text"
    result = {
        "flag": flag,
        "dest": action.dest,
        "label": action.dest.replace("_", " ").title(),
        "type": field_type,
        "group": group,
        "help": action.help if action.help is not argparse.SUPPRESS else None,
        "required": bool(action.required),
        "repeatable": isinstance(action, argparse._AppendAction),
        "multi": action.nargs in ("+", "*") or (isinstance(action.nargs, int) and action.nargs > 1),
    }
    default = (_plain_default(action.default) if isinstance(action, argparse.BooleanOptionalAction)
               else False if is_bool else _plain_default(action.default))
    if default is not None:
        result["default"] = default
    if choices:
        result["options"] = choices
    false_flag = (next((item for item in option_strings if item.startswith("--no-")), None)
                  if isinstance(action, argparse.BooleanOptionalAction) else None)
    if false_flag:
        result["falseFlag"] = false_flag
    return {key: value for key, value in result.items() if value is not None}


def command_schema(name: str, description: str) -> dict:
    module_name = MODULES.get(name, f"dvbfixer.{name}")
    parser = _capture_parser(module_name)
    groups: list[dict] = []
    inputs: list[dict] = []
    flags: list[dict] = []
    has_output = False
    for group in parser._action_groups:
        if group.title in {"Global logging", "Batch mode"}:
            continue
        group_fields: list[str] = []
        for action in group._group_actions:
            if isinstance(action, argparse._HelpAction):
                continue
            if action.dest == "output" or "--output" in action.option_strings:
                has_output = True
                continue
            field = _field(action, group.title)
            if action.option_strings:
                flags.append(field)
                group_fields.append(field["flag"])
            else:
                field["name"] = action.dest
                field["nargs"] = action.nargs
                inputs.append(field)
        if group_fields:
            groups.append({"name": group.title, "fields": group_fields})
    return {
        "name": name,
        "label": name.replace("_", " ").title(),
        "description": description,
        "category": CATEGORIES[name],
        "inputs": inputs,
        "flags": flags,
        "groups": groups,
        "outputExtension": OUTPUT_EXTENSIONS.get(name, ".pdb"),
        "outputMode": OUTPUT_MODES.get(name, "file"),
        "hasOutput": has_output,
        "outputKind": "report" if name in {"doctor", "diagnose"} else "artifact",
        "batch": name in BATCH,
        "successCodes": [0, 1] if name == "diagnose" else [0],
        "specialized": name == "homology",
    }


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Fail if the committed schema is stale")
    args = parser.parse_args(argv)
    schema = [command_schema(name, description) for name, description in COMMANDS.items()]
    payload = json.dumps(schema, indent=2, ensure_ascii=False)
    content = (
        "// Generated by scripts/gen_gui_spec.py; do not edit by hand.\n"
        "import type { CommandDef } from './dvbfixer-spec'\n\n"
        f"export const GENERATED_COMMANDS = {payload} as CommandDef[]\n"
    )
    if args.check:
        if not OUTPUT.exists() or OUTPUT.read_text() != content:
            raise SystemExit("GUI command schema is stale; run python scripts/gen_gui_spec.py")
        print("GUI command schema is in sync with argparse.")
        return
    OUTPUT.write_text(content)
    print(f"Wrote {OUTPUT.relative_to(ROOT)} ({len(schema)} commands)")


if __name__ == "__main__":
    main()
