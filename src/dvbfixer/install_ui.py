"""Dependency-free numbered menus for installer terminals."""

from __future__ import annotations

from dvbfixer.install_progress import waiting_for_input


def ask(prompt: str) -> str:
    with waiting_for_input():
        print("Waiting for your answer.", flush=True)
        return input(prompt)


def choose(title: str, options: list[tuple[str, str]], *, default: str | None = None) -> str:
    with waiting_for_input():
        return _choose(title, options, default=default)


def _choose(title: str, options: list[tuple[str, str]], *, default: str | None) -> str:
    print(f"\n{title}", flush=True)
    for number, (value, label) in enumerate(options, 1):
        suffix = " (default)" if value == default else ""
        print(f"  {number}. {label}{suffix}", flush=True)
    print("  0. Cancel", flush=True)
    print("Waiting for your choice — enter a number; no command is running.", flush=True)
    while True:
        answer = input("Choice: ").strip()
        if answer == "0":
            raise EOFError()
        if not answer and default is not None:
            return default
        if answer.isdecimal() and len(answer) < 8 and 1 <= int(answer) <= len(options):
            return options[int(answer) - 1][0]
        # Preserve existing scripted terminal answers, without requiring users
        # to type implementation identifiers in the normal numbered UI.
        if answer in {value for value, _ in options}:
            return answer
        print("Choose one of the displayed numbers.", flush=True)


def confirm(title: str, *, default: bool = False) -> bool:
    # yes/no aliases are accepted for older scripted terminals.
    return choose(title, [("yes", "Yes"), ("no", "No")], default="yes" if default else "no") == "yes"


def addons(options: list[tuple[str, str]]) -> list[str]:
    with waiting_for_input():
        return _addons(options)


def _addons(options: list[tuple[str, str]]) -> list[str]:
    print("\nOptional add-ons (comma-separated numbers; Enter for none):", flush=True)
    for number, (_, label) in enumerate(options, 1):
        print(f"  {number}. {label}", flush=True)
    print("  0. None", flush=True)
    print("Waiting for your selection (Enter for none).", flush=True)
    while True:
        answer = input("Add-on numbers: ").strip()
        if answer in {"", "0"}:
            return []
        values = []
        for item in answer.split(","):
            item = item.strip()
            if item.isdecimal() and len(item) < 8 and 1 <= int(item) <= len(options):
                values.append(options[int(item) - 1][0])
            elif item in {value for value, _ in options}:
                values.append(item)
            else:
                break
        else:
            return list(dict.fromkeys(values))
        print("Choose numbers from the list, or Enter for none.", flush=True)
