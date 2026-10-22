"""CLI boundary for alternate-state selection, including CIF normalization."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
from pathlib import Path

from dvbfixer.altloc_selection import AltlocError, select_altlocs


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog="dvbfixer altlocs", description="Analyze or select exactly one alternate state per residue; never modify the input.")
    io = parser.add_argument_group("Input / output")
    io.add_argument("input", help="PDB, ENT, CIF or mmCIF structure")
    io.add_argument("-o", "--output", help="New single-state PDB; omit to analyze only")
    io.add_argument("--report-json", help="JSON analysis and selection report")
    io.add_argument("--selection-plan", help="Version-1 JSON plan bound to source SHA-256")
    options = parser.add_argument_group("Selection")
    options.add_argument("--policy", choices=("manual", "occupancy", "first", "prefer"), default="occupancy", help="Residue-coherent automatic or manual state selection")
    options.add_argument("--prefer", help="Preferred altloc label (case-sensitive)")
    options.add_argument("--fallback", choices=("error", "occupancy", "first"), default="error", help="Fallback when preferred label is absent")
    options.add_argument("--select", action="append", default=[], metavar="CHAIN:RESID:LABEL", help="Manual override (repeatable); insertion codes supported, _ means blank chain; ambiguous MODEL/occurrence requires JSON plan")
    options.add_argument("--dry-run", action="store_true", help="Analyze without publishing a structure")
    options.add_argument("--include-structure", action="store_true", help="Include normalized PDB text in analysis JSON for graphical preview")
    from dvbfixer.batch import add_runtime_help

    add_runtime_help(parser)
    return parser.parse_args(argv)


def _write_new(path: Path, content: str) -> tuple[int, int]:
    """Publish no-replace from a private temporary file on the same filesystem."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=".altlocs-", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as stream:
            stream.write(content)
        identity = Path(temporary).stat()
        os.link(temporary, path)
        return identity.st_dev, identity.st_ino
    finally:
        Path(temporary).unlink(missing_ok=True)


def _unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise AltlocError(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    source = Path(args.input).resolve()
    destinations = [Path(value).resolve() for value in (args.output, args.report_json) if value]
    if len(set([source, *destinations])) != 1 + len(destinations):
        raise SystemExit("ERROR: input, output and report must be different files")
    if any(path.exists() or path.is_symlink() for path in destinations):
        raise SystemExit("ERROR: output/report already exists; choose new paths")
    created: list[tuple[Path, tuple[int, int], str]] = []
    try:
        for destination in destinations:
            destination.parent.mkdir(parents=True, exist_ok=True)
        with source.open("rb") as stream:
            source_bytes = stream.read(50 * 1024 * 1024 + 1)
        if len(source_bytes) > 50 * 1024 * 1024:
            raise AltlocError("input exceeds 50 MiB")
        digest = hashlib.sha256(source_bytes).hexdigest()
        with tempfile.TemporaryDirectory(prefix=".altlocs-", dir=None if not destinations else destinations[0].parent) as work:
            from dvbfixer.structure_input import is_cif_path, normalize_structure

            snapshot = Path(work) / ("source" + source.suffix)
            snapshot.write_bytes(source_bytes)
            normalized = normalize_structure(snapshot, Path(work) / "normalized.pdb") if is_cif_path(source) else None
            text = normalized.pdb_path.read_text(encoding="utf-8") if normalized else source_bytes.decode("utf-8")
            overrides = {}
            if args.selection_plan:
                with Path(args.selection_plan).open("rb") as stream:
                    plan_bytes = stream.read(2 * 1024 * 1024 + 1)
                if len(plan_bytes) > 2 * 1024 * 1024:
                    raise AltlocError("selection plan exceeds 2 MiB")
                plan = json.loads(plan_bytes, object_pairs_hook=_unique_object)
                if not isinstance(plan, dict) or set(plan) != {"schemaVersion", "sourceSha256", "selections"} or type(plan["schemaVersion"]) is not int or plan["schemaVersion"] != 1:
                    raise AltlocError("invalid selection plan schema")
                if plan["sourceSha256"] != digest:
                    raise AltlocError("selection plan source digest mismatch; analyze the current file again")
                overrides = plan["selections"]
                if not isinstance(overrides, dict) or any(not isinstance(k, str) or not isinstance(v, str) for k, v in overrides.items()):
                    raise AltlocError("plan selections must map residue keys to one label")
            # Resolve shorthand against the same Python inventory as GUI.
            _, inventory = select_altlocs(text, source_sha256=digest, policy="manual")
            for value in args.select:
                fields = value.split(":")
                if len(fields) != 3:
                    raise AltlocError("--select must be CHAIN:RESID:LABEL")
                chain, number, label = fields
                chain = " " if chain == "_" else chain
                matches = [row for row in inventory["residues"] if row["chain"] == chain
                           and row["residueNumber"] + row["insertionCode"] == number]
                if len(matches) != 1:
                    raise AltlocError("--select is missing or ambiguous; use an exact JSON selection plan")
                key = matches[0]["key"]
                if key in overrides:
                    raise AltlocError("duplicate selection override")
                overrides[key] = label
            output, report = select_altlocs(
                text, source_sha256=digest, policy=args.policy, prefer=args.prefer,
                fallback=args.fallback, overrides=overrides,
            )
            report["chainMap"] = normalized.chain_map if normalized else {}
            from dvbfixer import __version__

            report["toolVersion"] = __version__
            report["outputWritten"] = False
            if args.include_structure:
                report["normalizedPdb"] = text
            if args.output and not args.dry_run and output is None:
                raise AltlocError("unresolved states; select exactly one state per alternate residue")
            if args.output and not args.dry_run:
                target = Path(args.output)
                identity = _write_new(target, output or "")
                created.append((target, identity, hashlib.sha256((output or "").encode()).hexdigest()))
                report["outputWritten"] = True
            if args.report_json:
                _write_new(Path(args.report_json), json.dumps(report, indent=2, ensure_ascii=False) + "\n")
            for warning in report["warnings"]:
                print(f"WARNING: {warning}")
            print(f"Altlocs: {len(report['residues'])} residues; {len(report['unresolved'])} unresolved; {report['removedAtoms']} atoms removed")
            if not args.report_json:
                print(json.dumps(report, indent=2, ensure_ascii=False))
    except (AltlocError, OSError, ValueError) as exc:
        # A report publication failure must not leave an apparently successful
        # output. Remove only our own still-unmodified file, never a replacement.
        for target, identity, written_digest in created:
            try:
                stat = target.lstat()
                if not target.is_symlink() and (stat.st_dev, stat.st_ino) == identity and hashlib.sha256(target.read_bytes()).hexdigest() == written_digest:
                    target.unlink()
            except OSError:
                pass
        raise SystemExit(f"ERROR: {exc}") from exc


if __name__ == "__main__":
    main()
