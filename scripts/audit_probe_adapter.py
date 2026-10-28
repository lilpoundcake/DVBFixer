"""Validate a fail-closed Probe 2.18 research adapter, without changing production.

The adapter resolves every emitted atom against the exact PDB, preserves contact
classes, and reports requested/actual engine and fallback status even when clean.
Passing these controls is not comparative precision/recall qualification.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path


def atom_inventory(text):
    atoms = {}
    if sum(line.startswith("MODEL ") for line in text.splitlines()) > 1:
        raise ValueError("Multiple models need an explicit model-selection policy")
    for line in text.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        key = (line[21:22].strip(), line[22:26].strip(), line[26:27].strip(),
               line[12:16].strip(), line[16:17].strip())
        if key in atoms:
            raise ValueError(f"Ambiguous PDB atom identity: {key}")
        coordinates = [float(line[start:start + 8]) for start in [30, 38, 46]]
        if not all(math.isfinite(x) for x in coordinates):
            raise ValueError("Nonfinite PDB coordinates")
        atoms[key] = {"resname": line[17:20].strip(), "serial": line[6:11].strip(),
                      "coordinates": coordinates}
    if not atoms:
        raise ValueError("No atoms")
    return atoms


def normalize_fields(text):
    lines = []
    for line in text.splitlines(keepends=True):
        if line.startswith(("ATOM  ", "HETATM")):
            name, element = line[12:16].strip(), line[76:78].strip()
            if len(element) == 1 and len(name) < 4:
                line = line[:12] + f" {name:<3}" + line[16:]
        lines.append(line)
    result = "".join(lines)
    if atom_inventory(text) != atom_inventory(result):
        raise ValueError("Formatting changed atom identity/coordinates")
    return result


def parse_contacts(raw, inventory):
    """Strict supported-schema parser; never strip fixed-width atom fields."""
    # Probe 2.18 readPDBrecs.c uppercases atom/residue names, but preserves
    # chains, insertion codes and altloc. Resolve its lossy spelling only
    # through an injective mapping to original input identities.
    aliases = probe_aliases(inventory)
    contacts, overlaps = [], {}
    for line in raw.splitlines():
        if not line.strip():
            continue
        fields = line.split(":")
        if len(fields) != 20 or fields[0] or fields[1] != "1->1" or fields[2] not in {"wc", "cc", "so", "bo", "hb"}:
            raise ValueError("Unsupported Probe condensed record")
        if int(fields[5]) <= 0:
            raise ValueError("Invalid dot count")
        values = [float(fields[i]) for i in [*range(6, 13), *range(15, 20)]]
        if not all(math.isfinite(v) for v in values):
            raise ValueError("Nonfinite Probe record")
        identities = []
        for field in fields[3:5]:
            if len(field) != 16:
                raise ValueError("Unsupported Probe atom-field width")
            key = (field[:2].strip(), field[2:6].strip(), field[6:7].strip(),
                   field[11:15].strip(), field[15:16].strip())
            original = aliases.get(key)
            if original is None or inventory[original]["resname"].upper() != field[7:10].strip():
                raise ValueError(f"Unresolved Probe atom identity: {field!r}")
            identities.append(original)
        pair = tuple(sorted(identities))
        contact = {"atoms": pair, "type": fields[2], "dot_count": int(fields[5]), "minimum_gap_angstrom": values[0]}
        contacts.append(contact)
        if fields[2] in {"bo", "so"} and -values[0] >= .4 - 1e-12:
            entry = overlaps.setdefault(pair, {"atoms": pair, "overlap_angstrom": 0., "contact_types": set()})
            entry["overlap_angstrom"] = max(entry["overlap_angstrom"], -values[0])
    for contact in contacts:
        if contact["atoms"] in overlaps:
            overlaps[contact["atoms"]]["contact_types"].add(contact["type"])
    return contacts, [{**entry, "contact_types": sorted(entry["contact_types"])} for entry in overlaps.values()]


def probe_aliases(inventory):
    aliases = {}
    for key in inventory:
        chain, resid, icode, atom, altloc = key
        emitted = (chain, resid, icode, atom.upper(), altloc)
        if emitted in aliases:
            raise ValueError("Probe uppercasing makes atom identities ambiguous")
        aliases[emitted] = key
    return aliases


def run_adapter(binary, path, work, nuclear=False, timeout=60):
    report = {"requested_engine": "molprobity-probe", "actual_engine": None,
              "status": "not_started", "findings": [], "radii": "nuclear" if nuclear else "electron"}
    try:
        inventory = atom_inventory(path.read_text())
        probe_aliases(inventory)
        version = subprocess.run([str(binary), "-version"], capture_output=True, timeout=timeout)
        (work / "version.stdout").write_bytes(version.stdout)
        (work / "version.stderr").write_bytes(version.stderr)
        if version.returncode or version.stdout.strip() != b"probe.2.18.211005":
            raise ValueError("Unsupported or failed Probe version preflight")
        command = [str(binary), *(["-nuclear"] if nuclear else []), "-unformated", "-condensed", "-Q",
                   "-mc", "-het", "-nowaters", "-self", "ALL", str(path)]
        report["command"] = command
        run = subprocess.run(command, capture_output=True, timeout=timeout)
        (work / "stdout.raw").write_bytes(run.stdout)
        (work / "stderr.raw").write_bytes(run.stderr)
        report["returncode"] = run.returncode
        if run.returncode:
            raise ValueError("Probe returned nonzero")
        if b"ERROR" in run.stderr.upper():
            raise ValueError("Probe emitted an error")
        contacts, findings = parse_contacts(run.stdout.decode("utf-8", errors="strict"), inventory)
        report.update({"actual_engine": "molprobity-probe", "status": "complete", "findings": findings,
                       "contact_count": len(contacts), "ambiguous_hb_overlap_pairs": sum("hb" in f["contact_types"] for f in findings)})
    except (OSError, ValueError, UnicodeError, subprocess.TimeoutExpired) as exc:
        if isinstance(exc, subprocess.TimeoutExpired):
            (work / "stdout.raw").write_bytes(exc.stdout or b"")
            (work / "stderr.raw").write_bytes(exc.stderr or b"")
        report.update({"status": "failed", "failure": type(exc).__name__ + ": " + str(exc)})
    return report


def fallback(report, path):
    from openmm.app import PDBFile

    from dvbfixer.diagnose.steric import clashes_python

    if report["status"] != "complete":
        pdb = PDBFile(str(path))
        findings = clashes_python(pdb.topology, pdb.positions, clash_warn_a=.4, clash_error_a=.5)
        report.update({"actual_engine": "python", "fallback": True, "fallback_findings_count": len(findings)})
    else:
        report["fallback"] = False
    return report


def identity_control(source):
    """Relabel the established methane clash; geometry/bonds are untouched."""
    text = source.read_text()
    lines = []
    for line in text.splitlines(keepends=True):
        if line.startswith(("ATOM  ", "HETATM")):
            upper = line[21] == "A"
            chain, number, icode = ("H", 82, "A") if upper else ("h", 82, "B")
            line = line[:21] + chain + f"{number:4d}" + icode + line[27:]
        lines.append(line)
    return "".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True, help="Existing linux-fresh audit")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--input", type=Path, action="append", default=[])
    args = parser.parse_args()
    args.probe, args.audit, args.output = args.probe.resolve(), args.audit.resolve(), args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "audit-script.py").write_bytes(Path(__file__).read_bytes())
    prior = json.loads((args.audit / "probe-electron/summary.json").read_text())
    paths = [Path(r["source_input"]) for r in prior["runs"]] + [p.resolve() for p in args.input]
    synthetic = args.output / "case-sensitive-insertion-codes.pdb"
    synthetic.write_text(identity_control(next(p for p in paths if p.stem == "clashing-methanes")))
    blank = args.output / "blank-chain-negative-resid.pdb"
    text = synthetic.read_text()
    blank.write_text("".join(line[:21] + " " + f"{-1:4d}" + "A" + line[27:]
                            if line.startswith(("ATOM  ", "HETATM")) and line[21] == "H" else line
                            for line in text.splitlines(keepends=True)))
    paths += [synthetic, blank]
    lower_atom = args.output / "lowercase-atom-name.pdb"
    lower_atom.write_text("".join(line[:12] + " h1 " + line[16:]
                                 if line.startswith(("ATOM  ", "HETATM")) and line[12:16].strip() == "H1" else line
                                 for line in synthetic.read_text().splitlines(keepends=True)))
    paths.append(lower_atom)
    hinge = next(p for p in paths if p.stem == "hinge-openmm")
    misaligned = args.output / "misaligned-hinge.pdb"
    misaligned.write_text("".join(line[:12] + f"{line[12:16].strip():<4}" + line[16:]
                                 if line.startswith(("ATOM  ", "HETATM")) else line
                                 for line in hinge.read_text().splitlines(keepends=True)))
    assert atom_inventory(hinge.read_text()) == atom_inventory(misaligned.read_text())
    paths.append(misaligned)
    report = {"binary_sha256": hashlib.sha256(args.probe.read_bytes()).hexdigest(), "runs": [], "failures": []}
    for index, original in enumerate(paths):
        for formatted in [False, True]:
            for nuclear in [False, True]:
                work = args.output / f"{index}-{original.stem}-aligned{int(formatted)}-nuclear{int(nuclear)}"
                work.mkdir()
                path = work / "input.pdb"
                source = original.read_text()
                path.write_text(normalize_fields(source) if formatted else source)
                row = fallback(run_adapter(args.probe, path, work, nuclear), path)
                row.update({"source": str(original), "source_sha256": hashlib.sha256(original.read_bytes()).hexdigest(),
                            "aligned": formatted, "nuclear": nuclear})
                report["runs"].append(row)
    control = next(p for p in paths if p.stem == "clashing-methanes")
    for failure in ["missing", "timeout", "nonzero", "malformed", "nonfinite", "unsupported-version", "mixed-malformed"]:
        work = args.output / f"failure-{failure}"
        work.mkdir()
        binary = work / "probe"
        if failure != "missing":
            good_raw = next(args.output.glob("1-*-aligned1-nuclear0/stdout.raw")).read_text()
            body = "import sys,time\n"
            body += "if '-version' in sys.argv:\n    print('probe.9.99' if " + repr(failure == "unsupported-version") + " else 'probe.2.18.211005')\n    raise SystemExit(0)\n"
            if failure == "timeout":
                body += "time.sleep(10)\n"
            elif failure == "nonzero":
                body += "raise SystemExit(1)\n"
            elif failure == "nonfinite":
                fields = good_raw.splitlines()[0].split(":")
                fields[6] = "nan"
                body += f"print({':'.join(fields)!r})\n"
            else:
                body += f"print({(good_raw + 'unsupported' if failure == 'mixed-malformed' else 'unsupported')!r})\n"
            binary.write_text(f"#!{sys.executable}\n" + body)
            binary.chmod(0o700)
        row = fallback(run_adapter(binary, control, work, timeout=.2), control)
        row["case"] = failure
        assert row["status"] == "failed" and row["actual_engine"] == "python" and row["fallback_findings_count"] == 9
        report["failures"].append(row)
    collision = args.output / "ambiguous-atom-case.pdb"
    collision.write_text("".join(line[:12] + " h1 " + line[16:]
                                if line.startswith(("ATOM  ", "HETATM")) and line[12:16].strip() == "H2" else line
                                for line in synthetic.read_text().splitlines(keepends=True)))
    work = args.output / "failure-ambiguous-atom-case"
    work.mkdir()
    row = fallback(run_adapter(args.probe, collision, work), collision)
    row["case"] = "ambiguous-atom-case"
    assert row["status"] == "failed" and "ambiguous" in row["failure"] and row["actual_engine"] == "python"
    report["failures"].append(row)
    report["normalized_real_runs_complete"] = all(r["status"] == "complete" for r in report["runs"] if r["aligned"])
    report["format_control"] = []
    for nuclear in [False, True]:
        baseline = next(r for r in report["runs"] if r["source"] == str(hinge) and r["aligned"] and r["nuclear"] == nuclear)
        corrected = next(r for r in report["runs"] if r["source"] == str(misaligned) and r["aligned"] and r["nuclear"] == nuclear)
        raw = next(r for r in report["runs"] if r["source"] == str(misaligned) and not r["aligned"] and r["nuclear"] == nuclear)
        assert baseline["findings"] == corrected["findings"]
        report["format_control"].append({"nuclear": nuclear, "coordinates_unchanged": True,
                                         "misaligned_overlap_pairs": len(raw["findings"]),
                                         "normalized_overlap_pairs": len(corrected["findings"]),
                                         "matches_baseline_pairs_and_metrics": True})
    (args.output / "report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    manifest = {str(p.relative_to(args.output)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in sorted(args.output.rglob("*")) if p.is_file()}
    (args.output / "manifest.sha256.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
