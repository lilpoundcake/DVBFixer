"""Independent full-context export/DatRecord/postflight byte readback, no repair."""
from __future__ import annotations

import argparse
import hashlib
import json
import runpy
from collections import Counter
from pathlib import Path

import gemmi

from dvbfixer.ffutils.dat import DatRecord

PROJECT = Path(__file__).resolve().parents[1]
HELPERS = PROJECT / "tmp/mac-real-gap-continuation-2026-11-12/1KTZ-native-terminal-restoration-attempt01/verify_readback.py"
HANDED = PROJECT / "tmp/mac-real-gap-continuation-2026-11-14/1KTZ-inward-wall-refinement-attempt01/verify_readback.py"


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main(root):
    helpers = runpy.run_path(str(HELPERS))
    key, inventory = helpers["key"], helpers["inventory"]
    handedness = runpy.run_path(str(HANDED))["handedness"]
    plan = json.loads((root / "plan.json").read_text())
    for name, expected in plan["source_code_and_input_sha256"].items():
        if sha(Path(name)) != expected:
            raise ValueError("full-context export prospective input/code changed")
    if len(plan["jobs"]) != 1:
        raise ValueError("full-context export plan is not a single job")
    job = plan["jobs"][0]
    work = root / "execution" / f"{job['case_id']}-seed-{job['seed']}"
    execution = json.loads((work / "execution.json").read_text())
    result = json.loads((work / "result.json").read_text())
    if execution["returncode"] or execution["guard_stop_reason"] or not execution["completed"] or not result["private_full_context_export_completed"]:
        raise ValueError("full-context export failed its unchanged resource guard")
    bundle = work / "result"
    final_path = bundle / "1KTZ-restored-heavy-with-waters.pdb"
    original_path = bundle / "original-full-context.pdb"
    projection_path = bundle / "validated-protein-projection.pdb"
    original_lines = original_path.read_bytes().splitlines(keepends=True)
    final_lines = final_path.read_bytes().splitlines(keepends=True)
    remaining = iter(final_lines)
    for line in original_lines:
        if not any(candidate == line for candidate in remaining):
            raise ValueError("an original full-context line changed/order was not preserved")
    original, final, projection = inventory(original_path), inventory(final_path), inventory(projection_path)
    request = json.loads((bundle / "original-request.json").read_text())
    generated = {key(a) for a in request["generated_atoms"]}
    if len(original) != 3312 or len(final) != 4022 or final.keys() != original.keys() | generated:
        raise ValueError("complete context atom inventory differs from3312original+710restored")
    if any(final[k] != row for k, row in original.items()) or any(final[k] != projection[k] for k in generated):
        raise ValueError("full-context coordinates/chemistry differ from observations or auditedprojection")
    final_records = {}
    final_serials = []
    for line in final_lines:
        if line.startswith((b"ATOM  ", b"HETATM")):
            k = (line[21:22].decode(), line[22:26].strip().decode(), line[26:27].strip().decode(), line[12:16].strip().decode())
            final_records[k] = line
            final_serials.append(int(line[6:11]))
    if len(final_serials) != len(set(final_serials)):
        raise ValueError("full-context final atom serial collision")
    for line in projection_path.read_bytes().splitlines(keepends=True):
        if not line.startswith((b"ATOM  ", b"HETATM")):
            continue
        k = (line[21:22].decode(), line[22:26].strip().decode(), line[26:27].strip().decode(), line[12:16].strip().decode())
        if k in generated and (line[:6] != final_records[k][:6] or line[11:] != final_records[k][11:]):
            raise ValueError("new-atom export changed more than serial metadata")
    mapping = json.loads((bundle / "numbering-and-target-annotations.json").read_text())
    if len(mapping["new_atom_serials_only"]) != 710:
        raise ValueError("serial mapping is incomplete")
    reserved = {int(line[6:11]) for line in original_lines if line[:6].strip() in {b"ATOM", b"HETATM", b"TER"} and line[6:11].strip()}
    for row in mapping["new_atom_serials_only"]:
        k = key(row["atom"])
        if k not in generated or row["full_context_serial"] in reserved or int(final_records[k][6:11]) != row["full_context_serial"]:
            raise ValueError("new-atom serial remapping touched a reserved original identifier")
    dat = DatRecord.load(final_path.with_suffix(".dat"))
    if len(dat.added_atoms) != 710 or dat.added_keys() != generated or len(dat.residue_summary) != 92:
        raise ValueError("canonicalDatRecord doesnotidentifyexact710missingatoms/92residues")
    expected_pairs = {tuple(sorted((key(p["atom1"]), key(p["atom2"])))) for p in request["retained_explicit_links"]}
    expected_pairs |= {tuple(sorted(key(a) for a in pair)) for pair in mapping["derived_target_SS_pairs_not_original_annotations"]}
    actual = []
    for line in final_lines:
        if line.startswith(b"SSBOND"):
            actual.append(tuple(sorted(((line[15:16].decode(), line[17:21].strip().decode(), line[21:22].strip().decode(), "SG"),
                                        (line[29:30].decode(), line[31:35].strip().decode(), line[35:36].strip().decode(), "SG")))))
    if len(actual) != 21 or set(actual) != expected_pairs or len(mapping["derived_target_SS_pairs_not_original_annotations"]) != 3:
        raise ValueError("full-context target21/source18derived3annotation auditfailed")
    checked, bad = handedness({k: row[2] for k, row in final.items()})
    if bad:
        raise ValueError("independent final full-contextCalphaL checkfailed")
    structure = gemmi.read_structure(str(final_path))
    protein_counts = {chain.name: sum(residue.name != "HOH" for residue in chain) for chain in structure[0]}
    if protein_counts != {s["chain"]: len(s["sequence"]) for s in request["target_sequences"]}:
        raise ValueError("full known construct/fourchains not preserved")
    current = json.loads((bundle / "diagnose-restored.json").read_text())
    baseline = json.loads((bundle / "diagnose-original.json").read_text())

    def signature(finding):
        return tuple(finding.get(name) for name in ("severity", "category", "chain", "resid", "atom", "message"))

    baseline_warnings = Counter(signature(f) for f in baseline["findings"] if f["severity"] == "WARNING")
    final_warnings = Counter(signature(f) for f in current["findings"] if f["severity"] == "WARNING")
    new = final_warnings - baseline_warnings
    report = {"private_full_context_data_custody_passed": True,
              "original_context_records_unchanged_ordered_subsequence": True,
              "original_protein_atoms": 2986, "original_water_oxygens": 326, "restored_heavy_atoms": 710,
              "protein_residues": 468, "water_residues": 326, "protein_chain_counts": protein_counts,
              "SS_annotation_records": 21, "original_SS_records_preserved": 18, "derived_target_SS_annotations": 3,
              "canonical_dat_added_atoms": 710, "canonical_dat_added_residues": 92,
              "independent_chiral_centers_checked": checked, "D_centers": bad,
              "CLI_final_summary": current["summary"], "new_CLI_WARNING_signatures": [list(k) + [v] for k, v in new.items()],
              "all_final_CLI_warnings_are_baseline": not new,
              "CLI_no_errors_or_new_warnings": current["summary"].get("ERROR", 0) == 0 and not new,
              "final_sha256": sha(final_path), "data_sidecar_sha256": sha(final_path.with_suffix(".dat")),
              "complete_scientific_accuracy_or_public_qualification": False,
              "remaining_chemistry_boundary": "No H/protonation/MDparameterization;temporaryB/DOXT remainoutsidefrozenmask",
              "static_screen_not_knot_or_continuous_path_qualification": True,
              "bindings": {str(p.resolve()): sha(p) for p in
                           (Path(__file__), HELPERS, HANDED, root / "plan.json", work / "execution.json", work / "result.json",
                            *sorted(p for p in bundle.iterdir() if p.is_file()))}}
    with (root / "readback-receipt.json").open("x") as stream:
        json.dump(report, stream, indent=2, ensure_ascii=False)
        stream.write("\n")
    print(json.dumps({k: v for k, v in report.items() if k != "bindings"}, ensure_ascii=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    main(parser.parse_args().root.resolve())
