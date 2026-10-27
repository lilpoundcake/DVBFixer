"""Export existing 9CQD evidence for visual review; do not run new modeling."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from collections import Counter
from dataclasses import asdict
from pathlib import Path

import numpy as np

from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest, ResidueIdentity


def source_motif(text: str, request: DiffusionRequest) -> str:
    """Selected fixed atoms, original IDs/serials/B; TER at omitted peptide spans."""
    if request.protein_context is None:
        raise ValueError("missing explicit motif context")
    selected = {item.atom for item in request.protein_context.fixed_conditioning}
    positions = {r: i for p in request.sequence_placements
                 for i, r in zip(p.observed_target_indices, p.observed_residues, strict=True)}
    lines = ["REMARK 999 VISUAL REVIEW: SOURCE-IDENTITY FIXED CONDITIONING ONLY\n",
             "REMARK 999 GENERATED GAP ATOMS ARE ABSENT; OMITTED SPANS SEPARATED BY TER\n"]
    seen = set()
    previous = None
    for line in text.splitlines(keepends=True):
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atom = AtomIdentity(line[21], line[22:26].strip(), line[26].strip(), line[12:16].strip())
        if atom not in selected:
            continue
        if atom in seen:
            raise ValueError("duplicate source motif atom")
        r = ResidueIdentity(atom.chain, atom.residue_number, atom.insertion_code)
        if previous is not None and r != previous:
            if r.chain != previous.chain or positions[r] != positions[previous] + 1:
                lines.append("TER\n")
        lines.append(line if line.endswith("\n") else line + "\n")
        seen.add(atom)
        previous = r
    if seen != selected:
        raise ValueError("motif does not match selected conditioning atom identities")
    lines.extend(["TER\n", "END\n"])
    return "".join(lines)


def coordinate_records(text: str) -> Counter:
    """Compare motif content independently of generator chain/number relabeling."""
    return Counter((l[12:16], l[17:20], l[30:66]) for l in text.splitlines()
                   if l.startswith(("ATOM  ", "HETATM")))


def main() -> None:
    from scripts.run_cc89_lrd_ablation import atom_records

    from dvbfixer.align import _apply_transform_preserving_headers
    from dvbfixer.biological_assembly import _transform_anisou
    from dvbfixer.model.diffusion.geometry import weighted_kabsch
    from dvbfixer.model.diffusion.protein_context import independent_gap_requests

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--experiment", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    experiment, baseline, output = args.experiment, args.baseline, args.output
    report = json.loads((experiment / "comparison.json").read_text())
    completed = json.loads((experiment / "completed.json").read_text())
    if completed != {"rows": 12, "errors": 0} or len(report) != 12 or any("error" in row for row in report):
        raise ValueError("requires all twelve completed evidence rows")
    request = DiffusionRequest.from_json((experiment / "request.json").read_text())
    source_text = (experiment / "source.pdb").read_text()
    source = atom_records(source_text)
    fit = tuple(a for a in request.fixed_atoms if a.chain == "W" and a.atom_name in {"N", "CA", "C", "O"})
    if len(fit) != 844:
        raise ValueError("unexpected fitting population")
    output.mkdir(parents=True, exist_ok=False)
    for subdir in ("original-frame", "kabsch-W", "motifs", "evidence"):
        (output / subdir).mkdir()
    shutil.copy2(experiment / "source.pdb", output / "source.pdb")
    for name in ("request.json", "manifest.json", "comparison.json", "completed.json"):
        shutil.copy2(experiment / name, output / "evidence" / name)
    catalog = []
    transforms = {}

    def export_structure(path, name, kind, metrics=None):
        data = path.read_bytes()
        sha = hashlib.sha256(data).hexdigest()
        if metrics is not None and sha != metrics["candidate_sha256"]:
            raise ValueError("candidate differs from measured evidence")
        shutil.copy2(path, output / "original-frame" / name)
        records = atom_records(data.decode())
        transform = weighted_kabsch(np.array([records[a][0] for a in fit]), np.array([source[a][0] for a in fit]))
        destination = output / "kabsch-W" / name
        _apply_transform_preserving_headers(path, destination, transform.rotation, transform.translation)
        if any(l.startswith("ANISOU") for l in data.decode().splitlines()):
            destination.write_text("".join(_transform_anisou(l, transform.rotation) if l.startswith("ANISOU") else l
                                           for l in destination.read_text().splitlines(keepends=True)))
        saved = atom_records(destination.read_text())
        if set(saved) != set(records) or not all(np.allclose(saved[a][0], transform.apply(np.array([records[a][0]]))[0],
                                                             atol=0.000501, rtol=0) for a in records):
            raise ValueError("aligned copy did not preserve transformed coordinates")
        transforms[name] = {"rotation": transform.rotation.tolist(), "translation": transform.translation.tolist(),
                            "fit_chain": "W", "fit_atom_count": len(fit), "one_transform_for_whole_complex": True}
        catalog.append({"file": name, "kind": kind, "original_sha256": sha,
                        "geometry_pass": metrics["geometry_only_pass_not_acceptance"] if metrics else "not assessed here",
                        "strict_pass": metrics["strict_validation"]["passed"] if metrics else "not assessed here",
                        "anchor_backbone_rmsd_W_aligned_A": metrics["combined_anchor_backbone_rmsd_W_aligned_A"] if metrics else ""})

    for arm in ("reinjection-off", "reinjection-on"):
        export_structure(baseline / "initial" / arm / "candidate.pdb", f"{arm}__raw.pdb", "raw")
    for row in report:
        name = f"{row['arm']}__{row['policy']}__k{row['mean_k_kJ_mol_nm2']}.pdb"
        export_structure(experiment / row["candidate"], name, "minimized", row)

    motif_index = []
    for i, single in enumerate(independent_gap_requests(request)):
        gap = single.gaps[0]
        name = f"gap-{i}__{gap.chain}{gap.generated_residues[0].residue_number}-{gap.generated_residues[-1].residue_number}"
        motif = source_motif(source_text, single)
        (output / "motifs" / f"{name}__source-ids.pdb").write_text(motif)
        for arm in ("reinjection-off", "reinjection-on"):
            archived = baseline / "initial" / arm / f"gap-{i}" / "motif-1.pdb"
            if coordinate_records(motif) != coordinate_records(archived.read_text()):
                raise ValueError("source-identity motif differs from archived engine motif content")
            shutil.copy2(archived, output / "motifs" / f"{name}__{arm}__engine.pdb")
        motif_index.append({"gap": asdict(gap), "conditioning_atom_count": len(single.protein_context.fixed_conditioning),
                            "context_spans": [asdict(s) for s in single.protein_context.spans]})
    (output / "evidence" / "motifs.json").write_text(json.dumps(motif_index, indent=2) + "\n")
    (output / "evidence" / "transforms.json").write_text(json.dumps(transforms, indent=2) + "\n")
    with (output / "structures.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(catalog[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(catalog)
    (output / "README.md").write_text("""# 9CQD restrained-anchor visual review

- `source.pdb`: original experimental input, common reference for overlays.
- `original-frame/`: two raw predictions and all twelve minimization outputs, unchanged.
- `kabsch-W/`: the same fourteen complexes fitted to source using 844 observed W-chain N/CA/C/O atoms; includes W anchors, excludes generated atoms. One transform for ALL chains; X is not separately fitted. Coordinates rounded only when writing PDB.
- `motifs/*source-ids.pdb`: actual fixed conditioning atom subsets with source chain/residue IDs and B factors; TER separates omitted peptide spans. Not full-complex predictions. No stale full-structure sequence/link metadata is copied.
- `motifs/*engine.pdb`: byte-identical archived generator motif inputs, with generator proxy chain IDs and local numbering. Their atom/coordinate/B/occupancy content is checked against the source-ID motif, ignoring relabeling.
- `structures.tsv`: output labels and geometry/strict status. `k` units: kJ/mol/nm².
- `evidence/`: frozen request/manifest, full metrics, transforms and context spans.
- `SHA256SUMS`: checksums of every other exported file.

View source plus matching kabsch-W outputs for local comparison. Anchors: W25/W29 and X139/X142. Generated gaps: W26–28 and X140–141. For literal source preservation inspect original-frame instead. B factors on source atoms are inherited observations, not recalculated uncertainty of the relaxed output; values on generated atoms are not experimental B factors. No structure here passes strict production acceptance. Geometry passes are NOT evidence of agreement with experimental density. Do not treat viewer-inferred bonds across missing spans as experimental connectivity.
""")
    files = sorted(p for p in output.rglob("*") if p.is_file())
    (output / "SHA256SUMS").write_text("".join(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(output)}\n" for p in files))
    print(f"Exported {len(list(output.rglob('*.pdb')))} PDBs to {output}")


if __name__ == "__main__":
    main()
