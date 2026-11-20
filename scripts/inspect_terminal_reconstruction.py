"""Source-only metrics and matched-view plots for an already audited terminal run.

No generation, energy evaluation, coordinate fitting, repair or reference input.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
from scripts.audit_retained_full_context_controls import atoms
from scripts.hybrid_gap_reconstruction import key
from scripts.mac_reconstruction_comparison import digest, save, target_axis
from scripts.terminal_reconstruction_research import terminal_torsions

from dvbfixer.model.diffusion.contract import DiffusionRequest


def inspect(root):
    receipt = json.loads((root / "independent-readback.json").read_text())
    if receipt["plan_sha256"] != digest(root / "plan.json"):
        raise ValueError("requires an independent readback of this run")
    request = DiffusionRequest.from_json((root / "request.json").read_text())
    generated = list(map(key, request.generated_atoms))
    order = target_axis(request)[request.target_sequences[0].chain]
    length = len(request.gaps[0].generated_residues)
    report = {"new_diffusion_calls": 0, "new_minimization_calls": 0, "reference_coordinates_read": False,
              "plan_sha256": digest(root / "plan.json"), "arms": {}}
    for arm in ("idpforge", "modeller"):
        folder = root / arm
        packed, final = atoms(folder / "packed.pdb"), atoms(folder / "minimized.pdb")
        delta = np.array([np.array(final[a]["xyz"]) - packed[a]["xyz"] for a in generated])
        norms = np.linalg.norm(delta, axis=1)
        q0, p0 = terminal_torsions(order, packed, length)
        q1, p1 = terminal_torsions(order, final, length)
        difference = (np.array(q1)-q0+math.pi) % (2*math.pi)-math.pi
        findings = json.loads((folder / "minimized.diagnose.json").read_text())["findings"]
        warnings = [f for f in findings if f["severity"] == "WARNING"]
        report["arms"][arm] = {"generated_heavy_RMS_from_packed_A": float(np.sqrt(np.mean(norms**2))),
            "generated_heavy_max_from_packed_A": float(norms.max()), "phi_psi_change_degrees": np.degrees(difference).tolist(),
            "anchor_phi_change_degrees": math.degrees((p1-p0+math.pi)%(2*math.pi)-math.pi),
            "warning_categories": dict(Counter(f["category"] for f in warnings)), "warnings": warnings,
            "candidate_sha256": digest(folder / "minimized.pdb")}
    report["interpretation"] = "both error-level geometry gates pass; IDPForge-derived cis non-PRO and large proposal displacement remain visible; no superiority claim"
    save(root / "comparison-inspection.json", report)
    return report


def plot(root):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from openmm.app import PDBFile

    source = atoms(root / "source.pdb")
    panels = [("MODELLER raw (source fixed)", root / "modeller-source-restored-baseline.pdb", "#28704f"),
              ("MODELLER + common processing", root / "modeller/minimized.pdb", "#326ca8"),
              ("IDPForge torsions + common processing", root / "idpforge/minimized.pdb", "#a44d2d")]
    points = [atoms(path) for _, path, _ in panels]
    variable = [a for a in points[0] if a not in source]
    region = np.array([r["xyz"] for rows in points for a, r in rows.items() if a in variable])
    anchor = np.array(source["A", "5", "", "CA"]["xyz"])
    region = np.vstack((region, anchor))
    center = (region.min(axis=0) + region.max(axis=0))/2
    radius = float(np.ptp(region, axis=0).max())/2 + 3
    source_ca = np.array([r["xyz"] for a, r in source.items() if a[3] == "CA"])
    fig = plt.figure(figsize=(16, 10))
    for row, (elev, azim) in enumerate(((20, -65), (65, 35))):
        for column, ((title, path, color), coords) in enumerate(zip(panels, points, strict=True)):
            ax = fig.add_subplot(2, 3, row*3+column+1, projection="3d")
            ax.plot(*source_ca.T, color="#858585", alpha=.4, linewidth=1)
            structure = PDBFile(str(path))
            def identity(a):
                return a.residue.chain.id, str(a.residue.id), a.residue.insertionCode.strip(), a.name
            for a, b in structure.topology.bonds():
                ka, kb = identity(a), identity(b)
                if ka not in variable and kb not in variable:
                    continue
                segment = np.array([coords[ka]["xyz"], coords[kb]["xyz"]])
                ax.plot(*segment.T, color=color, linewidth=1.5)
            for a in variable:
                p = coords[a]["xyz"]
                element = coords[a]["element"]
                ax.scatter(*p, color={"N": "#3155ad", "O": "#b83333", "S": "#bc9a22"}.get(element, color), s=15)
                if a[3] == "CA":
                    ax.text(*p, f" {a[1]}", fontsize=9)
            ax.scatter(*anchor, color="black", s=35, marker="s")
            ax.text(*anchor, " A5 fixed", fontsize=9)
            if column == 2:
                p = np.array([coords["A", "1", "", "C"]["xyz"], coords["A", "2", "", "N"]["xyz"]])
                ax.plot(*p.T, color="red", linewidth=4)
                ax.text(*p.mean(axis=0), " cis M1–Q2", color="red", fontsize=9)
            for setter, c in zip((ax.set_xlim, ax.set_ylim, ax.set_zlim), center, strict=True):
                setter(c-radius, c+radius)
            ax.set_box_aspect((1, 1, 1))
            ax.view_init(elev=elev, azim=azim)
            ax.set_axis_off()
            if row == 0:
                ax.set_title(title, fontsize=11)
    fig.suptitle("1UBQ A1–4: identical observed scaffold; matched views, no reference structure", fontsize=14)
    fig.text(.5, .03, "Thin gray: fixed observed Cα context. Colored sticks: new heavy atoms. Red bond: retained cis non-PRO warning.\n"
             "All three have zero geometry ERROR findings; the pre-existing C-terminal OXT is missing. Not a model-accuracy ranking.",
             ha="center", fontsize=10)
    fig.subplots_adjust(left=.01, right=.99, top=.91, bottom=.09, wspace=.01, hspace=.01)
    fig.savefig(root / "matched-view-comparison.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    print(json.dumps(inspect(args.root.resolve())))
    plot(args.root.resolve())
