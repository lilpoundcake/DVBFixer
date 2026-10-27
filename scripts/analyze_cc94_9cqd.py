"""Score every cc94/cc89 ODE attempt and export structures for visual review."""

from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

import numpy as np
from scripts.cc94_9cqd_metrics import score
from scripts.cc94_9cqd_study import CELLS, ROOT, SEEDS, WORK, check_protocol, digest, write_json

from dvbfixer.model.diffusion.benchmark import load_pdb_coordinates
from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest, RunnerResult

LABELS = {
    "N0": "cc94_ode_no_reinjection",
    "N1": "cc94_ode_reinjection",
    "B0": "cc89_ode_no_reinjection",
    "B1": "cc89_ode_reinjection",
}


def rows():
    output = []
    for directory in sorted((WORK / "runs").glob("*/*/*")):
        cell, seed, attempt = directory.parts[-3:]
        request = DiffusionRequest.from_json((directory / "request.json").read_text())
        generated = set(request.generated_atoms)
        native = np.load(directory / "native-coordinates.npz")["coordinates"]
        axis = [
            AtomIdentity(**a)
            for a in json.loads((directory / "native-axis.json").read_text())["atoms"]
        ]
        gi = [i for i, a in enumerate(axis) if a in generated]
        xi = [i for i, a in enumerate(axis) if a.chain == "X" and a in set(request.fixed_atoms)]
        distances = np.linalg.norm(native[gi, None, :] - native[None, xi, :], axis=-1)
        for stage in ["raw", "refined"]:
            result_path = directory / f"{stage}-result.json"
            if not result_path.exists():
                continue
            result = RunnerResult.from_json(result_path.read_text())
            (candidate,) = result.candidates
            pdb_path = directory / candidate.coordinate_artifact.path
            assert digest(pdb_path) == candidate.coordinate_artifact.sha256
            geometry = score(pdb_path.read_text(), (directory / "source.pdb").read_text(), request)
            validation = json.loads((directory / f"{stage}-validation.json").read_text())[0][
                "summary"
            ]
            values = {m["name"]: m["value"] for m in validation["metrics"]}
            output.append(
                {
                    "cell": cell,
                    "model": CELLS[cell][0],
                    "reinjection": CELLS[cell][1],
                    "seed": int(seed[5:]),
                    "attempt": attempt,
                    "stage": stage,
                    "pdb": str(pdb_path.relative_to(WORK)),
                    "pdb_sha256": digest(pdb_path),
                    "native_MX_min_A": float(distances.min()),
                    "native_MX_pairs_below_1_9_A": int((distances < 1.9).sum()),
                    "min_MX_A": geometry["X"]["minimum_distance_angstrom"],
                    "MX_pairs_below_1_9_A": geometry["X"]["pairs_below_1_9_angstrom"],
                    "segment_min_A": geometry["X"]["backbone_segment_minimum_angstrom"],
                    "segment_alarm": geometry["X"]["near_crossing_alarm"],
                    "authoritative_severe_clashes": geometry["authoritative_severe_clash_count"],
                    "authoritative_severe_clash_pairs": geometry["authoritative_severe_clashes"],
                    "repository_severe_clashes": values.get("severe-steric-overlaps"),
                    "junction_CN_A": values.get("junction-cn-distance-max"),
                    "fixed_max_displacement_A": geometry["fixed_heavy_max_displacement_angstrom"],
                    "gate_failures": validation["hard_gate_failures"],
                    "validation_metrics": values,
                    "all_gates_pass": validation["passed"]
                    and geometry["authoritative_severe_clash_count"] == 0
                    and not candidate.postprocessing_failures,
                    "postprocessing_failures": list(candidate.postprocessing_failures),
                    "resource_metrics": {
                        "peak_ram_bytes": result.resource_metrics.peak_ram_bytes,
                        "peak_vram_bytes": result.resource_metrics.peak_vram_bytes,
                        "wall_time_seconds": result.resource_metrics.wall_time_seconds,
                    },
                }
            )
    return output


def export_visual(records):
    destination = ROOT / "gap_visual_review_9cqd_cc94"
    marker = destination / "study.json"
    if destination.exists() and (
        not marker.exists()
        or json.loads(marker.read_text())["protocol_sha256"] != digest(WORK / "protocol.json")
    ):
        raise ValueError("visual export directory is not owned by this study")
    destination.mkdir(exist_ok=True)
    write_json(
        marker, {"protocol_sha256": digest(WORK / "protocol.json"), "study": "9cqd-cc94-ODE"}
    )
    source = ROOT / "tests/fixtures/protein_context/9cqd_split.pdb"
    shutil.copyfile(source, destination / "reference_9CQD.pdb")
    exports = []
    for record in records:
        repeat = record["attempt"] == "repeat-1"
        relative = Path("repeats") / record["stage"] if repeat else Path(record["stage"])
        relative /= (
            f"{LABELS[record['cell']]}_seed{record['seed']}{'_repeat' if repeat else ''}.pdb"
        )
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(WORK / record["pdb"], target)
        assert digest(target) == record["pdb_sha256"]
        exports.append(
            {
                "file": relative.as_posix(),
                "cell": record["cell"],
                "seed": record["seed"],
                "stage": record["stage"],
                "repeat": repeat,
                "sha256": record["pdb_sha256"],
            }
        )
    write_json(destination / "files.json", exports)
    colors = {"N0": "cyan", "N1": "green", "B0": "red", "B1": "orange"}
    for stage in ["raw", "refined"]:
        entries = [r for r in records if r["stage"] == stage and r["attempt"] != "repeat-1"]
        if not entries:
            continue
        pml = [
            "reinitialize",
            "load reference_9CQD.pdb, reference",
            "hide everything",
            "show cartoon, reference",
            "color gray70, reference",
            "color yellow, reference and chain X and resi 28-34+51-54+72-79",
        ]
        for item in exports:
            if item["stage"] != stage or item["repeat"]:
                continue
            obj = f"{item['cell']}_seed{item['seed']}"
            pml.extend(
                [
                    f"load {item['file']}, {obj}",
                    f"hide everything, {obj}",
                    f"unbond ({obj} and chain M and resi 155-162), ({obj} and not chain M)",
                    f"show sticks, {obj} and chain M and resi 155-162",
                    f"color {colors[item['cell']]}, {obj}",
                    f"group {LABELS[item['cell']]}, {obj}",
                ]
            )
            if item["seed"] != 7:
                pml.append(f"disable {obj}")
        pml.extend(
            [
                "bg_color white",
                "set stick_radius, 0.12",
                "zoom reference and chain M and resi 163-170, 15",
            ]
        )
        (destination / f"inspect_{stage}.pml").write_text("\n".join(pml) + "\n")
        plot(entries, destination / f"{stage}-overlay.png")
    (destination / "README.md").write_text("""# 9CQD: cc94 и cc89 на одинаковом ODE-расписании

Сравнивается N-концевой хвост M155–162 рядом с цепью X. Все PDB содержат полный
комплекс в общей исходной системе координат; align/super не нужны.

- N0: cc94 без реинъекции; N1: cc94 с реинъекцией.
- B0: cc89 без реинъекции; B1: cc89 с реинъекцией.
- Во всех режимах: 500 шагов, step_scale=1.0, s_churn=0, фиксированная последовательность.
- raw/: до OpenMM; refined/: после одинакового generated-only OpenMM Reference уточнения.
- repeats/: независимые повторы seed 7, не дополнительные независимые seeds.

Запустить PyMOL из этой папки: `pymol inspect_raw.pml` или `pymol inspect_refined.pml`.
По умолчанию включён seed 7 каждой модели; остальные — в списке объектов.
Серый — исходный комплекс, жёлтый — участки X28–34/X51–54/X72–79.
Голубой/зелёный хвост — cc94, красный/оранжевый — cc89.
""")
    return destination


def plot(entries, path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from dvbfixer.model.diffusion.contract import ArtifactReference

    coordinates = {}
    for row in entries:
        coordinates[(row["cell"], row["seed"])] = load_pdb_coordinates(
            WORK, ArtifactReference(row["pdb"], row["pdb_sha256"])
        )
    original = next(iter(coordinates.values()))
    xca = np.asarray(
        [
            original[a]
            for a in sorted(original, key=lambda a: (a.chain, int(a.residue_number), a.atom_name))
            if a.chain == "X" and a.atom_name == "CA"
        ]
    )
    center = np.asarray(
        [original[AtomIdentity("M", str(r), "", "CA")] for r in range(163, 170)]
    ).mean(axis=0)
    fig = plt.figure(figsize=(14, 12))
    for panel, cell in enumerate(CELLS, 1):
        ax = fig.add_subplot(2, 2, panel, projection="3d")
        ax.plot(*xca.T, color="gray", alpha=0.6, lw=1)
        for seed in SEEDS:
            coords = coordinates.get((cell, seed))
            if coords is None:
                continue
            tail = np.asarray(
                [
                    coords[AtomIdentity("M", str(r), "", name)]
                    for r in range(155, 164)
                    for name in ["N", "CA", "C"]
                ]
            )
            ax.plot(*tail.T, lw=1.8, label=str(seed))
        ax.set(
            title=LABELS[cell],
            xlim=(center[0] - 20, center[0] + 20),
            ylim=(center[1] - 20, center[1] + 20),
            zlim=(center[2] - 20, center[2] + 20),
        )
        ax.set_box_aspect((1, 1, 1))
        ax.view_init(elev=25, azim=45)
        ax.legend(title="seed", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main():
    protocol = check_protocol()
    records = rows()
    write_json(WORK / "analysis/per-seed.json", records)
    summary = {}
    for stage in ["raw", "refined"]:
        summary[stage] = {}
        for cell in CELLS:
            items = [
                r
                for r in records
                if r["stage"] == stage and r["cell"] == cell and r["attempt"] != "repeat-1"
            ]
            summary[stage][cell] = {
                "n": len(items),
                "MX_penetration_cases": sum(r["MX_pairs_below_1_9_A"] > 0 for r in items),
                "native_MX_penetration_cases": sum(
                    r["native_MX_pairs_below_1_9_A"] > 0 for r in items
                ),
                "segment_alarm_cases": sum(r["segment_alarm"] for r in items),
                "severe_clash_cases": sum(r["authoritative_severe_clashes"] > 0 for r in items),
                "all_gate_pass_cases": sum(r["all_gates_pass"] for r in items),
            }
    write_json(WORK / "analysis/summary.json", summary)
    repeatability = {}
    for stage in ["raw", "refined"]:
        for cell in ["N1", "B1"]:
            items = [
                r for r in records if r["stage"] == stage and r["cell"] == cell and r["seed"] == 7
            ]
            if len(items) == 2:
                repeatability[f"{stage}-{cell}"] = {
                    "same_pdb_sha256": items[0]["pdb_sha256"] == items[1]["pdb_sha256"]
                }
    write_json(WORK / "analysis/repeatability.json", repeatability)
    columns = [
        "cell",
        "model",
        "reinjection",
        "seed",
        "attempt",
        "stage",
        "min_MX_A",
        "MX_pairs_below_1_9_A",
        "segment_min_A",
        "segment_alarm",
        "authoritative_severe_clashes",
        "junction_CN_A",
        "native_MX_min_A",
        "native_MX_pairs_below_1_9_A",
        "all_gates_pass",
        "pdb",
    ]
    with (WORK / "analysis/per-seed.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for r in records:
            writer.writerow({k: r[k] for k in columns})
    visual = export_visual(records)
    shutil.copyfile(WORK / "analysis/per-seed.csv", visual / "per-seed.csv")
    record = {
        "schema_version": 1,
        "protocol": protocol,
        "model_inventory": json.loads((WORK / "model-inventory.json").read_text()),
        "summary": summary,
        "repeatability": repeatability,
        "records": records,
    }
    write_json(WORK / "analysis/results.json", record)
    print(json.dumps(summary, indent=2))
    print("Visual structures:", visual)


if __name__ == "__main__":
    main()
