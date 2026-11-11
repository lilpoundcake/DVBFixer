#!/usr/bin/env python3
"""No-sampling topology/trajectory limitation controls, never a fatal policy.

An analytically linked pair of disjoint polygons is a legitimate mathematical
control, not proof of a native protein knot. A synthetic interpolated traversal
illustrates why clean endpoints cannot certify an unsaved optimizer trajectory.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scripts.mac_reconstruction_comparison import save

from dvbfixer.model.diffusion.backbone_intersections import BackboneSegment, screen_segments
from dvbfixer.model.diffusion.contract import AtomIdentity
from dvbfixer.model.diffusion.loop_surface_probe import probe_loop_surface
from dvbfixer.model.diffusion.scope import _parse_coordinate_records, _parse_explicit_links


def segments(points, chain):
    atoms = [AtomIdentity(chain, str(i+1), "", "CA") for i in range(len(points))]
    return tuple(BackboneSegment(atoms[i], atoms[(i+1) % len(points)], tuple(points[i]), tuple(points[(i+1) % len(points)]))
                 for i in range(len(points)))


def linked_control(rotation=None, translation=None):
    angles = np.linspace(0, 2*np.pi, 24, endpoint=False)+.17
    first = np.column_stack((5*np.cos(angles), 5*np.sin(angles), np.zeros_like(angles)))
    second = np.column_stack((5+5*np.cos(angles), np.zeros_like(angles), 5*np.sin(angles)))
    if rotation is not None:
        first, second = first @ rotation.T, second @ rotation.T
    if translation is not None:
        first, second = first+translation, second+translation
    first_segments, second_segments = segments(first, "A"), segments(second, "a")
    return {"segments": asdict(screen_segments(first_segments+second_segments, frozenset())),
        "surface": asdict(probe_loop_surface(first, second_segments)),
        "truth": "analytically linked disjoint polygonal circles; not an impossible configuration",
        "native_protein_topology_qualification": False}


def trajectory_control():
    angles = np.linspace(0, 2*np.pi, 16, endpoint=False)+.11
    polygon = np.column_stack((5*np.cos(angles), 5*np.sin(angles), np.zeros_like(angles)))
    first, second = AtomIdentity("B", "1", "", "N"), AtomIdentity("B", "1", "", "CA")
    snapshots = {}
    for label, z in (("start", (2., 3.)), ("intermediate", (-.5, .5)), ("end", (-3., -2.))):
        path = BackboneSegment(first, second, (-1., .23, z[0]), (1., .23, z[1]))
        snapshots[label] = asdict(probe_loop_surface(polygon, (path,)))
    return {"snapshots": snapshots, "trajectory": "declared synthetic linear interpolation, not a recorded optimizer trajectory",
        "endpoint_cleanliness_proves_no_passage": False, "public_gate": False}


def disulfide_control():
    lines = []
    for serial, chain, code, name, xyz in ((1, "A", "A", "SG", (0., 0., 0.)), (2, "a", "", "SG", (2.03, 0., 0.))):
        lines.append(f"ATOM  {serial:5d}  {name:<3} CYS {chain}{82:4d}{code or ' '}   "
                     f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{1:6.2f}{20:6.2f}           S\n")
    record = list(" "*80)
    record[:6] = "SSBOND"
    record[7:10], record[11:14], record[15], record[17:21], record[21] = "  1", "CYS", "A", "  82", "A"
    record[25:28], record[29], record[31:35] = "CYS", "a", "  82"
    lines.append("".join(record)+"\n")
    atoms, _, _ = _parse_coordinate_records(lines)
    links = _parse_explicit_links(lines, atoms)
    expected = tuple(sorted((AtomIdentity("A", "82", "A", "SG"), AtomIdentity("a", "82", "", "SG"))))
    if links != {expected}:
        raise ValueError("shared explicit-link parser lost chain case or insertion code")
    return {"SSBOND_preserved": True, "chain_case_and_insertion_code_preserved": True,
        "parser": "shared diffusion.scope explicit-link parser; no new distance inference",
        "disulfide_surface_or_knot_qualification": False}


def run(output):
    angle = .43
    rotation = np.array([[np.cos(angle), -np.sin(angle), 0.], [np.sin(angle), np.cos(angle), 0.], [0., 0., 1.]])
    original = linked_control()
    transformed = linked_control(rotation, np.array([17., -31., 6.]))
    trajectory, disulfide = trajectory_control(), disulfide_control()
    if original["segments"]["alarms"] or len(original["surface"]["piercings"]) != 1:
        raise ValueError("linked-circle control does not match its analytic geometry")
    if len(transformed["surface"]["piercings"]) != 1 or transformed["segments"]["alarms"]:
        raise ValueError("rigid transform changed diagnostic results")
    if any(trajectory["snapshots"][label]["piercings"] for label in ("start", "end")) or not trajectory["snapshots"]["intermediate"]["piercings"]:
        raise ValueError("synthetic endpoint/intermediate traversal control failed")
    output.mkdir(parents=True, exist_ok=False)
    save(output / "controls.json", {"linked_geometry": original, "rigid_transform": transformed,
        "trajectory_limit": trajectory, "explicit_disulfide": disulfide,
        "all_controls_passed": True, "public_fatal_gate_enabled": False, "native_topology_or_trajectory_qualified": False})
    (output / "entry-script.py").write_bytes(Path(__file__).read_bytes())
    print("Linked geometry, rigid invariance, endpoint limitation and explicit SSBOND transport controls passed; no topology qualification.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.output.resolve())
