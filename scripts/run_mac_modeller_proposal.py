#!/usr/bin/env python3
"""One seeded, private full LoopModel proposal; no held-out-reference input.

Matches the production chirality restraints and gap-only loop-MD mask. Initial
automodel optimization is NOT source-fixed; native drift is measured downstream.
Unlike the production wrapper this retains the unswept native proposal verbatim.
One initial model and one fast loop derivative, with no score-based selection.
"""
from __future__ import annotations

import argparse
import hashlib
import os
from pathlib import Path

from scripts.mac_reconstruction_comparison import (
    DEVELOPMENT_SEEDS,
    modeller_alignment,
    normalize_modeller,
    save,
)

from dvbfixer.model.diffusion.contract import DiffusionRequest


def run(request_path: Path, source_path: Path, output: Path, seed: int, *, validation_policy: Path | None = None):
    permission = None
    if validation_policy is not None:
        from scripts.mac_validation_policy import require_policy
        permission = require_policy("modeller", seed, validation_policy)
    if seed not in DEVELOPMENT_SEEDS and permission is None:
        raise ValueError("only declared development seeds are allowed")
    request = DiffusionRequest.from_json(request_path.read_text())
    data = source_path.read_bytes()
    if hashlib.sha256(data).hexdigest() != request.normalized_pdb.sha256:
        raise ValueError("source differs from frozen request")
    if request.retained_explicit_links or request.heterogen_contexts:
        raise ValueError("first comparator slice excludes crosslinks and heterogens")
    alignment = modeller_alignment(request)
    output.mkdir(parents=True, exist_ok=False)
    (output / "source.pdb").write_bytes(data)
    (output / "alignment.pir").write_text(alignment)
    save(output / "declaration.json", {"seed": seed, "modeller_rand_seed": -seed,
        "validation_permission": permission,
        "initial_models": 1, "loop_models": 1, "loop_md": "fast", "candidate_selection": "single loop output, no fallback",
        "initial_optimization": "all atoms as in production; not source-preserving",
        "loop_optimization": "all missing residues together; no observation extension",
        "native_chirality_sweep": False, "reference_coordinates_read": False, "qualification": False})
    previous = Path.cwd()
    try:
        os.chdir(output)
        from modeller import Environ, Selection, log
        from modeller.automodel import LoopModel, refine

        from dvbfixer.model.modeller_run import _add_chirality_restraints

        class Comparator(LoopModel):
            def special_restraints(self, aln):
                super().special_restraints(aln)
                _add_chirality_restraints(self)

            def select_loop_atoms(self):
                loops = self.loops(self.read_alignment(), minlength=1, maxlength=9999,
                                   insertion_ext=0, deletion_ext=0, include_termini=True)
                selected = Selection(loops).only_std_residues()
                if not len(selected):
                    raise ValueError("empty gap selection; refuse all-model fallback")
                return selected

        log.none()
        env = Environ(rand_seed=-seed)
        env.io.atom_files_directory = [str(output)]
        model = Comparator(env, alnfile="alignment.pir", knowns="source", sequence="target")
        model.starting_model = model.ending_model = 1
        model.loop.starting_model = model.loop.ending_model = 1
        model.loop.md_level = refine.fast
        model.make()
        def summary(rows):
            return [{key: (str(value) if key == "failure" and value is not None else value)
                     for key, value in row.items() if key in {"name", "num", "failure", "molpdf"}}
                    for row in rows]
        save(output / "native-outputs.json", {"modeller_version": __import__("modeller").__version__,
            "initial": summary(model.outputs), "loop": summary(model.loop.outputs),
            "pdfterms": "opaque MODELLER Values not serialized; full engine files retained"})
        if len(model.loop.outputs) != 1 or model.loop.outputs[0]["failure"] is not None:
            raise ValueError("single MODELLER loop proposal failed; no initial-model substitution")
        native = Path(model.loop.outputs[0]["name"]).read_bytes()
        (output / "native.pdb").write_bytes(native)
        normalized, report = normalize_modeller(request, native.decode())
        (output / "mapped-native.pdb").write_text(normalized)
        save(output / "mapping.json", report)
    finally:
        os.chdir(previous)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request", required=True, type=Path)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--validation-policy", type=Path)
    args = parser.parse_args()
    run(args.request.resolve(), args.source.resolve(), args.output.resolve(), args.seed,
        validation_policy=None if args.validation_policy is None else args.validation_policy.resolve())
