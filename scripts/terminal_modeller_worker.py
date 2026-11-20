"""One source-only generated-residue MODELLER proposal for a private comparison."""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path


def main(root):
    from modeller import Environ, Selection, log
    from modeller.automodel import LoopModel, refine
    from scripts.mac_reconstruction_comparison import (
        modeller_alignment,
        normalize_modeller,
        save,
        target_axis,
    )

    from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest
    from dvbfixer.model.modeller_run import _add_chirality_restraints

    request = DiffusionRequest.from_json((root / "request.json").read_text())
    order = target_axis(request)[request.target_sequences[0].chain]
    generated = set(request.generated_atoms)
    calls = []

    class FixedContext(LoopModel):
        def select_atoms(self):
            selected, identities = [], []
            for residue, identity in zip(self.residues, order, strict=True):
                for atom in residue.atoms:
                    key = AtomIdentity(identity.chain, identity.residue_number, identity.insertion_code, atom.name)
                    if key in generated:
                        selected.append(atom)
                        identities.append(key)
            if len(identities) != len(generated) or set(identities) != generated:
                raise ValueError("MODELLER generated-only selection differs from request")
            calls.append({"selected_atoms": len(selected), "observed_selected": 0})
            return Selection(*selected)

        def select_loop_atoms(self):
            return self.select_atoms()

        def special_restraints(self, aln):
            super().special_restraints(aln)
            _add_chirality_restraints(self)

        def to_iupac(self):
            # No full-model atom swaps or hidden source restoration.
            calls.append({"operation": "to_iupac-not-applied"})

    (root / "alignment.pir").write_text(modeller_alignment(request))
    os.chdir(root)
    log.none()
    env = Environ(rand_seed=-42)
    env.io.atom_files_directory = [str(root)]
    model = FixedContext(env, alnfile="alignment.pir", knowns="source", sequence="target")
    model.starting_model = model.ending_model = 1
    model.loop.starting_model = model.loop.ending_model = 1
    model.loop.md_level = refine.fast
    model.final_malign3d = False
    model.make()
    save(root / "selection-calls.json", calls)
    if len(model.loop.outputs) != 1 or model.loop.outputs[0]["failure"] is not None:
        raise ValueError("single MODELLER loop failed; no fallback")
    native = Path(model.loop.outputs[0]["name"]).read_bytes()
    (root / "native.pdb").write_bytes(native)
    mapped, mapping = normalize_modeller(request, native.decode())
    (root / "mapped-native.pdb").write_text(mapped)
    save(root / "mapping.json", mapping)
    save(root / "completed.json", {"version": __import__("modeller").__version__, "seed": 42,
         "initial_models": 1, "loop_models": 1, "candidate_selection": False, "coordinate_restoration": False,
         "selection_calls": len(calls), "python": sys.version})


if __name__ == "__main__":
    output = Path(sys.argv[1]).resolve()
    try:
        main(output)
    except Exception as exc:
        import traceback
        with (output / "failure.json").open("x") as stream:
            json.dump({"error": str(exc), "traceback": traceback.format_exc()}, stream, indent=2)
        raise
