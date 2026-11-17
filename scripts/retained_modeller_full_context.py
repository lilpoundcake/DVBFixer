"""Private single-protein/full-context MODELLER controls; never a public backend."""
from __future__ import annotations

import math
import os
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path

from scripts.mac_reconstruction_comparison import normalize_modeller, save, target_axis
from scripts.run_retained_6z38_restoration import compare_findings

from dvbfixer.ffutils.dat import DatRecord
from dvbfixer.model.diffusion.contract import AtomIdentity


@dataclass(frozen=True)
class Atom:
    identity: AtomIdentity
    name: str
    element: str
    xyz: tuple[float, float, float]
    line: str


def read_atoms(text, *, allow_fractional_heterogens=False):
    result = {}
    for line in text.splitlines(keepends=True):
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atom = AtomIdentity(line[21], line[22:26].strip(), line[26].strip(), line[12:16].strip())
        element = line[76:78].strip()
        occupancy = float(line[54:60])
        occupancy_valid = occupancy == 1 or (
            allow_fractional_heterogens and line.startswith("HETATM") and 0 < occupancy <= 1)
        xyz = tuple(float(line[i:i+8]) for i in (30, 38, 46))
        if (atom in result or line[16].strip() or not occupancy_valid
                or element in {"H", "D"} or not all(math.isfinite(x) for x in xyz)):
            raise ValueError("ambiguous/nonheavy full-context atom")
        result[atom] = Atom(atom, line[17:20].strip(), element, xyz, line)
    if not result:
        raise ValueError("empty full-context input")
    return result


def residue_groups(text):
    groups, seen = [], set()
    for atom in read_atoms(text, allow_fractional_heterogens=True).values():
        key = (atom.identity.chain, atom.identity.residue_number, atom.identity.insertion_code)
        if not groups or groups[-1][0] != key:
            if key in seen:
                raise ValueError("noncontiguous residue identity")
            seen.add(key)
            groups.append((key, []))
        groups[-1][1].append(atom)
    return groups


def merge_projection(projection, context):
    """Insert only new protein atoms; every original context record stays exact."""
    before, after = read_atoms(context, allow_fractional_heterogens=True), read_atoms(projection)
    protein = {key: atom for key, atom in before.items() if atom.line.startswith("ATOM  ")}
    if any(not atom.line.startswith("ATOM  ") for atom in after.values()):
        raise ValueError("projection must contain only canonical protein ATOM records")
    if (not protein.keys() <= after.keys() or any(
            after[key].line[:6] + after[key].line[11:] != atom.line[:6] + atom.line[11:]
            for key, atom in protein.items())
            or tuple(key for key in after if key in protein) != tuple(protein)):
        raise ValueError("projection changed or reordered original protein observations")
    if (set(after) - set(protein)) & (set(before) - set(protein)):
        raise ValueError("new protein atom collides with context identity")
    lines = context.splitlines(keepends=True)
    if any(line.startswith(("MODEL ", "ENDMDL", "CONECT", "LINK  ", "SSBOND")) for line in lines):
        raise ValueError("this private control does not qualify linked/multimodel contexts")
    for index, line in enumerate(lines):
        if not line.startswith("ANISOU"):
            continue
        if (index == 0 or not lines[index - 1].startswith(("ATOM  ", "HETATM"))
                or line[6:27] != lines[index - 1][6:27]
                or line[76:78] != lines[index - 1][76:78]):
            raise ValueError("ANISOU does not immediately match its original atom")
        for start in (28, 35, 42, 49, 56, 63):
            int(line[start:start + 7])
    serial = max(int(line[6:11]) for line in lines if line.startswith(("ATOM  ", "HETATM", "TER   ")) and line[6:11].strip()) + 1
    rendered = []
    for key, atom in after.items():
        if key in protein:
            # Input preparation may have allocated these serials before context
            # was reattached. Only serial aliases are allowed, never a new name,
            # coordinate, occupancy, B factor or other observed record field.
            rendered.append((key, protein[key].line))
        else:
            if serial > 99999:
                raise ValueError("new atom serial exceeds PDB capacity")
            rendered.append((key, atom.line[:6] + f"{serial:5d}" + atom.line[11:]))
            serial += 1
    output, cursor, observed = [], 0, 0
    for line in lines:
        if line.startswith("ATOM  "):
            key = AtomIdentity(line[21], line[22:26].strip(), line[26].strip(), line[12:16].strip())
            while cursor < len(rendered) and rendered[cursor][0] != key:
                if rendered[cursor][0] in protein:
                    raise ValueError("projection observation order changed during context insertion")
                output.append(rendered[cursor][1])
                cursor += 1
            if cursor == len(rendered):
                raise ValueError("projection lost an original observation")
            output.append(line)
            cursor += 1
            observed += 1
        else:
            # Leave each source atom/ANISOU pair adjacent. Trailing new protein
            # atoms are inserted only after the final original atom's tensor.
            if observed == len(protein) and not line.startswith("ANISOU"):
                while cursor < len(rendered):
                    output.append(rendered[cursor][1])
                    cursor += 1
            output.append(line)
    if cursor != len(rendered) or observed != len(protein):
        raise ValueError("full-context insertion did not consume the exact protein axis")
    merged = "".join(output)
    final = read_atoms(merged, allow_fractional_heterogens=True)
    if len(final) != len(before) + len(after) - len(protein):
        raise ValueError("context merge inventory changed")
    if any(final[key].line != atom.line for key, atom in before.items()):
        raise ValueError("context merge changed an observation")
    return merged


def alignment(request, context):
    if len(request.target_sequences) != 1:
        raise ValueError("private slice supports one protein chain only")
    target = request.target_sequences[0]
    observed = set(request.sequence_placements[0].observed_target_indices)
    template = "".join(letter if i in observed else "-" for i, letter in enumerate(target.sequence))
    groups = residue_groups(context)
    heterogens = [atoms for _, atoms in groups if atoms[0].line.startswith("HETATM")]
    if any(atom.identity.chain != target.chain for atom in read_atoms(context, allow_fractional_heterogens=True).values()):
        raise ValueError("private slice needs one explicit full-context chain")
    tail = "".join("w" if atoms[0].name == "HOH" else "." for atoms in heterogens)
    return (f">P1;source\nstructureX:source:FIRST:{target.chain}:LAST:{target.chain}::::\n{template}{tail}*\n"
            f">P1;target\nsequence:target::::::::\n{target.sequence}{tail}*\n")


def native_transport(request, context, native):
    """Ordinal identity transport only. Never fit, copy, reflect or rename coordinates."""
    groups = residue_groups(native)
    count = len(request.target_sequences[0].sequence)
    context_groups = [atoms for _, atoms in residue_groups(context) if atoms[0].line.startswith("HETATM")]
    if len(groups) != count + len(context_groups):
        raise ValueError("native full-context residue count changed")
    polymer = "".join(atom.line for _, atoms in groups[:count] for atom in atoms) + "END\n"
    projection, mapping = normalize_modeller(request, polymer)
    mapped = read_atoms(projection)
    context_name_transports = []
    for (_, atoms), expected in zip(groups[count:], context_groups, strict=True):
        by_name = {atom.identity.atom_name: atom for atom in atoms}
        # MODELLER's documented water residue can name its one heavy oxygen OH2.
        if (len(expected) == 1 and expected[0].name == "HOH" and expected[0].identity.atom_name == "O"
                and set(by_name) == {"OH2"}):
            by_name = {"O": by_name["OH2"]}
            context_name_transports.append({"source": asdict(expected[0].identity), "native_name": "OH2"})
        if set(by_name) != {atom.identity.atom_name for atom in expected}:
            raise ValueError("native context atom inventory changed")
        for reference in expected:
            atom = by_name[reference.identity.atom_name]
            if reference.identity in mapped:
                raise ValueError("context/protein identity collision")
            mapped[reference.identity] = Atom(reference.identity, reference.name, reference.element, atom.xyz, atom.line)
    return projection, mapped, {"protein": mapping, "water_name_transports": context_name_transports}


def fixed_audit(context, transported):
    source = read_atoms(context, allow_fractional_heterogens=True)
    if not source.keys() <= transported.keys():
        raise ValueError("native output lost context observations")
    changes = [{"atom": asdict(key), "distance_angstrom": math.dist(atom.xyz, transported[key].xyz)}
               for key, atom in source.items() if atom.xyz != transported[key].xyz]
    return {"fixed_atoms": len(source), "changed_atoms": changes,
            "exact_at_PDB_precision": not changes, "fitting_or_coordinate_copy": False}


def findings(text, path):
    from openmm.app import PDBFile

    from dvbfixer.diagnose.chemistry import run_all
    from dvbfixer.diagnose.report import Severity
    from dvbfixer.diagnose.steric import clashes_python
    from dvbfixer.model.diffusion.quality import check_ramachandran, check_sidechain_chi12

    path.write_text(text)
    pdb = PDBFile(str(path))
    rows = run_all(pdb.topology, pdb.positions)
    rows += check_ramachandran(pdb.topology, pdb.positions)
    rows += check_sidechain_chi12(pdb.topology, pdb.positions)
    rows += clashes_python(pdb.topology, pdb.positions, include_water=True, clash_warn_a=.6, clash_error_a=.9)
    return [row for row in rows if row.severity != Severity.INFO]


def ownership(request, projection, path):
    atoms = read_atoms(projection)
    dat = DatRecord(description="Private full-context MODELLER proposal; mask ownership, not acceptance",
                    target_chains={target.chain: len(target.sequence) for target in request.target_sequences})
    for identity in request.generated_atoms:
        atom = atoms[identity]
        dat.added_atoms.append({"chain": identity.chain, "resid": identity.residue_number, "icode": identity.insertion_code,
                               "resname": atom.name, "atom": identity.atom_name, "element": atom.element})
        row = dat.residue_summary.setdefault(f"{identity.chain}/{atom.name}{identity.residue_number}{identity.insertion_code}",
                                            {"heavy": 0, "hydrogen": 0})
        row["heavy"] += 1
    if path.exists():
        raise ValueError("refuse to overwrite ownership sidecar")
    dat.save(path)


def run(request, source, context, output, seed):
    from modeller import Environ, Selection, log
    from modeller.automodel import LoopModel, refine
    from modeller.optimizers import ConjugateGradients

    from dvbfixer.ffutils.geometry import ChiralityError, assert_all_l
    from dvbfixer.model.diffusion.anchor_acceptance import _structure
    from dvbfixer.model.diffusion.pdb_materialize import materialize_candidate_pdb
    from dvbfixer.model.diffusion.staged_runtime import assessment_adapter
    from dvbfixer.model.modeller_run import _add_chirality_restraints

    output.mkdir(exist_ok=False)
    (output / "source.pdb").write_text(context)
    (output / "alignment.pir").write_text(alignment(request, context))
    baseline = findings(context, output / "baseline.pdb")
    axis = target_axis(request)[request.target_sequences[0].chain]
    generated = set(request.generated_atoms)
    calls = []

    class FixedContext(LoopModel):
        def select_atoms(self):
            selected, names = [], []
            for residue, identity in zip(tuple(self.residues)[:len(axis)], axis, strict=True):
                for atom in residue.atoms:
                    key = AtomIdentity(identity.chain, identity.residue_number, identity.insertion_code, atom.name)
                    if key in generated:
                        selected.append(atom)
                        names.append(key)
            if len(set(names)) != len(names) or set(names) != generated:
                raise ValueError("native generated-only selection changed")
            calls.append({"operation": "selection", "atoms": [asdict(key) for key in names]})
            return Selection(*selected)

        def select_loop_atoms(self):
            return self.select_atoms()

        def special_restraints(self, aln):
            super().special_restraints(aln)
            _add_chirality_restraints(self)

        def to_iupac(self):
            # Prevent the documented full-model atom swaps rather than repairing
            # observed coordinates afterwards. Generated geometry is still audited.
            calls.append({"operation": "to_iupac", "policy": "not applied; preserve named observations",
                          "source_coordinates_copied": False})

    previous = Path.cwd()
    try:
        os.chdir(output)
        log.none()
        env = Environ(rand_seed=-seed)
        env.io.atom_files_directory = [str(output)]
        env.io.hetatm = env.io.water = True
        model = FixedContext(env, alnfile="alignment.pir", knowns="source", sequence="target")
        model.starting_model = model.ending_model = 1
        model.loop.starting_model = model.loop.ending_model = 1
        model.loop.md_level = refine.fast
        model.final_malign3d = False
        model.make()
        save(output / "native-interface-calls.json", calls)
        save(output / "native-outputs.json", {"modeller_version": __import__("modeller").__version__,
            "loop": [{key: str(value) if key == "failure" and value is not None else value
                      for key, value in row.items() if key in {"name", "failure", "num", "molpdf"}}
                     for row in model.loop.outputs]})
        if len(model.loop.outputs) != 1 or model.loop.outputs[0]["failure"] is not None:
            raise ValueError("single native loop failed; no initial-model substitution")
        native_path = output / "native.pdb"
        native_path.write_bytes(Path(model.loop.outputs[0]["name"]).read_bytes())
        memory_path = output / "post-make-memory.pdb"
        model.write(file=str(memory_path))
        _, memory, _ = native_transport(request, context, memory_path.read_text())
        _, raw_parent, _ = native_transport(request, context, native_path.read_text())
        cg_parent_exact = memory.keys() == raw_parent.keys() and all(memory[key].xyz == raw_parent[key].xyz for key in memory)
        stages, raw_coordinates = [], None
        for stage in ("raw", "generated-only-CG"):
            if stage != "raw":
                if (not cg_parent_exact or not stages[-1]["native_fixed"]["exact_at_PDB_precision"]
                        or not stages[-1]["global_all_l"]):
                    break
                if stages[-1]["accepted"]:
                    break
                # One deterministic continuation, never a random restart or anchor release.
                ConjugateGradients().optimize(model.select_atoms(), max_iterations=1000)
                native_path = output / "generated-only-CG-native.pdb"
                model.write(file=str(native_path))
            projection, transported, mapping = native_transport(request, context, native_path.read_text())
            native_fixed = fixed_audit(context, transported)
            (output / f"{stage}-projection.pdb").write_text(projection)
            save(output / f"{stage}-mapping.json", mapping)
            ownership(request, projection, output / f"{stage}.dat")
            coordinates = {key: read_atoms(projection)[key].xyz for key in generated}
            if raw_coordinates is None:
                raw_coordinates = coordinates
            max_motion = max(math.dist(raw_coordinates[key], xyz) for key, xyz in coordinates.items())
            structure = _structure(projection, "native-protein.pdb")
            all_l = True
            try:
                assert_all_l(structure.pdb.topology, structure.pdb.positions)
            except ChiralityError:
                all_l = False
            # This copy is explicitly diagnostic on a native-fixation failure.
            atoms = read_atoms(projection)
            identities = tuple(request.generated_atoms)
            restored = materialize_candidate_pdb(source, identities, [atoms[key].xyz for key in identities],
                tuple(atoms[key].name for key in identities), request, elements=tuple(atoms[key].element for key in identities))
            # This record-custody check belongs to the source-record-preserving
            # stage; actual native fixation is measured separately above.
            assessment = assessment_adapter(request, source)(restored)
            merged = merge_projection(restored, context)
            (output / f"{stage}-full-context.pdb").write_text(merged)
            comparison = compare_findings(baseline, findings(merged, output / f"{stage}-diagnostic-input.pdb"))
            accepted = (native_fixed["exact_at_PDB_precision"] and all_l and assessment.geometry_passed
                        and not assessment.admission_failures and not comparison["new_or_changed_findings"]
                        and max_motion <= 3.)
            stages.append({"stage": stage, "native_fixed": native_fixed, "global_all_l": all_l,
                           "assessment": asdict(assessment), "full_diagnostics": comparison,
                           "max_generated_motion_from_raw_angstrom": max_motion, "movement_limit_angstrom": 3.,
                           "accepted": accepted})
        save(output / "stages.json", {"source_findings": [asdict(row) for row in baseline], "stages": stages,
             "new_geometry_by_stage": {row["stage"]: dict(Counter(f["category"] for f in row["full_diagnostics"]["new_or_changed_findings"])) for row in stages},
             "accepted": stages[-1]["accepted"], "public_capability_or_MD_readiness": False,
             "CG_in_memory_parent_matches_retained_raw_at_PDB_precision": cg_parent_exact,
             "postprocessing_policy": "documented to_iupac skipped, no coordinate copying before native audit"})
        save(output / "native-interface-calls-final.json", calls)
    finally:
        os.chdir(previous)
