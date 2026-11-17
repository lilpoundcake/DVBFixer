"""Water-preserving export fixes NEW serial collisions, not original records."""
from dataclasses import replace

import pytest

from dvbfixer.ffutils.dat import DatRecord
from dvbfixer.model.diffusion.rfd3_custody_records import CustodyError, digest, read_records
from dvbfixer.model.diffusion.rfd3_full_context import (
    materialize_full_context,
    target_ssbond_record,
)
from tests.test_rfd3_refinement import case


def materialize(data):
    request, source, context, projection = data
    return materialize_full_context(request, source, context, projection,
        context_sha256=digest(context), projection_sha256=digest(projection),
        target_pairs=tuple((p.atom1, p.atom2) for p in request.retained_explicit_links))


def test_water_collision_reassigns_only_missing_atom_serials_and_dat_keys(tmp_path):
    request, source, context, projection = case()
    context = context.replace(b"HETATM  900", b"HETATM   13")
    result = materialize((request, source, context, projection))
    atoms, _, _, _ = read_records(result.text)
    assert len(atoms) == 17 and len({a.serial for a in atoms}) == 17
    assert result.original_protein_atoms == 12 and result.original_water_atoms == 1
    for line in context.decode().splitlines(keepends=True):
        assert line in result.text
    water = next(a for a in atoms if a.residue_name == "HOH")
    assert water.serial == 13 and water.record.encode() in context
    assert len(result.new_atom_serial_mapping) == 4
    assert min(new for _, _, new in result.new_atom_serial_mapping) == 14
    path = tmp_path / "restored.dat"
    result.dat.save(path, verbose=False)
    dat = DatRecord.load(path)
    assert dat.added_keys() == {(a.chain, a.residue_number, a.insertion_code, a.atom_name) for a in request.generated_atoms}
    assert len(dat.added_atoms) == 4 and not result.scientific_qualification


@pytest.mark.parametrize("problem", ["source", "hash", "conect", "target-graph", "capacity"])
def test_invalid_full_context_exports_fail_without_changing_input(problem):
    request, source, context, projection = case()
    if problem == "source":
        projection = projection.replace(b" 20.00", b" 30.00", 1)
    elif problem == "hash":
        request = replace(request, normalized_pdb=replace(request.normalized_pdb, sha256="0" * 64))
    elif problem == "conect":
        context = context.replace(b"END\n", b"CONECT    1    2\nEND\n")
    elif problem == "target-graph":
        request = replace(request, retained_explicit_links=())
    else:
        context = context.replace(b"HETATM  900", b"HETATM99999")
    with pytest.raises(CustodyError):
        materialize((request, source, context, projection))


def test_derived_target_ssbond_keeps_case_and_insertion_code():
    from dvbfixer.model.diffusion.contract import AtomIdentity

    first, last = AtomIdentity("a", "7", "B", "SG"), AtomIdentity("A", "16", "", "SG")
    line = target_ssbond_record(19, (first, last), 2.05)
    assert line[15] == "a" and line[29] == "A" and line[21] == "B"
    assert int(line[7:10]) == 19 and float(line[73:78]) == 2.05 and line.endswith("\n")
