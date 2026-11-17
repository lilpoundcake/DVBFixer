from types import SimpleNamespace

import pytest
from scripts.audit_retained_6z38_restoration import atoms
from scripts.audit_retained_full_context_controls import atoms as context_atoms
from scripts.audit_retained_full_context_controls import check_stage
from scripts.retained_modeller_full_context import (
    alignment,
    fixed_audit,
    merge_projection,
    read_atoms,
    residue_groups,
)
from scripts.run_retained_full_context_controls import declare

from dvbfixer.ffutils.dat import DatRecord
from dvbfixer.model.diffusion.contract import AtomIdentity, ResidueIdentity


def line(serial, name, resname, resid, x, *, record="ATOM", chain="a"):
    return (f"{record:6s}{serial:5d} {name:^4s} {resname:3s} {chain}{resid:4d}    "
            f"{x:8.3f}{0.:8.3f}{0.:8.3f}{1.:6.2f}{10.:6.2f}          {name[0]:>2s}  \n")


def test_full_context_merge_preserves_waters_ligands_headers_and_serials():
    protein = line(1, "CA", "ALA", 1, 0.)
    context = "HEADER    CONTROL\n" + protein + line(2, "C1", "PHB", 2, 4., record="HETATM") + line(3, "O", "HOH", 3, 8., record="HETATM") + "END\n"
    new = line(2, "N", "GLY", 4, 12.)
    merged = merge_projection(protein + new + "END\n", context)
    assert context.splitlines()[0] == merged.splitlines()[0]
    assert protein in merged
    assert all(atom.line in merged for atom in read_atoms(context).values())
    rows = list(read_atoms(merged).values())
    assert len({int(atom.line[6:11]) for atom in rows}) == 4
    assert rows[1].line[6:11].strip() == "4"


def test_full_context_merge_refuses_observation_change_or_collision():
    protein = line(1, "CA", "ALA", 1, 0.)
    water = line(2, "O", "HOH", 2, 4., record="HETATM")
    with pytest.raises(ValueError, match="observations"):
        merge_projection(line(1, "CA", "ALA", 1, .001), protein + water + "END\n")
    with pytest.raises(ValueError, match="collides"):
        merge_projection(protein + line(3, "O", "GLY", 2, 5.), protein + water + "END\n")


def test_prepared_protein_serial_alias_uses_original_full_context_record():
    prepared = line(2, "CA", "ALA", 1, 0.)
    full_record = line(4, "CA", "ALA", 1, 0.)
    context = full_record + line(2, "O", "HOH", 2, 4., record="HETATM") + "END\n"
    merged = merge_projection(prepared + "END\n", context)
    assert full_record in merged and prepared not in merged
    with pytest.raises(ValueError, match="observations"):
        merge_projection(prepared.replace("10.00", "11.00") + "END\n", context)


def test_native_context_motion_is_not_repaired_by_source_copy():
    context = line(1, "CA", "ALA", 1, 0.) + line(2, "O", "HOH", 2, 4., record="HETATM")
    moved = line(1, "CA", "ALA", 1, 0.) + line(2, "O", "HOH", 2, 4.001, record="HETATM")
    assert fixed_audit(context, read_atoms(context))["exact_at_PDB_precision"]
    report = fixed_audit(context, read_atoms(moved))
    assert not report["exact_at_PDB_precision"]
    assert len(report["changed_atoms"]) == 1


def test_alignment_includes_each_water_and_ligand_without_cropping():
    request = SimpleNamespace(target_sequences=(SimpleNamespace(chain="a", sequence="AAA"),),
                              sequence_placements=(SimpleNamespace(observed_target_indices=(0, 2)),))
    context = (line(1, "CA", "ALA", 1, 0.) + line(2, "CA", "ALA", 3, 3.)
               + line(3, "C1", "FAD", 4, 6., record="HETATM") + line(4, "O", "HOH", 5, 8., record="HETATM"))
    pir = alignment(request, context)
    assert "A-A.w*" in pir and "AAA.w*" in pir
    assert "FIRST:a:LAST:a" in pir


def test_noncontiguous_residue_identity_refuses_transport():
    with pytest.raises(ValueError, match="noncontiguous"):
        residue_groups(line(1, "N", "ALA", 1, 0.) + line(2, "CA", "ALA", 2, 1.) + line(3, "C", "ALA", 1, 2.))


def test_partial_occupancy_is_refused_not_normalized():
    text = line(1, "CA", "ALA", 1, 0.)
    with pytest.raises(ValueError, match="ambiguous"):
        read_atoms(text[:54] + "  0.50" + text[60:])


def test_fractional_unique_context_is_preserved_without_normalizing_protein(tmp_path):
    protein = line(1, "CA", "ALA", 1, 0.)
    ligand = line(2, "C1", "WNP", 3, 8., record="HETATM")
    ligand = ligand[:54] + "  0.48" + ligand[60:]
    context = protein + ligand + "END\n"
    merged = merge_projection(protein + line(3, "CA", "ALA", 2, 3.), context)
    assert ligand in merged
    path = tmp_path / "context.pdb"
    path.write_text(merged)
    assert len(context_atoms(path)) == 3
    fractional_protein = protein[:54] + "  0.48" + protein[60:]
    with pytest.raises(ValueError, match="ambiguous"):
        read_atoms(fractional_protein, allow_fractional_heterogens=True)
    path.write_text(fractional_protein + ligand)
    with pytest.raises(ValueError, match="unsupported"):
        context_atoms(path)


def test_zero_occupancy_or_alternate_context_is_still_refused():
    ligand = line(2, "C1", "WNP", 3, 8., record="HETATM")
    for text in (ligand[:54] + "  0.00" + ligand[60:], ligand[:16] + "A" + ligand[17:]):
        with pytest.raises(ValueError, match="ambiguous"):
            read_atoms(text, allow_fractional_heterogens=True)


def anisou(atom_line):
    return "ANISOU" + atom_line[6:27] + " " + "".join(f"{value:7d}" for value in (1000, 1200, 1400, -100, 200, -300)) + atom_line[70:]


def test_anisou_stays_adjacent_and_byte_identical_after_atom_insertion():
    first = line(1, "CA", "ALA", 1, 0.)
    last = line(2, "CA", "ALA", 3, 6.)
    ligand = line(3, "C1", "WNP", 5, 12., record="HETATM")
    source = "HEADER    TENSOR CONTROL\n" + first + anisou(first) + last + anisou(last) + "TER\n" + ligand + anisou(ligand) + "END\n"
    projection = first + line(4, "CA", "ALA", 2, 3.) + last + line(5, "CA", "ALA", 4, 9.) + "END\n"
    merged = merge_projection(projection, source)
    assert first + anisou(first) in merged
    assert last + anisou(last) in merged
    assert ligand + anisou(ligand) in merged
    assert [row for row in merged.splitlines() if row.startswith("ANISOU")] == [row for row in source.splitlines() if row.startswith("ANISOU")]
    stream = iter(merged.splitlines(keepends=True))
    assert all(any(row == original for row in stream) for original in source.splitlines(keepends=True))
    assert len(read_atoms(merged, allow_fractional_heterogens=True)) == 5


def test_anisou_dangling_or_wrong_atom_is_refused():
    atom = line(1, "CA", "ALA", 1, 0.)
    other = line(2, "CA", "ALA", 2, 3.)
    for source in (anisou(atom) + atom, atom + anisou(other), atom + anisou(atom) + anisou(atom)):
        with pytest.raises(ValueError, match="ANISOU"):
            merge_projection(atom, source + "END\n")


def independent_fixture(tmp_path):
    first, second = ResidueIdentity("a", "1", ""), ResidueIdentity("a", "2", "")
    generated = AtomIdentity("a", "2", "", "CA")
    request = SimpleNamespace(target_sequences=(SimpleNamespace(chain="a", sequence="AA"),),
        sequence_placements=(SimpleNamespace(observed_target_indices=(0,), observed_residues=(first,)),),
        gaps=(SimpleNamespace(target_interval=SimpleNamespace(start=1, stop=2), generated_residues=(second,), right_anchor=None),),
        generated_atoms=(generated,))
    protein = line(1, "CA", "ALA", 1, 0.)
    ligand = line(3, "C1", "FAD", 3, 6., record="HETATM")
    water = line(4, "O", "HOH", 4, 8., record="HETATM")
    native = protein + line(2, "CA", "ALA", 2, 3.) + ligand + line(4, "OH2", "HOH", 4, 8., record="HETATM")
    full = protein + line(5, "CA", "ALA", 2, 3.) + ligand + water
    paths = {}
    for name, text in (("source", protein + ligand + water), ("native", native), ("full", full)):
        path = tmp_path / f"{name}.pdb"
        path.write_text(text + "END\n")
        paths[name] = path
    dat_path = tmp_path / "generated.dat"
    record = DatRecord(added_atoms=[{"chain": "a", "resid": "2", "icode": "", "resname": "ALA", "atom": "CA", "element": "C"}])
    record.save(dat_path)
    return request, paths, dat_path


def test_independent_full_context_readback_verifies_native_water_and_dat(tmp_path):
    request, paths, dat_path = independent_fixture(tmp_path)
    report = check_stage(request, atoms(paths["source"]), atoms(paths["native"]), atoms(paths["full"]), dat_path)
    assert report["changed_native_observations"] == 0
    assert report["fixed_atoms"] == 3
    assert report["generated_atoms"] == 1
    assert report["canonical_dat_ownership_exact"]


def test_independent_full_context_readback_does_not_excuse_restored_native_water(tmp_path):
    request, paths, dat_path = independent_fixture(tmp_path)
    text = paths["native"].read_text().replace("   8.000", "   8.001")
    paths["native"].write_text(text)
    report = check_stage(request, atoms(paths["source"]), atoms(paths["native"]), atoms(paths["full"]), dat_path)
    assert report["changed_native_observations"] == 1
    assert report["max_native_observation_motion_angstrom"] == pytest.approx(.001)


def test_independent_full_context_readback_rejects_dat_that_owns_observations(tmp_path):
    request, paths, dat_path = independent_fixture(tmp_path)
    record = DatRecord.load(dat_path)
    record.added_atoms.append({"chain": "a", "resid": "1", "icode": "", "resname": "ALA", "atom": "CA", "element": "C"})
    record.save(dat_path)
    with pytest.raises(ValueError, match="sidecar"):
        check_stage(request, atoms(paths["source"]), atoms(paths["native"]), atoms(paths["full"]), dat_path)


def test_declaration_refuses_duplicate_cases_before_reserving_output(tmp_path):
    output = tmp_path / "run"
    with pytest.raises(ValueError, match="distinct"):
        declare(output, entries=("1CC6", "1CC6"))
    assert not output.exists()


@pytest.mark.parametrize("entries", [(), ("1CC6", "3C7D", "6ZH0")])
def test_declaration_stays_bounded_before_reserving_output(tmp_path, entries):
    output = tmp_path / "run"
    with pytest.raises(ValueError, match="one or two distinct"):
        declare(output, entries=entries)
    assert not output.exists()
