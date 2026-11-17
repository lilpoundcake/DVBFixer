"""Private fixed-water obstacle checks must not hide short pairs as bonds."""
from dataclasses import replace

import numpy as np
import pytest

from dvbfixer.diagnose.report import Severity
from dvbfixer.model.diffusion.rfd3_custody_records import CustodyError, digest
from dvbfixer.model.diffusion.rfd3_water_context import assess_water_context
from tests.test_rfd3_source_link_custody import assemble, fixture


def water(x=50., y=50., z=50., occupancy=1., altloc=" "):
    return (f"HETATM  900  O  {altloc}HOH W 101    {x:8.3f}{y:8.3f}{z:8.3f}"
            f"{occupancy:6.2f}{20.:6.2f}           O\n").encode()


def full_context(source, record):
    return source.replace(b"END\n", record + b"END\n")


def test_water_reference_is_immutable_assessment_context_not_native_conditioning():
    request, source, _, _ = fixture()
    context = full_context(source, water(occupancy=.5, altloc="A"))
    result = assess_water_context(request, source, context, context_sha256=digest(context))
    assert result.water_atoms == result.nonunit_water_occupancies == result.water_altlocs == 1
    assert result.source_context_clear
    assert result.restored_region_clear is None
    assert not result.native_water_conditioning and not result.scientific_acceptance


def test_severe_source_contact_is_not_removed_as_inferred_covalent_bond():
    request, source, _, _ = fixture()
    context = full_context(source, water(.1, 0, 0))
    result = assess_water_context(request, source, context, context_sha256=digest(context))
    assert not result.source_context_clear
    assert any(c.protein_atom.atom_name == "N" and c.distance_angstrom < .2 and c.severity == Severity.ERROR for c in result.baseline_contacts)


def test_canonical_heavy_hbond_envelope_remains_distinct_from_clashes():
    request, source, _, _ = fixture()
    # S-O at2.6A exceeds the warning overlap without the canonical H-bond
    # envelope,unlike N-O at this distance which is already below threshold.
    context = full_context(source, water(3.1, -1.5, 1.3))
    result = assess_water_context(request, source, context, context_sha256=digest(context))
    assert not any(c.protein_atom.residue_number == "1" and c.protein_atom.atom_name == "SG" for c in result.baseline_contacts)


def test_new_missing_atom_contact_is_checked_against_all_original_waters():
    case = fixture()
    request, source, _, raw = case
    context = full_context(source, water())
    coordinates = np.array(raw.coordinates)
    index = next(i for i, atom in enumerate(raw.atom_axis) if atom.residue_number == "2" and atom.atom_name == "CA")
    coordinates[index] = [50.01, 50., 50.]
    candidate = assemble((*case[:3], replace(raw, coordinates=tuple(map(tuple, coordinates)))))
    result = assess_water_context(request, source, context, context_sha256=digest(context), candidate=candidate.text.encode())
    assert result.source_context_clear
    assert result.restored_region_clear is False
    assert any(c.distance_angstrom < .02 and c.protein_atom.residue_number == "2" for c in result.restored_contacts)


@pytest.mark.parametrize("mutation", ["digest", "protein-coordinate", "extra-component", "duplicate-water", "water-hydrogen", "link"])
def test_unsupported_or_changed_full_context_is_not_silently_dropped(mutation):
    request, source, _, _ = fixture()
    context = full_context(source, water())
    expected = digest(context)
    if mutation == "digest":
        context += b"REMARK changed\n"
    elif mutation == "protein-coordinate":
        line = next(line for line in source.splitlines(keepends=True) if line.startswith(b"ATOM"))
        context = context.replace(line, line[:30] + b" 999.000" + line[38:])
    elif mutation == "extra-component":
        context = context.replace(b"HOH", b"LIG")
    elif mutation == "duplicate-water":
        context = full_context(source, water() + water().replace(b"  900", b"  901"))
    elif mutation == "water-hydrogen":
        context = context.replace(b" O  ", b" H1 ").replace(b"           O", b"           H")
    else:
        context = context.replace(b"END\n", b"CONECT    1  900\nEND\n")
    if mutation != "digest":
        expected = digest(context)
    with pytest.raises(CustodyError):
        assess_water_context(request, source, context, context_sha256=expected)


def test_changed_fixed_candidate_records_fail_even_without_water_contacts():
    case = fixture()
    request, source, _, _ = case
    context = full_context(source, water())
    candidate = assemble(case).text
    line = next(line for line in source.decode().splitlines(keepends=True) if line.startswith("ATOM"))
    candidate = candidate.replace(line, line[:30] + " 999.000" + line[38:])
    with pytest.raises(CustodyError, match="fixed source records"):
        assess_water_context(request, source, context, context_sha256=digest(context), candidate=candidate.encode())
