"""Derivative custody is distinct from native raw fidelity and geometry."""
import pytest

from dvbfixer.model.diffusion.rfd3_custody_records import CustodyError, digest, read_records
from dvbfixer.model.diffusion.rfd3_refinement import audit_refined_bytes, build_refinement_input
from tests.test_rfd3_source_link_custody import assemble, fixture
from tests.test_rfd3_water_context import full_context, water


def case():
    original = fixture()
    request, source, _, _ = original
    raw = assemble(original).text.encode()
    context = full_context(source, water())
    return request, source, context, raw


def audit(data, final, *, bound=3.):
    request, source, context, raw = data
    return audit_refined_bytes(request, source, context, raw, final, raw_sha256=digest(raw),
                               context_sha256=digest(context), max_displacement_angstrom=bound)


def move(raw, *, residue=2, distance=.2):
    lines = raw.decode().splitlines(keepends=True)
    index = next(i for i, line in enumerate(lines) if line.startswith("ATOM") and int(line[22:26]) == residue)
    line = lines[index]
    lines[index] = line[:30] + f"{float(line[30:38]) + distance:8.3f}" + line[38:]
    return "".join(lines).encode()


def test_only_missing_coordinates_can_change_without_scientific_acceptance():
    data = case()
    report = audit(data, move(data[3]))
    assert report.observations_byte_identical == 12
    assert report.retained_link_pairs == 1
    assert report.changed_generated_atoms == 1
    assert report.generated_max_movement_angstrom == pytest.approx(.2)
    assert not report.scientific_acceptance


@pytest.mark.parametrize("mutation", ["fixed", "link", "header", "metadata", "extra-line", "too-far"])
def test_undeclared_derivative_changes_are_refused(mutation):
    data = case()
    raw = data[3]
    if mutation == "fixed":
        final = move(raw, residue=1)
    elif mutation == "link":
        final = b"".join(line for line in raw.splitlines(keepends=True) if not line.startswith(b"SSBOND"))
    elif mutation == "header":
        final = raw.replace(b"TER\n", b"REMARK changed\n", 1)
    elif mutation == "metadata":
        final = raw.replace(b" 20.00", b" 30.00", 1)
    elif mutation == "extra-line":
        final = raw + b"REMARK extra\n"
    else:
        final = move(raw, distance=4.)
    with pytest.raises(CustodyError):
        audit(data, final)


def test_ff_input_adds_all_original_waters_with_unambiguous_temporary_serials():
    request, source, context, raw = case()
    result = build_refinement_input(request, source, context, raw, context_sha256=digest(context), raw_sha256=digest(raw))
    atoms, _, _, _ = read_records(result.text)
    assert result.water_atoms == 1 and len(atoms) == 17
    assert len({atom.serial for atom in atoms}) == 17
    oxygen = next(atom for atom in atoms if atom.residue_name == "HOH")
    assert oxygen.xyz == (50., 50., 50.) and oxygen.identity.chain == "W"
    assert result.water_serial_mapping == ((900, oxygen.serial),)
    assert raw[:-4].decode() in result.text


@pytest.mark.parametrize("bound", [0, -1, float("inf"), float("nan"), True])
def test_no_implicit_invalid_movement_bound(bound):
    data = case()
    with pytest.raises(CustodyError):
        audit(data, data[3], bound=bound)


def test_source_context_and_raw_hashes_are_bound():
    request, source, context, raw = case()
    with pytest.raises(CustodyError):
        audit_refined_bytes(request, source, context, raw, raw, raw_sha256="0" * 64,
                            context_sha256=digest(context), max_displacement_angstrom=3.)
    with pytest.raises(CustodyError):
        build_refinement_input(request, source, context + b"\n", raw,
                               context_sha256=digest(context), raw_sha256=digest(raw))


def test_existing_kernel_freezes_water_oxygen_not_disposable_hydrogens(monkeypatch):
    from pathlib import Path

    from openmm.unit import dalton

    from dvbfixer.model.diffusion import boundary_refinement as kernel
    from dvbfixer.model.diffusion.contract import AtomIdentity, ResidueIdentity

    source = Path(__file__).parent / "fixtures/8cz8/8cz8_a_u.pdb"
    lines = [line for line in source.read_text().splitlines(keepends=True)
             if line.startswith("ATOM  ") and line[21] == "C" and 6 <= int(line[22:26]) <= 8
             and line[76:78].strip() != "H"]
    mobile = tuple(AtomIdentity("C", "7", "", line[12:16].strip()) for line in lines if int(line[22:26]) == 7)
    text = "".join(lines) + "TER\n" + water(100., 100., 100.).decode() + "END\n"
    freeze = kernel._freeze_observed_heavy_atoms
    observations = []

    def inspect(system, topology, original, generated):
        freeze(system, topology, original, generated)
        for atom in topology.atoms():
            if atom.residue.name == "HOH":
                observations.append((atom.element.symbol, system.getParticleMass(atom.index).value_in_unit(dalton)))

    monkeypatch.setattr(kernel, "_freeze_observed_heavy_atoms", inspect)
    result = kernel.refine_generated_region(text, generated_residues=(ResidueIdentity("C", "7"),),
        generated_atoms=mobile, platform_name="CPU", max_iterations=10, restart_count=0, random_seed=101)
    assert set(result.coordinates_angstrom) == set(mobile)
    assert observations.count(("O", 0.)) == 1
    assert sum(element == "H" and mass > 0 for element, mass in observations) == 2
