"""Offline contract controls; invented coordinates, no RFD3/model imports."""
from dataclasses import replace

import numpy as np
import pytest

from dvbfixer.model.diffusion.contract import (
    DIFFUSION_SCHEMA_VERSION,
    ArtifactReference,
    AtomIdentity,
    DiffusionRequest,
    ExplicitLink,
    GapRegion,
    ResidueIdentity,
    SequencePlacement,
    TargetInterval,
    TargetSequence,
)
from dvbfixer.model.diffusion.rfd3_custody import (
    RawPrecision,
    RawProposal,
    materialize_private_candidate,
)
from dvbfixer.model.diffusion.rfd3_custody_records import (
    CustodyError,
    build_ledger,
    digest,
    expected_names,
)
from dvbfixer.model.diffusion.rfd3_custody_validation import audit_final_bytes, audit_geometry

PRECISION = RawPrecision(1e-8)  # Synthetic algebra only, NOT a sampling tolerance.


def ssbond(chain, first=1, last=3, code=""):
    chars = list(" " * 80)
    for start, value in [(0, "SSBOND"), (7, "  1"), (11, "CYS"), (15, chain),
                         (17, f"{first:4d}"), (25, "CYS"), (29, chain), (31, f"{last:4d}"), (35, code)]:
        chars[start:start + len(value)] = value
    return "".join(chars) + "\n"


def fixture(chains=("A",), *, code="", oxt=False):
    source, fixed, generated, placements, gaps, links, targets = [], [], [], [], [], [], []
    coords = {}
    serial = 1
    points = {
        1: {"N": (0, 0, 0), "CA": (1.45, 0, 0), "C": (2, 1.4, 0), "O": (1.7, 2.5, 0), "CB": (1.8, -.8, -1.2), "SG": (3.1, -1.5, -1.3)},
        2: {"N": (3.1, 1.4, 0), "CA": (4.2, 1.8, 0), "C": (4.9, .7, 0), "O": (4.5, -.4, 0)},
        3: {"N": (5.5, 1, 0), "CA": (6.9, 1, 0), "C": (7.5, 2.4, 0), "O": (7.3, 3.5, 0), "CB": (7.3, .2, -1.2), "SG": (5.1, -1.5, -1.3)},
    }
    for chain_index, chain in enumerate(chains):
        residues = tuple(ResidueIdentity(chain, str(i), code if i == 3 else "") for i in (1, 2, 3))
        source.append(ssbond(chain, code=code))
        for number, atoms in points.items():
            for name, point in atoms.items():
                atom = AtomIdentity(chain, str(number), code if number == 3 else "", name)
                xyz = tuple(float(v) for v in np.asarray(point) + [0, chain_index * 12, 0])
                coords[atom] = xyz
                if number == 2:
                    generated.append(atom)
                    continue
                fixed.append(atom)
                source.append(f"ATOM  {serial:5d}  {name:<3} CYS {chain}{number:4d}{atom.insertion_code or ' '}   "
                              f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{1.:6.2f}{20.:6.2f}          {name[0]:>2s}\n")
                serial += 1
        if oxt:
            atom = AtomIdentity(chain, "3", code, "OXT")
            fixed.append(atom)
            coords[atom] = (8., 2., 0.)
            source.append(f"ATOM  {serial:5d}  OXT CYS {chain}   3{code or ' '}   "
                          f"{8.:8.3f}{2.:8.3f}{0.:8.3f}{1.:6.2f}{20.:6.2f}           O\n")
            serial += 1
        source.append("TER\n")
        targets.append(TargetSequence(chain, "CGC"))
        placements.append(SequencePlacement(chain, 3, (0, 2), (residues[0], residues[2])))
        gaps.append(GapRegion(chain, TargetInterval(1, 2), residues[0], residues[2], (residues[1],), residues))
        links.append(ExplicitLink(AtomIdentity(chain, "1", "", "SG"), AtomIdentity(chain, "3", code, "SG"), "PDB"))
    source.append("END\n")
    data = "".join(source).encode()
    request = DiffusionRequest(DIFFUSION_SCHEMA_VERSION, ArtifactReference("source.pdb", digest(data)),
        tuple(targets), tuple(placements), tuple(gaps), tuple(fixed), tuple(generated), tuple(links), 1, (7,))
    ledger = build_ledger(request, data, data)
    names = expected_names(request)
    axis = tuple(atom for atom in (*fixed, *generated) if atom.atom_name != "OXT")
    raw = RawProposal(axis, tuple(names[a] for a in axis), tuple(coords[a] for a in axis),
                      tuple(a in fixed for a in axis), (True,) * len(axis))
    return request, data, ledger, raw


def assemble(case):
    request, source, ledger, raw = case
    return materialize_private_candidate(request, ledger, source, source, raw, precision=PRECISION)


def verify(case, text):
    request, source, ledger, raw = case
    return audit_final_bytes(request, ledger, source, source, raw, text.encode(), precision=PRECISION)


@pytest.mark.parametrize("chains,code,oxt", [(("A",), "", False), (("a",), "A", True), (("A", "a"), "A", True)])
def test_custody_preserves_exact_observations_links_and_OXT(chains, code, oxt):
    case = fixture(chains, code=code, oxt=oxt)
    result = assemble(case)
    audit = verify(case, result.text)
    assert audit.source_link_custody
    assert audit.retained_link_pairs == len(chains)
    assert audit.restored_atoms == 4 * len(chains)
    assert audit.geometry_passed is None
    assert not audit.scientific_acceptance and not result.scientific_acceptance
    assert audit.native_full_graph_transport == "not-assessed-by-this-auditor"
    for atom in case[2].atoms:
        assert atom.record in result.text


def test_only_generated_atoms_are_supplied_to_shared_materializer(monkeypatch):
    from dvbfixer.model.diffusion import rfd3_custody as adapter
    original = adapter.materialize_candidate_pdb
    case = fixture(oxt=True)
    seen = []

    def spy(source, axis, *args, **kwargs):
        seen.append(axis)
        return original(source, axis, *args, **kwargs)

    monkeypatch.setattr(adapter, "materialize_candidate_pdb", spy)
    assemble(case)
    assert seen == [case[0].generated_atoms]


def test_one_common_rotated_frame_is_resolved_without_moving_source_records():
    case = fixture(("A", "a"))
    raw = case[3]
    rotation = np.array([[0., -1, 0], [1, 0, 0], [0, 0, 1]])
    transformed = np.array(raw.coordinates) @ rotation + [13., -7., 3.]
    case = (*case[:3], replace(raw, coordinates=tuple(map(tuple, transformed))))
    assert verify(case, assemble(case).text).raw_max_residual_angstrom < 1e-10


@pytest.mark.parametrize("mutation", ["drift", "reflect", "degenerate", "nan", "mask", "sequence-mask", "name", "duplicate", "extra"])
def test_raw_mutations_fail_before_source_restoration(mutation):
    case = fixture()
    raw = case[3]
    xyz = np.array(raw.coordinates)
    if mutation == "drift":
        xyz[0, 0] += .2
    elif mutation == "reflect":
        xyz[:, 0] *= -1
    elif mutation == "degenerate":
        xyz[:, 2] = 0
    elif mutation == "nan":
        xyz[-1, 0] = np.nan
    elif mutation == "mask":
        raw = replace(raw, fixed_mask=(False, *raw.fixed_mask[1:]))
    elif mutation == "sequence-mask":
        raw = replace(raw, sequence_fixed_mask=(False,) * len(raw.atom_axis))
    elif mutation == "name":
        raw = replace(raw, residue_names=("ALA", *raw.residue_names[1:]))
    elif mutation == "duplicate":
        raw = replace(raw, atom_axis=(raw.atom_axis[1], *raw.atom_axis[1:]))
    else:
        raw = replace(raw, atom_axis=(*raw.atom_axis, AtomIdentity("B", "1", "", "N")))
    raw = replace(raw, coordinates=tuple(map(tuple, xyz)))
    with pytest.raises(CustodyError):
        assemble((*case[:3], raw))


@pytest.mark.parametrize("mutation", ["link-delete", "endpoint", "observation", "restored-coordinate", "extra-link"])
def test_final_bytes_are_independently_reaudited(mutation):
    case = fixture()
    text = assemble(case).text
    if mutation == "link-delete":
        text = "".join(line for line in text.splitlines(keepends=True) if not line.startswith("SSBOND"))
    elif mutation == "endpoint":
        text = text.replace(ssbond("A"), ssbond("A", last=2))
    elif mutation == "extra-link":
        text = text.replace("END\n", "CONECT    1    2\nEND\n")
    else:
        lines = text.splitlines(keepends=True)
        index = next(i for i, line in enumerate(lines) if line.startswith("ATOM") and
                     int(line[22:26]) == (1 if mutation == "observation" else 2))
        lines[index] = lines[index][:30] + f"{9.:8.3f}" + lines[index][38:]
        text = "".join(lines)
    with pytest.raises(CustodyError):
        verify(case, text)


@pytest.mark.parametrize("line", ["CONECT    1 9999\n", "CONECT    1    1\n", "CONECT    1    X\n"])
def test_bad_source_conect_is_not_silently_dropped(line):
    request, source, _, _ = fixture()
    changed = source.replace(b"END\n", line.encode() + b"END\n")
    request = replace(request, normalized_pdb=ArtifactReference("source.pdb", digest(changed)))
    with pytest.raises(CustodyError):
        build_ledger(request, changed, changed)


def test_valid_reciprocal_conect_keeps_original_serial_endpoint_identity():
    request, source, _, raw = fixture()
    records = b"CONECT    6   12\nCONECT   12    6\n"
    changed = source.replace(b"END\n", records + b"END\n")
    request = replace(request, normalized_pdb=ArtifactReference("source.pdb", digest(changed)))
    ledger = build_ledger(request, changed, changed)
    case = request, changed, ledger, raw
    result = assemble(case)
    assert records.decode() in result.text
    assert verify(case, result.text).retained_link_pairs == 1
    assert len(ledger.link_records) == 3  # SSBOND plus two reciprocal CONECTs.


def test_source_request_and_context_are_bound():
    case = fixture()
    request, source, ledger, raw = case
    for s, ctx, req in [(source + b"\n", source, request), (source, source + b"\n", request),
                        (source, source, replace(request, seeds=(11,)))]:
        with pytest.raises(CustodyError):
            materialize_private_candidate(req, ledger, s, ctx, raw, precision=PRECISION)


def test_image_and_duplicate_source_annotations_refused():
    request, source, _, _ = fixture()
    image = ssbond("A").rstrip("\n")
    image = image[:59] + "  2555" + image[65:] + "\n"
    for changed in (source.replace(ssbond("A").encode(), image.encode()), source.replace(b"END\n", ssbond("A").encode() + b"END\n")):
        req = replace(request, normalized_pdb=ArtifactReference("source.pdb", digest(changed)))
        with pytest.raises(CustodyError):
            build_ledger(req, changed, changed)


def test_custody_pass_is_not_a_geometry_pass_for_invented_bad_coordinates():
    case = fixture()
    result = assemble(case)
    assert verify(case, result.text).source_link_custody
    try:
        assessment = audit_geometry(case[0], case[1], result.text.encode())
    except Exception as error:
        # The real chirality checker can reject these invented source coordinates.
        from dvbfixer.ffutils.geometry import ChiralityError
        assert isinstance(error, ChiralityError)
    else:
        assert not assessment.geometry_passed or assessment.admission_failures


@pytest.mark.parametrize("value", [0., -1., float("nan"), float("inf"), True])
def test_precision_has_no_implicit_or_nonfinite_default(value):
    with pytest.raises(CustodyError):
        RawPrecision(value)


def test_independent_chain_fit_is_not_a_valid_common_frame():
    case = fixture(("A", "a"))
    raw = case[3]
    points = np.array(raw.coordinates)
    for i, atom in enumerate(raw.atom_axis):
        if atom.chain == "a":
            points[i, 0] += 1.
    with pytest.raises(CustodyError, match="drift before restoration"):
        assemble((*case[:3], replace(raw, coordinates=tuple(map(tuple, points)))))


def test_drift_refusal_precedes_materialization(monkeypatch):
    from dvbfixer.model.diffusion import rfd3_custody as adapter
    case = fixture()
    points = np.array(case[3].coordinates)
    points[0, 0] += 1.
    invoked = []
    monkeypatch.setattr(adapter, "materialize_candidate_pdb", lambda *a, **kw: invoked.append(True))
    with pytest.raises(CustodyError, match="drift before restoration"):
        assemble((*case[:3], replace(case[3], coordinates=tuple(map(tuple, points)))))
    assert not invoked


def test_materializer_accident_cannot_rewrite_fixed_coordinates(monkeypatch):
    from dvbfixer.model.diffusion import rfd3_custody as adapter
    case = fixture()
    original = adapter.materialize_candidate_pdb

    def bad_materializer(*args, **kwargs):
        text = original(*args, **kwargs)
        record = case[2].atoms[0].record
        return text.replace(record, record[:30] + f"{99.:8.3f}" + record[38:])

    monkeypatch.setattr(adapter, "materialize_candidate_pdb", bad_materializer)
    with pytest.raises(CustodyError, match="materializer changed"):
        assemble(case)


@pytest.mark.parametrize("mutation", ["serial", "occupancy", "altloc", "element"])
def test_source_record_corruption_refused(mutation):
    request, source, ledger, _ = fixture()
    original = ledger.atoms[1].record
    if mutation == "serial":
        changed = original[:6] + f"{1:5d}" + original[11:]
    elif mutation == "occupancy":
        changed = original[:54] + f"{.5:6.2f}" + original[60:]
    elif mutation == "altloc":
        changed = original[:16] + "A" + original[17:]
    else:
        changed = original[:76] + " S" + original[78:]
    data = source.replace(original.encode(), changed.encode())
    request = replace(request, normalized_pdb=ArtifactReference("source.pdb", digest(data)))
    with pytest.raises(CustodyError):
        build_ledger(request, data, data)


def test_anisou_is_bound_and_preserved_not_just_coordinates():
    request, source, ledger, raw = fixture()
    atom = ledger.atoms[0].record
    anisou = "ANISOU" + atom[6:28] + "  10000  11000  12000    100    200    300       N  \n"
    data = source.replace(atom.encode(), (atom + anisou).encode())
    request = replace(request, normalized_pdb=ArtifactReference("source.pdb", digest(data)))
    case = request, data, build_ledger(request, data, data), raw
    result = assemble(case)
    assert anisou in result.text
    assert verify(case, result.text).source_link_custody
    with pytest.raises(CustodyError):
        verify(case, result.text.replace("  10000", "  20000"))


def test_nonboolean_mask_cannot_coerce_into_observation_ownership():
    case = fixture()
    raw = replace(case[3], fixed_mask=tuple(int(value) for value in case[3].fixed_mask))
    with pytest.raises(CustodyError, match="exact boolean axes"):
        assemble((*case[:3], raw))


def test_extra_generated_element_change_fails_final_readback():
    case = fixture()
    text = assemble(case).text
    lines = text.splitlines(keepends=True)
    index = next(i for i, line in enumerate(lines) if line.startswith("ATOM") and int(line[22:26]) == 2)
    lines[index] = lines[index][:76] + " S" + lines[index][78:]
    with pytest.raises(CustodyError, match="element"):
        verify(case, "".join(lines))
