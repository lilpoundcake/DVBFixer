"""Source binding/readback controls; not production-backend qualification."""
from __future__ import annotations

import hashlib
from dataclasses import replace
from pathlib import Path

import pytest

from dvbfixer.domain.structure_identity import ComponentInstanceRef
from dvbfixer.pdb_component_reconstruction import _format_generated_atom
from dvbfixer.protein_completion.geometry import l_cb
from dvbfixer.protein_completion.service import complete_sidechain
from dvbfixer.protein_completion.source import (
    SourceCompletionError,
    audit_materialized_source,
    materialize_source,
    plan_source,
)

TARGET = ComponentInstanceRef(1, "a", "82", "A")
POINTS = {"N": (-1.46, 0., 0.), "CA": (0., 0., 0.), "C": (.55, 1.43, 0.), "O": (1.76, 1.66, 0.)}


def atom(serial, name, xyz, residue="ALA", chain="a", number="82", code="A"):
    line = _format_generated_atom(serial=serial, name=name, element=name[0], position=xyz,
                                  residue_name=residue, chain_id=chain, sequence_number=number,
                                  insertion_code=code)
    return ("ATOM  "+line[6:]).encode()


def source(residue="ALA"):
    return b"HEADER    PRIVATE SOURCE CONTROL\n"+b"".join(atom(i, n, p, residue) for i, (n, p) in enumerate(POINTS.items(), 1))+b"TER\nEND\n"


@pytest.mark.parametrize("residue", ("ALA", "SER", "CYS"))
def test_generated_only_pdb_readback_and_canonical_dat(residue, tmp_path):
    src = source(residue)
    req = plan_source(src, TARGET)
    candidate = complete_sidechain(req).selected
    output = materialize_source(src, req, candidate)
    assert output.source_sha256 == hashlib.sha256(src).hexdigest()
    assert output.output_sha256 == hashlib.sha256(output.pdb).hexdigest()
    assert len(output.generated) == len(req.missing_atoms)
    assert audit_materialized_source(src, req, candidate, output.pdb) == output.generated
    record = output.dat_record()
    assert {a["atom"] for a in record.added_atoms} == set(req.missing_atoms)
    assert all(a["chain"] == "a" and a["icode"] == "A" and a["resname"] == residue for a in record.added_atoms)
    path = tmp_path/"result.dat"
    record.save(path)
    assert record.load(path).added_atoms == record.added_atoms
    # Every old physical record is retained verbatim, including bare TER.
    old = iter(output.pdb.splitlines(keepends=True))
    assert all(any(found == line for found in old) for line in src.splitlines(keepends=True))


@pytest.mark.parametrize("mutation", (
    lambda r: replace(r, observations=r.observations[:-1]),
    lambda r: replace(r, observations=(replace(r.observations[0], xyz=(9., 9., 9.)), *r.observations[1:])),
    lambda r: replace(r, source_sha256="0"*64),
    lambda r: replace(r, missing_atoms=()),
))
def test_caller_cannot_supply_truncated_or_tampered_context(mutation):
    src = source()
    req = plan_source(src, TARGET)
    candidate = complete_sidechain(req).selected
    with pytest.raises(SourceCompletionError, match="binding"):
        materialize_source(src, mutation(req), candidate)


@pytest.mark.parametrize("tail", (b"MODEL        1\n", b"ENDMDL\n", b"LINK\n", b"SSBOND\n", b"MASTER\n",
                                  b"CONECT    1  999\n", b"CONECT    1\n", b"CONECTxxxxx    2\n",
                                  b"ANISOU    1\n", b"TER       1\n"))
def test_unsupported_layout_and_dangling_records_refuse(tail):
    with pytest.raises(ValueError):
        plan_source(source().replace(b"END\n", tail+b"END\n"), TARGET)


@pytest.mark.parametrize("change", (
    lambda s: s.replace(b"  1.00", b"  0.50", 1),
    lambda s: s.replace(b" ALA", b" VAL", 1),
    lambda s: s.replace(b" ALA", b" LIG"),
    lambda s: s.replace(b" ALA", b"AALA", 1),
    lambda s: s.replace(b"END\n", b"END\nREMARK after END\n"),
    lambda s: s.replace(b"\n", b"\r\n"),
    lambda s: s[:-1],
    lambda s: s.replace(b"TER\n", atom(9, "N", (8., 0., 0.), number="83")+atom(10, "CB", (9., 0., 0.))),
))
def test_source_ambiguity_occupancy_unknown_chemistry_and_global_undefined_refuse(change):
    with pytest.raises(ValueError):
        plan_source(change(source()), TARGET)


def test_anisou_and_sparse_serials_preserved_without_collision():
    src = source()
    first = src.splitlines(keepends=True)[1]
    anisou = b"ANISOU"+first[6:28]+b"   1000   1000   1000      0      0      0"+first[70:]
    src = src.replace(first, first+anisou).replace(b"TER\n", b"TER   99998\n")
    req = plan_source(src, TARGET)
    out = materialize_source(src, req, complete_sidechain(req).selected)
    assert anisou in out.pdb and b"ATOM  99999" in out.pdb
    exhausted = src.replace(b"99998", b"99999")
    req = plan_source(exhausted, TARGET)
    with pytest.raises(ValueError, match="limits"):
        materialize_source(exhausted, req, complete_sidechain(req).selected)


def test_chain_case_and_insertion_codes_remain_distinct():
    src = source()
    for chain, code, offset in (("A", "A", 30), ("a", "B", 60)):
        points = POINTS | {"CB": l_cb(POINTS["N"], POINTS["CA"], POINTS["C"])}
        extra = b"".join(atom(i+offset, n, tuple(x+offset for x in p), chain=chain, code=code)
                         for i, (n, p) in enumerate(points.items(), 1))
        src = src.replace(b"END\n", extra+b"TER\nEND\n")
    req = plan_source(src, TARGET)
    assert len({o.atom.component for o in req.observations}) == 3
    out = materialize_source(src, req, complete_sidechain(req).selected)
    assert all(a.atom.component == TARGET for a in out.generated)


def test_global_d_ca_elsewhere_is_rejected_not_reflected():
    cb = l_cb(POINTS["N"], POINTS["CA"], POINTS["C"])
    points = POINTS | {"CB": (cb[0], cb[1], -cb[2])}
    extra = b"".join(atom(i+10, n, tuple(x+30 for x in p), chain="B") for i, (n, p) in enumerate(points.items(), 1))
    with pytest.raises(SourceCompletionError, match="non-L"):
        plan_source(source().replace(b"END\n", extra+b"END\n"), TARGET)


@pytest.mark.parametrize("resid,missing", ((28, ("CB",)), (46, ("CB",)),
    (20, ("OG",)), (20, ("CB", "OG")), (57, ("OG",)), (57, ("CB", "OG")),
    (65, ("OG",)), (65, ("CB", "OG"))))
def test_real_explicit_ubiquitin_1_72_projection(resid, missing):
    original = (Path(__file__).resolve().parents[1]/"tests/fixtures/protein_context/1UBQ.pdb").read_bytes()
    # Explicit source projection, not occupancy rewriting or full-source admission.
    src = b"".join(l for l in original.splitlines(keepends=True) if l.startswith(b"ATOM  ")
                   and int(l[22:26]) <= 72 and not (int(l[22:26]) == resid and l[12:16].strip().decode() in missing))
    src += b"TER\nEND\n"
    target = ComponentInstanceRef(1, "A", str(resid))
    req = plan_source(src, target)
    assert set(req.missing_atoms) == set(missing)
    candidate = complete_sidechain(req).selected
    assert candidate is not None
    output = materialize_source(src, req, candidate)
    assert {a.atom.name for a in output.generated} == set(missing)
    # Mature output is a no-op request with no newly missing atoms.
    assert plan_source(output.pdb, target).missing_atoms == ()


def test_independent_output_audit_refuses_moved_observation_or_new_records():
    src = source()
    req = plan_source(src, TARGET)
    candidate = complete_sidechain(req).selected
    output = materialize_source(src, req, candidate).pdb
    for tampered in (output.replace(b"  -1.460", b"  -1.461"),
                     output.replace(b"END\n", b"REMARK injected\nEND\n"),
                     output.replace(b"  1.00  0.00", b"  1.00  1.00")):
        with pytest.raises(ValueError):
            audit_materialized_source(src, req, candidate, tampered)


def test_pdb_rounding_can_turn_a_raw_geometry_pass_into_a_refusal():
    # Fixed precision boundary control: raw CB--water distance 2.620063 A,
    # rounded CB--water distance 2.619994 A, crossing the -0.6 A clearance gate.
    water = atom(10, "O", (.589, -.600, 1.405), residue="HOH", chain="W", number="1", code="")
    src = source().replace(b"END\n", water+b"END\n")
    req = plan_source(src, TARGET)
    candidate = complete_sidechain(req).selected
    assert candidate is not None
    with pytest.raises(SourceCompletionError, match="rounded PDB rejected"):
        materialize_source(src, req, candidate)


def test_readback_refuses_shifted_generated_coordinates_even_if_geometry_passes():
    src = source()
    req = plan_source(src, TARGET)
    candidate = complete_sidechain(req).selected
    out = materialize_source(src, req, candidate).pdb
    changed = out.replace(b"   0.536", b"   0.537")
    assert changed != out
    with pytest.raises(SourceCompletionError, match="rounding bound"):
        audit_materialized_source(src, req, candidate, changed)


def test_conect_target_graph_is_not_silently_ignored():
    with pytest.raises(SourceCompletionError, match="noncanonical target/context bond"):
        plan_source(source().replace(b"END\n", b"CONECT    1    4\nEND\n"), TARGET)


def test_peptide_continuity_respects_ter_and_binding_includes_links():
    # Neighbour contains only a GLY backbone; global CA gate properly skips GLY.
    prev = atom(10, "C", (-2.10, -.70, .80), residue="GLY", number="81", code="")
    src = prev+source()
    req = plan_source(src, TARGET)
    assert len(req.links) == 1
    result = complete_sidechain(req)
    assert result.selected is not None
    with pytest.raises(SourceCompletionError, match="binding"):
        materialize_source(src, replace(req, links=()), result.selected)
    blocked = plan_source(prev+b"TER\n"+source(), TARGET)
    assert blocked.links == ()


def test_full_ubiquitin_source_is_not_silently_cropped_or_reoccupied():
    original = (Path(__file__).resolve().parents[1]/"tests/fixtures/protein_context/1UBQ.pdb").read_bytes()
    with pytest.raises(SourceCompletionError):
        plan_source(original, ComponentInstanceRef(1, "A", "28"))
    # Removing only the unsupported historical MASTER does not bypass occupancy.
    without_master = b"".join(l for l in original.splitlines(keepends=True) if not l.startswith(b"MASTER"))
    with pytest.raises(SourceCompletionError, match="occupancy"):
        plan_source(without_master, ComponentInstanceRef(1, "A", "28"))


@pytest.mark.parametrize("record", (b"ATOM\n", b"TERxxx\n"))
def test_malformed_short_structural_records_are_not_silently_ignored(record):
    with pytest.raises(SourceCompletionError):
        plan_source(source().replace(b"END\n", record+b"END\n"), TARGET)


def test_nonpeptide_conect_and_nearby_cysteine_are_refused():
    extra = atom(10, "C", (20., 0., 0.), residue="GLY", chain="B", number="1", code="")
    src = source().replace(b"END\n", extra+b"CONECT    1   10\nEND\n")
    with pytest.raises(SourceCompletionError, match="nonpeptide"):
        plan_source(src, TARGET)
    points = POINTS | {"CB": l_cb(POINTS["N"], POINTS["CA"], POINTS["C"]), "SG": (4., 0., 0.)}
    extra = b"".join(atom(i+10, n, tuple(x+3 for x in p), residue="CYS", chain="B", number="1", code="")
                     for i, (n, p) in enumerate(points.items(), 1))
    with pytest.raises(SourceCompletionError, match="disulfide"):
        plan_source(source("CYS").replace(b"END\n", extra+b"END\n"), TARGET)


def test_generated_record_type_and_serial_are_part_of_readback_contract():
    src = source()
    req = plan_source(src, TARGET)
    candidate = complete_sidechain(req).selected
    out = materialize_source(src, req, candidate).pdb
    for changed in (out.replace(b"ATOM      5", b"HETATM    5"), out.replace(b"ATOM      5", b"ATOM     15")):
        assert changed != out
        with pytest.raises(SourceCompletionError, match="generated record"):
            audit_materialized_source(src, req, candidate, changed)
