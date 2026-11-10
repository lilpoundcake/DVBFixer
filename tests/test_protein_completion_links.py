"""Independent fixed-column and chemical-admission checks for link metadata."""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

from dvbfixer.domain.structure_identity import ComponentInstanceRef
from dvbfixer.pdb_component_reconstruction import _format_generated_atom, _parse_atom_lines
from dvbfixer.protein_completion.geometry import l_cb
from dvbfixer.protein_completion.links import parse_declared_links
from dvbfixer.protein_completion.service import complete_sidechain
from dvbfixer.protein_completion.source import materialize_source, plan_source

TARGET = ComponentInstanceRef(1, "a", "82", "A")
BACKBONE = {"N": (-1.46, 0., 0.), "CA": (0., 0., 0.), "C": (.55, 1.43, 0.), "O": (1.76, 1.66, 0.)}


def atom(serial, name, point, *, chain="a", number="82", code="A", residue="ALA"):
    text = _format_generated_atom(serial=serial, name=name, element=name[0], position=point,
                                  residue_name=residue, chain_id=chain, sequence_number=number,
                                  insertion_code=code)
    return "ATOM  "+text[6:]


def source():
    return "".join(atom(i, n, p) for i, (n, p) in enumerate(BACKBONE.items(), 1))+"TER\nEND\n"


def link(left=("N", "ALA", "a", "82", "A"), right=("CA", "ALA", "a", "82", "A"),
         symmetry="1555", distance=""):
    fields = list(" "*80)
    fields[:6] = "LINK  "
    for data, starts in ((left, (12, 17, 21, 22, 26)), (right, (42, 47, 51, 52, 56))):
        for value, start, width in zip(data, starts, (4, 3, 1, 4, 1)):
            fields[start:start+width] = f"{value:>{width}}"
    fields[59:65] = f"{symmetry:>6}"
    fields[66:72] = f"{symmetry:>6}"
    fields[73:78] = f"{distance:>5}"
    return "".join(fields)+"\n"


def ssbond(left=("B", "1", ""), right=("b", "1", "A"), symmetry="1555", distance="2.03"):
    fields = list(" "*80)
    fields[:6] = "SSBOND"
    fields[7:10] = "  1"
    fields[11:14] = fields[25:28] = "CYS"
    for data, starts in ((left, (15, 17, 21)), (right, (29, 31, 35))):
        for value, start, width in zip(data, starts, (1, 4, 1)):
            fields[start:start+width] = f"{value:>{width}}"
    fields[59:65] = f"{symmetry:>6}"
    fields[66:72] = f"{symmetry:>6}"
    fields[73:78] = f"{distance:>5}"
    return "".join(fields)+"\n"


def context_pair():
    cb = l_cb(BACKBONE["N"], BACKBONE["CA"], BACKBONE["C"])
    points = BACKBONE | {"CB": cb, "SG": (2., -1., -2.)}
    result = ""
    for j, (chain, code, delta) in enumerate((("B", "", 0.), ("b", "A", 2.03))):
        for i, (name, xyz) in enumerate(points.items(), 1):
            result += atom(10+j*10+i, name, (xyz[0]+30+delta, xyz[1], xyz[2]),
                           chain=chain, number="1", code=code, residue="CYS")
        result += "TER\n"
    return result


@pytest.mark.parametrize("metadata", ("CONECT    1    2\n", link(), link(symmetry=""),
                                      "CONECT    1    2\nCONECT    2    1\n"))
def test_canonical_target_graph_metadata_is_preserved(metadata):
    src = source().replace("END\n", metadata+"END\n").encode()
    request = plan_source(src, TARGET)
    result = complete_sidechain(request)
    out = materialize_source(src, request, result.selected)
    assert metadata.encode() in out.pdb
    assert request.links == ()  # Internal links never become external exclusions.


@pytest.mark.parametrize("metadata,match", (
    (link(right=("O", "ALA", "a", "82", "A")), "noncanonical"),
    (link(right=("CB", "ALA", "a", "82", "A")), "endpoint"),
    (link(right=("CA", "ALA", "A", "82", "A")), "endpoint"),
    (link(right=("CA", "ALA", "a", "82", "B")), "endpoint"),
    (link(right=("CA", "SER", "a", "82", "A")), "endpoint"),
    (link(symmetry="2555"), "symmetry"),
    (link(distance="9.99"), "disagrees"),
    (link(distance="nan"), "distance"),
    (link(distance="-1"), "distance"),
    (link(distance="oops"), "distance"),
    ("CONECT    1    2    2\n", "multiplicity"),
    (link(right=("N", "ALA", "a", "82", "A")), "self"),
    ("LINKxx\n", "prefix"),
    ("SSBONDx\n", "prefix"),
    ("LINK\n", "truncated"),
))
def test_link_metadata_cannot_override_chemistry_or_identity(metadata, match):
    with pytest.raises(ValueError, match=match):
        plan_source(source().replace("END\n", metadata+"END\n").encode(), TARGET)


@pytest.mark.parametrize("kind", ("SSBOND", "LINK", "CONECT", "all"))
def test_valid_observed_context_disulfide_is_preserved_and_bound(kind):
    records = {"SSBOND": ssbond(),
               "LINK": link(("SG", "CYS", "B", "1", ""), ("SG", "CYS", "b", "1", "A"), distance="2.03"),
               "CONECT": "CONECT   16   26\nCONECT   26   16\n"}
    metadata = "".join(records.values()) if kind == "all" else records[kind]
    src = source().replace("END\n", context_pair()+metadata+"END\n").encode()
    req = plan_source(src, TARGET)
    assert len(req.links) == 1
    assert {a.component.chain_id for a in req.links[0]} == {"B", "b"}
    assert {a.component.insertion_code for a in req.links[0]} == {"", "A"}
    result = complete_sidechain(req)
    out = materialize_source(src, req, result.selected)
    assert metadata.encode() in out.pdb
    with pytest.raises(ValueError, match="binding"):
        materialize_source(src, replace(req, links=()), result.selected)


@pytest.mark.parametrize("metadata,match", ((ssbond(symmetry="2555"), "symmetry"),
    (ssbond(distance="1.0"), "disagrees"), (ssbond(right=("b", "1", "")), "endpoint"),
    (ssbond().replace("CYS", "SER", 1), "canonical CYS"),
    (ssbond().replace("  1", "xxx", 1), "serial")))
def test_ssbond_bad_symmetry_geometry_or_endpoints_refuse(metadata, match):
    with pytest.raises(ValueError, match=match):
        plan_source(source().replace("END\n", context_pair()+metadata+"END\n").encode(), TARGET)


def test_disulfide_with_two_partners_refuses_across_record_types():
    src = source().replace("END\n", context_pair()+"END\n")
    # Another complete CYS beside the same B/SG, without a repeated identity.
    extra = context_pair().split("TER\n")[0].replace("CYS B", "CYS C")
    lines = []
    for i, line in enumerate(extra.splitlines(keepends=True), 31):
        lines.append(line[:6]+f"{i:5d}"+line[11:])
    metadata = ssbond()+link(("SG", "CYS", "B", "1", ""), ("SG", "CYS", "C", "1", ""))
    # Move new C/SG to +2.00 A from B/SG, retaining L CA and other atom bytes.
    for i, line in enumerate(lines):
        if line[12:16].strip() == "SG":
            lines[i] = line[:30]+f"{34.:8.3f}"+line[38:]
    with pytest.raises(ValueError, match="multiple disulfide"):
        plan_source(src.replace("END\n", "".join(lines)+metadata+"END\n").encode(), TARGET)


def test_parser_does_not_silently_ignore_link_altloc():
    metadata = link()
    metadata = metadata[:16]+"A"+metadata[17:]
    lines = (source()+metadata).splitlines(keepends=True)
    with pytest.raises(ValueError, match="altloc"):
        parse_declared_links(lines, _parse_atom_lines(lines))


def test_link_peptide_must_agree_with_same_chain_file_order_and_ter():
    prev = atom(10, "C", (-2.10, -.70, .80), residue="GLY", number="81", code="")
    metadata = link(("C", "GLY", "a", "81", ""), ("N", "ALA", "a", "82", "A"))
    src = (prev+source()).replace("END\n", metadata+"END\n").encode()
    req = plan_source(src, TARGET)
    assert len(req.links) == 1
    out = materialize_source(src, req, complete_sidechain(req).selected)
    assert metadata.encode() in out.pdb
    blocked = (prev+"TER\n"+source()).replace("END\n", metadata+"END\n").encode()
    with pytest.raises(ValueError, match="nonpeptide"):
        plan_source(blocked, TARGET)


def test_disulfide_target_refuses_even_when_atoms_are_complete():
    src = source().replace("END\n", context_pair()+ssbond()+"END\n").encode()
    with pytest.raises(ValueError, match="disulfide-linked target"):
        plan_source(src, ComponentInstanceRef(1, "B", "1"))


def test_absent_reported_length_does_not_bypass_actual_disulfide_geometry():
    context = context_pair()
    lines = context.splitlines(keepends=True)
    for i, line in enumerate(lines):
        if line.startswith("ATOM  ") and int(line[6:11]) == 26:
            lines[i] = line[:30]+f"{45.:8.3f}"+line[38:]
    src = source().replace("END\n", "".join(lines)+ssbond(distance="")+"END\n").encode()
    with pytest.raises(ValueError, match="observed disulfide distance"):
        plan_source(src, TARGET)


def test_missing_disulfide_endpoint_is_not_generated_as_free_cysteine():
    context = "".join(l for l in context_pair().splitlines(keepends=True)
                      if not (l.startswith("ATOM  ") and int(l[6:11]) == 26))
    src = source().replace("END\n", context+ssbond()+"END\n").encode()
    with pytest.raises(ValueError, match="endpoint"):
        plan_source(src, ComponentInstanceRef(1, "b", "1", "A"))


def test_real_1dqj_explicit_ala9_cys23_cys88_projection():
    original = (Path(__file__).resolve().parents[1]/"tests/fixtures/regressions/1DQJ.pdb").read_text()
    lines = original.splitlines(keepends=True)
    metadata = next(l for l in lines if l.startswith("SSBOND   1"))
    parts = [metadata]
    for number in (9, 23, 88):
        parts.extend(l for l in lines if l.startswith("ATOM  ") and l[21] == "A" and int(l[22:26]) == number
                     and not (number == 9 and l[12:16].strip() == "CB"))
        parts.append("TER\n")
    src = ("".join(parts)+"END\n").encode()
    target = ComponentInstanceRef(1, "A", "9")
    req = plan_source(src, target)
    assert len(req.links) == 1 and {a.component.sequence_number for a in req.links[0]} == {"23", "88"}
    candidate = complete_sidechain(req).selected
    assert candidate is not None
    result = materialize_source(src, req, candidate)
    assert metadata.encode() in result.pdb
    assert [a.atom.name for a in result.generated] == ["CB"]
