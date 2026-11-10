"""Synthetic format/interpolation/transport controls, NOT Dunbrack data or sampling."""
from __future__ import annotations

import hashlib
import math
from dataclasses import replace

import pytest

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.protein_completion.backbone import backbone_angles
from dvbfixer.protein_completion.contracts import CompletionRequest, Observation
from dvbfixer.protein_completion.geometry import l_cb
from dvbfixer.protein_completion.rotamers import (
    LibraryAuthority,
    RotamerLibraryError,
    load_single_chi,
    parse_single_chi,
)
from dvbfixer.protein_completion.service import complete_sidechain
from dvbfixer.protein_completion.validation import audit_candidate


def authority(category="SER", kind="synthetic-control"):
    return LibraryAuthority("synthetic-format-control-not-dunbrack", category,
                            "authored numeric test rows", kind, "synthetic-test",
                            hashlib.sha256(b"synthetic-control-no-external-data").hexdigest())


def row(phi, psi, label, probability, mean, sigma=10, category="SER"):
    return f"{category} {phi} {psi} 19 {label} 0 0 0 {probability} {mean} 0 0 0 {sigma} 0 0 0\n"


def table(cells=None, category="SER", means=(62, 178, -66), probabilities=(0.7, 0.2, 0.1)):
    cells = [(p, q) for p in range(-180, 180, 10) for q in range(-180, 180, 10)] if cells is None else cells
    return ("# Synthetic test data; not from any rotamer database.\n"
            + "".join(row(p, q, i+1, prob, mean, category=category)
                      for p, q in cells for i, (prob, mean) in enumerate(zip(probabilities, means)))).encode()


def parse(data, category="SER", kind="synthetic-control"):
    return parse_single_chi(data, expected_sha256=hashlib.sha256(data).hexdigest(),
                            authority=authority(category, kind))


@pytest.fixture(scope="module")
def library():
    return parse(table())


def empirical_request(name="SER", with_cb=False):
    component = ComponentInstanceRef(1, "a", "82", "A", occurrence=2)
    points = {"N": (-1.46, 0., 0.), "CA": (0., 0., 0.),
              "C": (0.55, 1.43, 0.), "O": (1.76, 1.66, 0.)}
    if with_cb:
        points["CB"] = l_cb(points["N"], points["CA"], points["C"])
    observed = tuple(Observation(ExactAtomRef(component, k), k[0], v) for k, v in points.items())
    previous = Observation(ExactAtomRef(ComponentInstanceRef(1, "a", "81"), "C"), "C", (-2.1, -0.7, 0.8))
    following = Observation(ExactAtomRef(ComponentInstanceRef(1, "a", "83"), "N"), "N", (0.0, 2.15, -0.9))
    links = ((previous.atom, ExactAtomRef(component, "N")), (ExactAtomRef(component, "C"), following.atom))
    missing = (("OG",) if with_cb else ("CB", "OG")) if name == "SER" else (("SG",) if with_cb else ("CB", "SG"))
    return CompletionRequest("b"*64, component, name, observed+(previous, following), missing, links)


def test_full_synthetic_table_is_bound_and_exact_nodes_reproduce_rows(library):
    assert len(library.cells) == 36*36
    assert library.authority.kind == "synthetic-control"
    for phi, psi in ((-60., -40.), (0., 0.), (170., -180.)):
        choices = library.lookup(phi, psi)
        assert [r.label for r in choices] == [1, 2, 3]
        assert [r.probability for r in choices] == pytest.approx([0.7, 0.2, 0.1])
        assert [r.mean_degrees for r in choices] == pytest.approx([62, 178, -66])
        assert [r.sigma_degrees for r in choices] == pytest.approx([10, 10, 10])


def test_local_loader_checks_bytes_before_parsing(tmp_path):
    data = table([(-60, -40)])
    path = tmp_path/"synthetic.lib"
    path.write_bytes(data)
    loaded = load_single_chi(path, expected_sha256=hashlib.sha256(data).hexdigest(), authority=authority())
    assert loaded.data_sha256 == hashlib.sha256(data).hexdigest()
    with pytest.raises(RotamerLibraryError, match="checksum"):
        load_single_chi(path, expected_sha256="0"*64, authority=authority())
    path.write_bytes(b"<html>not a library</html>")
    with pytest.raises(RotamerLibraryError, match="checksum"):
        load_single_chi(path, expected_sha256=loaded.data_sha256, authority=authority())
    with pytest.raises(RotamerLibraryError, match="regular file"):
        load_single_chi(tmp_path/"absent", expected_sha256="0"*64, authority=authority())


def test_sparse_data_does_not_acquire_full_external_authority():
    with pytest.raises(RotamerLibraryError, match="complete periodic grid"):
        parse(table([(-60, -40)]), kind="external-local")


def test_table_metadata_changes_model_identity(library):
    altered = replace(library, authority=replace(library.authority, library_id="another-library"))
    assert altered.data_sha256 == library.data_sha256
    assert altered.model_digest != library.model_digest
    with pytest.raises(TypeError):
        library.cells[(-60, -40)] = ()


def test_parsed_view_is_bound_even_if_raw_digest_is_reused(library):
    copied = dict(library.cells)
    rows = copied[(-60, -40)]
    copied[(-60, -40)] = (replace(rows[0], mean_degrees=63), *rows[1:])
    changed = replace(library, cells=copied)
    copied.clear()
    assert changed.data_sha256 == library.data_sha256
    assert changed.view_sha256 != library.view_sha256
    assert changed.model_digest != library.model_digest
    assert len(changed.cells) == 36*36


def test_zero_probability_rotamer_is_not_fabricated_with_epsilon():
    loaded = parse(table([(-60, -40)], probabilities=(1, 0, 0)))
    choices = loaded.lookup(-60, -40)
    assert len(choices) == 1 and choices[0].label == 1 and choices[0].probability == 1


def test_probability_weighted_circular_moments_are_explicit():
    data = b""
    for phi, p, mean in ((-60, .8, 0), (-50, .2, 90)):
        data += "".join(row(phi, -40, label, probability, mu)
                         for label, probability, mu in ((1, p, mean), (2, 1-p, 178), (3, 0, -66))).encode()
    first = next(r for r in parse(data).lookup(-55, -40) if r.label == 1)
    assert first.probability == pytest.approx(.5)
    assert first.mean_degrees == pytest.approx(math.degrees(math.atan2(.2, .8)))


def test_antipodal_mean_is_refused_instead_of_picking_an_arbitrary_direction():
    data = b""
    for phi, mean in ((-60, 0), (-50, 180)):
        data += "".join(row(phi, -40, i+1, p, (mean, 178, -66)[i])
                         for i, p in enumerate((.7, .2, .1))).encode()
    with pytest.raises(RotamerLibraryError, match="ambiguous circular mean"):
        parse(data).lookup(-55, -40)


def test_backbone_dependent_probabilities_are_bilinear_not_nearest_bin():
    data = b""
    for phi, psi, probs in ((-60, -40, (.8, .1, .1)), (-50, -40, (.4, .4, .2)),
                            (-60, -30, (.2, .5, .3)), (-50, -30, (.6, .2, .2))):
        data += "".join(row(phi, psi, i+1, p, (60, 180, -60)[i]) for i, p in enumerate(probs)).encode()
    choices = {r.label: r for r in parse(data).lookup(-55, -35)}
    assert choices[1].probability == pytest.approx(.5)
    assert choices[2].probability == pytest.approx(.3)
    assert choices[3].probability == pytest.approx(.2)


def test_circular_interpolation_crosses_plus_minus_180_without_zero():
    data = b""
    for phi, mean in ((170, 179), (-180, -179)):
        data += "".join(row(phi, 0, i+1, p, (60, mean, -60)[i]) for i, p in enumerate((.6, .3, .1))).encode()
    choices = {r.label: r for r in parse(data).lookup(175, 0)}
    assert abs(choices[2].mean_degrees) == pytest.approx(180)
    assert choices[2].sigma_degrees == pytest.approx(math.sqrt(101))


def test_boundary_aliases_match_but_conflicts_and_duplicates_refuse():
    a = table([(-180, -180)])
    b = table([(180, 180)])
    assert len(parse(a+b).cells) == 1
    with pytest.raises(RotamerLibraryError, match="periodic boundary"):
        parse(a+table([(180, 180)], means=(63, 178, -66)))
    with pytest.raises(RotamerLibraryError, match="duplicate"):
        parse(a+a)


@pytest.mark.parametrize("phi,psi", ((180., 180.), (-180., -180.), (540., -540.), (-900., 900.)))
def test_lookup_is_periodic_on_both_axes(library, phi, psi):
    assert library.lookup(phi, psi) == library.lookup(-180, -180)


@pytest.mark.parametrize("phi", (math.nan, math.inf, True, "-60"))
def test_bad_query_angles_refuse(library, phi):
    with pytest.raises(RotamerLibraryError, match="finite torsion"):
        library.lookup(phi, 0)


def test_missing_interpolation_support_never_substitutes_another_bin():
    lib = parse(table([(-60, -40)]))
    assert lib.lookup(-60, -40)
    with pytest.raises(RotamerLibraryError, match="missing interpolation cell"):
        lib.lookup(-59, -40)


@pytest.mark.parametrize("transform", (
    lambda d: d.replace(b"SER", b"CYS"),
    lambda d: d.replace(b"0.7", b"nan"),
    lambda d: d.replace(b"0.7", b"1.7"),
    lambda d: d.replace(b"0.7", b"0.6"),
    lambda d: d.replace(b"-60 -40", b"-61 -40"),
    lambda d: d.replace(b"-60 -40", b"-190 -40"),
    lambda d: d.replace(b" 10 0 0 0", b" -1 0 0 0"),
    lambda d: d.replace(b" 10 0 0 0", b" 0 0 0 0"),
    lambda d: d.replace(b" 19 1 0 0 0", b" 19 1 2 0 0"),
    lambda d: d.replace(b" 19 1 0 0 0", b" 19 4 0 0 0"),
    lambda d: d+b"SER -60 -40 invalid\n",
    lambda d: b"# only comments\n",
    lambda d: d.splitlines(keepends=True)[0]+d.splitlines(keepends=True)[1],
))
def test_malformed_or_incompatible_tables_fail_closed(transform):
    with pytest.raises(RotamerLibraryError):
        parse(transform(table([(-60, -40)])))


def test_limits_and_plain_file_requirement(tmp_path):
    with pytest.raises(RotamerLibraryError, match="empty or oversized"):
        parse(b"")
    with pytest.raises(RotamerLibraryError, match="ASCII"):
        parse(b"\xff")
    with pytest.raises(RotamerLibraryError, match="oversized rotamer line"):
        parse(b"#"+b"a"*4096)
    path = tmp_path/"data"
    path.write_bytes(table([(-60, -40)]))
    link = tmp_path/"link"
    link.symlink_to(path)
    with pytest.raises(RotamerLibraryError, match="regular file"):
        load_single_chi(link, expected_sha256="0"*64, authority=authority())


@pytest.mark.parametrize("with_cb", (False, True))
def test_empirical_service_transports_synthetic_means_prior_and_exact_atoms(library, with_cb):
    req = empirical_request(with_cb=with_cb)
    baseline = complete_sidechain(req)
    result = complete_sidechain(req, library=library)
    assert result.status == "completed"
    assert result.model_digest == library.model_digest
    assert result.selected.rotamer_evidence is not None
    assert result.selected.rotamer_evidence.probability == pytest.approx(.7)
    assert {a.atom.name for a in result.selected.additions} == set(req.missing_atoms)
    assert result.selected.torsion_degrees == pytest.approx(62)
    assert all(a.candidate.rotamer_evidence is None for a in baseline.attempts)
    assert audit_candidate(req, result.selected, library=library).accepted
    assert not audit_candidate(req, result.selected).accepted


def test_free_cysteine_requires_explicit_free_category_not_combined_cys(library):
    req = empirical_request("CYS")
    result = complete_sidechain(req, library=library)
    assert result.status == "unsupported" and result.reason == "rotamer_category_mismatch"
    free = parse(table(category="CYH"), category="CYH")
    assert complete_sidechain(req, library=free).status == "completed"
    with pytest.raises(RotamerLibraryError):
        authority(category="CYS")


def test_phi_psi_require_actual_links_not_numeric_adjacency(library):
    req = empirical_request()
    assert all(math.isfinite(v) for v in backbone_angles(req))
    for links in ((), req.links[:1], req.links[1:]):
        result = complete_sidechain(replace(req, links=links), library=library)
        assert result.status == "unsupported" and "peptide neighbours" in result.reason
        assert result.selected is None and not result.attempts


def test_library_shortage_refuses_instead_of_using_staggered_control():
    req = empirical_request()
    sparse = parse(table([(-60, -40)]))
    result = complete_sidechain(req, library=sparse)
    assert result.status == "unsupported" and "missing interpolation cell" in result.reason


@pytest.mark.parametrize("field,value", (("probability", .999), ("mean_degrees", 14),
                                         ("sigma_degrees", 179), ("phi_degrees", 0),
                                         ("model_digest", "0"*64), ("category", "CYH")))
def test_audit_recomputes_rotamer_evidence(library, field, value):
    req = empirical_request()
    candidate = complete_sidechain(req, library=library).selected
    forged = replace(candidate, rotamer_evidence=replace(candidate.rotamer_evidence, **{field: value}))
    assert not audit_candidate(req, forged, library=library).accepted


def test_empirical_audit_does_not_accept_staggered_candidate_or_other_library(library):
    req = empirical_request()
    baseline = complete_sidechain(req).selected
    assert not audit_candidate(req, baseline, library=library).accepted
    candidate = complete_sidechain(req, library=library).selected
    other = replace(library, authority=replace(library.authority, library_id="different-authority"))
    assert not audit_candidate(req, candidate, library=other).accepted


def test_likely_but_clashing_rotamer_never_beats_feasible_lower_prior(library):
    req = empirical_request()
    first = complete_sidechain(req, library=library).selected
    cb = next(a.xyz for a in first.additions if a.atom.name == "CB")
    og = next(a.xyz for a in first.additions if a.atom.name == "OG")
    outside = tuple(o+1.3*(o-c)/math.dist(og, cb) for o, c in zip(og, cb))
    blocker = Observation(ExactAtomRef(ComponentInstanceRef(1, "L", "1"), "O1"), "O", outside)
    changed = replace(req, observations=req.observations+(blocker,))
    result = complete_sidechain(changed, library=library)
    assert result.status == "completed"
    assert not result.attempts[0].audit.accepted
    assert result.attempts[0].candidate.rotamer_evidence.probability == pytest.approx(.7)
    assert result.selected.rotamer_evidence.probability < .7
