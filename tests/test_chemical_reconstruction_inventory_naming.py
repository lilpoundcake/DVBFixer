"""Inventory naming proofs do not qualify stereo-specific H coordinates."""
from __future__ import annotations

from dataclasses import replace

import pytest
from scripts import chemical_reconstruction_naming_controls as controls
from scripts.chemical_reconstruction_peptide_controls import inventory, load_ccd

from dvbfixer.chemical_reconstruction.inventory_naming import transport_inventory
from dvbfixer.chemical_reconstruction.peptide_state import PeptideStateError, derive_local_peptide
from dvbfixer.domain.force_field_naming import ForceFieldTarget, NamingProfile
from dvbfixer.domain.structure_identity import ComponentInstanceRef


def transport(raw):
    return transport_inventory(raw, target=ForceFieldTarget.AMBER19SB, profile=NamingProfile.STANDARD)


@pytest.mark.parametrize("name", ["ALA", "THR", "ILE", "LYS", "PHE"])
def test_declared_raw_inventories_match_joined_chemistry_without_mutating_ccd(name) -> None:
    raw = controls.raw_inventory(name)
    before = repr(raw)
    result = transport(raw)
    parent = load_ccd(name)
    authority = derive_local_peptide(parent, load_ccd("ALA"), ComponentInstanceRef(1, "a", "82", "B"), result.canonical)
    assert repr(raw) == before
    assert authority.parent == parent
    assert result.raw == raw
    assert len(result.mapping) == len(raw.atoms) == len(result.canonical.atoms)
    names = dict((a, b) for a, b, _ in result.mapping)
    assert {frozenset(names[n] for n in b) for b in raw.bonds} == result.canonical.bonds
    assert result.canonical.source_sha256 == raw.source_sha256
    assert authority.identity.chain_id == "a" and authority.identity.insertion_code == "B"
    assert {a.name for a in authority.state.graph.atoms} == result.canonical.atoms


def test_ile_uses_one_simultaneous_policy_and_old_refusal_is_unchanged() -> None:
    result = transport(controls.raw_inventory("ILE"))
    names = {a: b for a, b, _ in result.mapping}
    assert {a: names[a] for a in ("CD", "HD1", "HD2", "HD3", "HG11")} == {
        "CD": "CD1", "HD1": "HD11", "HD2": "HD12", "HD3": "HD13", "HG11": "HG13"}
    assert names["HG12"] == "HG12"
    with pytest.raises(PeptideStateError, match="inventory-mismatch"):
        derive_local_peptide(load_ccd("ILE"), load_ccd("ALA"), ComponentInstanceRef(1, "A", "2"), inventory("ILE"))


@pytest.mark.parametrize("aliases", [(('N', 'CA'),), (('N', 'CA'), ('N', 'C'))])
def test_prealiased_or_duplicate_alias_input_is_not_raw(aliases) -> None:
    with pytest.raises(PeptideStateError, match="not-raw"):
        transport(replace(controls.raw_inventory("ALA"), aliases=aliases))


def test_implicit_identity_collision_refuses() -> None:
    raw = controls.raw_inventory("ILE")
    with pytest.raises(PeptideStateError, match="naming-collision"):
        transport(replace(raw, atoms=raw.atoms | {"CD1"}))


@pytest.mark.parametrize("name", ["N\u0000", "N CA", "-C", "Å", "N/CA"])
def test_invalid_atom_names_refuse(name) -> None:
    raw = controls.raw_inventory("ALA")
    with pytest.raises(PeptideStateError, match="atom-name"):
        transport(replace(raw, atoms=raw.atoms | {name}))


@pytest.mark.parametrize("token", ["+CD", "--C", "+", "+C ", "", "absent"])
def test_unknown_or_unqualified_external_references_refuse(token) -> None:
    raw = controls.raw_inventory("ILE")
    with pytest.raises(PeptideStateError):
        transport(replace(raw, external_tokens=(("N", token),)))


@pytest.mark.parametrize("rows", [[("-C", "N")], (("-C",),), (("N", "CA"),), (["-C", "N"],)])
def test_external_row_shape_and_immutability_refuse(rows) -> None:
    with pytest.raises(PeptideStateError):
        transport(replace(controls.raw_inventory("ALA"), external_tokens=rows))


def test_external_signs_order_and_repeated_rows_are_preserved() -> None:
    raw = controls.raw_inventory("ILE")
    rows = (("-C", "CA", "N", "HG11"), ("CA", "+N", "C", "O"), ("-C", "N"), ("-C", "N"))
    result = transport(replace(raw, external_tokens=rows))
    assert result.canonical.external_tokens == (
        ("-C", "CA", "N", "HG13"), ("CA", "+N", "C", "O"), ("-C", "N"), ("-C", "N"))


@pytest.mark.parametrize("field", ["mapping", "canonical", "raw"])
def test_transport_record_binds_mapping_inventory_and_source_pin(field) -> None:
    result = transport(controls.raw_inventory("ALA"))
    changes = {
        "mapping": result.mapping[:-1],
        "canonical": replace(result.canonical, source_sha256="1"*64),
        "raw": replace(result.raw, source_sha256="2"*64),
    }
    with pytest.raises(PeptideStateError, match="binding-mismatch"):
        replace(result, **{field: changes[field]})


def test_naming_cannot_hide_a_connectivity_mismatch_or_choose_a_variant() -> None:
    raw = controls.raw_inventory("ALA")
    bad = transport(replace(raw, bonds=raw.bonds - {frozenset(("CA", "CB"))}))
    with pytest.raises(PeptideStateError, match="inventory-mismatch"):
        derive_local_peptide(load_ccd("ALA"), load_ccd("ALA"), ComponentInstanceRef(1, "A", "2"), bad.canonical)
    with pytest.raises(PeptideStateError, match="variant-projection-unqualified"):
        derive_local_peptide(load_ccd("LYS"), load_ccd("ALA"), ComponentInstanceRef(1, "A", "2"),
                             transport(controls.raw_inventory("LYN")).canonical)


@pytest.mark.parametrize("field,value", [("source", " "), ("state_name", "ILE "), ("external_tokens", (("-C", None),))])
def test_invalid_metadata_refuses(field, value) -> None:
    with pytest.raises(PeptideStateError):
        transport(replace(controls.raw_inventory("ALA"), **{field: value}))


def test_no_implicit_profile_or_other_dialect_qualification() -> None:
    raw = controls.raw_inventory("ALA")
    with pytest.raises(PeptideStateError, match="policy-unqualified"):
        transport_inventory(raw, target=ForceFieldTarget.CHARMM, profile=NamingProfile.STANDARD)
    with pytest.raises(PeptideStateError, match="policy-unqualified"):
        transport_inventory(raw, target=ForceFieldTarget.AMBER19SB, profile=NamingProfile.GROMACS)


def test_pin_failure_precedes_output_reservation(tmp_path, monkeypatch) -> None:
    altered = tmp_path / "owner.py"
    altered.write_text("not the pinned naming owner")
    monkeypatch.setattr(controls, "OWNER", altered)
    output = tmp_path / "evidence"
    with pytest.raises(ValueError, match="pin mismatch"):
        controls.run_controls(output)
    assert not output.exists()


def test_raw_rtp_pin_failure(tmp_path, monkeypatch) -> None:
    altered = tmp_path / "input.rtp"
    altered.write_text("not the pinned RTP")
    monkeypatch.setattr(controls, "RTP", altered)
    with pytest.raises(ValueError, match="pin"):
        controls.raw_inventory("ALA")


def test_six_frozen_controls_and_no_overwrite(tmp_path) -> None:
    output = tmp_path / "controls"
    result = controls.run_controls(output)
    assert result["declared"] == result["executed_entries"] == 6
    assert [c["status"] for c in result["controls"]] == ["conditional-inventory-proof-passed"]*5 + ["refused"]
    assert result["controls"][-1]["reason"] == "explicit-variant-projection-unqualified"
    assert result["coordinate_constructor_executions"] == 0
    assert not result["h_geometry_qualified"] and not result["production_enabled"]
    with pytest.raises(FileExistsError):
        controls.run_controls(output)
