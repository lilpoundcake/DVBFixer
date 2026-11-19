"""Source policy/binding regressions; replays are not new scientific trials."""
from __future__ import annotations

import subprocess
import sys
from dataclasses import fields, replace

import pytest
from scripts.chemical_reconstruction_naming_controls import ROOT, raw_inventory
from scripts.chemical_reconstruction_peptide_controls import POLICY, load_ccd
from scripts.chemical_reconstruction_source_controls import run_controls

from dvbfixer.chemical_reconstruction.admission import admit
from dvbfixer.chemical_reconstruction.construction import propose
from dvbfixer.chemical_reconstruction.contracts import (
    ChemistryKind,
    Connection,
    Constraint,
    ConstraintKind,
    LinkKind,
    ReconstructionRequest,
    Unit,
    digest_record,
)
from dvbfixer.chemical_reconstruction.inventory_naming import transport_inventory
from dvbfixer.chemical_reconstruction.observed_peptide import (
    compile_peptide_plan,
    prepare_observed_probe,
)
from dvbfixer.chemical_reconstruction.peptide_state import derive_local_peptide
from dvbfixer.chemical_reconstruction.source_preparation import (
    audit_source_transport,
    prepare_source,
    prepare_source_probe,
)
from dvbfixer.domain.force_field_naming import ForceFieldTarget, NamingProfile
from dvbfixer.domain.nonprotein_reconstruction import ExternalLink
from dvbfixer.domain.structure_identity import ComponentInstanceRef

MASK = frozenset({"CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"})


def source():
    return (ROOT / "tests/fixtures/protein_context/1UBQ.pdb").read_bytes()


def authority():
    inventory = transport_inventory(raw_inventory("PHE"), target=ForceFieldTarget.AMBER19SB, profile=NamingProfile.STANDARD)
    return derive_local_peptide(load_ccd("PHE"), load_ccd("ALA"), ComponentInstanceRef(1, "A", "45"), inventory.canonical)


def probe():
    return prepare_source_probe(prepare_source(source()), authority(), withheld=MASK)


def add_records(content, records):
    lines = content.decode().splitlines(keepends=True)
    index = next(i for i, line in enumerate(lines) if line[:6].strip() == "END")
    lines[index:index] = records
    return "".join(lines).encode()


def test_mixed_conect_projection_drops_whole_dangling_record_and_preserves_every_other_byte() -> None:
    valid = "CONECT    1    2\n"
    bad = "CONECT    1    299999\n"
    original = add_records(source(), [valid, bad, "CONECT\n", "CONECT    1\n"])
    prepared = prepare_source(original)
    assert len(prepared.removed_metadata) == 3
    assert prepared.prepared_text.encode() == add_records(source(), [valid])
    assert prepared.original_text.encode() == original
    assert [r.reason for r in prepared.removed_metadata] == ["dangling-CONECT", "empty-CONECT", "empty-CONECT"]
    assert prepared.removed_metadata[0].original_line == bad
    for r in prepared.removed_metadata:
        assert original.decode().splitlines(keepends=True)[r.line_number-1] == r.original_line


@pytest.mark.parametrize("record", ["CONECT    1    X\n", "CONECT    1   2\n", "CONECT    1    1\n"])
def test_malformed_or_self_conect_refuses_instead_of_sanitizing(record) -> None:
    with pytest.raises(ValueError):
        prepare_source(add_records(source(), [record]))


@pytest.mark.parametrize("field", ["prepared_text", "removed_metadata"])
def test_preparation_record_recomputes_exact_projection(field) -> None:
    prepared = prepare_source(add_records(source(), ["CONECT    199999\n"]))
    changes = {"prepared_text": prepared.prepared_text.replace("HEADER", "REMARK", 1), "removed_metadata": ()}
    with pytest.raises(ValueError, match="projection-binding"):
        replace(prepared, **{field: changes[field]})


def test_fractional_context_is_exact_retained_and_conservatively_admitted() -> None:
    prepared = prepare_source(source())
    assert prepared.original_text == prepared.prepared_text and not prepared.removed_metadata
    p = prepare_source_probe(prepared, authority(), withheld=MASK)
    r = p.request
    assert len(r.context_policy.fractional_atoms) == 57
    assert len(r.observations) == 653 and p.source_atom_count == 660
    assert all(o.occupancy < 1 for o in r.observations if o.identity in r.context_policy.fractional_atoms)
    result = propose(r, compile_peptide_plan(p, POLICY))
    assert result.status == "declared-geometry-passed"
    assert result.raw_audit.passed and result.rounded_audit.passed
    assert not result.geometry_approved and not result.source_coverage_verified and not result.md_ready
    assert admit(r).unassessed_context


def test_original_request_and_observed_adapter_keep_their_strict_occupancy_refusal() -> None:
    a = authority()
    old = prepare_observed_probe(source(), a, withheld=MASK)
    assert "partial-occupancy-or-altloc-requires-explicit-selection" in {f.code for f in admit(old.request).findings}
    r = probe().request
    plain = ReconstructionRequest(**{f.name: getattr(r, f.name) for f in fields(ReconstructionRequest)})
    assert plain.digest != r.digest
    assert "partial-occupancy-or-altloc-requires-explicit-selection" in {f.code for f in admit(plain).findings}


@pytest.mark.parametrize("field", ["source", "evidence", "scope"])
def test_prepared_context_policy_binds_exact_source_evidence_and_scope(field) -> None:
    r = probe().request
    policy = r.context_policy
    changes = {"source": {"prepared_source_sha256": "0"*64}, "evidence": {"preparation_digest": "0"*64},
               "scope": {"fractional_atoms": frozenset()}}
    with pytest.raises(ValueError):
        replace(r, context_policy=replace(policy, **changes[field]))


@pytest.mark.parametrize("role", ["constraint", "connection", "mobile", "resolved"])
def test_fractional_context_must_not_become_a_geometry_or_chemistry_authority(role) -> None:
    r = probe().request
    ref = next(iter(r.context_policy.fractional_atoms))
    target_n = next(o.identity for o in r.observations if o.identity.component in r.targets and o.identity.name == "N")
    with pytest.raises(ValueError):
        if role == "constraint":
            q = Constraint(ConstraintKind.DISTANCE, (ref, target_n), Unit.ANGSTROM, 0., 1000.)
            replace(r, constraints=(*r.constraints, q))
        elif role == "connection":
            q = Connection(ExternalLink(ref, target_n, 1.), LinkKind.COVALENT)
            replace(r, connections=(*r.connections, q))
        elif role == "mobile":
            replace(r, ownership=replace(r.ownership, mobile=frozenset({ref})))
        else:
            state = next(c.resolved for c in r.components if c.resolved is not None)
            components = tuple(replace(c, kind=ChemistryKind.PROTEIN, states=(state,)) if c.identity == ref.component else c
                               for c in r.components)
            replace(r, components=components)


@pytest.mark.parametrize("occupancy", ["  0.50", "  0.00"])
def test_withholding_cannot_hide_fractional_target_occupancy(occupancy) -> None:
    lines = source().decode().splitlines(keepends=True)
    i = next(i for i, line in enumerate(lines) if line.startswith("ATOM") and line[22:26].strip() == "45" and line[12:16].strip() == "CB")
    lines[i] = lines[i][:54]+occupancy+lines[i][60:]
    prepared = prepare_source("".join(lines).encode())
    with pytest.raises(ValueError, match="target-partial-occupancy"):
        prepare_source_probe(prepared, authority(), withheld=MASK)


def test_fractional_context_is_not_discounted_out_of_clash_measurement() -> None:
    p = probe()
    r = p.request
    h = propose(r, compile_peptide_plan(p, POLICY)).raw.additions[0].position
    ref = next(iter(r.context_policy.fractional_atoms))
    lines = source().decode().splitlines(keepends=True)
    i = next(i for i, line in enumerate(lines) if line.startswith(("ATOM  ", "HETATM"))
             and line[21].strip() == ref.component.chain_id and line[22:26].strip() == ref.component.sequence_number
             and line[12:16].strip() == ref.name)
    lines[i] = lines[i][:30]+f"{h[0]+.2:8.3f}{h[1]:8.3f}{h[2]:8.3f}"+lines[i][54:]
    prepared = prepare_source("".join(lines).encode())
    changed = prepare_source_probe(prepared, authority(), withheld=MASK)
    assert audit_source_transport(prepared, changed).passed
    r = changed.request
    result = propose(r, compile_peptide_plan(changed, POLICY))
    assert result.status == "rejected"
    assert "generated-context-minimum-separation" in result.raw_audit.findings
    assert next(o.occupancy for o in r.observations if o.identity == ref) < 1


def test_observation_digest_and_independent_source_audit_reject_stale_or_forged_transport() -> None:
    prepared = prepare_source(source())
    p = prepare_source_probe(prepared, authority(), withheld=MASK)
    r = p.request
    observations = (replace(r.observations[0], position=(10., 20., 30.)), *r.observations[1:])
    with pytest.raises(ValueError, match="observation-binding"):
        replace(r, observations=observations)
    # Even an explicitly rehashed DTO is not proof of faithful source transport.
    forged = replace(r, observations=observations,
                      context_policy=replace(r.context_policy, observations_digest=digest_record(observations)))
    report = audit_source_transport(prepared, replace(p, request=forged))
    assert not report.passed and "source-observation-transport" in report.findings


def test_withheld_reference_position_cannot_drive_the_constructor() -> None:
    p = probe()
    result = propose(p.request, compile_peptide_plan(p, POLICY))
    lines = source().decode().splitlines(keepends=True)
    for i, line in enumerate(lines):
        if line.startswith("ATOM") and line[22:26].strip() == "45" and line[12:16].strip() in MASK:
            lines[i] = line[:30]+f"{900.:8.3f}{900.:8.3f}{900.:8.3f}"+line[54:]
    other = prepare_source_probe(prepare_source("".join(lines).encode()), authority(), withheld=MASK)
    assert other.request.source_sha256 != p.request.source_sha256
    observed = propose(other.request, compile_peptide_plan(other, POLICY))
    assert observed.status == "declared-geometry-passed"
    assert observed.raw.additions == result.raw.additions


def test_altloc_and_multiple_models_require_explicit_selection() -> None:
    lines = source().decode().splitlines(keepends=True)
    i = next(i for i, line in enumerate(lines) if line.startswith("ATOM"))
    lines[i] = lines[i][:16]+"A"+lines[i][17:]
    with pytest.raises(ValueError, match="altloc"):
        prepare_source("".join(lines).encode())
    atom = next(line for line in lines if line.startswith("ATOM"))
    with pytest.raises(ValueError):
        prepare_source(("MODEL        1\n"+atom+"ENDMDL\nMODEL        2\n"+atom+"ENDMDL\n").encode())


def test_three_frozen_controls_preserve_files_and_nonoverwriting_evidence(tmp_path) -> None:
    output = tmp_path / "controls"
    result = run_controls(output)
    assert result["declared"] == result["executed_entries"] == 3
    assert [r["status"] for r in result["controls"]] == ["declared-geometry-passed"]*2+["refused"]
    assert result["controls"][-1]["reason"] == "explicit-altloc-selection-required"
    assert not result["full_complex_or_h_qualification"]
    for c in result["controls"]:
        assert (output / c["id"] / "original.pdb").read_bytes() == (ROOT / c["source"]).read_bytes()
    with pytest.raises(FileExistsError):
        run_controls(output)


def test_conect_owner_move_preserves_legacy_alias_without_importing_cli_for_source_preparation() -> None:
    from dvbfixer.pdbutils.io import update_conect
    from dvbfixer.renumber import update_conect as legacy
    assert legacy is update_conect
    assert legacy("CONECT    1    299999\n", {1: 10, 2: 20}) is None
    assert legacy("CONECT    1    2\n", {1: 10, 2: 20}) == "CONECT   10   20\n"
    code = "import sys; import dvbfixer.chemical_reconstruction.source_preparation; assert 'dvbfixer.renumber' not in sys.modules"
    subprocess.run([sys.executable, "-c", code], check=True)
