from __future__ import annotations

import json
from dataclasses import asdict

import pytest

from dvbfixer.domain.force_field_naming import (
    ForceFieldTarget,
    NamingConversionError,
    NamingErrorCode,
    NamingProfile,
    NamingRuleId,
    ResidueIdentity,
    VariantOverride,
)
from dvbfixer.force_field_naming import (
    NamingConversionRequest,
    convert_force_field_naming,
)


def _atom(
    serial: int,
    atom: str,
    residue: str,
    chain: str = "A",
    number: int = 1,
    *,
    icode: str = "",
    altloc: str = "",
    record: str = "ATOM",
    suffix: str = "      10.123  11.234  12.345  0.75 19.25           C  ",
) -> str:
    atom_field = atom[:4] if len(atom) == 4 else f" {atom:<3}"
    residue_field = residue if len(residue) == 4 else f"{residue:<3} "
    return (
        f"{record:<6}{serial:5d} {atom_field}{altloc or ' '}{residue_field}"
        f"{chain}{number:4d}{icode or ' '}{suffix}\n"
    )


def _convert(
    text: str,
    *,
    target: ForceFieldTarget = ForceFieldTarget.CHARMM,
    profile: NamingProfile = NamingProfile.STANDARD,
    overrides: dict[ResidueIdentity, str] | None = None,
):
    return convert_force_field_naming(
        NamingConversionRequest(
            pdb_text=text,
            target=target,
            profile=profile,
            variant_overrides=tuple(
                VariantOverride(identity, variant)
                for identity, variant in (overrides or {}).items()
            ),
        )
    )


def test_noop_has_empty_structured_report_and_model_one() -> None:
    text = "HEADER unchanged\n" + _atom(1, "CA", "HIS") + "END\n"
    result = _convert(text)
    assert result.pdb_text == text
    assert result.model == 1
    assert result.changes == ()
    assert result.diagnostics == ()
    assert result.summary.changed_atoms == 0
    assert result.summary.changed_residues == 0
    assert result.summary.variant_changes == 0
    assert result.summary.atom_name_changes == 0


def test_structured_result_is_json_serializable_for_future_adapters() -> None:
    result = _convert(
        _atom(1, "N", "HIS"),
        overrides={ResidueIdentity("A", "1", ""): "HIE"},
    )
    payload = json.loads(json.dumps(asdict(result)))
    assert payload["changes"][0]["rule_ids"] == ["explicit-variant"]
    assert payload["summary"]["variant_changes"] == 1


def test_explicit_variant_precedes_source_and_his_is_not_inferred() -> None:
    text = _atom(1, "N", "HIE", number=1) + _atom(2, "N", "HIS", number=2)
    result = _convert(
        text,
        overrides={ResidueIdentity("A", "1", ""): "HID"},
    )
    assert " HSD A   1" in result.pdb_text
    assert " HIS A   2" in result.pdb_text
    assert result.changes[0].rule_ids == (NamingRuleId.EXPLICIT_VARIANT,)


def test_insertion_code_siblings_and_chain_case_are_isolated() -> None:
    text = (
        _atom(1, "N", "HIS", "A", 7)
        + _atom(2, "N", "HIS", "A", 7, icode="B")
        + _atom(3, "N", "HIS", "a", 7)
    )
    result = _convert(
        text,
        overrides={ResidueIdentity("A", "7", ""): "HIE"},
    )
    lines = result.pdb_text.splitlines()
    assert _resname(lines[0]) == "HSE"
    assert _resname(lines[1]) == "HIS"
    assert _resname(lines[2]) == "HIS"


def _resname(line: str) -> str:
    return line[17:21].strip() if line[20:21].strip() else line[17:20].strip()


def test_preflight_rejects_target_atom_collision() -> None:
    text = _atom(1, "HB1", "LYS") + _atom(2, "HB3", "LYS")
    with pytest.raises(NamingConversionError) as caught:
        _convert(text, target=ForceFieldTarget.AMBER, profile=NamingProfile.GROMACS)
    assert caught.value.code is NamingErrorCode.NAMING_COLLISION
    assert caught.value.details["source_atom_names"] == ["HB1", "HB3"]


def test_preflight_rejects_duplicate_source_atom_identity() -> None:
    text = _atom(1, "HB3", "SER") + _atom(2, "HB3", "SER")
    with pytest.raises(NamingConversionError) as caught:
        _convert(text, target=ForceFieldTarget.AMBER, profile=NamingProfile.GROMACS)
    assert caught.value.code is NamingErrorCode.NAMING_COLLISION
    assert caught.value.details["atom_name"] == "HB3"


def test_duplicate_coordinate_serial_is_malformed() -> None:
    text = _atom(1, "N", "ALA") + _atom(1, "CA", "ALA")
    with pytest.raises(NamingConversionError) as caught:
        _convert(text)
    assert caught.value.code is NamingErrorCode.MALFORMED_PDB


def test_multiple_models_are_rejected_before_conversion() -> None:
    text = "MODEL        1\nENDMDL\nMODEL        2\nENDMDL\n"
    with pytest.raises(NamingConversionError) as caught:
        _convert(text)
    assert caught.value.code is NamingErrorCode.MULTI_MODEL_UNSUPPORTED


def test_single_model_number_is_reported() -> None:
    result = _convert("MODEL        7\n" + _atom(1, "N", "ALA") + "ENDMDL\n")
    assert result.model == 7


@pytest.mark.parametrize(("source", "target"), [("ASH", "ASPP"), ("GLH", "GLUP")])
def test_four_character_charmm_names_parse_and_render_idempotently(
    source: str, target: str
) -> None:
    overrides = {ResidueIdentity("A", "1", ""): source}
    first = _convert(_atom(1, "N", source), overrides=overrides)
    assert _resname(first.pdb_text) == target
    second = _convert(first.pdb_text, overrides=overrides)
    assert second.pdb_text == first.pdb_text
    assert second.changes == ()


def test_lyn_to_lsn_pair_is_unique_and_repeated_conversion_is_noop() -> None:
    text = _atom(1, "HZ2", "LYN") + _atom(2, "HZ3", "LYN")
    overrides = {ResidueIdentity("A", "1", ""): "LYN"}
    first = _convert(text, overrides=overrides)
    names = [line[12:16].strip() for line in first.pdb_text.splitlines()]
    assert names == ["HZ1", "HZ2"]
    second = _convert(first.pdb_text, overrides=overrides)
    assert second.pdb_text == first.pdb_text
    assert second.changes == ()


def test_already_target_lyn_pair_is_not_shifted_again() -> None:
    text = _atom(1, "HZ1", "LYN") + _atom(2, "HZ2", "LYN")
    result = _convert(text)
    assert [line[12:16].strip() for line in result.pdb_text.splitlines()] == ["HZ1", "HZ2"]


@pytest.mark.parametrize("middle", ["", "HZ2"])
def test_mixed_lyn_pair_fails_instead_of_duplicating_names(middle: str) -> None:
    text = _atom(1, "HZ1", "LYN")
    if middle:
        text += _atom(2, middle, "LYN")
    text += _atom(3, "HZ3", "LYN")
    with pytest.raises(NamingConversionError) as caught:
        _convert(text)
    assert caught.value.code is NamingErrorCode.NAMING_COLLISION


def test_anisou_and_full_ter_follow_residue_and_bare_ter_is_preserved() -> None:
    atom = _atom(9, "HZ2", "LYN")
    anisou = _atom(
        9,
        "HZ2",
        "LYN",
        record="ANISOU",
        suffix="    1000   1000   1000      0      0      0       H  ",
    )
    full_ter = f"TER   {10:5d}      {'LYN':<3} A{1:4d} \n"
    text = atom + anisou + full_ter + "TER\n"
    result = _convert(text, overrides={ResidueIdentity("A", "1", ""): "LYN"})
    lines = result.pdb_text.splitlines(keepends=True)
    assert lines[0][12:16].strip() == "HZ1"
    assert lines[1][12:16].strip() == "HZ1"
    assert _resname(lines[1]) == "LSN"
    assert _resname(lines[2]) == "LSN"
    assert lines[3] == "TER\n"
    assert result.summary.changed_atoms == 1


def test_anisou_labels_must_match_coordinate_record() -> None:
    text = _atom(9, "HZ2", "LYN") + _atom(
        9,
        "HZ3",
        "LYN",
        record="ANISOU",
        suffix="    1000   1000   1000      0      0      0       H  ",
    )
    with pytest.raises(NamingConversionError) as caught:
        _convert(text)
    assert caught.value.code is NamingErrorCode.MALFORMED_PDB


def test_bytes_outside_supported_name_fields_are_preserved() -> None:
    source = _atom(123, "HB3", "SER", "z", 42, icode="Q", altloc="B")
    result = _convert(
        source,
        target=ForceFieldTarget.AMBER,
        profile=NamingProfile.GROMACS,
    )
    changed = result.pdb_text
    assert changed[:12] == source[:12]
    assert changed[16:17] == source[16:17]
    assert changed[17:] == source[17:]


@pytest.mark.parametrize("record", ["ATOM", "HETATM", "ANISOU"])
def test_malformed_supported_records_fail_with_stable_code(record: str) -> None:
    with pytest.raises(NamingConversionError) as caught:
        _convert(f"{record}\n")
    assert caught.value.code is NamingErrorCode.MALFORMED_PDB


def test_pure_api_rejects_ambiguous_tuple_identity() -> None:
    with pytest.raises(TypeError, match="VariantOverride"):
        NamingConversionRequest(
            pdb_text="END\n",
            target=ForceFieldTarget.AMBER,
            variant_overrides=(("A", "1"),),  # type: ignore[arg-type]
        )


def test_invalid_variant_has_stable_error_code() -> None:
    with pytest.raises(NamingConversionError) as caught:
        VariantOverride(ResidueIdentity("A", "1"), "HIS")
    assert caught.value.code is NamingErrorCode.INVALID_VARIANT


def test_duplicate_override_has_stable_error_code() -> None:
    override = VariantOverride(ResidueIdentity("A", "1"), "HIE")
    with pytest.raises(NamingConversionError) as caught:
        NamingConversionRequest(
            pdb_text="END\n",
            target=ForceFieldTarget.AMBER,
            variant_overrides=(override, override),
        )
    assert caught.value.code is NamingErrorCode.AMBIGUOUS_VARIANT


def test_unsupported_target_has_stable_error_code() -> None:
    with pytest.raises(NamingConversionError) as caught:
        NamingConversionRequest(pdb_text="END\n", target="opls")  # type: ignore[arg-type]
    assert caught.value.code is NamingErrorCode.UNSUPPORTED_NAMING_DIRECTION
