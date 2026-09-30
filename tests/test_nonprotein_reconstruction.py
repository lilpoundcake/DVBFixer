from __future__ import annotations

import hashlib
import json
from dataclasses import replace
from pathlib import Path

import pytest

from dvbfixer.domain.nonprotein_reconstruction import (
    ChemicalAtom,
    ChemicalBond,
    ChemicalGraph,
    ComponentClass,
    ExternalLink,
    LocalPinnedCcdAuthority,
    ObservedAtom,
    OnlineEnrichmentAuthority,
    ParameterizationRequest,
    ReconstructionRequest,
    ReconstructionStatus,
    UserMappedAuthority,
)
from dvbfixer.domain.parameterization import ParameterizationRoute
from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.nonprotein_reconstruction import (
    CCD_FORMAT,
    chemical_graph_digest,
    decide_parameterization,
    publish_reconstruction,
    reconstruct_and_publish,
    reconstruct_nonprotein,
)


def _component(chain: str = "d", *, icode: str = "A", altloc: str = "B") -> ComponentInstanceRef:
    return ComponentInstanceRef(2, chain, "82", icode, altloc, 3)


def _graph(*, component_id: str = "LIG", element: str = "O", state_id: str = "anion") -> ChemicalGraph:
    return ChemicalGraph(
        component_id,
        state_id,
        (
            ChemicalAtom("C1", "C", 0, (0.0, 0.0, 0.0), stereo="R"),
            ChemicalAtom("C2", "C", 0, (1.0, 0.0, 0.0), aromatic=True),
            ChemicalAtom("N1", "N", 0, (0.0, 1.0, 0.0), aromatic=True),
            ChemicalAtom("X1", element, -1 if element == "O" else 2, (0.0, 0.0, 1.0)),
        ),
        (
            ChemicalBond("C1", "C2", 1.5, 1.0, aromatic=True),
            ChemicalBond("C1", "N1", 1.5, 1.0, aromatic=True),
            ChemicalBond("C1", "X1", 1.0, 1.0),
        ),
        -1 if element == "O" else 2,
        6.0,
        14.0,
    )


def _authority(graph: ChemicalGraph | None = None) -> UserMappedAuthority:
    graph = graph or _graph()
    return UserMappedAuthority(graph, "user:ligand-map", chemical_graph_digest(graph))


def _request(
    *,
    component: ComponentInstanceRef | None = None,
    authority: object | None = None,
) -> ReconstructionRequest:
    component = component or _component()
    atoms = (
        ObservedAtom(ExactAtomRef(component, "C1"), "C", (10.0, 20.0, 30.0)),
        ObservedAtom(ExactAtomRef(component, "C2"), "C", (11.0, 20.0, 30.0)),
        ObservedAtom(ExactAtomRef(component, "N1"), "N", (10.0, 21.0, 30.0)),
    )
    return ReconstructionRequest(component, atoms, authority or _authority(), 7.4)  # type: ignore[arg-type]


def _snapshot(path: Path, states: list[dict[str, object]]) -> str:
    data = {
        "format": CCD_FORMAT,
        "version": "ccd-test-2026-09",
        "components": {"LIG": {"states": states}},
    }
    content = json.dumps(data, sort_keys=True).encode()
    path.write_bytes(content)
    return hashlib.sha256(content).hexdigest()


def _state(graph: ChemicalGraph, **updates: object) -> dict[str, object]:
    result: dict[str, object] = {
        "id": graph.state_id,
        "net_charge": graph.net_charge,
        "ph_min": graph.ph_min,
        "ph_max": graph.ph_max,
        "priority": graph.priority,
        "atoms": [
            {
                "name": atom.name,
                "element": atom.element,
                "formal_charge": atom.formal_charge,
                "ideal_position": list(atom.ideal_position),
                "aromatic": atom.aromatic,
                "stereo": atom.stereo,
            }
            for atom in graph.atoms
        ],
        "bonds": [
            {
                "atom1": bond.atom1,
                "atom2": bond.atom2,
                "order": bond.order,
                "ideal_length": bond.ideal_length,
                "aromatic": bond.aromatic,
            }
            for bond in graph.bonds
        ],
    }
    result.update(updates)
    return result


def test_user_mapped_class_a_reconstructs_only_missing_atoms() -> None:
    request = _request()
    result = reconstruct_nonprotein(request)

    assert result.status is ReconstructionStatus.SUCCEEDED
    assert result.component_class is ComponentClass.A
    assert result.component == _component()
    assert [atom.identity.name for atom in result.coordinates] == ["C1", "C2", "N1", "X1"]
    assert result.coordinates[-1].position == pytest.approx((10.0, 20.0, 31.0))
    assert result.coordinates[-1].source.value == "source-template"
    assert all(movement.distance == 0.0 for movement in result.movements)
    assert result.microstate is not None and result.microstate.reason == "user graph is a locked microspecies"
    assert result.output_digest and not result.geometry_approved and not result.md_ready


def test_exact_identity_preserves_case_model_icode_altloc_and_occurrence() -> None:
    lower = reconstruct_nonprotein(_request(component=_component("d")))
    upper = reconstruct_nonprotein(_request(component=_component("D")))

    assert lower.component != upper.component
    assert lower.coordinates[0].identity.component == ComponentInstanceRef(2, "d", "82", "A", "B", 3)
    assert lower.output_digest != upper.output_digest


def test_pinned_local_snapshot_selects_explicit_ph_state_deterministically(tmp_path: Path) -> None:
    graph = _graph()
    acidic = replace(graph, state_id="acid", ph_min=0.0, ph_max=5.9)
    basic = replace(graph, state_id="anion", ph_min=6.0, ph_max=14.0)
    snapshot = tmp_path / "ccd.json"
    digest = _snapshot(snapshot, [_state(acidic), _state(basic)])
    request = _request(authority=LocalPinnedCcdAuthority(snapshot, digest, "LIG"))

    first = reconstruct_nonprotein(request)
    second = reconstruct_nonprotein(request)

    assert first == second
    assert first.microstate is not None and first.microstate.selected_state == "anion"
    assert first.authority is not None and first.authority.source_digest == digest
    assert first.authority.source_version == "ccd-test-2026-09"


def test_overlapping_equal_priority_microstates_are_ambiguous(tmp_path: Path) -> None:
    graph = _graph()
    snapshot = tmp_path / "ccd.json"
    digest = _snapshot(
        snapshot,
        [
            _state(replace(graph, state_id="one", ph_min=0.0, ph_max=10.0)),
            _state(replace(graph, state_id="two", ph_min=5.0, ph_max=14.0)),
        ],
    )

    result = reconstruct_nonprotein(
        _request(authority=LocalPinnedCcdAuthority(snapshot, digest, "LIG"))
    )

    assert result.status is ReconstructionStatus.AMBIGUOUS
    assert result.component_class is ComponentClass.D
    assert result.findings[0].code == "microstate-ambiguous"
    assert result.microstate is not None
    assert result.microstate.candidate_scores == (("one", 0), ("two", 0))
    assert not result.coordinates and result.output_digest is None


@pytest.mark.parametrize(
    ("reconstruction_request", "expected_class", "code"),
    [
        (
            replace(
                _request(),
                external_links=(
                    ExternalLink(
                        ExactAtomRef(_component(), "C1"),
                        ExactAtomRef(ComponentInstanceRef(2, "A", "10"), "SG"),
                    ),
                ),
            ),
            ComponentClass.B,
            "class-b-unsupported",
        ),
        (_request(authority=_authority(_graph(element="ZN"))), ComponentClass.C, "class-c-unsupported"),
        (_request(authority=_authority(_graph(component_id="HEM"))), ComponentClass.C, "class-c-unsupported"),
        (_request(authority=_authority(_graph(element="XE"))), ComponentClass.D, "class-d-unsupported"),
    ],
)
def test_linked_and_metal_classes_have_explicit_refusal_contracts(
    reconstruction_request: ReconstructionRequest, expected_class: ComponentClass, code: str
) -> None:
    result = reconstruct_nonprotein(reconstruction_request)

    assert result.status is ReconstructionStatus.UNSUPPORTED
    assert result.component_class is expected_class
    assert result.findings[0].code == code
    assert not result.coordinates and result.output_digest is None


def test_online_authority_is_explicitly_refused() -> None:
    result = reconstruct_nonprotein(_request(authority=OnlineEnrichmentAuthority("example")))

    assert result.status is ReconstructionStatus.UNSUPPORTED
    assert result.findings[0].code == "online-enrichment-disabled"


def test_digest_mismatch_fails_before_geometry_or_publication(tmp_path: Path) -> None:
    snapshot = tmp_path / "ccd.json"
    _snapshot(snapshot, [_state(_graph())])
    request = _request(authority=LocalPinnedCcdAuthority(snapshot, "0" * 64, "LIG"))

    result, published = reconstruct_and_publish(request, tmp_path, "result")

    assert result.status is ReconstructionStatus.INVALID_INPUT
    assert result.findings[0].code == "authority-digest-mismatch"
    assert published is None
    assert not (tmp_path / "result").exists()


def test_insufficient_anchor_geometry_fails_closed() -> None:
    request = _request()
    result = reconstruct_nonprotein(replace(request, observed_atoms=request.observed_atoms[:2]))

    assert result.status is ReconstructionStatus.UNSUPPORTED
    assert result.findings[0].code == "insufficient-anchors"
    assert not result.coordinates


def test_inverted_observed_stereochemistry_fails_validation() -> None:
    request = _request()
    inverted = request.observed_atoms + (
        ObservedAtom(ExactAtomRef(request.component, "X1"), "O", (10.0, 20.0, 29.0)),
    )

    result = reconstruct_nonprotein(replace(request, observed_atoms=inverted))

    assert result.status is ReconstructionStatus.FAILED
    assert result.findings[0].code == "stereochemistry"
    assert not result.coordinates and result.output_digest is None


def test_user_graph_security_bound_is_structured_invalid_input() -> None:
    graph = _graph()
    oversized = replace(
        graph,
        atoms=tuple(
            ChemicalAtom(f"C{index}", "C", 0, (float(index), 0.0, 0.0))
            for index in range(513)
        ),
        bonds=(),
        net_charge=0,
    )

    result = reconstruct_nonprotein(_request(authority=_authority(oversized)))

    assert result.status is ReconstructionStatus.INVALID_INPUT
    assert result.findings[0].code == "authority-bounds"
    assert result.input_digest is None


def test_parameterization_requires_approved_success_and_never_claims_md_ready() -> None:
    geometry = reconstruct_nonprotein(_request())
    unapproved = decide_parameterization(
        ParameterizationRequest(geometry, False, "amber", "explicit-v1", allow_isolated_ligand_candidate=True)
    )
    approved = decide_parameterization(
        ParameterizationRequest(geometry, True, "amber", "explicit-v1", allow_isolated_ligand_candidate=True)
    )
    user_template = decide_parameterization(
        ParameterizationRequest(geometry, True, "amber", "explicit-v1", validated_user_template=True)
    )

    assert unapproved.status is ReconstructionStatus.UNSUPPORTED
    assert approved.route is ParameterizationRoute.GAFF_CANDIDATE
    assert user_template.route is ParameterizationRoute.USER_TEMPLATE
    assert not approved.md_ready and "have not run" in approved.reason


def test_parameterization_rejects_tampered_geometry_digest() -> None:
    geometry = reconstruct_nonprotein(_request())
    changed_atom = replace(geometry.coordinates[-1], position=(99.0, 99.0, 99.0))
    tampered = replace(geometry, coordinates=geometry.coordinates[:-1] + (changed_atom,))

    decision = decide_parameterization(
        ParameterizationRequest(tampered, True, "amber", "explicit-v1", exact_native_match=True)
    )

    assert decision.status is ReconstructionStatus.INVALID_INPUT
    assert "digest" in decision.reason


def test_success_bundle_is_atomic_deterministic_and_contains_provenance(tmp_path: Path) -> None:
    result, published = reconstruct_and_publish(_request(), tmp_path, "ligand-geometry")

    assert result.status is ReconstructionStatus.SUCCEEDED
    assert published == tmp_path / "ligand-geometry"
    assert sorted(path.name for path in published.iterdir()) == ["geometry.json", "provenance.json"]
    provenance = json.loads((published / "provenance.json").read_text())
    assert provenance["md_ready"] is False
    assert provenance["component"]["alternate_location"] == "B"
    assert provenance["coordinates"][-1]["source"] == "source-template"


def test_publication_rolls_back_staging_on_rename_failure(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    result = reconstruct_nonprotein(_request())

    def fail_rename(_source: Path, _destination: Path) -> None:
        raise OSError("simulated publication failure")

    monkeypatch.setattr("dvbfixer.nonprotein_reconstruction.os.replace", fail_rename)
    with pytest.raises(OSError, match="simulated"):
        publish_reconstruction(result, tmp_path, "result")

    assert not (tmp_path / "result").exists()
    assert list(tmp_path.iterdir()) == []


def test_publication_rejects_path_traversal(tmp_path: Path) -> None:
    result = reconstruct_nonprotein(_request())

    with pytest.raises(ValueError, match="portable name"):
        publish_reconstruction(result, tmp_path, "../escape")


def test_publication_rejects_symlink_root(tmp_path: Path) -> None:
    result = reconstruct_nonprotein(_request())
    real_root = tmp_path / "real"
    real_root.mkdir()
    linked_root = tmp_path / "linked"
    linked_root.symlink_to(real_root, target_is_directory=True)

    with pytest.raises(ValueError, match="non-symlink"):
        publish_reconstruction(result, linked_root, "result")
