from dataclasses import asdict

import numpy as np
import pytest
from scripts.audit_1non_diffusion_full import native_coordinates

from dvbfixer.model.diffusion.contract import (
    DIFFUSION_SCHEMA_VERSION,
    ArtifactReference,
    AtomIdentity,
    DiffusionRequest,
    GapRegion,
    ResidueIdentity,
    SequencePlacement,
    TargetInterval,
    TargetSequence,
)
from dvbfixer.model.diffusion.protein_atoms import ATOM37_ORDER, CANONICAL_HEAVY_ATOMS


def transport_only_view():
    # Synthetic CA-only transport axis, NOT a chemically complete modeling input.
    targets, placements, gaps, fixed, generated = [], [], [], [], []
    for chain in ("A", "a", "B", "b"):
        targets.append(TargetSequence(chain, "A"*179))
        observed = tuple(i for i in range(179) if not 72 <= i < 90)
        residues = tuple(ResidueIdentity(chain, str(i+1)) for i in range(179))
        placements.append(SequencePlacement(chain, 179, observed, tuple(residues[i] for i in observed)))
        gaps.append(GapRegion(chain, TargetInterval(72, 90), residues[71], residues[90], residues[72:90],
                              (residues[71], *residues[72:90], residues[90])))
        fixed += [AtomIdentity(chain, str(i+1), "", "CA") for i in observed]
        generated += [AtomIdentity(chain, str(i+1), "", "CA") for i in range(72, 90)]
    request = DiffusionRequest(schema_version=DIFFUSION_SCHEMA_VERSION,
        normalized_pdb=ArtifactReference("synthetic.pdb", "0"*64),
        target_sequences=tuple(targets), sequence_placements=tuple(placements), gaps=tuple(gaps),
        fixed_atoms=tuple(fixed), generated_atoms=tuple(generated), retained_explicit_links=(),
        candidate_count=1, seeds=(23,), backend_options=())
    mapping = [[asdict(a), asdict(a)] for a in (*fixed, *generated)]
    return {"engine_view": {"request": asdict(request), "atom_to_source": mapping}}


def test_native_decoder_keeps_all_four_chains_and_generated_observed_mapping():
    view = transport_only_view()
    native = np.zeros((716, 37, 3))
    native[:, ATOM37_ORDER.index("CA"), 0] = np.arange(716)
    decoded = native_coordinates(view, native)
    assert len(decoded) == 716
    assert decoded[AtomIdentity("A", "73", "", "CA")][0] == 72
    assert decoded[AtomIdentity("a", "73", "", "CA")][0] == 179+72
    assert decoded[AtomIdentity("b", "179", "", "CA")][0] == 715


def test_native_decoder_refuses_local_axis_and_nonfinite_endpoint():
    view = transport_only_view()
    with pytest.raises(ValueError, match="complete tetramer"):
        native_coordinates(view, np.zeros((354, 37, 3)))
    native = np.zeros((716, 37, 3))
    native[0, 0, 0] = np.nan
    with pytest.raises(ValueError, match="finite complete"):
        native_coordinates(view, native)


def test_oxt_array_slot_is_not_a_native_output_atom():
    view = transport_only_view()
    oxt = asdict(AtomIdentity("A", "179", "", "OXT"))
    view["engine_view"]["atom_to_source"].append([oxt, oxt])
    native = np.ones((716, 37, 3))
    decoded = native_coordinates(view, native)
    assert AtomIdentity("A", "179", "", "OXT") not in decoded
    assert len(decoded) == 716


def test_declared_input_passthrough_mask_keeps_oxt_but_not_unowned_output_atoms():
    view = transport_only_view()
    oxt = asdict(AtomIdentity("A", "179", "", "OXT"))
    view["engine_view"]["request"]["fixed_atoms"] += (oxt,)
    view["engine_view"]["atom_to_source"].append([oxt, oxt])
    native = np.ones((716, 37, 3))
    mask = np.zeros((716, 37), dtype=bool)
    for name in CANONICAL_HEAVY_ATOMS["ALA"]:
        mask[:, ATOM37_ORDER.index(name)] = True
    mask[178, ATOM37_ORDER.index("OXT")] = True
    decoded = native_coordinates(view, native, mask)
    assert AtomIdentity("A", "179", "", "OXT") in decoded
    mask[179, ATOM37_ORDER.index("OXT")] = True
    with pytest.raises(ValueError, match="exact fixed input atoms"):
        native_coordinates(view, native, mask)
