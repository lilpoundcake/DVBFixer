from __future__ import annotations

import pytest

from dvbfixer.model.diffusion.heterogen_context import (
    HeterogenContextError,
    build_smiles_heterogen_contexts,
)


def _hetatm(serial: int, atom: str, element: str, x: float) -> str:
    return (
        f"HETATM{serial:5d} {atom:>4} LIG X 501    "
        f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00          {element:>2}\n"
    )


def test_smiles_context_records_authoritative_graph_and_atom_mapping() -> None:
    text = (
        _hetatm(1, "C1", "C", 0.0)
        + _hetatm(2, "O1", "O", 1.4)
        + "CONECT    1    2\nEND\n"
    )

    contexts = build_smiles_heterogen_contexts(text, {"LIG": "CO"})

    assert len(contexts) == 1
    assert contexts[0].canonical_smiles == "CO"
    assert [atom.atom_name for atom in contexts[0].atoms] == ["C1", "O1"]
    assert len(contexts[0].authoritative_graph_sha256) == 64


def test_smiles_context_rejects_unmapped_and_covalently_attached_heterogens() -> None:
    text = (
        _hetatm(1, "C1", "C", 0.0)
        + _hetatm(2, "O1", "O", 1.4)
        + "ATOM      3  CA  ALA A   1       2.800   0.000   0.000  1.00  0.00           C  \n"
        + "CONECT    1    2    3\nEND\n"
    )

    with pytest.raises(HeterogenContextError, match="authoritative SMILES"):
        build_smiles_heterogen_contexts(text, {})
    with pytest.raises(HeterogenContextError, match="covalently attached"):
        build_smiles_heterogen_contexts(text, {"LIG": "CO"})
