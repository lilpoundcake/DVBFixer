"""Private final-byte auditor, not the adapter's self-reported acceptance.

Shares canonical source/link readers and the raw-fit primitive explicitly; does
not call materialization or trust a worker's ledger, coordinates or verdicts.
"""
from __future__ import annotations

from dataclasses import dataclass

from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.rfd3_custody import RawPrecision, RawProposal, audit_raw
from dvbfixer.model.diffusion.rfd3_custody_records import (
    CustodyError,
    SourceLedger,
    build_ledger,
    digest,
    expected_names,
    read_records,
)


@dataclass(frozen=True)
class CustodyAudit:
    final_sha256: str
    observed_atoms: int
    restored_atoms: int
    retained_link_pairs: int
    raw_max_residual_angstrom: float
    raw_sha256: str
    precision: RawPrecision
    canonical_request_sha256: str
    source_sha256: str
    context_sha256: str
    source_link_custody: bool = True
    # A different criterion; never retroactively turn the old failure into pass.
    native_full_graph_transport: str = "not-assessed-by-this-auditor"
    geometry_passed: bool | None = None
    scientific_acceptance: bool = False


def audit_final_bytes(
    request: DiffusionRequest, ledger: SourceLedger, source: bytes, context: bytes,
    raw: RawProposal, final: bytes, *, precision: RawPrecision,
) -> CustodyAudit:
    rebuilt = build_ledger(request, source, context)
    if ledger != rebuilt:
        raise CustodyError("audit source/request/context differs from frozen ledger")
    raw_audit, aligned = audit_raw(request, rebuilt, raw, precision)
    atoms, records, links, auxiliary = read_records(final.decode("utf-8"))
    names = expected_names(request)
    current = {atom.identity: atom for atom in atoms}
    if set(current) != set(names) or any(atom.residue_name != names[atom.identity] or atom.element != atom.identity.atom_name[0] for atom in atoms):
        raise CustodyError("final inventory, element or known sequence changed")
    if any(current[atom.identity].record != atom.record for atom in rebuilt.atoms):
        raise CustodyError("final observation record changed")
    if records != rebuilt.link_records or links != rebuilt.links or auxiliary != rebuilt.auxiliary_records:
        raise CustodyError("final explicit links or auxiliary records changed")
    # Check actual emitted coordinates against the raw proposal, at PDB rounding
    # precision, rather than trusting the adapter's generated-only claim.
    for atom in request.generated_atoms:
        wanted = tuple(float(f"{value:8.3f}") for value in aligned[atom])
        if current[atom].xyz != wanted:
            raise CustodyError("final restored coordinates differ from audited raw frame")
    return CustodyAudit(digest(final), len(rebuilt.atoms), len(request.generated_atoms), len(links),
                        raw_audit.max_observed_residual_angstrom, raw_audit.raw_sha256, precision,
                        rebuilt.canonical_request_sha256, rebuilt.source_sha256, rebuilt.context_sha256)


def audit_geometry(request: DiffusionRequest, source: bytes, final: bytes):
    """Separate scientific geometry finding; never grants context/admission.

    Heavy scientific dependencies are lazy. Existing versioned geometry and
    global chirality checks are reused, not reimplemented as distance heuristics.
    """
    if digest(source) != request.normalized_pdb.sha256:
        raise CustodyError("geometry source differs from frozen request")
    from dvbfixer.ffutils.geometry import assert_all_l
    from dvbfixer.model.diffusion.anchor_acceptance import _structure
    from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

    text = final.decode("utf-8")
    structure = _structure(text, "private-custody-final.pdb")
    assert_all_l(structure.pdb.topology, structure.pdb.positions)
    return assessment_adapter(request, source.decode("utf-8"))(text)
