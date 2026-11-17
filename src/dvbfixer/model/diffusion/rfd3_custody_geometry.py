"""Private geometry of declared links, including links between fixed atoms.

Legacy generated/junction diagnostics alone do not establish the geometry of an
immutable source link. Here only request-declared pairs enter a diagnostic-only
topology. No native/model topology or coordinates are changed, and no link is
inferred from proximity. Canonical diagnostic thresholds remain authoritative.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from dvbfixer.diagnose.report import Finding
from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.rfd3_custody_records import CustodyError, read_records

if TYPE_CHECKING:
    from dvbfixer.model.diffusion.refinement_stages import RepairAssessment


@dataclass(frozen=True)
class PrivateGeometryAudit:
    generated_and_junction: RepairAssessment
    declared_link_findings: tuple[Finding, ...]
    scientific_acceptance: bool = False

    @property
    def geometry_passed(self) -> bool:
        # Canonical warnings remain warnings in evidence, but require review for
        # an immutable declared link rather than becoming an automatic pass.
        return (self.generated_and_junction.geometry_passed
                and not self.generated_and_junction.admission_failures and not self.declared_link_findings)


def declared_link_findings(request: DiffusionRequest, final: bytes) -> tuple[Finding, ...]:
    from openmm import Vec3
    from openmm.app import Topology, element
    from openmm.unit import angstrom

    from dvbfixer.diagnose.chemistry import check_bond_lengths, check_disulfides

    atoms, _, pairs, _ = read_records(final.decode("utf-8"))
    expected = {tuple(sorted((link.atom1, link.atom2))) for link in request.retained_explicit_links}
    if set(pairs) != expected:
        raise CustodyError("declared geometry requires exact request link endpoints")
    rows = {atom.identity: atom for atom in atoms}
    findings = []
    for pair in sorted(expected):
        topology = Topology()
        chains, residues = {}, {}
        endpoints, points = [], []
        for identity in pair:
            row = rows[identity]
            chain = chains.setdefault(identity.chain, None)
            if chain is None:
                chain = topology.addChain(identity.chain)
                chains[identity.chain] = chain
            key = identity.chain, identity.residue_number, identity.insertion_code
            residue = residues.get(key)
            if residue is None:
                residue = topology.addResidue(row.residue_name, chain, id=identity.residue_number,
                                               insertionCode=identity.insertion_code)
                residues[key] = residue
            endpoints.append(topology.addAtom(identity.atom_name, element.get_by_symbol(row.element), residue))
            points.append(Vec3(*row.xyz))
        topology.addBond(*endpoints)
        positions = points * angstrom
        # Canonical bond diagnostics catch stretched declared links even when a
        # distance-based disulfide finder would omit them. For SG-SG, evaluate
        # canonical close-pair geometry too; only this declared pair is present.
        findings.extend(check_bond_lengths(topology, positions))
        if all(identity.atom_name == "SG" for identity in pair):
            findings.extend(check_disulfides(topology, positions))
    return tuple(findings)


def audit_geometry_with_declared_links(request: DiffusionRequest, source: bytes, final: bytes) -> PrivateGeometryAudit:
    """Separate reports; warnings on immutable links require review, not auto-pass."""
    from dvbfixer.model.diffusion.rfd3_custody_validation import audit_geometry

    assessment = audit_geometry(request, source, final)
    findings = declared_link_findings(request, final)
    return PrivateGeometryAudit(assessment, findings)
