"""Declared fixed S-S geometry is not implied by byte-preserving custody."""
from dataclasses import replace

import pytest

from dvbfixer.diagnose.report import Severity
from dvbfixer.model.diffusion.contract import ArtifactReference
from dvbfixer.model.diffusion.refinement_stages import RepairAssessment
from dvbfixer.model.diffusion.rfd3_custody_geometry import (
    PrivateGeometryAudit,
    declared_link_findings,
)
from dvbfixer.model.diffusion.rfd3_custody_records import build_ledger, digest
from tests.test_rfd3_source_link_custody import assemble, fixture, verify


@pytest.mark.parametrize("distance,clear", [(2.05, True), (7.0, False), (.1, False), (2.25, False)])
def test_declared_link_metric_checks_fixed_endpoints_without_proximity_inference(distance, clear):
    request, source, ledger, raw = fixture()
    first = next(atom for atom in ledger.atoms if atom.identity.residue_number == "1" and atom.identity.atom_name == "SG")
    second = next(atom for atom in ledger.atoms if atom.identity.residue_number == "3" and atom.identity.atom_name == "SG")
    xyz = (first.xyz[0] + distance, first.xyz[1], first.xyz[2])
    line = second.record[:30] + "".join(f"{v:8.3f}" for v in xyz) + second.record[54:]
    changed = source.replace(second.record.encode(), line.encode())
    request = replace(request, normalized_pdb=ArtifactReference("source.pdb", digest(changed)))
    coordinates = list(raw.coordinates)
    coordinates[raw.atom_axis.index(second.identity)] = xyz
    raw = replace(raw, coordinates=tuple(coordinates))
    case = request, changed, build_ledger(request, changed, changed), raw
    result = assemble(case)
    assert verify(case, result.text).source_link_custody
    findings = declared_link_findings(request, result.text.encode())
    assert (not findings) == clear
    report = PrivateGeometryAudit(RepairAssessment((), ()), findings)
    assert report.geometry_passed == clear
    assert not report.scientific_acceptance
    if not clear and distance != 2.25:
        assert any(finding.severity == Severity.ERROR for finding in findings)
    elif distance == 2.25:
        assert any(finding.severity == Severity.WARNING for finding in findings)
