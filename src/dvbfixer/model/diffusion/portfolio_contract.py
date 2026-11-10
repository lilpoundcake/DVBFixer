"""Version-one private multi-parent publication; disjoint from sampling bundles."""
from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

from dvbfixer.model.diffusion.contract import (
    ArtifactReference,
    DiffusionContractError,
    _JsonContract,
)
from dvbfixer.model.diffusion.gap_portfolio import MAX_PARENTS, GapChoice


@dataclass(frozen=True, slots=True)
class PortfolioParentArtifacts:
    candidate_id: str
    seed: int
    coordinates: ArtifactReference
    sampler_trace: ArtifactReference
    staged_evidence: ArtifactReference
    raw: ArtifactReference
    runner_result: ArtifactReference

    def __post_init__(self):
        if not isinstance(self.candidate_id, str) or not self.candidate_id or len(self.candidate_id) > 256:
            raise DiffusionContractError("portfolio parent needs a bounded candidate ID")
        if type(self.seed) is not int or self.seed < 0:
            raise DiffusionContractError("portfolio parent seed must be a nonnegative integer")
        if any(not isinstance(a, ArtifactReference) for a in self.artifacts().values()):
            raise DiffusionContractError("portfolio parent needs typed artifact references")

    def artifacts(self):
        return {name: getattr(self, name) for name in
                ("coordinates", "sampler_trace", "staged_evidence", "raw", "runner_result")}


@dataclass(frozen=True, slots=True)
class PortfolioPreparationArtifacts:
    original: ArtifactReference
    prepared: ArtifactReference
    report: ArtifactReference
    request_snapshot: ArtifactReference

    def artifacts(self):
        return {"input-original": self.original, "input-prepared": self.prepared,
                "input-repair-report": self.report, "input-repair-request": self.request_snapshot}


@dataclass(frozen=True, slots=True)
class PortfolioEvidence(_JsonContract):
    request: ArtifactReference
    source: ArtifactReference
    selected: ArtifactReference
    candidate: ArtifactReference
    search_protocol: ArtifactReference
    search_report: ArtifactReference
    selection_source: ArtifactReference
    parents: tuple[PortfolioParentArtifacts, ...]
    choices: tuple[GapChoice, ...]
    repair_report: ArtifactReference | None = None
    repair_selector_source: ArtifactReference | None = None
    repair_kernel_source: ArtifactReference | None = None
    repair_conditioning_source: ArtifactReference | None = None
    preparation: PortfolioPreparationArtifacts | None = None
    kind: str = "experimental-gap-portfolio-evidence"
    scientific_qualification: bool = False
    public_profile_enabled: bool = False

    _schema_version: ClassVar[int] = 1

    def __post_init__(self):
        super(PortfolioEvidence, self).__post_init__()
        if (type(self.schema_version) is not int or self.kind != "experimental-gap-portfolio-evidence"
                or self.scientific_qualification is not False or self.public_profile_enabled is not False):
            raise DiffusionContractError("portfolio evidence cannot claim public qualification")
        if not 1 <= len(self.parents) <= MAX_PARENTS or len({(p.candidate_id, p.seed) for p in self.parents}) != len(self.parents):
            raise DiffusionContractError("portfolio needs bounded distinct parents")
        if not self.choices or len(self.choices) > 4096 or tuple(c.gap_index for c in self.choices) != tuple(range(len(self.choices))):
            raise DiffusionContractError("portfolio choices must cover ordered whole gaps")
        if any(type(c.gap_index) is not int or type(c.parent_index) is not int or not 0 <= c.parent_index < len(self.parents)
               for c in self.choices):
            raise DiffusionContractError("portfolio choice references an absent parent")
        repair = (self.repair_report, self.repair_selector_source, self.repair_kernel_source, self.repair_conditioning_source)
        if any(a is not None for a in repair) and not all(a is not None for a in repair):
            raise DiffusionContractError("portfolio local repair requires its report and source snapshots")
        artifacts = self.artifacts()
        if any(not isinstance(a, ArtifactReference) for a in artifacts.values()):
            raise DiffusionContractError("portfolio needs typed artifact references")
        if len({a.path for a in artifacts.values()}) != len(artifacts):
            raise DiffusionContractError("portfolio artifact paths must be distinct")

    def artifacts(self):
        values = {name: getattr(self, name) for name in ("request", "source", "selected", "candidate",
            "search_protocol", "search_report", "selection_source", "repair_report", "repair_selector_source",
            "repair_kernel_source", "repair_conditioning_source") if getattr(self, name) is not None}
        for i, parent in enumerate(self.parents):
            values.update({f"parent-{i:02d}-{name}": a for name, a in parent.artifacts().items()})
        if self.preparation is not None:
            values.update(self.preparation.artifacts())
        return values


@dataclass(frozen=True, slots=True)
class PortfolioBundleIndex(_JsonContract):
    evidence: ArtifactReference
    dat: ArtifactReference
    validation: ArtifactReference
    resume_key: str
    kind: str = "experimental-gap-portfolio-bundle"
    scientific_qualification: bool = False
    public_profile_enabled: bool = False

    _schema_version: ClassVar[int] = 1

    def __post_init__(self):
        super(PortfolioBundleIndex, self).__post_init__()
        if (type(self.schema_version) is not int or self.kind != "experimental-gap-portfolio-bundle"
                or self.scientific_qualification is not False or self.public_profile_enabled is not False):
            raise DiffusionContractError("portfolio bundle cannot claim public qualification")
        ArtifactReference("resume-key", self.resume_key)
        if (any(not isinstance(a, ArtifactReference) for a in (self.evidence, self.dat, self.validation))
                or len({a.path for a in (self.evidence, self.dat, self.validation)}) != 3):
            raise DiffusionContractError("portfolio index needs distinct typed references")
