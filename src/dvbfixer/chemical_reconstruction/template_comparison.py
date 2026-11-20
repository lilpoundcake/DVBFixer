"""Private retained-context seeded heavy baseline, not a graph strategy or writer."""
from __future__ import annotations

import hashlib
import math
from contextlib import redirect_stderr, redirect_stdout
from dataclasses import dataclass
from importlib.metadata import version
from importlib.util import find_spec
from io import StringIO
from pathlib import Path
from typing import Literal

from dvbfixer.domain.structure_identity import ExactAtomRef
from dvbfixer.pdb_component_reconstruction import _component_ref
from dvbfixer.protein_completion.links import parse_declared_links

from .contracts import Observation, ReconstructionRequest, Unit, _sha, _tuple, digest_record
from .observed_peptide import _parse_source


@dataclass(frozen=True)
class SeededHeavyBaseline:
    request_digest: str
    masked_source_sha256: str
    status: Literal["coordinates-proposed", "refused", "unavailable"]
    reason: str
    raw: tuple[Observation, ...] = ()
    rounded: tuple[Observation, ...] = ()
    temporary_atoms: tuple[Observation, ...] = ()
    known_maximum_movement_angstrom: float | None = None
    global_l_passed: bool = False
    environment: tuple[tuple[str, str], ...] = ()
    log: str = ""
    method: str = "seeded-PDBFixer-heavy-helper-1-to-5-Reference-no-H-no-reflect-v1"

    def __post_init__(self) -> None:
        _sha(self.request_digest)
        _sha(self.masked_source_sha256)
        for values in (self.raw, self.rounded, self.temporary_atoms):
            _tuple(values, Observation, "template comparison coordinates", maximum=100000)
        if self.method != "seeded-PDBFixer-heavy-helper-1-to-5-Reference-no-H-no-reflect-v1":
            raise ValueError("template-baseline-method")


def masked_probe_source(prepared: bytes, request: ReconstructionRequest) -> tuple[bytes, tuple[str, ...]]:
    """I/O-only exact mask adapter; constructor receives returned bytes only."""
    lines, atoms = _parse_source(prepared)
    refs = {a.line_index: ExactAtomRef(_component_ref(a), a.name) for a in atoms}
    serials = {a.serial: refs[a.line_index] for a in atoms}
    removed, output = [], []
    for index, line in enumerate(lines):
        remove = refs.get(index) in request.ownership.generated
        if line.startswith("CONECT"):
            ids = [int(line[i:i+5]) for i in range(6, len(line.rstrip("\n")), 5) if line[i:i+5].strip()]
            remove = any(serials.get(serial) in request.ownership.generated for serial in ids)
        if line[:6].strip() in {"LINK", "SSBOND"}:
            remove = any(link.left in request.ownership.generated or link.right in request.ownership.generated
                         for link in parse_declared_links([line], atoms))
        if remove:
            removed.append(line)
        else:
            output.append(line)
    masked = "".join(output).encode("ascii")
    validate_masked_source(masked, request)
    return masked, tuple(removed)


def validate_masked_source(source: bytes, request: ReconstructionRequest) -> None:
    lines, atoms = _parse_source(source)
    observations = tuple(Observation(ExactAtomRef(_component_ref(a), a.name), a.element, a.position, Unit.ANGSTROM,
                                     float(lines[a.line_index][54:60])) for a in atoms)
    if digest_record(observations) != digest_record(request.observations):
        raise ValueError("template-masked-source-exact-observation-transport")
    if (request.ownership.mobile or request.ownership.temporary or not request.ownership.generated
            or len(request.targets) != 1 or any(a.component not in request.targets for a in request.ownership.generated)):
        raise ValueError("template-baseline-exact-fixed-heavy-mask")


def inspect_temporary(request: ReconstructionRequest, raw: tuple[Observation, ...],
                        temporary: tuple[Observation, ...]) -> float:
    known = {o.identity: o for o in request.observations}
    all_atoms = {o.identity: o for o in temporary}
    if (len(all_atoms) != len(temporary) or set(all_atoms) != set(known) | request.ownership.generated
            or {o.identity for o in raw} != request.ownership.generated or len(raw) != len(request.ownership.generated)
            or any(all_atoms[o.identity] != o for o in raw)
            or any(all_atoms[a].element != o.element for a, o in known.items())):
        raise ValueError("template-temporary-extra-missing-or-ambiguous-atoms")
    return max(math.dist(o.position, all_atoms[a].position) for a, o in known.items())


def template_environment() -> tuple[tuple[str, str], ...]:
    spec = find_spec("pdbfixer")
    if spec is None or spec.origin is None:
        raise ModuleNotFoundError("pdbfixer")
    root = Path(spec.origin).parent
    return (("PDBFixer", version("pdbfixer")), ("OpenMM", version("openmm")), ("platform", "Reference"),
            ("soft.xml.sha256", hashlib.sha256((root/"soft.xml").read_bytes()).hexdigest()),
            ("PRO.pdb.sha256", hashlib.sha256((root/"templates/PRO.pdb").read_bytes()).hexdigest()))


def audit_seeded_baseline(request: ReconstructionRequest, masked: bytes, result: SeededHeavyBaseline) -> None:
    """Saved snapshot audit; metadata-only Topology, no PDBFixer/coordinate build."""
    validate_masked_source(masked, request)
    if result.request_digest != request.digest or result.masked_source_sha256 != hashlib.sha256(masked).hexdigest():
        raise ValueError("template-saved-source-binding")
    if result.status == "unavailable":
        if result.raw or result.rounded or result.temporary_atoms or result.global_l_passed:
            raise ValueError("template-unavailable-cannot-have-coordinates")
        return
    if result.environment != template_environment():
        raise ValueError("template-saved-environment-binding")
    if not result.raw:
        if result.status != "refused" or result.rounded or result.temporary_atoms or result.global_l_passed:
            raise ValueError("template-empty-proposal-not-a-pass")
        return
    drift = inspect_temporary(request, result.raw, result.temporary_atoms)
    rounded = tuple(Observation(o.identity, o.element, (float(f"{o.position[0]:.3f}"), float(f"{o.position[1]:.3f}"),
                               float(f"{o.position[2]:.3f}")), o.unit, o.occupancy) for o in result.raw)
    if result.known_maximum_movement_angstrom != drift or result.rounded != rounded:
        raise ValueError("template-saved-drift-or-rounding")
    from openmm import Vec3
    from openmm.app import Topology, element
    from openmm.unit import angstrom

    from dvbfixer.ffutils.geometry import ChiralityError, assert_all_l

    _, atoms = _parse_source(masked)
    names = {_component_ref(a): a.residue_name for a in atoms}
    topology = Topology()
    chains, residues = {}, {}
    for o in result.temporary_atoms:
        component = o.identity.component
        if component.chain_id not in chains:
            chains[component.chain_id] = topology.addChain(component.chain_id)
        if component not in residues:
            residues[component] = topology.addResidue(names[component], chains[component.chain_id],
                id=component.sequence_number, insertionCode=component.insertion_code)
        topology.addAtom(o.identity.name, element.get_by_symbol(o.element), residues[component])
    positions = [Vec3(*o.position)*angstrom for o in result.temporary_atoms]
    passed = True
    try:
        assert_all_l(topology, positions)
    except ChiralityError:
        passed = False
    if result.global_l_passed != passed or (result.status == "coordinates-proposed") != (passed and drift <= 1e-4):
        raise ValueError("template-saved-global-L-or-status")


def seeded_heavy_baseline(request: ReconstructionRequest, masked: bytes) -> SeededHeavyBaseline:
    """One scoped helper call; not five independent trials or a full legacy prep."""
    validate_masked_source(masked, request)
    wanted = request.ownership.generated
    lookup = {(a.component.chain_id, a.component.sequence_number, a.component.insertion_code, a.name): a
              for a in ({o.identity for o in request.observations} | wanted)}
    if len(lookup) != len(request.observations)+len(wanted):
        raise ValueError("template-source-identity-not-representable-in-OpenMM")
    pin = hashlib.sha256(masked).hexdigest()
    try:
        from openmm import Platform
        from openmm.unit import angstrom
        from pdbfixer import PDBFixer
    except ModuleNotFoundError as exc:
        return SeededHeavyBaseline(request.digest, pin, "unavailable", str(exc))
    from dvbfixer.ffutils.geometry import (
        ChiralityError,
        assert_all_l,
        rebuild_missing_atoms_with_retry,
    )
    from dvbfixer.profile_readiness import pdbfixer_resources_ready

    if not pdbfixer_resources_ready():
        return SeededHeavyBaseline(request.digest, pin, "unavailable", "pdbfixer-resources-not-ready")
    environment = template_environment()
    fixer = PDBFixer(pdbfile=StringIO(masked.decode("ascii")))
    fixer.platform = Platform.getPlatformByName("Reference")
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    target = next(iter(request.targets))
    residues = [r for r in fixer.topology.residues() if
                (r.chain.id, r.id, r.insertionCode.strip()) == (target.chain_id, target.sequence_number, target.insertion_code)]
    if len(residues) != 1:
        return SeededHeavyBaseline(request.digest, pin, "refused", "template-ambiguous-target-residue", environment=environment)
    residue = residues[0]
    missing = [a for a in fixer.missingAtoms.get(residue, []) if a.name in {r.name for r in wanted}]
    if {a.name for a in missing} != {r.name for r in wanted}:
        return SeededHeavyBaseline(request.digest, pin, "refused", "template-requested-mask-not-exactly-rebuildable", environment=environment)
    fixer.missingAtoms = {residue: missing}
    fixer.missingTerminals = {}
    log = StringIO()
    with redirect_stdout(log), redirect_stderr(log):
        rebuild_missing_atoms_with_retry(fixer, verbose=True, log_prefix="WP2P: ")
    positions = fixer.positions.value_in_unit(angstrom)
    # Map only through exact existing identities and the explicit generated mask.
    occupancies = {o.identity: o.occupancy for o in request.observations}
    atoms = []
    for atom in fixer.topology.atoms():
        r = atom.residue
        key = (r.chain.id, r.id, r.insertionCode.strip(), atom.name)
        if key not in lookup:
            return SeededHeavyBaseline(request.digest, pin, "refused", "template-unrequested-output-atom", environment=environment, log=log.getvalue())
        point = positions[atom.index]
        identity = lookup[key]
        occupancy = occupancies.get(identity, 1.)
        atoms.append(Observation(identity, atom.element.symbol, (float(point[0]), float(point[1]), float(point[2])), Unit.ANGSTROM, occupancy))
    temporary = tuple(atoms)
    raw = tuple(sorted((o for o in temporary if o.identity in wanted), key=lambda o: o.identity))
    drift = inspect_temporary(request, raw, temporary)
    rounded = tuple(Observation(o.identity, o.element, (float(f"{o.position[0]:.3f}"), float(f"{o.position[1]:.3f}"),
                               float(f"{o.position[2]:.3f}")), o.unit, o.occupancy) for o in raw)
    reason = "temporary-known-motion" if drift > 1e-4 else "seeded-heavy-coordinates-only-not-geometry-approved"
    passed = False
    try:
        assert_all_l(fixer.topology, fixer.positions)
        passed = True
    except ChiralityError as exc:
        reason = "global-L-refused:"+str(exc)
    status: Literal["coordinates-proposed", "refused"] = "coordinates-proposed" if passed and drift <= 1e-4 else "refused"
    return SeededHeavyBaseline(request.digest, pin, status, reason, raw, rounded, temporary, drift, passed, environment, log.getvalue())
