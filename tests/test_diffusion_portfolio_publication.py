"""Multi-parent publication is not a relabelled schema-4 candidate or a cache flag."""
from dataclasses import asdict, replace

import pytest

from dvbfixer.ffutils.dat import DatRecord
from dvbfixer.model.diffusion.contract import (
    ArtifactReference,
    DiffusionContractError,
    RunnerCandidate,
    SamplerTrace,
)
from dvbfixer.model.diffusion.gap_portfolio import GapChoice, digest
from dvbfixer.model.diffusion.portfolio_contract import PortfolioEvidence, PortfolioParentArtifacts
from dvbfixer.model.diffusion.portfolio_evidence import (
    json_bytes,
    portfolio_resume_key,
    strict_json,
)
from dvbfixer.model.diffusion.portfolio_publication import (
    publish_portfolio_bundle,
    read_portfolio_bundle,
    resume_portfolio_bundle,
)
from dvbfixer.model.diffusion.refinement_evidence import StagedRefinementEvidence, bind_evidence
from dvbfixer.model.diffusion.refinement_policy import POLICY_REVISION, RefinementPolicy
from tests.test_diffusion_anchor_acceptance import shift
from tests.test_diffusion_gap_portfolio import two_gaps
from tests.test_diffusion_validate import _runner_result, _trace_bytes


def setup_evidence(root, *, repair=False):
    request, source, good = two_gaps()
    request = replace(request, backend_options=RefinementPolicy(revision=POLICY_REVISION).backend_options())
    root.mkdir()

    def write(name, value):
        data = value.encode() if isinstance(value, str) else value
        (root / name).write_bytes(data)
        import hashlib

        return ArtifactReference(name, hashlib.sha256(data).hexdigest())

    parent_entries = []
    for i, seed in enumerate(request.seeds):
        # Complementary imperfect donors. Tests use known fixture coordinates,
        # not checkpoint/scientific execution evidence.
        text = shift(good, "65" if i == 0 else "69", 1.5, atom_name="C")
        coordinates = write(f"parent-{i}.pdb", text)
        raw = write(f"raw-{i}.pdb", text)
        stage = StagedRefinementEvidence(digest(request.to_json()), raw, coordinates,
            "generated-only", stage1=coordinates, stage1_input_sha256=raw.sha256)
        staged = write(f"stage-{i}.json", stage.to_json())
        base = SamplerTrace.from_json(_trace_bytes(request.fixed_atoms, request.generated_atoms).decode())
        trace = write(f"trace-{i}.json", bind_evidence(base, staged).to_json())
        candidate = RunnerCandidate(f"parent-{i}", seed, raw, trace, request.generated_atoms,
            tuple(r for g in request.gaps for r in g.generated_residues), None, "unit-stub")
        result = write(f"runner-{i}.json", _runner_result((candidate,)).to_json())
        parent_entries.append(PortfolioParentArtifacts(candidate.candidate_id, seed, coordinates, trace, staged, raw, result))
    choices = (GapChoice(0, 0 if repair else 1), GapChoice(1, 0))
    selected = shift(good, "65", 1.5, atom_name="C") if repair else good
    selected_ref = write("selected.pdb", selected)
    protocol = {"kind": "experimental-multiparent-gap-portfolio", "stage": "generated-only", "anchor_mobility": False,
        "seeds": list(request.seeds), "max_evaluations": 8, "beam_width": 2}
    selected_row = {"choices": [asdict(c) for c in choices], "coordinate_sha256": selected_ref.sha256}
    report = {"kind": "experimental-multiparent-gap-portfolio", "revision": "generated-only-gap-portfolio-v1",
        "publication_qualified": False, "selected": selected_row, "trials": [selected_row] * 3}
    values = {"request": write("request.json", request.to_json()), "source": write("source.pdb", source),
        "selected": selected_ref, "candidate": write("candidate.pdb", good),
        "search_protocol": write("protocol.json", json_bytes(protocol)),
        "search_report": write("report.json", json_bytes(report)),
        "selection_source": write("selection.py", "# unit selector snapshot\n"),
        "parents": tuple(parent_entries), "choices": choices}
    if repair:
        from dvbfixer.model.diffusion.gap_portfolio import plan_gap_repairs
        from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

        scope = plan_gap_repairs(request, assessment_adapter(request, source)(selected))
        operation = {"operation": "generated-only-issue-local-refinement-v1", "publication_qualified": False,
            "input_sha256": selected_ref.sha256, "output_sha256": digest(good),
            "scope": asdict(scope), "mobile_atoms": [asdict(a) for a in request.generated_atoms if a.residue_number == "65"],
            "refinement_seed": 7, "max_iterations": 250, "restart_count": 4}
        values.update(repair_report=write("repair.json", json_bytes(operation)),
            repair_selector_source=write("repair-selector.py", "# unit repair snapshot\n"),
            repair_kernel_source=write("kernel.py", "# unit kernel snapshot\n"),
            repair_conditioning_source=write("conditioning.py", "# unit conditioning snapshot\n"))
    return PortfolioEvidence(1, **values), request


@pytest.mark.parametrize("repair", [False, True])
def test_publication_readback_and_exact_resume_preserve_original_artifacts(tmp_path, repair):
    workspace = tmp_path / "run"
    evidence, request = setup_evidence(workspace, repair=repair)
    original = {p.name: p.read_bytes() for p in workspace.iterdir()}
    bundle = tmp_path / "bundle"
    publish_portfolio_bundle(evidence, workspace, bundle)
    published = read_portfolio_bundle(bundle)
    assert published.index.kind == "experimental-gap-portfolio-bundle"
    assert not published.index.public_profile_enabled
    assert published.dat.added_keys() == {(a.chain, a.residue_number, a.insertion_code, a.atom_name)
        for a in request.generated_atoms}
    assert not published.dat.added_keys() & {(a.chain, a.residue_number, a.insertion_code, a.atom_name)
        for a in request.fixed_atoms}
    assert resume_portfolio_bundle(bundle, evidence).candidate_bytes == (workspace / "candidate.pdb").read_bytes()
    assert original == {p.name: p.read_bytes() for p in workspace.iterdir()}
    changed = replace(evidence, selection_source=replace(evidence.selection_source, sha256="1" * 64))
    with pytest.raises(DiffusionContractError, match="resume input"):
        resume_portfolio_bundle(bundle, changed)
    with pytest.raises(FileExistsError):
        publish_portfolio_bundle(evidence, workspace, bundle)


def test_changed_receipt_cannot_override_geometry_and_generated_ownership(tmp_path):
    workspace = tmp_path / "run"
    evidence, _ = setup_evidence(workspace)
    bundle = tmp_path / "bundle"
    publish_portfolio_bundle(evidence, workspace, bundle)
    # Rehashing an incorrect .dat into the manifest does not grant ownership.
    dat = DatRecord.load(bundle / "candidate.dat")
    dat.added_atoms[0]["resid"] = "64"
    data = json_bytes(dat.to_dict())
    (bundle / "candidate.dat").write_bytes(data)
    index = strict_json((bundle / "bundle.json").read_bytes())
    index["dat"]["sha256"] = digest(data.decode())
    (bundle / "bundle.json").write_bytes(json_bytes(index))
    with pytest.raises(DiffusionContractError, match="exact generated-only"):
        read_portfolio_bundle(bundle)


@pytest.mark.parametrize("kind", ["symlink", "hardlink", "extra", "orphan"])
def test_bundle_requires_exact_private_payloads_and_committed_readiness(tmp_path, kind):
    workspace = tmp_path / "run"
    evidence, _ = setup_evidence(workspace)
    bundle = tmp_path / "bundle"
    if kind == "orphan":
        bundle.mkdir()
        (bundle / "diagnostic.txt").write_text("crash-left evidence")
        with pytest.raises(FileExistsError):
            publish_portfolio_bundle(evidence, workspace, bundle)
        with pytest.raises(FileNotFoundError):
            resume_portfolio_bundle(bundle, evidence)
        assert (bundle / "diagnostic.txt").is_file()
        return
    publish_portfolio_bundle(evidence, workspace, bundle)
    if kind == "extra":
        (bundle / "undeclared").write_text("not indexed")
    else:
        candidate = bundle / "candidate.pdb"
        kept = tmp_path / "kept.pdb"
        candidate.rename(kept)
        if kind == "symlink":
            candidate.symlink_to(kept)
        else:
            candidate.hardlink_to(kept)
    with pytest.raises((ValueError, DiffusionContractError)):
        read_portfolio_bundle(bundle)


def test_malformed_or_failed_results_cannot_reserve_a_bundle(tmp_path):
    workspace = tmp_path / "run"
    evidence, _ = setup_evidence(workspace)
    bad = shift((workspace / "candidate.pdb").read_text(), "69", 1.5, atom_name="C")
    (workspace / "candidate.pdb").write_text(bad)
    evidence = replace(evidence, candidate=replace(evidence.candidate, sha256=digest(bad)))
    with pytest.raises(DiffusionContractError, match="explicit repair"):
        publish_portfolio_bundle(evidence, workspace, tmp_path / "bundle")
    assert not (tmp_path / "bundle").exists()
    with pytest.raises(DiffusionContractError):
        replace(evidence, public_profile_enabled=True)
    assert portfolio_resume_key(evidence) != portfolio_resume_key(replace(evidence,
        choices=(GapChoice(0, 0), GapChoice(1, 0))))


def test_public_sampling_and_completion_readers_do_not_adopt_portfolio_kind(tmp_path):
    from dvbfixer.model.diffusion.atom_completion_publication import read_completion_bundle
    from dvbfixer.model.diffusion.publication import PublicationError, read_published_bundle

    workspace = tmp_path / "run"
    evidence, request = setup_evidence(workspace)
    bundle = tmp_path / "bundle"
    publish_portfolio_bundle(evidence, workspace, bundle)
    with pytest.raises(PublicationError):
        read_published_bundle(bundle, request=request, source_root=workspace, work_parent=tmp_path)
    with pytest.raises(DiffusionContractError):
        read_completion_bundle(bundle)


def test_read_only_resume_never_calls_a_selector_or_minimizer(tmp_path, monkeypatch):
    workspace = tmp_path / "run"
    evidence, _ = setup_evidence(workspace)
    bundle = tmp_path / "bundle"
    publish_portfolio_bundle(evidence, workspace, bundle)
    before = {p.name: (p.read_bytes(), p.stat().st_mtime_ns) for p in bundle.iterdir()}

    def forbidden(*args, **kwargs):
        pytest.fail("resume attempted scientific execution")

    monkeypatch.setattr("dvbfixer.model.diffusion.gap_portfolio.search_gap_portfolio", forbidden)
    monkeypatch.setattr("dvbfixer.model.diffusion.boundary_refinement.refine_generated_region", forbidden)
    resume_portfolio_bundle(bundle, evidence)
    assert before == {p.name: (p.read_bytes(), p.stat().st_mtime_ns) for p in bundle.iterdir()}


def test_rehashed_forged_validation_receipt_fails_independent_readback(tmp_path):
    workspace = tmp_path / "run"
    evidence, _ = setup_evidence(workspace)
    bundle = tmp_path / "bundle"
    publish_portfolio_bundle(evidence, workspace, bundle)
    receipt = strict_json((bundle / "validation.json").read_bytes())
    receipt["preparation_verified"] = True
    data = json_bytes(receipt)
    (bundle / "validation.json").write_bytes(data)
    index = strict_json((bundle / "bundle.json").read_bytes())
    index["validation"]["sha256"] = digest(data.decode())
    (bundle / "bundle.json").write_bytes(json_bytes(index))
    with pytest.raises(DiffusionContractError, match="receipt differs"):
        read_portfolio_bundle(bundle)


def test_preparation_binding_cannot_be_omitted_or_invented(tmp_path):
    from dvbfixer.model.diffusion.contract import BackendOption
    from dvbfixer.model.diffusion.portfolio_contract import PortfolioPreparationArtifacts

    workspace = tmp_path / "run"
    evidence, request = setup_evidence(workspace)
    dummy = []
    for i in range(4):
        name = f"fake-preparation-{i}"
        (workspace / name).write_bytes(b"dummy")
        dummy.append(ArtifactReference(name, digest("dummy")))
    with pytest.raises(DiffusionContractError, match="source preparation"):
        publish_portfolio_bundle(replace(evidence, preparation=PortfolioPreparationArtifacts(*dummy)),
            workspace, tmp_path / "invented")
    request = replace(request, backend_options=(*request.backend_options,
        BackendOption("input_preparation.revision", "bounded-observed-input-repair-v1"),
        BackendOption("input_preparation.original_sha256", evidence.source.sha256),
        BackendOption("input_preparation.report_sha256", "0" * 64)))
    (workspace / "request.json").write_text(request.to_json())
    evidence = replace(evidence, request=replace(evidence.request, sha256=digest(request.to_json())))
    with pytest.raises(DiffusionContractError, match="source preparation"):
        publish_portfolio_bundle(evidence, workspace, tmp_path / "omitted")


def test_manifest_commit_failure_preserves_inputs_and_never_establishes_readiness(tmp_path, monkeypatch):
    import dvbfixer.model.diffusion.bundle_commit as commit

    workspace = tmp_path / "run"
    evidence, _ = setup_evidence(workspace)
    original = {p.name: p.read_bytes() for p in workspace.iterdir()}
    transfer = commit._transfer

    def fault(staging, name, descriptor, linked):
        if name == "bundle.json":
            raise OSError("injected failure before readiness commit")
        return transfer(staging, name, descriptor, linked)

    monkeypatch.setattr(commit, "_transfer", fault)
    with pytest.raises(OSError, match="injected"):
        publish_portfolio_bundle(evidence, workspace, tmp_path / "bundle")
    assert not (tmp_path / "bundle/bundle.json").exists()
    assert list(tmp_path.glob(".portfolio-staging-*"))
    assert original == {p.name: p.read_bytes() for p in workspace.iterdir()}


def test_generated_dat_is_resolvable_by_downstream_identity_after_topology_load(tmp_path):
    from io import StringIO

    from openmm.app import PDBFile

    from dvbfixer.minimize.pipeline import resolve_new_atom_indices

    workspace = tmp_path / "run"
    evidence, request = setup_evidence(workspace)
    bundle = tmp_path / "bundle"
    publish_portfolio_bundle(evidence, workspace, bundle)
    published = read_portfolio_bundle(bundle)
    pdb = PDBFile(StringIO(published.candidate_bytes.decode()))
    indices = resolve_new_atom_indices(pdb.topology, published.dat.added_keys())
    expected = {a.index for a in pdb.topology.atoms()
        if (a.residue.chain.id, str(a.residue.id), a.residue.insertionCode.strip(), a.name) in
           {(x.chain, x.residue_number, x.insertion_code, x.atom_name) for x in request.generated_atoms}}
    assert indices == expected


@pytest.mark.parametrize("mutation", ["scope", "atoms", "seed-bool", "iterations", "restarts", "input"])
def test_rehashed_repair_claim_cannot_expand_mobility_or_change_execution_budget(tmp_path, mutation):
    workspace = tmp_path / "run"
    evidence, request = setup_evidence(workspace, repair=True)
    operation = strict_json((workspace / "repair.json").read_bytes())
    if mutation == "scope":
        operation["scope"]["gap_indices"] = [0, 1]
    elif mutation == "atoms":
        operation["mobile_atoms"].append(asdict(next(a for a in request.generated_atoms if a.residue_number == "69")))
    elif mutation == "seed-bool":
        operation["refinement_seed"] = True
    elif mutation == "iterations":
        operation["max_iterations"] = 0
    elif mutation == "restarts":
        operation["restart_count"] = 9
    else:
        operation["input_sha256"] = "a" * 64
    data = json_bytes(operation)
    (workspace / "repair.json").write_bytes(data)
    evidence = replace(evidence, repair_report=replace(evidence.repair_report, sha256=digest(data.decode())))
    with pytest.raises(DiffusionContractError, match="independent mobility/budget"):
        publish_portfolio_bundle(evidence, workspace, tmp_path / "bundle")
    assert not (tmp_path / "bundle").exists()


def test_dat_and_downstream_resolution_keep_chain_case_and_insertion_codes():
    from io import StringIO
    from types import SimpleNamespace

    from openmm.app import PDBFile

    from dvbfixer.minimize.pipeline import resolve_new_atom_indices
    from dvbfixer.model.diffusion.contract import AtomIdentity
    from dvbfixer.model.diffusion.portfolio_publication import portfolio_dat_record

    atoms = [AtomIdentity(chain, "82", code, "CA") for chain in ("D", "d") for code in ("", "A")]
    lines = [f"ATOM  {i + 1:5d}  CA  ALA {a.chain}{82:4d}{a.insertion_code or ' '}   "
        f"{float(i):8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 10.00           C  \n" for i, a in enumerate(atoms)]
    text = "".join(lines) + "END\n"
    request = SimpleNamespace(generated_atoms=(atoms[1], atoms[2]))
    dat = portfolio_dat_record(request, text.encode())
    assert dat.added_keys() == {("D", "82", "A", "CA"), ("d", "82", "", "CA")}
    pdb = PDBFile(StringIO(text))
    assert resolve_new_atom_indices(pdb.topology, dat.added_keys()) == {1, 2}


def test_portfolio_json_and_versions_fail_closed():
    from dvbfixer.model.diffusion.portfolio_contract import PortfolioBundleIndex

    with pytest.raises(DiffusionContractError, match="duplicate"):
        strict_json(b'{"kind": "a", "kind": "b"}')
    with pytest.raises(DiffusionContractError, match="nonfinite"):
        strict_json(b'{"value": NaN}')
    with pytest.raises(DiffusionContractError):
        PortfolioBundleIndex.from_dict({"schema_version": 4, "kind": "experimental-gap-portfolio-bundle"})


def test_bound_reader_rejects_oversized_payload_before_any_scientific_audit(tmp_path, monkeypatch):
    import dvbfixer.model.diffusion.portfolio_evidence as module

    workspace = tmp_path / "run"
    evidence, _ = setup_evidence(workspace)
    monkeypatch.setattr(module, "MAX_PORTFOLIO_TOTAL_BYTES", 1)
    with pytest.raises(ValueError, match="bounded"):
        publish_portfolio_bundle(evidence, workspace, tmp_path / "bundle")
    assert not (tmp_path / "bundle").exists()


def rewrite_parent_stage(root, evidence, index, text):
    """Adversarially rehash a claimed stage, including its trace/evidence chain."""
    parent = evidence.parents[index]
    (root / parent.coordinates.path).write_text(text)
    coordinates = replace(parent.coordinates, sha256=digest(text))
    staged = StagedRefinementEvidence.from_json((root / parent.staged_evidence.path).read_text())
    staged = replace(staged, selected=coordinates, stage1=coordinates)
    staged_text = staged.to_json()
    (root / parent.staged_evidence.path).write_text(staged_text)
    staged_reference = replace(parent.staged_evidence, sha256=digest(staged_text))
    trace = SamplerTrace.from_json((root / parent.sampler_trace.path).read_text())
    trace_text = bind_evidence(trace, staged_reference).to_json()
    (root / parent.sampler_trace.path).write_text(trace_text)
    changed = replace(parent, coordinates=coordinates, staged_evidence=staged_reference,
        sampler_trace=replace(parent.sampler_trace, sha256=digest(trace_text)))
    return replace(evidence, parents=evidence.parents[:index] + (changed,) + evidence.parents[index + 1:])


def test_unused_donor_cannot_hide_noncoordinate_changes_in_its_stage(tmp_path):
    workspace = tmp_path / "run"
    evidence, _ = setup_evidence(workspace)
    text = (workspace / evidence.parents[1].coordinates.path).read_text()
    # Gap 69 from donor 1 is NOT selected; that does not authorize altered B.
    changed = "".join(line[:60] + "999.00" + line[66:]
        if line.startswith("ATOM  ") and line[22:26].strip() == "69" else line
        for line in text.splitlines(keepends=True))
    evidence = rewrite_parent_stage(workspace, evidence, 1, changed)
    with pytest.raises(DiffusionContractError, match="stage-1.*projection"):
        publish_portfolio_bundle(evidence, workspace, tmp_path / "bundle")


def test_unused_donor_must_still_have_a_globally_all_l_stage(tmp_path):
    import numpy as np

    from dvbfixer.ffutils.geometry import ChiralityError
    from dvbfixer.model.diffusion.contract import AtomIdentity
    from dvbfixer.model.diffusion.gap_portfolio import _coordinates
    from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates

    workspace = tmp_path / "run"
    evidence, _ = setup_evidence(workspace)
    text = (workspace / evidence.parents[0].coordinates.path).read_text()
    xyz, _ = _coordinates(text)
    n, ca, c, cb = (xyz[AtomIdentity("C", "65", "", name)] for name in ("N", "CA", "C", "CB"))
    normal = np.cross(n - ca, c - ca)
    reflected = cb - 2 * np.dot(cb - ca, normal) / np.dot(normal, normal) * normal
    text = rewrite_generated_coordinates(text, {AtomIdentity("C", "65", "", "CB"): reflected})
    evidence = rewrite_parent_stage(workspace, evidence, 0, text)
    with pytest.raises(ChiralityError):
        publish_portfolio_bundle(evidence, workspace, tmp_path / "bundle")
