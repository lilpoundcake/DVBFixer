# Implementation Checklist

This is the canonical completion tracker for the plans in this directory.
Check an item only after its implementation and required verification are
complete. Detailed requirements remain in these authoritative sources:

- [API roadmap](api-roadmap.md) for API implementation and deployment.
- [DDD documentation roadmap](ddd-documentation-roadmap.md) for agent maps and domain documentation.
- [Continuation context](api-ddd-continuation-context.md) for historical evidence only.

Status is binary: `[x]` means complete and verified; `[ ]` means remaining,
partially complete, optional but not accepted, or not yet verified.

## API Phase 0: Decisions And Baseline

### Accepted Decisions

- [ ] Accept an ADR making Node/TypeScript the HTTP, workspace, artifact, and job shell.
- [ ] Accept an ADR making Python the sole owner of scientific naming policy.
- [ ] Accept an ADR for host-neutral route composition rather than a duplicate Python web service.
- [ ] Record `dvbfixer atom-names` as the dedicated naming command.
- [ ] Record naming conversion as synchronous and long-running workflows as asynchronous jobs.
- [ ] Record workspace artifact IDs as the only V1 server-side input reference.
- [ ] Record the V1 policy rejecting multiple `MODEL` blocks.

### Open Maintainer Decisions

- [x] Select an authentication model.
- [x] Define principal-to-workspace authorization semantics.
- [ ] Decide whether V1 needs direct upload/download conversion.
- [ ] Decide whether HTTP should accept a JSON override artifact as well as inline overrides.
- [ ] Decide whether `tleap-reduce` must honor the force-field naming option.

### Baseline Evidence

- [x] Add regression coverage for naming blockers 1-6 and 9 from the API roadmap.
- [ ] Complete end-to-end nucleic-acid naming coverage.
- [ ] Capture and document any still-needed representative golden inputs.
- [ ] Complete ADR and test evidence distinguishing current behavior from intended V1 behavior.

## API Phase 1: Naming Application Service

### Application Boundary

- [x] Add typed naming request, result, diagnostic, target, profile, override, identity, and report vocabulary.
- [x] Implement `convert_force_field_naming` as PDB text in and transformed text plus report out.
- [x] Keep file I/O outside the pure service.
- [x] Retain the atomic in-place compatibility adapter for existing pipeline callers.
- [x] Report changed atoms, residues, variants, and atom names separately.
- [x] Keep naming conversion distinct from canonical `rename` and chemical parameterization.

### V1 Semantics And Correctness

- [x] Accept PDB text and support AMBER and CHARMM with the GROMACS profile.
- [x] Preserve exact, case-sensitive chain IDs and insertion codes.
- [x] Use typed three-field override identities and restrict legacy two-field identities to blank insertion codes.
- [x] Apply explicit override, explicit source variant, then no inferred change precedence.
- [x] Do not infer HID/HIE/HIP from canonical HIS alone.
- [x] Validate before returning or committing output and leave sources untouched on failure.
- [x] Reject multiple `MODEL` blocks.
- [x] Detect source and target atom-name collisions per residue.
- [x] Update correlated `ANISOU` and parseable full `TER` records.
- [x] Preserve unrelated records and non-name atom fields.
- [x] Support four-character CHARMM names and LYN/LSN hydrogen pairs idempotently.
- [x] Remove ambiguous two-field override fallback from the public service.
- [ ] Add end-to-end DNA and RNA naming tests.
- [ ] Document unvalidated terminal 5-prime/3-prime hydroxyl cases as unsupported.
- [ ] Resolve `tleap-reduce` naming parity.
- [x] Document `standard` as suppressing shifts, not reversing GROMACS names.
- [x] Do not claim naming compatibility proves template or bond compatibility.

### Test Matrix

- [ ] Verify every `GROMACS_AMBER_ATOM_RENAMES` mapping.
- [ ] Verify AMBER and CHARMM caps and termini.
- [ ] Verify exact coexistence of residues such as `H:82` and `H:82A`.
- [ ] Verify one-residue and multi-chain proteins.
- [ ] Verify byte-identical repeat conversion across all supported mappings.
- [ ] Complete malformed-line coverage.
- [ ] Complete `ATOM`, `HETATM`, `ANISOU`, `TER`, and untouched-record matrix coverage.
- [ ] Verify empty/no-op conversion and report counts.

## API Phase 2: Dedicated CLI Adapter

- [x] Add `dvbfixer atom-names` with one input and one explicit output.
- [x] Add `--target-ff {amber,charmm}` and `--profile gromacs`.
- [x] Support variant overrides from JSON, `--dry-run`, and `--report-json`.
- [x] Call the pure naming service and keep the command distinct from `rename` and `convert`.
- [x] Register the command and regenerate CLI and GUI command references.
- [x] Add subprocess-level tests and deterministic conversion.
- [ ] Explicitly verify source preservation at CLI level.
- [ ] Explicitly validate the report schema at CLI level.
- [ ] Explicitly verify non-zero exit and no destination on unsafe input.
- [ ] Explicitly verify the selected CIF boundary behavior.

## API Phase 3: Versioned Naming API

### Core Vertical Slice

- [x] Add a V1 naming-conversion route with TypeBox runtime schemas and OpenAPI 3.1.
- [x] Resolve inputs by workspace artifact ID.
- [x] Invoke the Python command instead of duplicating naming rules in Node.
- [x] Validate bounded command reports and output digests.
- [x] Reload the current manifest before publication.
- [x] Register exactly one output with reproducible provenance and none for dry runs or failures.

### Request And Artifact Transaction

- [x] Authenticate the caller.
- [x] Authorize the principal before workspace path resolution.
- [x] Validate the request body and resolve `inputArtifactId` through the manifest.
- [ ] Explicitly verify accepted source artifact types.
- [ ] Enforce source artifact size limits before execution.
- [ ] Verify operation-directory privacy and cleanup behavior.
- [x] Parse the report, validate output existence/digest, and publish against the current manifest revision.
- [x] Return the output artifact and report.
- [ ] Define and verify retained failure-log behavior.

### Provenance, Observability, And Tests

- [x] Record source, target, profile, overrides, report, and digest provenance.
- [ ] Explicitly record command/service and DVBFixer package versions.
- [x] Add server-owned request IDs and structured API access logs.
- [x] Add bounded process-local service metrics.
- [ ] Add durable audit events and retention.
- [ ] Generate a client from OpenAPI.
- [ ] Verify stable validation and error envelopes across V1 routes.
- [x] Add authentication and authorization tests after implementation.
- [ ] Explicitly test path containment, source/request size limits, and naming timeouts.
- [x] Test artifact lookup, dry runs, one-artifact publication, concurrent manifests, and failures.
- [ ] Validate OpenAPI examples against runtime schemas.

## API Phase 4: Standalone Server

### Local Standalone Host

- [x] Extract host-neutral route composition and keep Vite as a thin adapter.
- [x] Add a bundled Node entry point serving the built GUI and API.
- [x] Bind to loopback by default and require explicit insecure remote acknowledgment.
- [x] Validate host, port, roots, and shutdown configuration.
- [x] Reject lexical and symlink-resolved static/data root overlap.
- [x] Return JSON for unknown API paths and use SPA fallback only for HTML requests.
- [x] Stop intake, reject new work, terminate children, drain HTTP, then close PostgreSQL.
- [x] Bound the emitted-server startup and shutdown smoke test.
- [x] Complete final review and verification of the standalone slice.
- [x] Commit the standalone slice.

### Public Deployment Hardening

- [x] Add authentication and principal-to-workspace authorization.
- [x] Add restrictive CORS.
- [x] Add upload and workspace limits plus process-wide child concurrency control.
- [x] Bound the process-wide child admission queue and reject overload before spawning.
- [x] Add bounded process-local request-rate limits before authentication.
- [ ] Add OS-level CPU, memory, and filesystem quotas.
  - [x] Add a fail-closed systemd/cgroup-v2 deployment profile and bounded-storage preflight.
  - [ ] Run privileged kernel-enforcement acceptance tests on the target deployment host.
- [x] Document current local deployment and security boundaries.
- [ ] Document a supported public deployment after security controls exist.
- [x] Add storage-level manifest locking or compare-and-swap.
- [ ] Replace process-local locks and event subscribers for multi-instance operation.
- [ ] Complete public-deployment smoke tests.

## API Phase 5: General DVBFixer API

### Route And Contract Policy

- [ ] Version managed-job routes under `/api/v1`.
- [ ] Retire the duplicate synchronous generic command route.
- [ ] Classify every operation as a deterministic transform or managed workflow.
- [ ] Keep minimization, modeling, and parameterization out of synchronous APIs.
- [ ] Preserve cancellation for long-running operations.
- [ ] Generate clients and contract tests from OpenAPI.
- [ ] Define API compatibility, versioning, and deprecation policy.
- [ ] Add durable scheduling, state, and event delivery before horizontal scaling.

### Operation Rollout

- [ ] Expose `diagnose`.
- [ ] Expose `conect`.
- [ ] Expose `renumber`.
- [ ] Version or migrate managed `prepare`, `minimize`, `model`, and `zbs`.

## API-Wide Quality Gates

- [ ] Use one runtime-schema source for every V1 route and OpenAPI description.
- [ ] Ensure no public route accepts unrestricted filesystem paths.
- [x] Keep scientific naming policy implemented once in Python.
- [ ] Record reproducible inputs, options, and versions for every output artifact.
- [x] Enforce workspace ownership before path resolution.
- [ ] Return structured reports for every exposed operation.
- [ ] Require durable workflow state before horizontal scaling.
- [ ] Adopt one common V1 error envelope and explicit compatibility policy.
- [ ] Test status mappings for 400, 401, 403, 404, 409, 413, 415, 422, and 500.
- [ ] Stabilize documented naming error codes.
- [ ] Run release verification: GUI checks, focused Python tests, generated-file checks, agent-doc checks, and the non-slow Python suite.

## DDD Documentation Roadmap

### Phases 0-2: Foundation And Domain Contexts

- [x] Define source precedence, statuses, references, delegation, and ownership rules.
- [x] Use dependency-free TOML maps and validate syntax, links, IDs, and evidence in CI.
- [x] Map force-field naming callers, boundaries, contracts, invariants, and known gaps.
- [x] Add Force-Field Naming, Structure Identity, and Parameterization contexts.
- [x] Add their task, contract, invariant, and change-group records.
- [x] Give every `dvbfixer.domain` export an active meaning or explicit vocabulary-only status.

### Phases 3-4: Artifacts And Pipelines

- [x] Document `.dat`, structure identity, connectivity, workspace, Homology, and topology artifacts.
- [x] Record owners, producers, consumers, keys, mutation policy, atomicity, cleanup, and failure behavior.
- [x] Add Structure Preparation, Topology Generation, Diagnostics, and Workflow Execution contexts.
- [x] Document boundaries without copying algorithms and reconcile cross-context contradictions.

### Phases 5-6: Tasks And Validation

- [x] Index high-risk `.dat`, rebuild, CONECT, topology, workspace, API, and fixture changes.
- [ ] Add new task records only when repeated high-risk work justifies them.
- [x] Validate TOML, metadata, paths, Markdown, linked records, evidence, sections, and indexed pages.
- [ ] Add symbol checks, stale-commit reporting, generated views, or PR checks only when justified.
- [x] Keep a general Python/TypeScript AST index out of the initial rollout.

### Ongoing Maintenance

- [x] Use `docs/agent/README.md` as the canonical agent-map entry point.
- [x] Maintain shared maps, contexts, validator, and tests through one coordinating editor.
- [x] Separate implemented behavior from vocabulary, proposals, research, and gaps.
- [x] Surface unresolved contradictions instead of choosing silently.
- [ ] Add context pages only when structured records cannot explain the concern.
- [x] Run agent-doc validation and `git diff --check` after every map change.
