# API and DDD continuation context

Last updated: 2026-09-18.

Purpose: give a future coding session enough verified context to continue the
API and DDD work without repeating repository discovery. This is a handoff
record, not a public contract. Proposed decisions remain proposals until an ADR
is accepted and code exists.

## Requested outcome

The initial request covered planning, followed by implementation of the first
LLM-oriented documentation slice:

1. Plan a DVBFixer API.
2. Prioritize an API for GROMACS-compatible atom naming.
3. Plan adoption of Domain-Driven Design documentation.
4. Preserve enough context for later implementation.
5. Add a validated task/contract/invariant map for coding agents.

Documents created in this pass:

- [DVBFixer API roadmap](api-roadmap.md)
- [DDD documentation roadmap](ddd-documentation-roadmap.md)
- this continuation context
- [`../agent/`](../agent/) knowledge maps and Force-field Naming context

No scientific production code, CLI arguments, or generated references were
changed. Documentation validation, its focused test, and the lightweight CI
check were added.

## Verified current state

### Naming implementation

- `src/dvbfixer/ffutils/ff_names.py::apply_variants_to_pdb_text` is the shared
  implementation relevant to GROMACS-compatible names.
- It mutates a PDB file in place and returns the number of changed atom lines.
- It accepts `target_ff="amber"` or `"charmm"` and an
  `include_gromacs_shifts` switch.
- It restores explicit residue variants and applies atom-name mappings for
  AMBER methylenes, termini, caps, nucleic acids, and CHARMM backbone/caps.
- It is called by legacy paths in prepare, minimize, and protonate.
- `zbs` forwards atom-naming options to prepare and minimize.
- The `tleap-reduce` early-return paths in prepare/protonate do not apply the
  final naming helper. A later minimize can mask this in some ZBS workflows.

Important distinction:

- `src/dvbfixer/rename.py` canonicalizes residue names such as `HIE` to `HIS`.
- It does not perform GROMACS-compatible atom naming.
- Do not expose or rename that operation as the planned naming conversion.

### Known naming risks found during discovery

- Terminal detection is chain-aware but not `MODEL`-aware.
- Public two-tuple variant keys can match insertion-code siblings through a
  fallback; new contracts must require insertion code explicitly.
- Target atom-name collisions are not rejected.
- `ANISOU` and `TER` labels are not rewritten with changed atom/residue labels.
- Four-character CHARMM residue names are written through a PDB compromise but
  later parsing reads the three-character slice, threatening idempotence.
- CHARMM `LYN` to `LSN` atom-pair conversion is not safely repeatable.
- End-to-end nucleic-acid naming coverage is incomplete; terminal hydroxyl
  cases are explicitly unverified.
- `standard` suppresses shifts; it is not a reverse conversion from existing
  GROMACS names.
- Naming compatibility is not proof of force-field template/bond compatibility.

Tests already worth extending:

- `tests/test_ff_names.py`
- `tests/test_variants_gromacs_lyn.py`
- `tests/test_terminal_caps.py`
- `tests/test_rename.py`

### Existing HTTP/backend implementation

- `gui/server/api-plugin.ts` is the current composition root and Vite middleware
  host.
- `gui/server/workspace-api.ts` owns versioned workspace manifests, contained
  paths, atomic writes, imports/downloads, and recoverable trash.
- `gui/server/managed-jobs.ts` owns persisted job records, one active run per
  workspace, cancellation, polling, and SSE.
- `gui/server/dvbfixer-runner.ts` invokes the CLI with timeout, abort,
  output-size bounds, and process-group termination.
- `gui/server/command-args.ts` allowlists command arguments from generated
  definitions.
- `gui/server/api-plugin.ts` still contains a duplicate synchronous generic
  `/api/dvbfixer/:command` path. The GUI primarily uses managed jobs.
- The server has no formal OpenAPI document, API version, shared runtime schema
  system, authentication, or workspace authorization.
- Active locks and SSE subscribers are process-local, so the job system is not
  ready for multi-instance deployment.
- The production GUI build is viewer-only; the backend is not shipped outside
  Vite development middleware.

Relevant backend tests:

- `gui/server/workspace-api.test.ts`
- `gui/server/managed-jobs.test.ts`
- `gui/server/dvbfixer-runner.test.ts`
- `gui/server/request-body.test.ts`
- `gui/server/generated-spec.test.ts`
- `gui/server/homology-api.test.ts`

### Existing domain model

- `src/dvbfixer/domain/structure_identity.py` owns immutable identity vocabulary
  and chain-ID allocation policy.
- `src/dvbfixer/domain/parameterization.py` owns parameterization route policy
  and the complex-cofactor guard.
- `src/dvbfixer/domain/__init__.py` exports those concepts.
- Active adapters include `src/dvbfixer/molecule_chains.py` and
  `src/dvbfixer/lig_params.py`.
- `docs/domain-model.md` and the first section of `ARCHITECTURE.md` explicitly
  call DDD partial and focused.
- `tests/test_scientific_domain.py` covers the existing domain policies.

Do not describe parsing, topology construction, atom matching, minimization, or
all five parameterization route labels as fully migrated domain implementations.

## Proposed decisions in the roadmaps

These are recommendations, not accepted ADRs:

1. Keep Node/TypeScript as the HTTP, workspace, artifact, and job shell.
2. Keep Python as the sole owner of scientific naming policy.
3. Extract route composition from Vite rather than create a duplicate Python
   web service.
4. Introduce a pure Python naming service: text in, transformed text plus a
   structured report out.
5. Keep file I/O in adapters and make API output non-destructive and atomic.
6. Add a dedicated CLI adapter rather than overload `dvbfixer rename`.
7. Make the first endpoint synchronous because naming conversion is a bounded
   text transformation; retain managed jobs for long-running workflows.
8. Make V1 workspace-scoped and reference input by artifact ID, never by a raw
   server path.
9. Publish `/api/v1` plus OpenAPI generated from the runtime validation schema.
10. Document Force-field Naming as the first new bounded context.

Decisions still requiring maintainer agreement:

- Final command name (`atom-names` is only a candidate).
- Runtime schema/OpenAPI library for the TypeScript server.
- Authentication model and principal-to-workspace authorization.
- Whether multi-model PDB is supported in V1 or rejected.
- Exact convention for four-character CHARMM residue names.
- Whether `ANISOU` and `TER` are rewritten or unsupported when affected.
- Whether direct upload/download conversion is needed in V1 or only after the
  workspace-artifact vertical slice.
- Whether explicit variant overrides are accepted inline, by JSON artifact, or
  both at the CLI boundary.

## Recommended next implementation session

Start with Phase 0 and Phase 1 of the API roadmap, not the HTTP route.

1. Read `docs/plans/api-roadmap.md` and accept or revise the open decisions.
2. Write ADRs for the naming boundary, command name, Node/Python split, and
   multi-model policy.
3. Add failing tests for model awareness, collisions, insertion codes,
   four-character CHARMM idempotence, and record consistency.
4. Extract `convert_pdb_naming` as a pure function with typed request/result,
   per-change records, and diagnostics.
5. Keep `apply_variants_to_pdb_text` as a compatibility wrapper while migrating
   prepare/minimize/protonate incrementally.
6. Run focused tests, then the non-slow suite and strict checks.

Do not begin by adding an HTTP handler around
`apply_variants_to_pdb_text`. Its in-place behavior, weak result contract, and
known edge cases are unsuitable as a public API boundary.

## Likely code touch points

Python core and CLI:

- `src/dvbfixer/ffutils/ff_names.py`
- a new naming application module, location to be decided after the bounded
  context ADR
- `src/dvbfixer/command_registry.py`
- `src/dvbfixer/cli.py`
- `src/dvbfixer/prepare/pipeline.py`
- `src/dvbfixer/minimize/pipeline.py`
- `src/dvbfixer/protonate.py`
- `src/dvbfixer/zbs.py`
- `tests/test_ff_names.py`
- new CLI/application tests

Generated/documentation surfaces after argparse changes:

- `scripts/gen_cli_reference.py`
- `scripts/gen_gui_spec.py`
- `docs/reference/`
- `gui/server/generated-dvbfixer-spec.ts`

HTTP server:

- `gui/server/api-plugin.ts`
- a new host-neutral API composition module
- a new production server entry point
- `gui/server/workspace-api.ts`
- `gui/server/managed-jobs.ts`
- `gui/server/dvbfixer-runner.ts`
- request/response schema and OpenAPI generation modules
- route and contract tests

DDD documentation:

- `docs/domain-model.md`
- `ARCHITECTURE.md`
- `docs/agent/` knowledge maps and context pages
- API ADRs and Force-field Naming context page first

## Repository rules that matter

- Preserve `keepIds=True` on `PDBFile.writeFile` calls.
- Never collapse `(chain, resid, icode)` identities to two-tuples.
- Chain IDs are case-sensitive.
- Do not hand-edit `docs/reference/*.md`; regenerate them.
- Scientific materialization belongs in Python, not React/Node route handlers.
- Use `DatRecord` for `.dat` files rather than hand-written JSON loading.
- Preserve workspace isolation and revision-conflict rebasing semantics.
- Do not modify unrelated untracked files. At discovery time the worktree
  already contained untracked `.claude/`, force-field directories, and
  `mini.yml`; they were not touched by this planning work.

## Verification commands for future code work

Focused Python tests should be run first. Candidate sequence:

```bash
pytest -q tests/test_ff_names.py tests/test_variants_gromacs_lyn.py tests/test_terminal_caps.py
python scripts/gen_cli_reference.py --check
python scripts/gen_gui_spec.py --check
ruff check src/dvbfixer
mypy src/dvbfixer/cli.py src/dvbfixer/ffutils src/dvbfixer/pdbutils src/dvbfixer/align.py
pytest -m 'not slow' -q
```

GUI/backend verification from `gui/`:

```bash
npm run typecheck
npm test -- --run
```

For documentation-only changes, validate relative Markdown links and inspect
`git diff --check` plus the final diff.

## Completion marker for this pass

- Repository discovery: complete.
- API plan: complete as a proposal.
- DDD documentation plan: revised for LLM navigation.
- Agent-map foundation, seven contexts, core artifact contracts, high-risk task
  groups, and cross-cutting invariants: implemented.
- Continuation context: complete.
- API implementation: not started.
- Accepted API ADRs: not started; create only after decisions are approved.
