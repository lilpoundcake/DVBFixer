# Remove the Tarantino Mutation Cluster

## Status

Completed on branch `worktree-agent-ab3b721599cbb67b1`. The removal began in
`e0d68e7` and was completed after explicit approval of the compatibility policy:
historical `mutationIds` and `mutationsResolved` values remain readable as
immutable provenance, but no database lookup, new write, editor, or execution
path remains.

## Boundary

Removed:

- PostgreSQL runtime, lifecycle, dependencies, local Compose service, launcher,
  environment configuration, health fields, and deployment claims.
- The mutation library editor and `/api/mutations` CRUD routes.
- Antibody Engineer UI, route, orchestration, numbering/reference helpers, and
  active workspace writes.

Retained:

- The React/Mol* shell, viewers, workspace manifests and authorization,
  imports/trash, managed DVBFixer jobs, naming, Homology, alignment,
  diagnostics, security middleware, observability, quotas, and standalone host.
- Non-destructive import of retired per-workspace indexes. Existing artifacts
  with command `antibody-engineer` remain generic workspace artifacts.
- Valid historical `mutationIds`, `mutationsResolved`, and related engineer
  provenance. IDs are opaque and are never resolved against a database.

## Completion Audit

- `pg`, `@types/pg`, `@mui/x-data-grid`, and mutation-only drag/drop packages
  are absent from `gui/package.json` and its lockfile.
- `gui/docker-compose.yml`, the database-aware development launcher, both
  panels, the antibody pipeline, and antibody numbering/reference modules are
  deleted.
- Shared API composition has no database initialization or close operation and
  does not inspect `DATABASE_URL` or `DVBFIXER_MUTATIONS_BACKUP_FILE`.
- Retired routes reach the normal JSON API 404 boundary for every supported
  HTTP method. Retained health, OpenAPI, workspaces, Homology, naming, and
  managed-job routes remain composed.
- The active panel registry retains DVBFixer, Homology, viewers, Library,
  Workspace, and Alignment and excludes both retired panels.
- Artifact metadata patching does not accept either legacy mutation field;
  unrelated edits preserve values already present in historical manifests.
- Current documentation and the tracked `.dsp` inventory no longer claim an
  active PostgreSQL, mutation-library, or Antibody Engineer capability.

## Remaining References

Every retained active-tree occurrence is compatibility or test evidence:

- `gui/server/workspace-api.ts` validates and imports historical provenance.
- `gui/server/workspace-api.test.ts` proves migration, preservation, and the
  absence of a write path.
- `gui/server/standalone.test.ts` names retired routes only to prove 404s and
  names retired environment variables only to prove they are ignored.
- `gui/src/app-panels.test.ts` names retired panel IDs only to prove they are
  absent while retained workflow panels remain registered.
- `docs/agent/contexts/workflow-execution.md` documents the compatibility
  boundary.
- `docs/pipelines.md` uses the generic domain phrase "antibody engineering" in
  a retained Homology modeling recipe; it does not refer to the retired panel.
- The generated GUI command specification and Python modules retain DVBFixer's
  independent `prepare --mutate` and antibody-numbering functionality. They do
  not use the retired mutation library, routes, database, or engineer workflow.
- This status document records the retired design and audit vocabulary.

No external database or backup is modified or deleted by this change.

## Verification

Completed checks:

- `cd gui && npm ci`
- `cd gui && npm run typecheck`
- `cd gui && npm test -- --run` (30 files, 234 tests)
- `cd gui && npm run lint`
- `cd gui && npm run check:openapi`
- `cd gui && npm run build` (includes the standalone smoke)
- `cd gui && npm run smoke:server`
- `micromamba run -n dvbfixer python scripts/check_agent_docs.py`
- `micromamba run -n dvbfixer python scripts/check_versions.py`
- `micromamba run -n dvbfixer python scripts/gen_gui_spec.py --check`
- `git diff --check`

The unrelated CLI-reference check reports pre-existing drift for the generated
command pages at base `e0d68e7`. This retirement changes no argparse code or CLI
reference content, so those pages were not regenerated in this branch.
