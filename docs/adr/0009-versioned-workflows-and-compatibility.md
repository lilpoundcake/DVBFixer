# ADR 0009: Versioned workflow API and compatibility

Date: 2026-09-21

Status: Accepted.

## Context

The original generic synchronous command route duplicated the managed-job
executor. Scientific commands can run long enough to outlive a request or its
serving process. The naming transform has a bounded synchronous contract;
general commands require persisted progress, cancellation, and artifact records.

## Decision

- `/api/v1/workspaces/{workspaceId}/naming-conversions` is the synchronous
  transform. All other current CLI operations use versioned managed jobs. The
  former `/api/jobs` and `/api/dvbfixer/:command` routes are retired; clients
  migrate to `/api/v1/workspaces/{workspaceId}/jobs`.
- Job requests, records, errors, and naming operations use runtime TypeBox
  schemas that produce OpenAPI 3.1 and a checked-in generated TypeScript client.
  Scientific reports are preserved in the job record, registered output
  artifacts, and bounded logs; stdout-only commands publish a text or JSON
  artifact. Each output records the request, content digests, and versions.
- V1 error responses contain `error.code`, `error.message`, and
  `error.requestId`; `error.details` is optional. Existing unversioned routes
  retain their response contracts. The common mappings are 400 invalid body,
  401 unauthenticated, 403 forbidden, 404 missing, 409 conflict, 413 oversized,
  415 unsupported media, 422 invalid scientific operation, and 500 unexpected
  server failure. Rate limits return 429, capacity 503, and timeout 504.
- Additive optional fields and new operations may ship within V1. A removal,
  required-field change, enum narrowing, status/code reassignment, or altered
  meaning requires V2. Publish migration guidance and keep V1 available for at
  least one minor release after announcing a V2 deprecation. Never repurpose an
  existing error code; reserve new codes for distinct conditions.
- Managed jobs persist their queued/running/terminal record before responding.
  A workspace-scoped durable run lock coordinates processes on **one host** with
  one local atomic filesystem; cancellation markers and SSE polling bridge
  processes. Stale work after a process crash becomes an explicit failed job,
  not an implicit retry. Multiple hosts, network filesystems, and distributed
  admission require a future external scheduler and shared state store.

## Consequences

The GUI uses generated path/request/response types for managed operations and
retains the authenticated SSE transport. The CLI spec remains available to
render command forms. Multi-process execution is limited to the documented
single-host storage boundary; process-local request rate and child admission
limits are still enforced independently by each instance.
