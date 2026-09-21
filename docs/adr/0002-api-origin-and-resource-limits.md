# ADR 0002: API Origin And Resource Limits

Status: accepted

Date: 2026-09-21

## Context

Authentication and workspace authorization do not bound browser origins,
request storage, retained workspace data, or simultaneous scientific child
processes. The Vite and standalone hosts must apply the same policy without
moving scientific behavior into Node.

## Decision

Install one host-neutral CORS middleware before authentication. Requests without
`Origin` remain available to non-browser clients. Browser requests are allowed
from the exact same origin or an exact canonical HTTP(S) origin listed in
`DVBFIXER_CORS_ALLOWED_ORIGINS`. Reject wildcard, `null`, malformed, duplicate,
private-network, and otherwise unapproved requests. Do not enable credentialed
CORS.

Same-origin derivation uses the backend connection and `Host` header. A TLS
terminating or Host-rewriting reverse proxy must explicitly allowlist its public
origin; forwarded headers are intentionally not trusted.

Limit workspace import bodies with `DVBFIXER_GUI_MAX_UPLOAD_BYTES` (256 MiB by
default). Before manifest or import publication, count all regular-file logical
bytes in the workspace, including hidden, failed, and trash data, and enforce
`DVBFIXER_GUI_WORKSPACE_QUOTA_BYTES` (5 GiB by default; `0` disables it). Reject
symlinks during accounting. Return `413` for oversized request bodies and `507`
for exhausted workspace storage.

Place a process-wide FIFO semaphore at the shared child-spawn boundary.
`DVBFIXER_MAX_CONCURRENT_PROCESSES` defaults to one and is capped at 64. A queued
request can be cancelled without spawning; shutdown rejects queued requests,
terminates active children, and waits for permit release.

## Consequences

Vite and standalone behavior is consistent and invalid configuration fails at
startup. Every DVBFixer invocation shares one concurrency bound regardless of
the route that requested it.

These controls are intentionally single-process and application-level. They do
not provide TLS, request-rate limiting, OS CPU or memory isolation, hard
filesystem quotas, or coordination across server instances. Child processes can
create files before a publication-time workspace check; retained failed and
trash files consume quota until an operator or future retention policy removes
them.
