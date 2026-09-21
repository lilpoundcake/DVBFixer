# ADR 0005: API Request Observability

## Status

Accepted.

## Decision

Install one request-observability middleware before CORS, rate limiting, and
authentication at the shared `/api` boundary. It generates a server-owned
`req_<UUID>` identifier, sets `X-Request-Id`, and makes V1 error envelopes reuse
that identifier. Client-provided request IDs are ignored.

When `DVBFIXER_ACCESS_LOG=json`, emit exactly one JSON line on response finish or
premature close. Records contain only a bounded method, coarse route template,
status, monotonic duration, outcome, and successfully authenticated principal.
They never contain headers, credentials, addresses, forwarded metadata, query
strings, bodies, raw paths, workspace/artifact identifiers, filenames,
subprocess output, or exception messages. Logging failures do not affect HTTP
handling.

## Scope

The log covers requests reaching the shared Vite/standalone API composition,
including CORS, rate-limit, authentication, route, and unknown-API outcomes. The
standalone loopback Host-header guard runs before this boundary and is excluded.
Static and SPA requests are excluded.

Access records go to stdout. The deployment manager owns forwarding, access,
rotation, and retention. These records are operational logs, not durable audit
events; metrics and an authorization audit model remain separate work.
