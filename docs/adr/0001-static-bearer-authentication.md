# ADR 0001: Static Bearer Principals And Workspace Ownership

Status: Accepted

Date: 2026-09-18

## Context

The standalone Node host exposes workspace files and scientific execution APIs.
Loopback binding and Host-header checks reduce local exposure, but they do not
identify callers or isolate workspaces when remote access is explicitly enabled.
The first authentication mechanism must remain small, auditable, host-neutral,
and compatible with a later OIDC or reverse-proxy adapter.

## Decision

The V1 server authenticates opaque bearer tokens against statically configured
SHA-256 digests. `DVBFIXER_AUTH_PRINCIPALS` is a versioned JSON object mapping
unique, case-sensitive principal IDs to unique token digests. Raw tokens are
never stored in manifests or server configuration.

Authentication runs once at the shared `api-routes.ts` composition boundary, so
Vite and the standalone host enforce the same policy. Only `GET /api/health` and
`GET /api/v1/openapi.json` are public. All other API routes require a bearer
credential. The browser keeps its token in `sessionStorage`, sends it through a
shared API client, and never places it in a URL.

Workspace manifests use schema version 2:

- `ownerPrincipalId` identifies the immutable owner.
- `acl` grants a principal either `reader` or `writer` access.
- Owners can read, write, delete, and replace the ACL.
- Writers can read and mutate workspace content and run workflows.
- Readers can inspect and download workspace content but cannot mutate it.
- Unlisted callers receive `404`, preventing workspace enumeration.
- Authenticated readers receive `403` for attempted writes.

Existing V1 manifests are assigned to `DVBFIXER_LEGACY_WORKSPACE_OWNER`. A
single configured principal is inferred; multiple principals require an
explicit owner. The zero-configuration loopback development mode uses a local
principal named `local`. When authentication is later enabled, that provisional
owner is transferred deterministically to the configured legacy owner.

Non-loopback standalone binding requires both a valid principal configuration
and `DVBFIXER_ALLOW_INSECURE_REMOTE=1`.

## Consequences

- Scientific policy remains in Python; authentication and workspace policy stay
  in the Node application shell.
- Principal identity is request-scoped and is never captured globally by a
  route registration.
- Native `EventSource` and direct API download links are replaced by
  authenticated fetch-based streams and Blob downloads.
- Removing or renaming an owner principal requires an explicit ownership
  migration; the server does not silently let the first caller claim data.
- Static bearer authentication does not provide TLS, restrictive CORS, quotas,
  audit retention, token rotation endpoints, or multi-instance coordination.
  Remote hosting remains unsupported on untrusted networks until those controls
  are implemented.
