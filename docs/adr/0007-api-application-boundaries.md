# ADR 0007: API Application Boundaries

## Status

Accepted.

## Decision

Node and TypeScript own the HTTP transport, authentication, workspace and
artifact registry, job lifecycle, and host integration. Python remains the sole
owner of scientific naming and structure-transformation policy. HTTP handlers
invoke public `dvbfixer` commands and validate their structured outputs rather
than reproducing scientific rules in Node.

API routes compose through the host-neutral `gui/server/api-routes.ts` boundary.
The Vite development plugin and standalone Node server are adapters over that
same composition root; there is no duplicate Python web service.

`dvbfixer atom-names` is the dedicated force-field naming command. The V1 naming
operation is a synchronous deterministic transform whose input is a workspace
artifact ID, never a caller-supplied server path. Long-running scientific work
uses managed asynchronous jobs. Naming rejects multiple PDB `MODEL` blocks
rather than selecting or combining models implicitly.

## Consequences

- Scientific behavior is tested and versioned in the Python package.
- Node owns authorization, quotas, persistence, publication, and transport
  concerns without becoming a second scientific implementation.
- Vite and standalone hosting share middleware order and route behavior.
- New V1 operations must classify themselves as deterministic transforms or
  managed workflows and resolve server-side inputs through artifact identity.
