# DVBFixer API

DVBFixer exposes an initial versioned HTTP operation through a host-neutral Node
route composition shared by the Vite development adapter and the bundled
standalone server. The standalone host is a local-workspace boundary, not yet an
internet-facing production service: authentication, principal-to-workspace
authorization, CORS policy, quotas, and multi-instance locking remain planned
work.

The OpenAPI 3.1 document is available at:

```http
GET /api/v1/openapi.json
```

Runtime validation and OpenAPI are generated from the same TypeBox schemas in
`gui/server/naming-api-schema.ts`.

## Naming Conversion

```http
POST /api/v1/workspaces/{workspaceId}/naming-conversions
Content-Type: application/json
```

The operation converts PDB residue and atom names through `dvbfixer atom-names`.
It accepts an artifact ID, never a filesystem path, and resolves the artifact
inside the named workspace with traversal and symlink containment checks.

```json
{
  "inputArtifactId": "2a4d4d4e-2c82-4ea0-98ef-92df4cdd57f0",
  "target": {
    "forceField": "amber",
    "profile": "gromacs"
  },
  "variantOverrides": [
    {
      "chainId": "H",
      "residueNumber": "82",
      "insertionCode": "A",
      "variant": "HIE"
    }
  ],
  "outputName": "complex_gromacs.pdb",
  "dryRun": false
}
```

`variantOverrides`, `outputName`, and `dryRun` are optional. Overrides use exact,
case-sensitive chain, residue-number, and insertion-code identity. V1 accepts
PDB artifacts (`.pdb` or `.ent`) only. The server generates an output name when
one is omitted.

A successful conversion returns `201` and registers exactly one new structure
artifact. Its manifest metadata records the source artifact ID, request options,
DVBFixer version, operation/report schema versions, and source/report/output
SHA-256 digests. Helper reports, overrides, and logs are removed before
successful publication. A successful dry run returns `200`, includes the same conversion
result and candidate digest, and does not change the manifest or retain output.

All failures use a versioned error envelope with a request ID:

```json
{
  "error": {
    "code": "NAMING_COLLISION",
    "message": "Target atom name already exists",
    "details": {},
    "requestId": "req_..."
  }
}
```

Scientific conversion failures use their stable Python error codes and HTTP
`422`. Transport and workspace failures use API-specific codes. Unexpected
subprocess and filesystem details are not returned to clients.

## Limits

The JSON body uses the shared 2 MiB request limit. Naming-specific defaults are:

| Environment variable | Default | Purpose |
|---|---:|---|
| `DVBFIXER_NAMING_MAX_SOURCE_BYTES` | 50 MiB | Maximum source artifact size |
| `DVBFIXER_NAMING_MAX_REPORT_BYTES` | 20 MiB | Maximum CLI report size |
| `DVBFIXER_NAMING_TIMEOUT_MS` | 60,000 ms | Subprocess timeout |

The subprocess receives `SIGTERM` on timeout and escalates to `SIGKILL` after a
short grace period. Failed operation directories move below `runs/_failed` and
are never registered as visible artifacts.

## Standalone Host

From `gui/`, build and start the browser client and API together:

```bash
npm run build
npm start
```

The server defaults to `127.0.0.1:5173`, serves only the built `dist/` tree as
static content, and keeps workspace files behind contained API routes. Unknown
`/api/*` paths return JSON `404` responses rather than the SPA shell. Shutdown
stops accepting requests, cancels tracked DVBFixer subprocesses with signal
escalation, drains HTTP work, and then closes PostgreSQL.

Remote binding requires `DVBFIXER_ALLOW_INSECURE_REMOTE=1` as an explicit
acknowledgment. It remains unsupported for untrusted networks until the planned
security controls are implemented.
