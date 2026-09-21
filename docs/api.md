# DVBFixer API

DVBFixer exposes an initial versioned HTTP operation through a host-neutral Node
route composition shared by the Vite development adapter and the bundled
standalone server. Static bearer authentication and manifest-backed workspace
authorization, restrictive browser-origin handling, request-rate and workspace/
upload limits, and process-wide child concurrency limits are available. The host
is not yet an internet-facing production service: TLS, target-host resource
enforcement acceptance, audit retention, and distributed scheduling/event
delivery remain planned work.

Workspace manifest mutations use same-host cross-process filesystem locks and
compare revisions against the locked on-disk manifest. This prevents concurrent
Node processes sharing one local data root from losing manifest updates. The
lock is not a distributed lease and does not support multiple hosts or
filesystems without local atomic rename and reliable process identity.

A Linux systemd/cgroup-v2 profile for aggregate OS resource containment is
available in [`deployment.md`](deployment.md). Public deployment remains
unsupported until the target host passes its privileged enforcement checks and
the remaining TLS, audit, and durable-scheduling controls are supplied.

The OpenAPI 3.1 document is available at:

```http
GET /api/v1/openapi.json
```

Runtime validation and OpenAPI are generated from the same TypeBox schemas in
`gui/server/naming-api-schema.ts`.

## Authentication

`GET /api/health` and `GET /api/v1/openapi.json` are public. Every other API
route requires:

```http
Authorization: Bearer <43-character base64url token>
```

Missing or invalid credentials return `401` with `WWW-Authenticate: Bearer`.
The protected `GET /api/session` endpoint returns the authenticated principal.
Tokens are stored only in browser `sessionStorage` and are never sent in URLs.

Configure token digests, not raw tokens:

```bash
node -e "const c=require('node:crypto');const t=c.randomBytes(32).toString('base64url');console.log('token:',t);console.log('sha256:',c.createHash('sha256').update(t).digest('hex'))"
export DVBFIXER_AUTH_PRINCIPALS='{"version":1,"principals":[{"id":"alice","tokenSha256":"<sha256>"}]}'
```

Keep the printed token for the client and place only its digest in server
configuration. With multiple principals, set
`DVBFIXER_LEGACY_WORKSPACE_OWNER` to the principal that should own existing
workspaces.

Workspace manifests contain one owner and optional `reader`/`writer` ACL
entries. Unlisted principals receive `404`; readers receive `403` for writes.
Owners update sharing through the revisioned
`PATCH /api/workspaces/{workspaceId}/acl` endpoint.

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
PDB artifacts (`.pdb` or `.ent`) only. The artifact must be a visible workspace
`structure` whose contained source resolves to a regular file. Extension
matching is case-insensitive. The server generates an output name when one is
omitted.

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

Every request reaching the shared `/api` composition receives a server-generated
`X-Request-Id`; client-supplied values are ignored. V1 error envelopes reuse the
same ID. Allowed cross-origin responses expose this header.

## Access Logs

Set `DVBFIXER_ACCESS_LOG=json` to emit one compact JSON record to stdout when an
API response finishes or aborts. The default is `off`. Records use schema
`dvbfixer.http_access.v1` and include the request ID, method, coarse route label,
status, monotonic duration, completion outcome, and authenticated principal ID.

Logs intentionally omit headers, bearer credentials, client addresses,
forwarded headers, query strings, request/response bodies, workspace and
artifact IDs, filenames, subprocess output, exception messages, and raw paths.
Client-supplied request IDs are never trusted. Route labels are bounded templates
or coarse API groups to avoid sensitive high-cardinality data.

This is an API access log, not a durable authorization audit trail. It covers
the shared API composition, including CORS, rate-limit, authentication, public
endpoint, route, and unknown-API responses. A standalone loopback Host-header
rejection occurs before that boundary and is not included.

## Metrics

Set `DVBFIXER_METRICS=prometheus` to enable the protected `GET /api/metrics`
endpoint. The default is `off`. It returns Prometheus text format with bounded
HTTP request counters, route-only duration histograms, active-request gauges,
and DVBFixer child-process active/queued/limit/admission gauges.

Metrics never label by principal, request ID, client address, workspace,
artifact, job, filename, raw path, query, command arguments, exception text, or
subprocess output. The scrape excludes itself from HTTP metrics. Metrics are
process-local and reset on restart; scrape each instance independently and
aggregate externally. Scrapes remain authenticated and rate-limited.

Scientific conversion failures use their stable Python error codes and HTTP
`422`. Transport and workspace failures use API-specific codes. Unexpected
subprocess and filesystem details are not returned to clients.

## Limits

The JSON body uses the shared 2 MiB request limit. General server defaults are:

| Environment variable | Default | Purpose |
|---|---:|---|
| `DVBFIXER_GUI_MAX_UPLOAD_BYTES` | 256 MiB | Maximum workspace import request body |
| `DVBFIXER_GUI_WORKSPACE_QUOTA_BYTES` | 5 GiB | Logical bytes per workspace; `0` disables |
| `DVBFIXER_MAX_CONCURRENT_PROCESSES` | 1 | FIFO child-process limit, from 1 through 64 |
| `DVBFIXER_MAX_QUEUED_PROCESSES` | 16 | Maximum waiting child-process requests, from 1 through 256 |
| `DVBFIXER_RATE_LIMIT_REQUESTS` | 120 | Requests per direct client address/window; `0` disables |
| `DVBFIXER_RATE_LIMIT_WINDOW_MS` | 60,000 ms | Fixed rate-limit window, 1–3,600 seconds |
| `DVBFIXER_RATE_LIMIT_MAX_KEYS` | 10,000 | Maximum tracked client addresses |
| `DVBFIXER_ACCESS_LOG` | `off` | `off` or one-line `json` API access records on stdout |
| `DVBFIXER_METRICS` | `off` | `off` or protected Prometheus metrics via `/api/metrics` |
| `DVBFIXER_MUTATIONS_BACKUP_FILE` | `<gui>/mutations.json` | Mutable PostgreSQL backup; hardened deployments place it on bounded state |
| `DVBFIXER_OS_RESOURCE_LIMITS_REQUIRED` | `0` | Require the Linux cgroup/filesystem preflight when set to `1` |
| `DVBFIXER_OS_DATA_FILESYSTEM_MAX_BYTES` | unset | Maximum dedicated data-filesystem capacity; required by the OS preflight |
| `DVBFIXER_OS_TEMP_FILESYSTEM_MAX_BYTES` | unset | Maximum capacity of each private `/tmp` and `/var/tmp`; required by the OS preflight |
| `DVBFIXER_OS_CPU_QUOTA_PERCENT` | unset | Maximum accepted cgroup CPU quota percentage; required by the OS preflight |
| `DVBFIXER_OS_MEMORY_MAX_BYTES` | unset | Maximum accepted cgroup memory limit; required by the OS preflight |
| `DVBFIXER_OS_MEMORY_SWAP_MAX_BYTES` | unset | Maximum accepted cgroup swap limit; `0` disables swap; required by the OS preflight |
| `DVBFIXER_OS_TASKS_MAX` | unset | Maximum accepted cgroup task limit; required by the OS preflight |

Workspace accounting includes regular files in hidden, failed-run, and trash
directories and refuses symlinks. Request-size failures return `413`; a
workspace publication that exceeds its quota returns `507`. The storage limit
is an application-level, single-process steady-state check, not an OS filesystem
quota. A child process can temporarily create data before publication checks,
and multiple server instances do not coordinate limits.

The opt-in OS preflight verifies cgroup CPU, memory, swap, and task limits
against configured ceilings, a read-only root, a dedicated capacity-bounded
data filesystem, bounded private temporary filesystems, and containment of
mutable home/cache/backup paths. These limits
apply to the whole service and all descendants, not independently per workspace
or job. See [`deployment.md`](deployment.md).

The request limiter runs before authentication. CORS invokes the same limiter
before rejecting an Origin, so rejected Origins, invalid-token attempts, and
public endpoints consume the same direct-address allowance. Valid CORS
preflights terminate before admission; malformed and plain `OPTIONS` requests do
not bypass it. Exhausted clients receive
`429` with `Retry-After` and `RateLimit-*` headers. If the bounded address
tracker is full, an untracked client receives `503`; active entries are never
evicted to admit a new address. The server intentionally ignores forwarded-
address headers. Allowed cross-origin responses expose the rate-limit headers.
A reverse proxy therefore appears as one shared client and
should enforce its own edge limit; do not disable the backend limit without an
equivalent trusted control.

Child-process admission starts work immediately when a permit is available and
otherwise waits FIFO up to `DVBFIXER_MAX_QUEUED_PROCESSES`. Further requests are
rejected before spawning with an overload result; HTTP command and naming routes
map that condition to `503`. Command timeouts begin after admission, not while
waiting in the queue. Queue and concurrency limits are process-local rather than
distributed scheduling controls.

Naming-specific defaults are:

| Environment variable | Default | Purpose |
|---|---:|---|
| `DVBFIXER_NAMING_MAX_SOURCE_BYTES` | 50 MiB | Maximum source artifact size |
| `DVBFIXER_NAMING_MAX_REPORT_BYTES` | 20 MiB | Maximum CLI report size |
| `DVBFIXER_NAMING_TIMEOUT_MS` | 60,000 ms | Subprocess timeout |

These settings must be positive safe integers; the timeout may not exceed
2,147,483,647 ms, Node's timer maximum. Source size is checked before the
operation directory is created or a child process is admitted. The subprocess
receives `SIGTERM` on timeout and escalates to `SIGKILL` after a short grace
period; the API returns `504 NAMING_TIMEOUT`.

Operation directories and the shared `runs/_failed` directory use mode `0700`
on POSIX hosts. Server-owned logs, validated reports, and variant-override files
use mode `0600`. Successful publication removes those helpers, and a dry run
removes the whole operation directory. Other failures retain the bounded stdout,
stderr, report when available, overrides, and any partial output below
`runs/_failed`; publication-stage failures restore diagnostics removed during
the commit attempt. Failed runs are never registered as visible artifacts, are
counted against workspace quota, and currently remain until an operator removes
them. A quota-exhaustion failure deletes its operation data instead of retaining
more bytes.

## Browser Origins

Requests without an `Origin` header remain available to non-browser clients.
Browser requests are allowed from the request's exact same origin. Configure
additional origins as a bounded JSON array of exact canonical HTTP(S) origins:

```bash
export DVBFIXER_CORS_ALLOWED_ORIGINS='["https://structures.example.org"]'
```

When TLS terminates or `Host` is rewritten at a reverse proxy, list the public
HTTPS origin explicitly; the backend does not trust forwarded scheme or host
headers.

Wildcard origins, `null`, credentials mode, private-network preflights, malformed
or duplicate Origin headers, and unlisted cross-origin requests are rejected.
Preflights permit only the API methods and the `Authorization`, `Content-Type`,
`Accept`, and `X-File-Name` request headers.

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

Remote binding requires configured authentication and
`DVBFIXER_ALLOW_INSECURE_REMOTE=1` as an explicit acknowledgment. It remains
unsupported for untrusted networks until TLS, OS-level resource
enforcement is verified on the target host, audit retention, and distributed
scheduling/event controls are implemented.
