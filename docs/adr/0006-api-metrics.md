# ADR 0006: Process-Local API Metrics

## Status

Accepted.

## Decision

Provide an opt-in protected Prometheus endpoint at `GET /api/metrics`, enabled by
`DVBFIXER_METRICS=prometheus`. It reports bounded HTTP request totals, route-only
duration histograms, active requests, and child-process admission gauges. The
scrape itself is excluded from HTTP accounting.

Labels come only from fixed method, coarse route, status, and completion-outcome
taxonomies. Principal, request/client identity, forwarded headers, workspace,
artifact, job, filename, raw path, query, command arguments, errors, and process
output are prohibited.

## Consequences

- Metrics are process-local and reset on restart; fleet aggregation belongs in
  the monitoring system.
- The endpoint remains behind CORS, rate limiting, and static bearer
  authentication. Every authenticated principal can read it because operator
  roles do not yet exist.
- Scrapes have no database side effects and perform no recursive workspace scan.
- Metrics are operational telemetry, not durable audit events.
