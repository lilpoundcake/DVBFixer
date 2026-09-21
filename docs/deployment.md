# Resource-Bounded Deployment

DVBFixer's supported resource-containment profile is a single Linux systemd
service on a cgroup-v2 host. It limits the aggregate Node server, DVBFixer
processes, and their scientific descendants. It is not per-job or per-workspace
OS isolation.

This profile does not by itself make the API safe for the public internet. A
production operator must still provide TLS, audit retention, secret management,
database limits, backups, monitoring, and durable scheduling before horizontal
scaling.

## Storage Prerequisite

Mount a dedicated finite-capacity filesystem at `/var/lib/dvbfixer`. A normal
directory on the host root filesystem does not satisfy this requirement. Use a
fixed-size logical volume and filesystem. The example environment caps accepted
filesystem capacity at 64 GiB; the startup
preflight fails if the actual filesystem is larger.

The mount contains workspaces, the mutation backup, `HOME`, `TMPDIR`, and cache
state. The service also mounts private 8 GiB tmpfs filesystems on `/tmp` and
`/var/tmp`, preventing tools that ignore `TMPDIR` from filling the host root
filesystem. PostgreSQL is external to this boundary and needs independent CPU,
memory, task, and storage limits.

## Installation

1. Build and install the repository at `/opt/dvbfixer`, including the GUI's
   `dist/` and `dist-server/` outputs and the scientific environment.
2. Create an unprivileged `dvbfixer` user and group.
3. Mount the quota-backed filesystem at `/var/lib/dvbfixer`, owned by that user.
4. Create `/var/lib/dvbfixer/{tmp,home,cache}` with mode `0700`.
5. Copy `deploy/systemd/dvbfixer.env.example` to
   `/etc/dvbfixer/dvbfixer.env`, mode `0600`, and configure authentication,
   database, CORS, executable paths, and licensed software.
6. Install `deploy/systemd/dvbfixer.service` under `/etc/systemd/system/`.
7. Adjust CPU, memory, task, and temporary-storage limits with a systemd
   drop-in based on measured scientific workloads.
8. Run `systemd-analyze verify dvbfixer.service`, then enable and start it.

The supplied service defaults to `CPUQuota=200%`, `MemoryHigh=12G`,
`MemoryMax=16G`, `MemorySwapMax=0`, `TasksMax=512`, and two private 8 GiB
temporary filesystems. These are deployment examples, not universal scientific
sizing recommendations. GPU memory is not controlled by the memory cgroup.

## Fail-Closed Preflight

`DVBFIXER_OS_RESOURCE_LIMITS_REQUIRED=1` makes standalone startup verify:

- Linux and cgroup v2 are active;
- `CPUQuota`, `MemoryMax`, `MemorySwapMax`, and `TasksMax` do not exceed the
  configured ceilings for the current cgroup;
- the service sees its root filesystem as read-only;
- no additional persistent filesystem is writable inside the service namespace;
- `DVBFIXER_GUI_DATA_DIR` is the root of a dedicated filesystem whose capacity
  does not exceed `DVBFIXER_OS_DATA_FILESYSTEM_MAX_BYTES`;
- `HOME`, `TMPDIR`, `XDG_CACHE_HOME`, and the mutation backup are inside that
  filesystem; and
- `/tmp` and `/var/tmp` are dedicated tmpfs mounts no larger than
  `DVBFIXER_OS_TEMP_FILESYSTEM_MAX_BYTES`.

The flag is disabled by default so local macOS, Windows, Vite, and loopback
development remain unchanged. Do not disable it in the resource-bounded systemd
profile.

## Access Log Retention

The systemd example enables `DVBFIXER_ACCESS_LOG=json`; records go to stdout and
therefore to the service journal. Configure journald forwarding, access control,
rotation, retention, and deletion for the deployment. These redacted request
records are operational access logs, not durable audit events. Do not grant log
readers broader access than API operators merely because bodies and credentials
are omitted.

## Metrics

The systemd example enables `DVBFIXER_METRICS=prometheus`. Scrape
`GET /api/metrics` with a configured bearer credential in the `Authorization`
header; never place credentials in the URL. The endpoint is rate-limited and
available to every authenticated API principal because the current static-token
model has no operator role.

Metrics reset on restart. Scrape each server process independently and aggregate
externally. Fixed labels omit principals, workspace identity, filenames, client
addresses, and other scientific or tenant data.

## Acceptance

Before treating a deployment as supported, test on the target host that CPU is
throttled, combined child memory triggers the cgroup limit, `TasksMax` blocks a
fork-heavy descendant, persistent and temporary writes stop at their storage
boundaries, and service shutdown removes session-detached grandchildren. Include
at least one external MSA tool and one AmberTools/Reduce workflow in this smoke
test. The repository's unprivileged tests validate configuration and fail-closed
preflight logic but cannot prove host-kernel enforcement.
