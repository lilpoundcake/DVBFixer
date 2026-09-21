# ADR 0004: Systemd Resource Boundary

## Status

Accepted.

## Decision

The supported OS resource-containment profile is one Linux systemd service on a
cgroup-v2 host. `CPUQuota`, `MemoryMax`, `MemorySwapMax`, and `TasksMax` limit
the aggregate Node server and every inherited scientific descendant. A dedicated
fixed-capacity filesystem contains workspace data and mutable home/cache/backup
state; private size-limited tmpfs mounts contain `/tmp` and `/var/tmp`.

Standalone startup fails closed when
`DVBFIXER_OS_RESOURCE_LIMITS_REQUIRED=1` unless it observes all of those bounds,
their configured ceilings, and a read-only root filesystem.
Persistent writable mounts other than the dedicated data filesystem are
rejected; kernel and memory-backed filesystems remain bounded by permissions and
the service memory cgroup.
Local Vite and loopback development leave the check disabled.

## Rationale

Per-process `rlimit` controls do not bound an entire process tree: DVBFixer can
launch AmberTools, Reduce, Modeller, MSA engines, Open Babel, xTB, and other
children. `RLIMIT_FSIZE` limits one file rather than aggregate storage. A service
cgroup and bounded filesystems cover descendants without reproducing scientific
launch logic in Node.

The existing Docker Compose file remains development-only PostgreSQL. Container
CPU and memory limits would not portably limit bind-mounted or named-volume
storage, and a production scientific image would be a separate deployment
project.

## Consequences

- Limits are aggregate per service, not per tenant, workspace, or job.
- GPU memory and externally managed PostgreSQL remain outside this boundary.
- A memory-limit event may terminate the API along with a child; systemd restarts
  the service and persisted managed-job recovery records the interruption.
- Operators must size limits for their workloads and run privileged enforcement
  tests on the target host before claiming the deployment is supported.
