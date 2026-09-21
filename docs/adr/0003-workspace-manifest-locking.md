# ADR 0003: Workspace Manifest Locking

## Status

Accepted.

## Decision

Serialize every workspace manifest mutation with an atomic directory lock below
`<dataRoot>/.workspace-locks`. Lock ownership records include a random token,
hostname, boot identity, PID, and process-start identity. A lock is reclaimed
only when its same-host owner is definitely dead; malformed or unverifiable
ownership fails closed. Release and stale-owner recovery rename the owned lock
directory before removal so they cannot delete a replacement owner's lock.

`saveWorkspace` compares its revision with the manifest read while locked.
Server-generated additive publications use `updateWorkspace` to mutate the
locked latest manifest. Workspace and artifact deletion use the same lock,
recoverable trash moves, rollback where possible, and durable deletion markers
so a stale revision-zero save cannot recreate a deleted workspace. One-time
migrations use a separate global lock namespace.

## Scope

This supports multiple Node processes on one host sharing a local filesystem.
It is not a distributed lease and does not claim safety across hosts, NFS-like
filesystems, or storage without atomic rename and reliable process identity.
Job scheduling, active-run locks, and SSE subscribers remain process-local.

## Consequences

- Concurrent manifest updates no longer silently overwrite each other.
- Client stale revisions return `409`; temporary lock contention returns `503`.
- Crashed same-host owners and deletion transactions are recoverable.
- Horizontal deployment still requires durable scheduling, event delivery, and
  distributed coordination.
