import crypto from 'node:crypto'
import fs from 'node:fs'
import os from 'node:os'
import path from 'node:path'
import { spawnSync } from 'node:child_process'

interface LockOwner {
  version: 1
  token: string
  pid: number
  hostname: string
  bootId: string | null
  processStartTicks: string | null
  acquiredAt: string
}

export class WorkspaceManifestLockError extends Error {
  readonly statusCode = 503
  readonly code = 'WORKSPACE_MANIFEST_LOCKED'
  readonly workspaceId: string

  constructor(workspaceId: string) {
    super(`workspace manifest is locked: ${workspaceId}`)
    this.name = 'WorkspaceManifestLockError'
    this.workspaceId = workspaceId
  }
}

export class WorkspaceRunLockError extends Error {
  readonly statusCode = 409
  readonly code = 'WORKSPACE_BUSY'
  readonly workspaceId: string

  constructor(workspaceId: string) {
    super(`another DVBFixer operation is already running in workspace: ${workspaceId}`)
    this.name = 'WorkspaceRunLockError'
    this.workspaceId = workspaceId
  }
}

function safeId(id: string): string {
  if (!/^[a-zA-Z0-9_-]+$/.test(id)) throw new Error('invalid workspace id')
  return id
}

function linuxBootId(): string | null {
  try { return fs.readFileSync('/proc/sys/kernel/random/boot_id', 'utf8').trim() || null } catch { return null }
}

function linuxProcessIdentity(pid: number): { state: string; startTicks: string } | null {
  try {
    const text = fs.readFileSync(`/proc/${pid}/stat`, 'utf8')
    const end = text.lastIndexOf(') ')
    if (end < 0) return null
    const fields = text.slice(end + 2).trim().split(/\s+/)
    if (fields.length < 20) return null
    return { state: fields[0], startTicks: fields[19] }
  } catch (error) {
    if ((error as NodeJS.ErrnoException).code === 'ENOENT') return null
    throw error
  }
}

const currentBootId = linuxBootId()

function processStartIdentity(pid: number): string | null {
  if (process.platform === 'linux') return linuxProcessIdentity(pid)?.startTicks ?? null
  const command = process.platform === 'win32'
    ? {
        executable: 'powershell.exe',
        args: ['-NoProfile', '-NonInteractive', '-Command',
          '$p=Get-Process -Id $args[0] -ErrorAction Stop; $p.StartTime.ToUniversalTime().Ticks', String(pid)],
      }
    : { executable: '/bin/ps', args: ['-o', 'lstart=', '-p', String(pid)] }
  const result = spawnSync(command.executable, command.args, {
    encoding: 'utf8', timeout: 1_000, windowsHide: true,
  })
  return result.status === 0 && result.stdout.trim() ? result.stdout.trim() : null
}

const currentStartTicks = processStartIdentity(process.pid)

function ownerRecord(): LockOwner {
  return {
    version: 1,
    token: crypto.randomUUID(),
    pid: process.pid,
    hostname: os.hostname(),
    bootId: currentBootId,
    processStartTicks: currentStartTicks,
    acquiredAt: new Date().toISOString(),
  }
}

function validOwner(value: unknown): value is LockOwner {
  if (!value || typeof value !== 'object') return false
  const owner = value as Record<string, unknown>
  return owner.version === 1 && typeof owner.token === 'string' &&
    Number.isSafeInteger(owner.pid) && Number(owner.pid) > 0 &&
    typeof owner.hostname === 'string' &&
    (owner.bootId === null || typeof owner.bootId === 'string') &&
    (owner.processStartTicks === null || typeof owner.processStartTicks === 'string') &&
    typeof owner.acquiredAt === 'string'
}

function readOwner(directory: string): LockOwner | null {
  try {
    const parsed: unknown = JSON.parse(fs.readFileSync(path.join(directory, 'owner.json'), 'utf8'))
    return validOwner(parsed) ? parsed : null
  } catch {
    return null
  }
}

function ownerDefinitelyDead(owner: LockOwner): boolean {
  if (owner.hostname !== os.hostname()) return false
  if (owner.bootId && currentBootId && owner.bootId !== currentBootId) return true
  if (process.platform === 'linux') {
    const identity = linuxProcessIdentity(owner.pid)
    if (!identity || identity.state === 'Z') return true
    if (owner.processStartTicks && identity.startTicks !== owner.processStartTicks) return true
    return false
  }
  const identity = processStartIdentity(owner.pid)
  if (owner.processStartTicks && identity) return identity !== owner.processStartTicks
  try {
    process.kill(owner.pid, 0)
    return false
  } catch (error) {
    return (error as NodeJS.ErrnoException).code === 'ESRCH'
  }
}

function writeOwner(directory: string, owner: LockOwner): void {
  fs.mkdirSync(directory, { mode: 0o700 })
  const file = path.join(directory, 'owner.json')
  const descriptor = fs.openSync(file, 'wx', 0o600)
  try {
    fs.writeFileSync(descriptor, `${JSON.stringify(owner)}\n`, 'utf8')
    fs.fsyncSync(descriptor)
  } finally {
    fs.closeSync(descriptor)
  }
  if (process.platform !== 'win32') {
    const directoryDescriptor = fs.openSync(directory, 'r')
    try { fs.fsyncSync(directoryDescriptor) } finally { fs.closeSync(directoryDescriptor) }
  }
}

function lockRoot(dataRoot: string): string {
  const root = path.join(path.resolve(dataRoot), '.workspace-locks')
  fs.mkdirSync(root, { recursive: true, mode: 0o700 })
  const stat = fs.lstatSync(root)
  if (!stat.isDirectory() || stat.isSymbolicLink()) {
    throw new Error('workspace manifest lock root must be a real directory')
  }
  return root
}

function acquireReapMarker(directory: string, observedToken: string, workspaceId: string): void {
  const marker = path.join(directory, `.reap-${observedToken}`)
  for (let attempt = 0; attempt < 3; attempt += 1) {
    const reaper = ownerRecord()
    const claim = path.join(directory, `.reap-${observedToken}.${reaper.token}.claim`)
    try {
      writeOwner(claim, reaper)
      fs.renameSync(claim, marker)
      return
    } catch (error) {
      fs.rmSync(claim, { recursive: true, force: true })
      const code = (error as NodeJS.ErrnoException).code
      if (code === 'ENOENT') throw new WorkspaceManifestLockError(workspaceId)
      if (code !== 'EEXIST' && code !== 'ENOTEMPTY') throw error
    }
    const markerOwner = readOwner(marker)
    if (!markerOwner || !ownerDefinitelyDead(markerOwner)) {
      throw new WorkspaceManifestLockError(workspaceId)
    }
    const staleMarker = `${marker}.${markerOwner.token}.stale`
    try {
      fs.renameSync(marker, staleMarker)
      fs.rmSync(staleMarker, { recursive: true, force: true })
    } catch (error) {
      if ((error as NodeJS.ErrnoException).code !== 'ENOENT') throw error
    }
  }
  throw new WorkspaceManifestLockError(workspaceId)
}

function acquire(
  dataRoot: string,
  workspaceId: string,
  lockName: string,
): { directory: string; owner: LockOwner } {
  const root = lockRoot(dataRoot)
  const directory = path.join(root, `${safeId(lockName)}.lock`)
  for (let attempt = 0; attempt < 3; attempt += 1) {
    const owner = ownerRecord()
    const claim = path.join(root, `.${workspaceId}.${owner.token}.claim`)
    try {
      writeOwner(claim, owner)
      fs.renameSync(claim, directory)
      return { directory, owner }
    } catch (error) {
      fs.rmSync(claim, { recursive: true, force: true })
      const code = (error as NodeJS.ErrnoException).code
      if (code !== 'EEXIST' && code !== 'ENOTEMPTY' && code !== 'ENOENT') throw error
    }

    const observed = readOwner(directory)
    if (!observed || !ownerDefinitelyDead(observed)) {
      throw new WorkspaceManifestLockError(workspaceId)
    }
    acquireReapMarker(directory, observed.token, workspaceId)
    const confirmed = readOwner(directory)
    if (confirmed?.token !== observed.token || !ownerDefinitelyDead(confirmed)) {
      fs.rmSync(path.join(directory, `.reap-${observed.token}`), { recursive: true, force: true })
      throw new WorkspaceManifestLockError(workspaceId)
    }
    const quarantine = path.join(root, `.${workspaceId}.${observed.token}.stale`)
    try {
      fs.renameSync(directory, quarantine)
    } catch (error) {
      if ((error as NodeJS.ErrnoException).code === 'ENOENT') continue
      throw error
    }
    fs.rmSync(quarantine, { recursive: true, force: true })
  }
  throw new WorkspaceManifestLockError(workspaceId)
}

function release(directory: string, owner: LockOwner): void {
  const marker = path.join(directory, `.release-${owner.token}`)
  fs.mkdirSync(marker)
  const current = readOwner(directory)
  if (current?.token !== owner.token) {
    fs.rmSync(marker, { recursive: true, force: true })
    throw new Error('workspace manifest lock ownership changed')
  }
  const released = `${directory}.${owner.token}.released`
  fs.renameSync(directory, released)
  fs.rmSync(released, { recursive: true })
}

export function withWorkspaceManifestLock<T>(
  dataRoot: string,
  workspaceId: string,
  callback: () => T,
): T {
  const lock = acquire(dataRoot, workspaceId, `workspace-${workspaceId}`)
  try {
    return callback()
  } finally {
    release(lock.directory, lock.owner)
  }
}

export function withGlobalMigrationLock<T>(dataRoot: string, callback: () => T): T {
  const lock = acquire(dataRoot, 'global migrations', 'system-migrations')
  try {
    return callback()
  } finally {
    release(lock.directory, lock.owner)
  }
}

/**
 * Acquire the durable per-workspace scientific-run lease.
 *
 * The lock record includes host, boot, PID, and process-start identity, so a
 * dead same-host owner can be reaped without stealing a live or unverifiable
 * lock. The returned closure releases only the caller's token.
 */
export function acquireWorkspaceRunLock(dataRoot: string, workspaceId: string): () => void {
  let lock: ReturnType<typeof acquire>
  try {
    lock = acquire(dataRoot, workspaceId, `run-${workspaceId}`)
  } catch (error) {
    if (error instanceof WorkspaceManifestLockError) throw new WorkspaceRunLockError(workspaceId)
    throw error
  }
  return () => release(lock.directory, lock.owner)
}
