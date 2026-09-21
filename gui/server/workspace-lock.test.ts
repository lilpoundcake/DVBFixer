import fs from 'node:fs'
import os from 'node:os'
import path from 'node:path'
import { spawn } from 'node:child_process'
import { afterEach, describe, expect, it } from 'vitest'
import {
  acquireWorkspaceRunLock, WorkspaceManifestLockError, WorkspaceRunLockError,
  withGlobalMigrationLock, withWorkspaceManifestLock,
} from './workspace-lock'

const directories: string[] = []

function temp(): string {
  const directory = fs.mkdtempSync(path.join(os.tmpdir(), 'dvbfixer-workspace-lock-'))
  directories.push(directory)
  return directory
}

afterEach(() => {
  for (const directory of directories.splice(0)) fs.rmSync(directory, { recursive: true, force: true })
})

describe('workspace manifest lock', () => {
  it('serializes scientific runs independently of manifest writes', () => {
    const root = temp()
    const release = acquireWorkspaceRunLock(root, 'workspace-a')
    try {
      expect(() => acquireWorkspaceRunLock(root, 'workspace-a')).toThrow(WorkspaceRunLockError)
      expect(withWorkspaceManifestLock(root, 'workspace-a', () => 'manifest')).toBe('manifest')
    } finally {
      release()
    }
    const releaseAgain = acquireWorkspaceRunLock(root, 'workspace-a')
    releaseAgain()
  })

  it('excludes a live holder and releases for the next transaction', () => {
    const root = temp()
    withWorkspaceManifestLock(root, 'workspace-a', () => {
      expect(() => withWorkspaceManifestLock(root, 'workspace-a', () => {}))
        .toThrow(WorkspaceManifestLockError)
    })
    expect(withWorkspaceManifestLock(root, 'workspace-a', () => 'released')).toBe('released')
  })

  it('keeps the global migration lock separate from workspace IDs', () => {
    const root = temp()
    withGlobalMigrationLock(root, () => {
      expect(withWorkspaceManifestLock(root, 'system-migrations', () => 'workspace')).toBe('workspace')
    })
  })

  it('recovers a lock whose same-host owner is definitely dead', () => {
    const root = temp()
    const lock = path.join(root, '.workspace-locks', 'workspace-workspace-a.lock')
    fs.mkdirSync(lock, { recursive: true })
    fs.writeFileSync(path.join(lock, 'owner.json'), JSON.stringify({
      version: 1,
      token: 'dead-owner',
      pid: 2_147_483_647,
      hostname: os.hostname(),
      bootId: null,
      processStartTicks: null,
      acquiredAt: new Date(0).toISOString(),
    }))

    expect(withWorkspaceManifestLock(root, 'workspace-a', () => 'recovered')).toBe('recovered')
    expect(fs.existsSync(lock)).toBe(false)
  })

  it('fails closed for malformed owner metadata', () => {
    const root = temp()
    const lock = path.join(root, '.workspace-locks', 'workspace-workspace-a.lock')
    fs.mkdirSync(lock, { recursive: true })
    fs.writeFileSync(path.join(lock, 'owner.json'), '{broken')
    expect(() => withWorkspaceManifestLock(root, 'workspace-a', () => {}))
      .toThrow(WorkspaceManifestLockError)
    expect(fs.existsSync(lock)).toBe(true)
  })

  it('rejects a symlinked lock root', () => {
    const root = temp()
    const outside = temp()
    fs.symlinkSync(outside, path.join(root, '.workspace-locks'))
    expect(() => withWorkspaceManifestLock(root, 'workspace-a', () => {}))
      .toThrow(/must be a real directory/)
  })

  it('does not steal a live cross-process lock and recovers it after owner death', async () => {
    const root = temp()
    const script = String.raw`
      const fs = require('node:fs');
      const os = require('node:os');
      const path = require('node:path');
      const root = process.argv[1];
      const lock = path.join(root, '.workspace-locks', 'workspace-workspace-a.lock');
      fs.mkdirSync(lock, { recursive: true });
      let bootId = null;
      let processStartTicks = null;
      try { bootId = fs.readFileSync('/proc/sys/kernel/random/boot_id', 'utf8').trim(); } catch {}
      try {
        const stat = fs.readFileSync('/proc/' + process.pid + '/stat', 'utf8');
        const fields = stat.slice(stat.lastIndexOf(') ') + 2).trim().split(/\s+/);
        processStartTicks = fields[19];
      } catch {}
      const owner = {
        version: 1, token: 'child-owner', pid: process.pid, hostname: os.hostname(),
        bootId, processStartTicks, acquiredAt: new Date().toISOString(),
      };
      fs.writeFileSync(path.join(lock, 'owner.json'), JSON.stringify(owner));
      process.stdout.write('ready\n');
      setInterval(() => {}, 1000);
    `
    const child = spawn(process.execPath, ['-e', script, root], { stdio: ['ignore', 'pipe', 'inherit'] })
    await new Promise<void>((resolve, reject) => {
      child.once('error', reject)
      child.stdout.once('data', () => resolve())
    })
    expect(() => withWorkspaceManifestLock(root, 'workspace-a', () => {}))
      .toThrow(WorkspaceManifestLockError)

    child.kill('SIGKILL')
    await new Promise<void>(resolve => child.once('exit', () => resolve()))
    expect(withWorkspaceManifestLock(root, 'workspace-a', () => 'recovered')).toBe('recovered')
  })

  it('recovers a reaper marker whose owner died before quarantining the stale lock', () => {
    const root = temp()
    const lock = path.join(root, '.workspace-locks', 'workspace-workspace-a.lock')
    const deadOwner = {
      version: 1,
      token: 'dead-owner',
      pid: 2_147_483_647,
      hostname: os.hostname(),
      bootId: null,
      processStartTicks: null,
      acquiredAt: new Date(0).toISOString(),
    }
    fs.mkdirSync(lock, { recursive: true })
    fs.writeFileSync(path.join(lock, 'owner.json'), JSON.stringify(deadOwner))
    const marker = path.join(lock, '.reap-dead-owner')
    fs.mkdirSync(marker)
    fs.writeFileSync(path.join(marker, 'owner.json'), JSON.stringify({ ...deadOwner, token: 'dead-reaper' }))

    expect(withWorkspaceManifestLock(root, 'workspace-a', () => 'recovered')).toBe('recovered')
    expect(fs.existsSync(lock)).toBe(false)
  })
})
