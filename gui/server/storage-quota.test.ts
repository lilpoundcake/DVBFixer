import { afterEach, describe, expect, it } from 'vitest'
import fs from 'node:fs'
import os from 'node:os'
import path from 'node:path'
import {
  assertWorkspaceQuota, availableWorkspaceBytes, DEFAULT_MAX_UPLOAD_BYTES,
  DEFAULT_WORKSPACE_QUOTA_BYTES, initializeStorageQuota, parseStorageQuotaSettings,
  workspaceLogicalBytes, WorkspaceQuotaExceededError,
} from './storage-quota'

const temporaryDirectories: string[] = []

function temp(): string {
  const directory = fs.mkdtempSync(path.join(os.tmpdir(), 'dvbfixer-quota-test-'))
  temporaryDirectories.push(directory)
  return directory
}

afterEach(() => {
  initializeStorageQuota({
    maxUploadBytes: DEFAULT_MAX_UPLOAD_BYTES,
    workspaceQuotaBytes: DEFAULT_WORKSPACE_QUOTA_BYTES,
  })
  temporaryDirectories.splice(0).forEach(directory => fs.rmSync(directory, { recursive: true, force: true }))
})

describe('storage quota settings', () => {
  it('uses safe defaults and parses explicit limits', () => {
    expect(parseStorageQuotaSettings({})).toEqual({
      maxUploadBytes: 256 * 1024 * 1024,
      workspaceQuotaBytes: 5 * 1024 * 1024 * 1024,
    })
    expect(parseStorageQuotaSettings({
      DVBFIXER_GUI_MAX_UPLOAD_BYTES: '123',
      DVBFIXER_GUI_WORKSPACE_QUOTA_BYTES: '0',
    })).toEqual({ maxUploadBytes: 123, workspaceQuotaBytes: 0 })
  })

  it('rejects unsafe, fractional, negative, and zero upload limits', () => {
    for (const value of ['0', '-1', '1.5', String(Number.MAX_SAFE_INTEGER + 1), 'nope']) {
      expect(() => parseStorageQuotaSettings({ DVBFIXER_GUI_MAX_UPLOAD_BYTES: value })).toThrow(/positive safe integer/)
    }
    expect(() => parseStorageQuotaSettings({ DVBFIXER_GUI_WORKSPACE_QUOTA_BYTES: '-1' })).toThrow(/non-negative safe integer/)
  })
})

describe('workspace logical byte accounting', () => {
  it('counts hidden, trash, and failed-run regular files recursively', () => {
    const root = temp()
    fs.mkdirSync(path.join(root, '.trash'), { recursive: true })
    fs.mkdirSync(path.join(root, 'runs', 'failed'), { recursive: true })
    fs.writeFileSync(path.join(root, '.hidden'), '12')
    fs.writeFileSync(path.join(root, '.trash', 'old.pdb'), '345')
    fs.writeFileSync(path.join(root, 'runs', 'failed', 'log.txt'), '6789')
    expect(workspaceLogicalBytes(root)).toBe(9)
  })

  it('does not follow symlinks and fails when one is present', () => {
    const root = temp()
    const outside = temp()
    fs.writeFileSync(path.join(outside, 'large.bin'), 'outside')
    fs.symlinkSync(path.join(outside, 'large.bin'), path.join(root, 'link'))
    expect(() => workspaceLogicalBytes(root)).toThrow(/must not contain symlinks/)
  })

  it('reports available bytes and typed 507 details', () => {
    const root = temp()
    fs.writeFileSync(path.join(root, 'data'), '1234')
    initializeStorageQuota({ maxUploadBytes: 10, workspaceQuotaBytes: 6 })
    expect(availableWorkspaceBytes(root)).toBe(2)
    expect(() => assertWorkspaceQuota(root, 3)).toThrow(WorkspaceQuotaExceededError)
    try { assertWorkspaceQuota(root, 3) } catch (error) {
      expect(error).toMatchObject({
        statusCode: 507,
        code: 'WORKSPACE_QUOTA_EXCEEDED',
        details: { quotaBytes: 6, usedBytes: 4, additionalBytes: 3, projectedBytes: 7 },
      })
    }
    initializeStorageQuota({ maxUploadBytes: 10, workspaceQuotaBytes: 0 })
    expect(availableWorkspaceBytes(root)).toBeNull()
    expect(() => assertWorkspaceQuota(root, Number.MAX_SAFE_INTEGER)).not.toThrow()
  })
})
