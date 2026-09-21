import { afterEach, describe, expect, it } from 'vitest'
import fs from 'node:fs'
import os from 'node:os'
import path from 'node:path'
import { createAuditSink, parseAuditLogConfig, pruneAuditLogs, type AuditEvent } from './audit-log'

const roots: string[] = []
function root(): string {
  const value = fs.mkdtempSync(path.join(os.tmpdir(), 'dvbfixer-audit-'))
  roots.push(value)
  return value
}
const event: AuditEvent = {
  schema: 'dvbfixer.audit.v1', timestamp: '2026-09-21T12:00:00.000Z',
  requestId: 'req_test', principalId: 'operator', method: 'POST',
  route: '/api/v1/workspaces/:workspaceId/jobs', status: 202, outcome: 'completed',
}

afterEach(() => roots.splice(0).forEach(value => fs.rmSync(value, { recursive: true, force: true })))

describe('durable audit log', () => {
  it('appends fs-backed JSONL with private permissions', () => {
    const dataRoot = root()
    const sink = createAuditSink(dataRoot, { enabled: true, retentionDays: 90 }, () => new Date(event.timestamp))!
    sink(event)
    sink({ ...event, requestId: 'req_second', status: 422 })
    const file = path.join(dataRoot, '_audit', 'audit-2026-09-21.jsonl')
    expect(fs.readFileSync(file, 'utf8').trim().split('\n').map(line => JSON.parse(line).requestId))
      .toEqual(['req_test', 'req_second'])
    if (process.platform !== 'win32') {
      expect(fs.statSync(path.dirname(file)).mode & 0o777).toBe(0o700)
      expect(fs.statSync(file).mode & 0o777).toBe(0o600)
    }
  })

  it('prunes only expired audit files', () => {
    const dataRoot = root()
    const directory = path.join(dataRoot, '_audit')
    fs.mkdirSync(directory)
    const old = path.join(directory, 'audit-2020-01-01.jsonl')
    const recent = path.join(directory, 'audit-2026-09-21.jsonl')
    const unrelated = path.join(directory, 'notes.txt')
    for (const file of [old, recent, unrelated]) fs.writeFileSync(file, 'x')
    fs.utimesSync(old, new Date(0), new Date(0))
    pruneAuditLogs(directory, 90, Date.parse('2026-09-21T12:00:00Z'))
    expect(fs.existsSync(old)).toBe(false)
    expect(fs.existsSync(recent)).toBe(true)
    expect(fs.existsSync(unrelated)).toBe(true)
  })

  it('validates mode and retention bounds', () => {
    expect(parseAuditLogConfig({ DVBFIXER_AUDIT_LOG: 'jsonl', DVBFIXER_AUDIT_RETENTION_DAYS: '30' }))
      .toEqual({ enabled: true, retentionDays: 30 })
    expect(() => parseAuditLogConfig({ DVBFIXER_AUDIT_LOG: 'text' })).toThrow('DVBFIXER_AUDIT_LOG')
    expect(() => parseAuditLogConfig({ DVBFIXER_AUDIT_RETENTION_DAYS: '0' })).toThrow('positive integer')
    expect(() => parseAuditLogConfig({ DVBFIXER_AUDIT_RETENTION_DAYS: '3651' })).toThrow('at most 3650')
  })

  it('rejects a symlinked audit root', () => {
    const dataRoot = root()
    const outside = root()
    fs.symlinkSync(outside, path.join(dataRoot, '_audit'), 'dir')
    expect(() => createAuditSink(dataRoot, { enabled: true, retentionDays: 90 })).toThrow('real directory')
  })
})
