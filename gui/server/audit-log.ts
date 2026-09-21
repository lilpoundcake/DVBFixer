import fs from 'node:fs'
import path from 'node:path'

export interface AuditEvent {
  schema: 'dvbfixer.audit.v1'
  timestamp: string
  requestId: string
  principalId: string | null
  method: string
  route: string
  status: number | null
  outcome: 'completed' | 'aborted'
}

export interface AuditLogConfig {
  enabled: boolean
  retentionDays: number
}

export function parseAuditLogConfig(environment: NodeJS.ProcessEnv = process.env): AuditLogConfig {
  const mode = environment.DVBFIXER_AUDIT_LOG ?? 'off'
  if (mode !== 'off' && mode !== 'jsonl') throw new Error('DVBFIXER_AUDIT_LOG must be off or jsonl')
  const rawDays = environment.DVBFIXER_AUDIT_RETENTION_DAYS ?? '90'
  if (!/^[1-9]\d*$/.test(rawDays)) throw new Error('DVBFIXER_AUDIT_RETENTION_DAYS must be a positive integer')
  const retentionDays = Number(rawDays)
  if (!Number.isSafeInteger(retentionDays) || retentionDays > 3650) {
    throw new Error('DVBFIXER_AUDIT_RETENTION_DAYS must be at most 3650')
  }
  return { enabled: mode === 'jsonl', retentionDays }
}

function auditRoot(dataRoot: string): string {
  const root = path.join(path.resolve(dataRoot), '_audit')
  fs.mkdirSync(root, { recursive: true, mode: 0o700 })
  const stat = fs.lstatSync(root)
  if (!stat.isDirectory() || stat.isSymbolicLink()) throw new Error('audit root must be a real directory')
  fs.chmodSync(root, 0o700)
  return root
}

export function pruneAuditLogs(root: string, retentionDays: number, now = Date.now()): void {
  const cutoff = now - retentionDays * 24 * 60 * 60 * 1000
  for (const entry of fs.readdirSync(root, { withFileTypes: true })) {
    if (!entry.isFile() || !/^audit-\d{4}-\d{2}-\d{2}\.jsonl$/.test(entry.name)) continue
    const file = path.join(root, entry.name)
    const stat = fs.lstatSync(file)
    if (stat.isFile() && !stat.isSymbolicLink() && stat.mtimeMs <= cutoff) fs.rmSync(file)
  }
}

export function createAuditSink(
  dataRoot: string,
  config: AuditLogConfig,
  now: () => Date = () => new Date(),
): ((event: AuditEvent) => void) | undefined {
  if (!config.enabled) return undefined
  const root = auditRoot(dataRoot)
  let lastPrunedDate = ''
  return event => {
    const current = now()
    const date = current.toISOString().slice(0, 10)
    if (date !== lastPrunedDate) {
      pruneAuditLogs(root, config.retentionDays, current.getTime())
      lastPrunedDate = date
    }
    const file = path.join(root, `audit-${date}.jsonl`)
    const descriptor = fs.openSync(file, fs.constants.O_WRONLY | fs.constants.O_CREAT | fs.constants.O_APPEND, 0o600)
    try {
      fs.writeSync(descriptor, `${JSON.stringify(event)}\n`, undefined, 'utf8')
      fs.fsyncSync(descriptor)
    } finally {
      fs.closeSync(descriptor)
    }
    fs.chmodSync(file, 0o600)
  }
}
