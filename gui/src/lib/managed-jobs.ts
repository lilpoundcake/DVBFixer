export type ManagedJobStatus = 'queued' | 'running' | 'succeeded' | 'failed' | 'cancelled'

/** Human-readable local time; retain original ISO timestamps in persisted data. */
export function formatJobTimestamp(value: string | null, locale?: string, timeZone?: string): string {
  if (!value) return 'Not started'
  const date = new Date(value)
  if (!Number.isFinite(date.getTime())) return 'Unknown time'
  return new Intl.DateTimeFormat(locale, {
    year: 'numeric', month: '2-digit', day: '2-digit',
    hour: '2-digit', minute: '2-digit', second: '2-digit', hour12: false,
    ...(timeZone ? { timeZone } : {}),
  }).format(date)
}

export interface ManagedJobRecord {
  version: 1
  id: string
  workspaceId: string
  status: ManagedJobStatus
  createdAt: string
  startedAt: string | null
  finishedAt: string | null
  command: string
  args: string[]
  exitCode: number | null
  stdoutLog: string
  stderrLog: string
  outputDir: string
  outputFile: string | null
  provenance?: { outputs: Array<{ artifactId: string; file: string; sha256: string }> }
  error?: string
}

export interface ManagedJobRequestPayload {
  workspaceId: string
  command: string
  inputFile: string
  inputs: Record<string, string | string[]>
  values: Record<string, unknown>
  fastaContent: string
}

/** Build the managed-job wire payload with the selected command included. */
export function createManagedJobRequest(
  command: string,
  request: Omit<ManagedJobRequestPayload, 'command'>,
): ManagedJobRequestPayload {
  return { ...request, command }
}

export function isActiveManagedJob(job: ManagedJobRecord | null | undefined): boolean {
  return job?.status === 'queued' || job?.status === 'running'
}

/** Prefer the active run when restoring a panel; otherwise show the newest recorded run. */
export function selectRestoredManagedJob(jobs: readonly ManagedJobRecord[]): ManagedJobRecord | null {
  return jobs.find(isActiveManagedJob) || jobs[0] || null
}

export function mergeManagedJobHistory(jobs: readonly ManagedJobRecord[], latest: ManagedJobRecord | null): ManagedJobRecord[] {
  const merged = latest ? [latest, ...jobs.filter(job => job.id !== latest.id || job.workspaceId !== latest.workspaceId)] : [...jobs]
  return merged.sort((a, b) => b.createdAt.localeCompare(a.createdAt))
}

/** Select within BOTH the active tool and workspace; never borrow another tool's output. */
export function selectToolManagedJob(jobs: readonly ManagedJobRecord[], workspaceId: string | undefined, command: string | undefined): ManagedJobRecord | null {
  if (!workspaceId || !command) return null
  return selectRestoredManagedJob(jobs.filter(job => job.workspaceId === workspaceId && job.command === command)
    .sort((a, b) => b.createdAt.localeCompare(a.createdAt)))
}

export function outputMatchesManagedJob(output: { jobId: string; workspaceId: string } | null, job: ManagedJobRecord | null, workspaceId: string | undefined): boolean {
  return Boolean(output && job && job.workspaceId === workspaceId && output.jobId === job.id && output.workspaceId === workspaceId)
}

export function managedJobStatusLabel(job: ManagedJobRecord): string {
  if (job.status === 'queued') return `Queued · ${job.command}`
  if (job.status === 'running') return `Running · ${job.command}`
  if (job.status === 'succeeded') return `Completed · ${job.outputFile || job.command}`
  if (job.status === 'cancelled') return `Cancelled · ${job.command}`
  return `Failed${job.exitCode === null ? '' : ` (exit ${job.exitCode})`} · ${job.command}`
}
