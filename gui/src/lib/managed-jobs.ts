import type { components } from '../generated/dvbfixer-api'

export type DiffusionOutcomes = NonNullable<components['schemas']['ManagedJobRecord']['diffusion']>
export type JobProgress = NonNullable<components['schemas']['ManagedJobRecord']['progress']>
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
  revision?: number
  cancellationRequestedAt?: string
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
  diffusion?: DiffusionOutcomes
  workflow?: components['schemas']['ManagedJobRecord']['workflow']
  modelStatus?: components['schemas']['ManagedJobRecord']['modelStatus']
  progress?: JobProgress
}

/** Counts are stage-local, not whole-job readiness or scientific acceptance. */
export function managedJobProgress(job: ManagedJobRecord): { label: string; percent?: number } {
  const progress = job.progress
  if (!progress || progress.version !== 1) return { label: 'Waiting for stage events' }
  const labels = {
    preparation: 'Input preparation', runner: 'Runner (sampling / refinement)', sampling: 'Sampling',
    refinement: 'Refinement', validation: 'Independent validation', publication: 'Publication',
  }
  const stage = labels[progress.stage]
  if (!stage) return { label: 'Waiting for stage events' }
  const label = stage + (typeof progress.seed === 'number' && Number.isSafeInteger(progress.seed) && progress.seed >= 0 ? ` · Seed ${progress.seed}` : '')
  const { completed, total } = 'completed' in progress ? progress : { completed: undefined, total: undefined }
  if (typeof completed !== 'number' || typeof total !== 'number' || !Number.isSafeInteger(completed) ||
    !Number.isSafeInteger(total) || total <= 0 || completed < 0 || completed > total) return { label }
  return { label: `${label} · ${completed}/${total}`, percent: 100 * completed / total }
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

/** Do not let an older GET snapshot overwrite a newer SSE update. */
export function acceptManagedJobUpdate(current: ManagedJobRecord | null, next: ManagedJobRecord): boolean {
  if (!current || current.id !== next.id || current.workspaceId !== next.workspaceId) return true
  if (current.revision !== undefined && next.revision !== undefined) return next.revision >= current.revision
  return isActiveManagedJob(current) || !isActiveManagedJob(next)
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
