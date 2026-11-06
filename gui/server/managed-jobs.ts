import type { IncomingMessage, ServerResponse } from 'node:http'
import crypto from 'node:crypto'
import fs from 'node:fs'
import path from 'node:path'
import os from 'node:os'
import { Value } from '@sinclair/typebox/value'
import { buildArgs } from './command-args'
import { readDiffusionPublication } from './diffusion-bundle'
import { COMMANDS } from './dvbfixer-spec'
import { runDvbfixerArgs, type DvbfixerRunOptions } from './dvbfixer-runner'
import { errorStatus, readRequestBody } from './request-body'
import { assertWorkspaceQuota, WorkspaceQuotaExceededError } from './storage-quota'
import type { ApiRouteHost } from './http-types'
import { ensureRequestId } from './request-observability'
import { acquireWorkspaceRunLock, WorkspaceRunLockError } from './workspace-lock'
import { ManagedJobRequestSchema } from './managed-jobs-schema'
import type { JobProgress } from './job-progress'
import { relocateDiffusionOutcomes, type DiffusionOutcomes } from './diffusion-outcomes'
import { executionPolicy } from './execution-policy'
import { compactRun, copyNativeRun, outputFiles, RUN_ARCHIVE, RUN_SUMMARY, type RunStatus } from './run-storage'
import {
  assertWorkspaceAccess, listWorkspaces, loadWorkspace, resolveWorkspaceFile, updateWorkspace, workspaceRoot, writeJsonAtomic,
} from './workspace-api'

export type ManagedJobStatus = 'queued' | 'running' | 'succeeded' | 'failed' | 'cancelled'

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
  requestedByPrincipalId: string
  error?: string
  progress?: JobProgress
  diffusion?: DiffusionOutcomes
  workflow?: { kind: 'homology'; projectId: string; projectName: string; snapshotSha256?: string }
  modelStatus?: 'success' | 'failed' | 'validation_failed'
  provenance?: ManagedJobProvenance
}

export interface ManagedJobProvenance {
  apiVersion: 1
  command: string
  serviceVersion: string
  dvbfixerVersion: string
  request: {
    inputs: Record<string, string | string[]>
    values: Record<string, unknown>
    fastaSha256: string | null
  }
  resolvedInputs: Array<{ artifactId: string | null; file: string; sha256: string }>
  outputs: Array<{ artifactId: string; file: string; sha256: string }>
}

export interface ManagedJobRequest {
  workspaceId: string
  command: string
  inputFile?: string
  inputs?: Record<string, string | string[]>
  values?: Record<string, any>
  fastaContent?: string
  workflow?: ManagedJobRecord['workflow']
}

interface ActiveRun { owner: string; controller?: AbortController }
const activeRuns = new Map<string, ActiveRun>()
const activeJobWork = new Set<Promise<void>>()
interface JobSubscriber { response: ServerResponse; authorize(): void }
const subscribers = new Map<string, Set<JobSubscriber>>()
let acceptingManagedJobs = true

function trackManagedJob(work: Promise<void>, jobId: string): void {
  activeJobWork.add(work)
  void work.then(() => activeJobWork.delete(work), error => {
    activeJobWork.delete(work)
    console.error(`ERROR: Could not finalize managed job ${jobId}: ${error instanceof Error ? error.message : String(error)}`)
  })
}

export async function drainManagedJobs(): Promise<void> {
  await Promise.allSettled([...activeJobWork])
}

export function resetManagedJobAdmission(): void {
  acceptingManagedJobs = true
}

export function shutdownManagedJobs(): void {
  acceptingManagedJobs = false
  for (const run of activeRuns.values()) run.controller?.abort()
  for (const responses of subscribers.values()) {
    for (const { response } of responses) response.end()
  }
  subscribers.clear()
}

function httpError(statusCode: number, message: string, code = 'INVALID_REQUEST'): Error & { statusCode: number; code: string } {
  const error = new Error(message) as Error & { statusCode: number; code: string }
  error.statusCode = statusCode
  error.code = code
  return error
}

export function acquireWorkspaceRun(
  workspaceId: string,
  owner: string,
  controller?: AbortController,
): (() => void) | null {
  if (activeRuns.has(workspaceId)) return null
  activeRuns.set(workspaceId, { owner, controller })
  return () => {
    if (activeRuns.get(workspaceId)?.owner === owner) activeRuns.delete(workspaceId)
  }
}

function jobDirectory(dataRoot: string, workspaceId: string, jobId: string): string {
  if (!/^[0-9a-f-]{36}$/.test(jobId)) throw httpError(400, 'invalid job id')
  return path.join(workspaceRoot(dataRoot, workspaceId), 'runs', `job_${jobId}`)
}

function recordPath(dataRoot: string, workspaceId: string, jobId: string): string {
  return path.join(jobDirectory(dataRoot, workspaceId, jobId), 'job.json')
}

function cancellationPath(dataRoot: string, workspaceId: string, jobId: string): string {
  return path.join(jobDirectory(dataRoot, workspaceId, jobId), 'cancel.requested')
}

function readCancellationRequest(file: string): string | undefined {
  let descriptor: number
  try { descriptor = fs.openSync(file, fs.constants.O_RDONLY | fs.constants.O_NOFOLLOW) } catch (error) {
    if ((error as NodeJS.ErrnoException).code === 'ENOENT') return undefined
    throw error
  }
  try {
    const stat = fs.fstatSync(descriptor)
    if (!stat.isFile() || stat.nlink !== 1 || stat.size > 64) throw new Error('invalid job cancellation marker')
    const value = fs.readFileSync(descriptor, 'utf8').trim()
    if (!/^\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d\.\d{3}Z$/.test(value) || !Number.isFinite(Date.parse(value))) throw new Error('invalid cancellation timestamp')
    return value
  } finally { fs.closeSync(descriptor) }
}

function persistJob(dataRoot: string, record: ManagedJobRecord): void {
  const cancelledAt = readCancellationRequest(cancellationPath(dataRoot, record.workspaceId, record.id))
  if (cancelledAt) record.cancellationRequestedAt = cancelledAt
  record.revision = (record.revision || 0) + 1
  writeJsonAtomic(recordPath(dataRoot, record.workspaceId, record.id), record)
  const key = `${record.workspaceId}:${record.id}`
  const event = `data: ${JSON.stringify(record)}\n\n`
  const listeners = subscribers.get(key)
  for (const listener of listeners || []) {
    try { listener.authorize(); listener.response.write(event) } catch {
      listener.response.end()
      listeners?.delete(listener)
    }
  }
  if (['succeeded', 'failed', 'cancelled'].includes(record.status)) {
    for (const { response } of listeners || []) response.end()
    subscribers.delete(key)
  }
}

function loadJob(dataRoot: string, workspaceId: string, jobId: string): ManagedJobRecord {
  const file = recordPath(dataRoot, workspaceId, jobId)
  if (!fs.existsSync(file)) throw httpError(404, `job not found: ${jobId}`)
  const record = JSON.parse(fs.readFileSync(file, 'utf8')) as ManagedJobRecord
  if (record.status === 'queued' || record.status === 'running') {
    const cancelledAt = readCancellationRequest(cancellationPath(dataRoot, workspaceId, jobId))
    if (cancelledAt) record.cancellationRequestedAt = cancelledAt
  }
  if ((record.status === 'queued' || record.status === 'running') && !activeRuns.has(workspaceId)) {
    // Another server process may own the durable workspace lease. Only mark a
    // persisted active record orphaned when the lease can actually be acquired.
    try {
      const release = acquireWorkspaceRunLock(dataRoot, workspaceId)
      try {
        record.status = 'failed'
        record.finishedAt = new Date().toISOString()
        record.error = 'DVBFixer server stopped before the job completed'
        persistJob(dataRoot, record)
      } finally {
        release()
      }
    } catch (error) {
      if (!(error instanceof WorkspaceRunLockError)) throw error
    }
  }
  return record
}

export function getManagedJob(dataRoot: string, workspaceId: string, jobId: string): ManagedJobRecord {
  return loadJob(dataRoot, workspaceId, jobId)
}

export function listManagedJobs(dataRoot: string, workspaceId: string): ManagedJobRecord[] {
  const runs = path.join(workspaceRoot(dataRoot, workspaceId), 'runs')
  if (!fs.existsSync(runs)) return []
  return fs.readdirSync(runs, { withFileTypes: true })
    .filter(entry => entry.isDirectory() && entry.name.startsWith('job_'))
    .flatMap(entry => {
      try { return [loadJob(dataRoot, workspaceId, entry.name.slice(4))] } catch { return [] }
    })
    .sort((a, b) => b.createdAt.localeCompare(a.createdAt))
}

/** ACL-filtered monitor view: active work is never hidden behind history pages. */
export function listJobMonitorPage(dataRoot: string, principalId: string, legacyOwnerPrincipalId: string, query: URLSearchParams) {
  const rawLimit = query.get('limit') || '30'
  if (!/^\d+$/.test(rawLimit) || Number(rawLimit) < 1 || Number(rawLimit) > 100) throw httpError(400, 'limit must be an integer from 1 to 100')
  const limit = Number(rawLimit)
  const cursor = query.get('before')
  let before: { time: string; key: string } | null = null
  if (cursor) {
    try {
      if (!/^[a-zA-Z0-9_-]{1,512}$/.test(cursor)) throw new Error('invalid cursor')
      const value = JSON.parse(Buffer.from(cursor, 'base64url').toString('utf8'))
      if (typeof value.time !== 'string' || value.time.length > 64 || !Number.isFinite(Date.parse(value.time)) ||
          typeof value.key !== 'string' || !/^[a-zA-Z0-9_-]+:[0-9a-f-]{36}$/.test(value.key)) throw new Error('invalid cursor')
      before = value
    } catch { throw httpError(400, 'invalid monitor history cursor') }
  }
  const readableWorkspaces = listWorkspaces(dataRoot, principalId, legacyOwnerPrincipalId)
  const all = readableWorkspaces
    .flatMap(workspace => listManagedJobs(dataRoot, workspace.id))
  const active = (job: ManagedJobRecord) => job.status === 'queued' || job.status === 'running'
  const key = (job: ManagedJobRecord) => `${job.workspaceId}:${job.id}`
  const time = (job: ManagedJobRecord) => job.finishedAt || job.createdAt
  const completed = all.filter(job => !active(job))
    .sort((a, b) => time(b).localeCompare(time(a)) || key(b).localeCompare(key(a)))
    .filter(job => !before || time(job) < before.time || (time(job) === before.time && key(job).localeCompare(before.key) < 0))
  const jobs = completed.slice(0, limit)
  const last = jobs.at(-1)
  return {
    workspaceIds: readableWorkspaces.map(workspace => workspace.id),
    activeJobs: all.filter(active).sort((a, b) => b.createdAt.localeCompare(a.createdAt)),
    jobs,
    nextCursor: completed.length > limit && last
      ? Buffer.from(JSON.stringify({ time: time(last), key: key(last) })).toString('base64url') : null,
  }
}

function listFiles(root: string): string[] {
  const files: string[] = []
  const walk = (directory: string) => {
    for (const entry of fs.readdirSync(directory, { withFileTypes: true })) {
      const full = path.join(directory, entry.name)
      if (entry.isDirectory()) walk(full)
      else if (entry.isFile()) files.push(path.relative(root, full).replace(/\\/g, '/'))
    }
  }
  walk(root)
  return files.sort()
}

function prepareArgs(
  dataRoot: string,
  request: ManagedJobRequest,
  directory: string,
  publicationRoot = directory,
): { args: string[]; inputBase: string; resolvedInputs: ManagedJobProvenance['resolvedInputs'] } {
  const def = COMMANDS.find(command => command.name === request.command)
  if (!def) throw httpError(400, `unknown command: ${request.command}`)
  const suppliedInputs = { ...(request.inputs || {}) }
  if (request.inputFile && def.inputs[0] && suppliedInputs[def.inputs[0].dest] === undefined) {
    suppliedInputs[def.inputs[0].dest] = request.inputFile
  }
  const positional: string[] = []
  const resolvedInputs: ManagedJobProvenance['resolvedInputs'] = []
  const workspace = loadWorkspace(dataRoot, request.workspaceId)
  for (const input of def.inputs) {
    const raw = suppliedInputs[input.dest]
    const items = Array.isArray(raw) ? raw : (typeof raw === 'string' && raw ? [raw] : [])
    if (input.required && items.length === 0) throw httpError(400, `${input.label} is required`)
    for (const item of items) {
      const resolved = resolveWorkspaceFile(dataRoot, request.workspaceId, item)
      let commandInput = resolved
      const bundleIndex = path.join(path.dirname(resolved), 'bundle.json')
      if (input.dest === 'input' && fs.existsSync(bundleIndex) && ['prepare', 'minimize', 'zbs'].includes(request.command)) {
        const publication = readDiffusionPublication(bundleIndex)
        if (!publication.candidates.includes(resolved)) throw new Error('input is not a declared diffusion candidate')
        commandInput = path.join(directory, `_input_${path.basename(resolved)}`)
        fs.copyFileSync(resolved, commandInput, fs.constants.COPYFILE_EXCL)
        const sidecar = publication.sidecars[resolved]
        fs.copyFileSync(sidecar, commandInput.replace(/\.pdb$/i, '.dat'), fs.constants.COPYFILE_EXCL)
      }
      positional.push(commandInput)
      const bytes = fs.readFileSync(resolved)
      resolvedInputs.push({
        artifactId: workspace.artifacts.find(artifact => artifact.file === item)?.id || null,
        file: item,
        sha256: crypto.createHash('sha256').update(bytes).digest('hex'),
      })
    }
  }
  const originalInput = resolvedInputs[0]?.file || positional[0]
  const inputBase = originalInput
    ? path.basename(originalInput, path.extname(originalInput)) : request.command
  const values = { ...(request.values || {}) }
  for (const field of def.flags.filter(field => field.type === 'artifact')) {
    const raw = values[field.flag]
    if (Array.isArray(raw)) {
      values[field.flag] = raw.map(item => resolveWorkspaceFile(dataRoot, request.workspaceId, String(item)))
    } else if (typeof raw === 'string' && raw.trim()) {
      const items = field.repeatable ? raw.split(',').map(item => item.trim()).filter(Boolean) : [raw]
      values[field.flag] = field.repeatable
        ? items.map(item => resolveWorkspaceFile(dataRoot, request.workspaceId, item))
        : resolveWorkspaceFile(dataRoot, request.workspaceId, raw)
    }
  }
  if (typeof request.fastaContent === 'string' && request.fastaContent.trim()) {
    const fasta = path.join(directory, `${inputBase}.fasta`)
    fs.writeFileSync(fasta, request.fastaContent)
    values['--fasta'] = fasta
  }
  const outputStem = path.join(publicationRoot, `${inputBase}_${request.command}`)
  const backendMode = def.backendOutputModes?.[String(values['--backend'])]
  const diffusionBundle = backendMode === 'directory'
  const outputTarget = (def.outputMode === 'directory' || diffusionBundle) ? path.join(publicationRoot, 'result')
    : def.outputMode === 'prefix' ? outputStem : `${outputStem}${def.outputExtension || ''}`
  const args = [...positional]
  if (def.hasOutput && def.outputMode !== 'stdout') args.push('-o', outputTarget)
  args.push(...buildArgs(request.command, values))
  if (publicationRoot !== directory) args.push('--diffusion-work-parent', publicationRoot)
  return { args, inputBase, resolvedInputs }
}

function registerOutputs(dataRoot: string, record: ManagedJobRecord, inputBase: string, publishedFiles?: string[], runStatus?: RunStatus): ManagedJobRecord {
  const def = COMMANDS.find(command => command.name === record.command)!
  const directory = jobDirectory(dataRoot, record.workspaceId, record.id)
  const control = new Set(['job.json', 'stdout.log', 'stderr.log', 'cancel.requested'])
  const files = (publishedFiles || listFiles(directory)).filter(file => !control.has(file))
  const modelFlag = record.command === 'zbs' ? '--model-backend' : '--backend'
  const usesDiffusion = record.args[record.args.indexOf(modelFlag) + 1] === 'diffusion'
  const manifest = usesDiffusion ? files.find(file => path.basename(file) === 'bundle.json') : undefined
  const diffusion = manifest ? readDiffusionPublication(path.join(directory, manifest)) : null
  const declared = record.args?.includes('-o') ? record.args[record.args.indexOf('-o') + 1] : null
  const declaredRelative = declared ? (declared.startsWith(`${directory}${path.sep}`)
    ? path.relative(directory, declared).replace(/\\/g, '/') : path.basename(declared)) : null
  const results = files.filter(file => file !== RUN_SUMMARY && file !== RUN_ARCHIVE)
  const primary = record.status !== 'succeeded' ? null : declaredRelative && results.includes(declaredRelative) ? declaredRelative
    : diffusion ? (diffusion.primary ? path.relative(directory, diffusion.primary).replace(/\\/g, '/') : null)
      : results.find(file => /\.(pdb|cif|mmcif)$/i.test(file)) || results[0] || null
  if (files.length === 0) return record
  const folder = path.relative(workspaceRoot(dataRoot, record.workspaceId), directory).replace(/\\/g, '/')
  const outputs: ManagedJobProvenance['outputs'] = []
  updateWorkspace(dataRoot, record.workspaceId, workspace => {
    assertWorkspaceAccess(workspace, record.requestedByPrincipalId, 'writer')
    for (const file of files) {
      const rejected = diffusion?.rejectedCoordinates?.includes(path.join(directory, file)) ||
        (diffusion?.status === 'validation_failed' && file.endsWith('.pdb'))
      const relative = `${folder}/${file}`
      if (workspace.artifacts.some(artifact => artifact.file === relative)) continue
      const artifactId = crypto.randomUUID()
      const digest = crypto.createHash('sha256').update(fs.readFileSync(path.join(directory, file))).digest('hex')
      outputs.push({ artifactId, file: relative, sha256: digest })
      workspace.artifacts.push({
        id: artifactId, file: relative,
        name: file === primary ? `${inputBase} → ${record.command}`
          : `${rejected ? '[rejected] ' : ''}${path.basename(file)}`,
        kind: !rejected && /\.(pdb|cif|mmcif)$/i.test(file) ? 'structure' : 'artifact',
         artifactType: def.outputKind, command: record.command, folder,
        run: { command: record.command, status: runStatus || record.status as RunStatus },
        hidden: path.basename(file).startsWith('_'),
        workflowProvenance: record.provenance && {
          ...record.provenance,
          outputs: [{ artifactId, file: relative, sha256: digest }],
        },
      })
    }
    return workspace
  })
  return {
    ...record,
    outputFile: primary ? `${folder}/${primary}` : null,
    ...(diffusion?.seedOutcomes ? { diffusion: relocateDiffusionOutcomes(diffusion.seedOutcomes,
      path.relative(workspaceRoot(dataRoot, record.workspaceId), diffusion.seedOutcomes.report).split(path.sep).join('/')) } : {}),
    provenance: record.provenance && { ...record.provenance, outputs },
  }
}

async function executeJob(
  dataRoot: string,
  record: ManagedJobRecord,
  request: ManagedJobRequest,
  controller: AbortController,
  release: () => void,
): Promise<void> {
  let current = record
  let nativeRoot: string | undefined
  let nativeCopied = false
  let finalized = false
  const cancelFile = cancellationPath(dataRoot, record.workspaceId, record.id)
  const cancellationPoll = setInterval(() => {
    if (fs.existsSync(cancelFile)) controller.abort()
  }, 250)
  cancellationPoll.unref()
  try {
    const directory = jobDirectory(dataRoot, record.workspaceId, record.id)
    assertWorkspaceAccess(
      loadWorkspace(dataRoot, record.workspaceId), record.requestedByPrincipalId, 'writer',
    )
    const usesDiffusion = (request.command === 'model' && request.values?.['--backend'] === 'diffusion') ||
      (request.command === 'zbs' && request.values?.['--model-backend'] === 'diffusion')
    if (usesDiffusion) nativeRoot = fs.mkdtempSync(path.join(os.tmpdir(), 'dvbfixer-workflow-'))
    const prepared = prepareArgs(dataRoot, request, directory, nativeRoot || directory)
    current = {
      ...current,
      args: prepared.args,
      provenance: current.provenance && {
        ...current.provenance,
        resolvedInputs: prepared.resolvedInputs,
      },
    }
    persistJob(dataRoot, current)
    const result = await runDvbfixerArgs(record.command, prepared.args, directory, {
      signal: controller.signal,
      onStart: () => {
        current = { ...current, status: 'running', startedAt: new Date().toISOString() }
        persistJob(dataRoot, current)
      },
      onProgress: progress => {
        current = { ...current, progress }
        persistJob(dataRoot, current)
      },
    })
    if (nativeRoot) {
      assertWorkspaceAccess(loadWorkspace(dataRoot, record.workspaceId), record.requestedByPrincipalId, 'writer')
      if (result.code === 0 && !controller.signal.aborted) {
        const indexes = listFiles(nativeRoot).filter(file => path.basename(file) === 'bundle.json')
        if (indexes.length !== 1) throw new Error('diffusion workflow did not publish exactly one complete bundle')
        for (const index of indexes) readDiffusionPublication(path.join(nativeRoot, index))
      }
      copyNativeRun(nativeRoot, directory)
      nativeCopied = true
      // Logs must point to retained workspace artifacts, not deleted staging.
      result.stdout = result.stdout.split(nativeRoot).join(directory)
      result.stderr = result.stderr.split(nativeRoot).join(directory)
    }
    assertWorkspaceQuota(workspaceRoot(dataRoot, record.workspaceId))
    fs.writeFileSync(path.join(directory, 'stdout.log'), result.stdout)
    fs.writeFileSync(path.join(directory, 'stderr.log'), result.stderr)
    const definition = COMMANDS.find(command => command.name === record.command)!
    const succeeded = definition.successCodes.includes(result.code)
    if (succeeded && definition.outputMode === 'stdout' && result.stdout) {
      const format = request.values?.['--format'] === 'json' ? 'json' : 'txt'
      fs.writeFileSync(path.join(directory, `${record.command}.${format}`), result.stdout, { mode: 0o600 })
    }
    const cancelled = controller.signal.aborted
    current = {
      ...current,
      status: cancelled ? 'cancelled' : (succeeded ? 'succeeded' : 'failed'),
      exitCode: result.code,
      finishedAt: new Date().toISOString(),
      error: cancelled ? 'DVBfixer run cancelled' : (succeeded ? undefined : result.stderr || `exit code ${result.code}`),
    }
    const files = listFiles(directory)
    const manifests = succeeded ? files.filter(file => path.basename(file) === 'bundle.json') : []
    const publications = manifests.map(file => readDiffusionPublication(path.join(directory, file)))
    let status: RunStatus = publications.some(item => item.status === 'validation_failed') ? 'validation_failed' : current.status as RunStatus
    const declaredPath = current.args.includes('-o') ? current.args[current.args.indexOf('-o') + 1] : null
    const declared = definition.outputMode === 'stdout' ? `${record.command}.${request.values?.['--format'] === 'json' ? 'json' : 'txt'}`
      : declaredPath ? path.relative(nativeRoot || directory, declaredPath).split(path.sep).join('/') : null
    let outputs = current.status === 'succeeded' ? outputFiles(files, declared, manifests.map(file => path.dirname(file))) : []
    const summary = await compactRun(directory, current.command, status, outputs)
    finalized = true
    if (fs.existsSync(cancelFile)) controller.abort()
    if (controller.signal.aborted && current.status !== 'cancelled') {
      current = { ...current, status: 'cancelled', outputFile: null, error: 'DVBFixer run cancelled' }
      status = 'cancelled'
      outputs = []
      summary.status = 'cancelled'
      summary.outputs = []
      summary.warning = 'Cancellation arrived during packaging; unpublished outputs remain in the run directory for inspection.'
      writeJsonAtomic(path.join(directory, RUN_SUMMARY), summary)
    }
    assertWorkspaceQuota(workspaceRoot(dataRoot, record.workspaceId))
    current = registerOutputs(dataRoot, current, prepared.inputBase,
      [...outputs, RUN_SUMMARY, ...(summary.diagnostics ? [RUN_ARCHIVE] : [])], status)
  } catch (error: any) {
    if (error instanceof WorkspaceQuotaExceededError) {
      const directory = jobDirectory(dataRoot, record.workspaceId, record.id)
      fs.rmSync(directory, { recursive: true, force: true })
      fs.mkdirSync(directory, { recursive: true })
    }
    current = {
      ...current, status: controller.signal.aborted ? 'cancelled' : 'failed',
      finishedAt: new Date().toISOString(), error: error?.message || String(error),
    }
  } finally {
    if (nativeRoot && !nativeCopied) {
      try {
        copyNativeRun(nativeRoot, jobDirectory(dataRoot, record.workspaceId, record.id))
        nativeCopied = true
      } catch {
        current.error = `${current.error || 'Run failed'}. Native diagnostics retained at ${nativeRoot}`
      }
    }
    if (nativeRoot && nativeCopied) fs.rmSync(nativeRoot, { recursive: true, force: true })
    clearInterval(cancellationPoll)
    fs.rmSync(cancelFile, { force: true })
    const directory = jobDirectory(dataRoot, record.workspaceId, record.id)
    const stdout = path.join(directory, 'stdout.log')
    const stderr = path.join(directory, 'stderr.log')
    if (!fs.existsSync(stdout)) fs.writeFileSync(stdout, '')
    if (!fs.existsSync(stderr)) fs.writeFileSync(stderr, current.error || '')
    if (!finalized) {
      try {
        const summary = await compactRun(directory, current.command, current.status as RunStatus, [])
        assertWorkspaceQuota(workspaceRoot(dataRoot, record.workspaceId))
        current = registerOutputs(dataRoot, current, current.command,
          [RUN_SUMMARY, ...(summary.diagnostics ? [RUN_ARCHIVE] : [])], current.status as RunStatus)
      } catch { /* Keep original failure and all evidence, even if publication fails. */ }
    }
    try { persistJob(dataRoot, current) } finally { release() }
  }
}

export function createManagedJob(
  dataRoot: string,
  request: ManagedJobRequest,
  requestedByPrincipalId = 'local',
  serviceVersion = 'unknown',
): ManagedJobRecord {
  if (!acceptingManagedJobs) throw httpError(503, 'server is shutting down')
  if (!request.workspaceId) throw httpError(400, 'workspaceId is required')
  if (typeof request.command !== 'string' || !request.command.trim()) {
    throw httpError(400, 'command is required')
  }
  loadWorkspace(dataRoot, request.workspaceId)
  const policy = executionPolicy(request.command)
  if (!policy) {
    throw httpError(422, `unknown command: ${request.command}`, 'UNSUPPORTED_OPERATION')
  }
  if (policy !== 'managed-workflow') {
    throw httpError(
      422,
      `${request.command} uses its dedicated synchronous V1 endpoint`,
      'EXECUTION_POLICY_MISMATCH',
    )
  }
  const reserved = reserveManagedJob(dataRoot, request, requestedByPrincipalId, serviceVersion)
  trackManagedJob(executeJob(dataRoot, reserved.record, request, reserved.controller, reserved.release), reserved.record.id)
  return reserved.record
}

function reserveManagedJob(dataRoot: string, request: ManagedJobRequest, requestedByPrincipalId: string, serviceVersion: string) {
  if (!acceptingManagedJobs) throw httpError(503, 'server is shutting down')
  assertWorkspaceQuota(workspaceRoot(dataRoot, request.workspaceId), 64 * 1024)
  const id = crypto.randomUUID()
  const controller = new AbortController()
  const releaseDurable = acquireWorkspaceRunLock(dataRoot, request.workspaceId)
  const releaseLocal = acquireWorkspaceRun(request.workspaceId, id, controller)
  if (!releaseLocal) {
    releaseDurable()
    throw httpError(409, 'another DVBfixer job is already running in this workspace')
  }
  const release = () => {
    releaseLocal()
    releaseDurable()
  }
  const directory = jobDirectory(dataRoot, request.workspaceId, id)
  try {
    fs.mkdirSync(directory, { recursive: false })
    const now = new Date().toISOString()
    const folder = path.relative(workspaceRoot(dataRoot, request.workspaceId), directory).replace(/\\/g, '/')
    const record: ManagedJobRecord = {
      version: 1, id, workspaceId: request.workspaceId, status: 'queued',
      createdAt: now, startedAt: null, finishedAt: null, command: request.command, args: [],
      exitCode: null, stdoutLog: `${folder}/stdout.log`, stderrLog: `${folder}/stderr.log`,
      outputDir: folder, outputFile: null, requestedByPrincipalId,
      ...(request.workflow ? { workflow: request.workflow } : {}),
      provenance: {
        apiVersion: 1,
        command: request.command,
        serviceVersion,
        dvbfixerVersion: serviceVersion,
        request: {
          inputs: {
            ...(request.inputFile ? { inputFile: request.inputFile } : {}),
            ...(request.inputs || {}),
          },
          values: { ...(request.values || {}) },
          fastaSha256: typeof request.fastaContent === 'string' && request.fastaContent.length
            ? crypto.createHash('sha256').update(request.fastaContent).digest('hex') : null,
        },
        resolvedInputs: [],
        outputs: [],
      },
    }
    for (const name of ['stdout.log', 'stderr.log']) fs.writeFileSync(path.join(directory, name), '')
    persistJob(dataRoot, record)
    return { record, controller, release, directory }
  } catch (error) {
    release()
    throw error
  }
}

export interface ManagedWorkflowContext {
  directory: string
  signal: AbortSignal
  isCancelled(): boolean
  processOptions: Pick<DvbfixerRunOptions, 'onStart' | 'onProgress'>
  setArguments(args: string[], resolvedInputs: ManagedJobProvenance['resolvedInputs']): void
}

export interface ManagedWorkflowResult {
  exitCode: number
  outputFile: string
  modelStatus: 'success' | 'failed' | 'validation_failed'
  diffusion?: DiffusionOutcomes
  stderr: string
}

/** Attach a composed workflow to the SAME durable jobs, leases and cancellation.
 * The operation owns scientific execution and normal packaging, not scheduling.
 * HTTP disconnects are deliberately not a cancellation signal for this path.
 */
export function createManagedWorkflowJob(
  dataRoot: string, request: ManagedJobRequest,
  operation: (context: ManagedWorkflowContext) => Promise<ManagedWorkflowResult>,
  requestedByPrincipalId = 'local', serviceVersion = 'unknown',
): ManagedJobRecord {
  assertWorkspaceAccess(loadWorkspace(dataRoot, request.workspaceId), requestedByPrincipalId, 'writer')
  const { record, controller, release, directory } = reserveManagedJob(dataRoot, request, requestedByPrincipalId, serviceVersion)
  let current = record
  const cancelFile = cancellationPath(dataRoot, record.workspaceId, record.id)
  const cancellationPoll = setInterval(() => { if (fs.existsSync(cancelFile)) controller.abort() }, 250)
  cancellationPoll.unref()
  const update = (patch: Partial<ManagedJobRecord>) => { current = { ...current, ...patch }; persistJob(dataRoot, current) }
  const work = Promise.resolve().then(async () => {
    try {
      if (controller.signal.aborted) throw new Error('DVBFixer run cancelled before execution')
      const result = await operation({ directory, signal: controller.signal,
        isCancelled: () => {
          if (fs.existsSync(cancelFile)) controller.abort()
          return controller.signal.aborted
        },
        processOptions: {
          onStart: () => {
            assertWorkspaceAccess(loadWorkspace(dataRoot, record.workspaceId), requestedByPrincipalId, 'writer')
            update({ status: 'running', startedAt: new Date().toISOString() })
          },
          onProgress: progress => update({ progress }),
        },
        setArguments: (args, resolvedInputs) => update({ args, provenance: current.provenance && { ...current.provenance, resolvedInputs } }),
      })
      const cancelled = controller.signal.aborted
      const succeeded = result.exitCode === 0 && !cancelled
      current = { ...current, status: cancelled ? 'cancelled' : succeeded ? 'succeeded' : 'failed',
        exitCode: result.exitCode, outputFile: succeeded ? result.outputFile || null : null,
        modelStatus: result.modelStatus, ...(succeeded && result.diffusion ? { diffusion: result.diffusion } : {}),
        error: cancelled ? 'DVBFixer run cancelled' : succeeded ? undefined : result.stderr || `exit code ${result.exitCode}` }
      const workspace = loadWorkspace(dataRoot, record.workspaceId)
      assertWorkspaceAccess(workspace, requestedByPrincipalId, 'writer')
      const outputs = workspace.artifacts.filter(artifact => artifact.file.startsWith(`${record.outputDir}/`))
        .map(artifact => ({ artifactId: artifact.id, file: artifact.file,
          sha256: crypto.createHash('sha256').update(fs.readFileSync(resolveWorkspaceFile(dataRoot, record.workspaceId, artifact.file))).digest('hex') }))
      if (current.provenance) {
        const provenance = { ...current.provenance, outputs }
        current.provenance = provenance
        const byId = new Map(outputs.map(output => [output.artifactId, output]))
        updateWorkspace(dataRoot, record.workspaceId, latest => {
          assertWorkspaceAccess(latest, requestedByPrincipalId, 'writer')
          for (const artifact of latest.artifacts) {
            const output = byId.get(artifact.id)
            if (output && output.file === artifact.file) artifact.workflowProvenance = { ...provenance, outputs: [output] }
          }
          return latest
        })
      }
    } catch (error) {
      current = { ...current, status: controller.signal.aborted ? 'cancelled' : 'failed', outputFile: null,
        error: error instanceof Error ? error.message : String(error) }
    } finally {
      clearInterval(cancellationPoll)
      try {
        fs.rmSync(cancelFile, { force: true })
        current.finishedAt = new Date().toISOString()
        try {
          for (const [name, fallback] of [['stdout.log', ''], ['stderr.log', current.error || '']]) {
            if (!fs.existsSync(path.join(directory, name)) || fs.statSync(path.join(directory, name)).size === 0) {
              fs.writeFileSync(path.join(directory, name), fallback)
            }
          }
          // Preparation errors can precede the composed workflow's own packaging.
          if (!fs.existsSync(path.join(directory, RUN_SUMMARY))) {
            await compactRun(directory, current.command, current.status as RunStatus, [])
          }
        } catch { /* Preserve the original failure and all unarchived evidence. */ }
        persistJob(dataRoot, current)
      } finally { release() }
    }
  })
  trackManagedJob(work, record.id)
  return record
}

export function cancelManagedJob(dataRoot: string, workspaceId: string, jobId: string): ManagedJobRecord {
  const record = loadJob(dataRoot, workspaceId, jobId)
  if (record.status !== 'queued' && record.status !== 'running') return record
  const active = activeRuns.get(workspaceId)
  try { fs.writeFileSync(cancellationPath(dataRoot, workspaceId, jobId), `${new Date().toISOString()}\n`, { mode: 0o600, flag: 'wx' }) } catch (error) {
    if ((error as NodeJS.ErrnoException).code !== 'EEXIST') throw error
  }
  if (active?.owner === jobId && active.controller) active.controller.abort()
  return loadJob(dataRoot, workspaceId, jobId)
}

function sendJson(response: ServerResponse, status: number, body: unknown): void {
  response.statusCode = status
  response.setHeader('Content-Type', 'application/json')
  response.setHeader('Cache-Control', 'no-store')
  response.end(JSON.stringify(body))
}

export function sendV1Error(
  response: ServerResponse,
  error: unknown,
  requestId: string,
): void {
  const candidate = error as { statusCode?: number; code?: string; message?: string }
  const status = errorStatus(error)
  const defaultCodes: Record<number, string> = {
    400: 'INVALID_REQUEST', 401: 'AUTHENTICATION_REQUIRED',
    403: 'WORKSPACE_ACCESS_DENIED', 404: 'NOT_FOUND',
    405: 'METHOD_NOT_ALLOWED',
    409: 'WORKSPACE_BUSY', 413: 'PAYLOAD_TOO_LARGE',
    415: 'UNSUPPORTED_MEDIA_TYPE', 422: 'UNPROCESSABLE_OPERATION',
    500: 'INTERNAL_ERROR', 503: 'SERVICE_UNAVAILABLE', 507: 'WORKSPACE_QUOTA_EXCEEDED',
  }
  sendJson(response, status, {
    error: {
      code: candidate.code || defaultCodes[status] || 'INTERNAL_ERROR',
      message: status >= 500 && status !== 503 && status !== 507
        ? 'unexpected server error' : candidate.message || 'request failed',
      requestId,
    },
  })
}

export function registerManagedJobApi(
  server: ApiRouteHost,
  dataRoot: string,
  principalSource: (request: IncomingMessage) => string = () => 'local',
  legacyOwnerPrincipalId = 'local',
  serviceVersion = 'unknown',
): void {
  server.middlewares.use('/api/v1', async (request: IncomingMessage, response: ServerResponse, next) => {
    const route = (request.url || '').split('?')[0]
    const match = route.match(/^\/workspaces\/([^/]+)\/jobs(?:\/([^/]+)(?:\/(events))?)?$/)
    if (!match && route !== '/jobs') return next()
    const requestId = ensureRequestId(request, response)
    try {
      const principalId = principalSource(request)
      if (route === '/jobs') {
        if (request.method !== 'GET') throw httpError(405, 'method not allowed', 'METHOD_NOT_ALLOWED')
        return sendJson(response, 200, listJobMonitorPage(dataRoot, principalId, legacyOwnerPrincipalId,
          new URL(request.url || '', 'http://localhost').searchParams))
      }
      if (!match) return next()
      let workspaceId: string
      try { workspaceId = decodeURIComponent(match[1]) } catch {
        throw httpError(400, 'workspace ID is not valid URL encoding', 'INVALID_WORKSPACE_ID')
      }
      if (!/^[a-zA-Z0-9_-]+$/.test(workspaceId)) {
        throw httpError(400, 'workspace ID contains unsupported characters', 'INVALID_WORKSPACE_ID')
      }
      const jobId = match[2] ? decodeURIComponent(match[2]) : ''
      const events = match[3] === 'events'
      if (request.method === 'POST' && !jobId) {
        if (!String(request.headers['content-type'] || '').toLowerCase().startsWith('application/json')) {
          throw httpError(415, 'Content-Type must be application/json', 'UNSUPPORTED_MEDIA_TYPE')
        }
        let body: ManagedJobRequest
        try {
          body = JSON.parse((await readRequestBody(request)).toString('utf8') || '{}') as ManagedJobRequest
        } catch (error) {
          if (errorStatus(error) === 413) throw error
          throw httpError(400, 'request body is not valid JSON', 'INVALID_JSON')
        }
        if (!Value.Check(ManagedJobRequestSchema, body)) {
          throw httpError(400, 'request body does not match ManagedJobRequest', 'INVALID_REQUEST')
        }
        if (body.workspaceId !== undefined && body.workspaceId !== workspaceId) {
          throw httpError(400, 'body workspaceId must match the route', 'WORKSPACE_ID_MISMATCH')
        }
        body.workspaceId = workspaceId
        assertWorkspaceAccess(loadWorkspace(dataRoot, workspaceId, legacyOwnerPrincipalId), principalId, 'writer')
        return sendJson(response, 202, createManagedJob(dataRoot, body, principalId, serviceVersion))
      }
      if (request.method === 'GET' && !jobId) {
        assertWorkspaceAccess(loadWorkspace(dataRoot, workspaceId, legacyOwnerPrincipalId), principalId)
        return sendJson(response, 200, listManagedJobs(dataRoot, workspaceId))
      }
      const workspace = loadWorkspace(dataRoot, workspaceId, legacyOwnerPrincipalId)
      assertWorkspaceAccess(workspace, principalId, request.method === 'DELETE' ? 'writer' : 'reader')
      if (request.method === 'GET' && events) {
        let record = getManagedJob(dataRoot, workspaceId, jobId)
        response.writeHead(200, {
          'Content-Type': 'text/event-stream; charset=utf-8', 'Cache-Control': 'no-cache, no-transform',
          Connection: 'keep-alive', 'X-Accel-Buffering': 'no',
        })
        response.write(`data: ${JSON.stringify(record)}\n\n`)
        if (['succeeded', 'failed', 'cancelled'].includes(record.status)) return response.end()
        const key = `${workspaceId}:${record.id}`
        const listeners = subscribers.get(key) || new Set<JobSubscriber>()
        const listener: JobSubscriber = { response, authorize: () => {
          assertWorkspaceAccess(loadWorkspace(dataRoot, workspaceId, legacyOwnerPrincipalId), principalId)
        } }
        listeners.add(listener)
        subscribers.set(key, listeners)
        let last = JSON.stringify(record)
        const poll = setInterval(() => {
          try {
            listener.authorize()
            record = getManagedJob(dataRoot, workspaceId, jobId)
            const serialized = JSON.stringify(record)
            if (serialized !== last) {
              last = serialized
              response.write(`data: ${serialized}\n\n`)
            }
            if (['succeeded', 'failed', 'cancelled'].includes(record.status)) response.end()
          } catch { response.end() }
        }, 500)
        poll.unref()
        const close = () => {
          clearInterval(poll)
          listeners.delete(listener)
          if (!listeners.size) subscribers.delete(key)
        }
        request.on('close', close)
        response.on('close', close)
        return
      }
      if (request.method === 'GET' && jobId && !events) {
        return sendJson(response, 200, getManagedJob(dataRoot, workspaceId, jobId))
      }
      if (request.method === 'DELETE' && jobId && !events) {
        return sendJson(response, 202, cancelManagedJob(dataRoot, workspaceId, jobId))
      }
      throw httpError(405, 'method not allowed', 'METHOD_NOT_ALLOWED')
    } catch (error: unknown) {
      return sendV1Error(response, error, requestId)
    }
  })
}
