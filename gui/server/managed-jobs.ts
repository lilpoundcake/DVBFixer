import type { IncomingMessage, ServerResponse } from 'node:http'
import type { ViteDevServer } from 'vite'
import crypto from 'node:crypto'
import fs from 'node:fs'
import path from 'node:path'
import { buildArgs } from './command-args'
import { COMMANDS } from './dvbfixer-spec'
import { runDvbfixerArgs } from './dvbfixer-runner'
import { errorStatus, readRequestBody } from './request-body'
import {
  loadWorkspace, resolveWorkspaceFile, saveWorkspace, workspaceRoot, writeJsonAtomic,
} from './workspace-api'

export type ManagedJobStatus = 'queued' | 'running' | 'succeeded' | 'failed' | 'cancelled'

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
  error?: string
}

export interface ManagedJobRequest {
  workspaceId: string
  command: string
  inputFile?: string
  inputs?: Record<string, string | string[]>
  values?: Record<string, any>
  fastaContent?: string
}

interface ActiveRun { owner: string; controller?: AbortController }
const activeRuns = new Map<string, ActiveRun>()
const subscribers = new Map<string, Set<ServerResponse>>()

function httpError(statusCode: number, message: string): Error & { statusCode: number } {
  const error = new Error(message) as Error & { statusCode: number }
  error.statusCode = statusCode
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

function persistJob(dataRoot: string, record: ManagedJobRecord): void {
  writeJsonAtomic(recordPath(dataRoot, record.workspaceId, record.id), record)
  const key = `${record.workspaceId}:${record.id}`
  const event = `data: ${JSON.stringify(record)}\n\n`
  for (const response of subscribers.get(key) || []) response.write(event)
  if (['succeeded', 'failed', 'cancelled'].includes(record.status)) {
    for (const response of subscribers.get(key) || []) response.end()
    subscribers.delete(key)
  }
}

function loadJob(dataRoot: string, workspaceId: string, jobId: string): ManagedJobRecord {
  const file = recordPath(dataRoot, workspaceId, jobId)
  if (!fs.existsSync(file)) throw httpError(404, `job not found: ${jobId}`)
  const record = JSON.parse(fs.readFileSync(file, 'utf8')) as ManagedJobRecord
  if ((record.status === 'queued' || record.status === 'running') && !activeRuns.has(workspaceId)) {
    record.status = 'failed'
    record.finishedAt = new Date().toISOString()
    record.error = 'GUI server stopped before the job completed'
    persistJob(dataRoot, record)
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
): { args: string[]; inputBase: string } {
  const def = COMMANDS.find(command => command.name === request.command)
  if (!def) throw httpError(400, `unknown command: ${request.command}`)
  const suppliedInputs = { ...(request.inputs || {}) }
  if (request.inputFile && def.inputs[0] && suppliedInputs[def.inputs[0].dest] === undefined) {
    suppliedInputs[def.inputs[0].dest] = request.inputFile
  }
  const positional: string[] = []
  for (const input of def.inputs) {
    const raw = suppliedInputs[input.dest]
    const items = Array.isArray(raw) ? raw : (typeof raw === 'string' && raw ? [raw] : [])
    if (input.required && items.length === 0) throw httpError(400, `${input.label} is required`)
    for (const item of items) positional.push(resolveWorkspaceFile(dataRoot, request.workspaceId, item))
  }
  const inputBase = positional[0]
    ? path.basename(positional[0], path.extname(positional[0])) : request.command
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
  const outputStem = path.join(directory, `${inputBase}_${request.command}`)
  const outputTarget = def.outputMode === 'directory' ? path.join(directory, 'result')
    : def.outputMode === 'prefix' ? outputStem : `${outputStem}${def.outputExtension || ''}`
  const args = [...positional]
  if (def.hasOutput && def.outputMode !== 'stdout') args.push('-o', outputTarget)
  args.push(...buildArgs(request.command, values))
  return { args, inputBase }
}

function registerOutputs(dataRoot: string, record: ManagedJobRecord, inputBase: string): ManagedJobRecord {
  const def = COMMANDS.find(command => command.name === record.command)!
  const directory = jobDirectory(dataRoot, record.workspaceId, record.id)
  const control = new Set(['job.json', 'stdout.log', 'stderr.log'])
  const files = listFiles(directory).filter(file => !control.has(file))
  const primary = files.find(file => /\.(pdb|cif|mmcif)$/i.test(file)) || files[0] || null
  if (!primary) return record
  const workspace = loadWorkspace(dataRoot, record.workspaceId)
  const folder = path.relative(workspaceRoot(dataRoot, record.workspaceId), directory).replace(/\\/g, '/')
  for (const file of files) {
    const relative = `${folder}/${file}`
    if (workspace.artifacts.some(artifact => artifact.file === relative)) continue
    workspace.artifacts.push({
      id: crypto.randomUUID(), file: relative,
      name: file === primary ? `${inputBase} → ${record.command}` : path.basename(file),
      kind: /\.(pdb|cif|mmcif)$/i.test(file) ? 'structure' : 'artifact',
      artifactType: def.outputKind, command: record.command, folder,
      hidden: path.basename(file).startsWith('_'),
    })
  }
  saveWorkspace(dataRoot, workspace)
  return { ...record, outputFile: `${folder}/${primary}` }
}

async function executeJob(
  dataRoot: string,
  record: ManagedJobRecord,
  request: ManagedJobRequest,
  controller: AbortController,
  release: () => void,
): Promise<void> {
  let current = record
  try {
    const directory = jobDirectory(dataRoot, record.workspaceId, record.id)
    const prepared = prepareArgs(dataRoot, request, directory)
    current = { ...current, args: prepared.args, status: 'running', startedAt: new Date().toISOString() }
    persistJob(dataRoot, current)
    const result = await runDvbfixerArgs(record.command, prepared.args, directory, { signal: controller.signal })
    fs.writeFileSync(path.join(directory, 'stdout.log'), result.stdout)
    fs.writeFileSync(path.join(directory, 'stderr.log'), result.stderr)
    const cancelled = controller.signal.aborted
    current = {
      ...current,
      status: cancelled ? 'cancelled' : (result.code === 0 ? 'succeeded' : 'failed'),
      exitCode: result.code,
      finishedAt: new Date().toISOString(),
      error: cancelled ? 'DVBfixer run cancelled' : (result.code === 0 ? undefined : result.stderr || `exit code ${result.code}`),
    }
    if (current.status === 'succeeded') current = registerOutputs(dataRoot, current, prepared.inputBase)
  } catch (error: any) {
    current = {
      ...current, status: controller.signal.aborted ? 'cancelled' : 'failed',
      finishedAt: new Date().toISOString(), error: error?.message || String(error),
    }
  } finally {
    const directory = jobDirectory(dataRoot, record.workspaceId, record.id)
    const stdout = path.join(directory, 'stdout.log')
    const stderr = path.join(directory, 'stderr.log')
    if (!fs.existsSync(stdout)) fs.writeFileSync(stdout, '')
    if (!fs.existsSync(stderr)) fs.writeFileSync(stderr, current.error || '')
    persistJob(dataRoot, current)
    release()
  }
}

export function createManagedJob(dataRoot: string, request: ManagedJobRequest): ManagedJobRecord {
  if (!request.workspaceId) throw httpError(400, 'workspaceId is required')
  if (typeof request.command !== 'string' || !request.command.trim()) {
    throw httpError(400, 'command is required')
  }
  loadWorkspace(dataRoot, request.workspaceId)
  if (!COMMANDS.some(command => command.name === request.command)) throw httpError(400, `unknown command: ${request.command}`)
  const id = crypto.randomUUID()
  const controller = new AbortController()
  const release = acquireWorkspaceRun(request.workspaceId, id, controller)
  if (!release) throw httpError(409, 'another DVBfixer job is already running in this workspace')
  const directory = jobDirectory(dataRoot, request.workspaceId, id)
  fs.mkdirSync(directory, { recursive: false })
  const now = new Date().toISOString()
  const folder = path.relative(workspaceRoot(dataRoot, request.workspaceId), directory).replace(/\\/g, '/')
  const record: ManagedJobRecord = {
    version: 1, id, workspaceId: request.workspaceId, status: 'queued',
    createdAt: now, startedAt: null, finishedAt: null, command: request.command, args: [],
    exitCode: null, stdoutLog: `${folder}/stdout.log`, stderrLog: `${folder}/stderr.log`,
    outputDir: folder, outputFile: null,
  }
  persistJob(dataRoot, record)
  void executeJob(dataRoot, record, request, controller, release)
  return record
}

export function cancelManagedJob(dataRoot: string, workspaceId: string, jobId: string): ManagedJobRecord {
  const record = loadJob(dataRoot, workspaceId, jobId)
  if (record.status !== 'queued' && record.status !== 'running') return record
  const active = activeRuns.get(workspaceId)
  if (!active || active.owner !== jobId || !active.controller) throw httpError(409, 'job is not active in this server')
  active.controller.abort()
  return record
}

function sendJson(response: ServerResponse, status: number, body: unknown): void {
  response.statusCode = status
  response.setHeader('Content-Type', 'application/json')
  response.setHeader('Cache-Control', 'no-store')
  response.end(JSON.stringify(body))
}

export function registerManagedJobApi(server: ViteDevServer, dataRoot: string): void {
  server.middlewares.use('/api/jobs', async (request: IncomingMessage, response: ServerResponse, next) => {
    try {
      const url = new URL(request.url || '/', 'http://localhost')
      const parts = url.pathname.split('/').filter(Boolean)
      const workspaceId = url.searchParams.get('workspaceId') || ''
      if (request.method === 'POST' && parts.length === 0) {
        const body = JSON.parse((await readRequestBody(request)).toString('utf8') || '{}') as ManagedJobRequest
        return sendJson(response, 202, createManagedJob(dataRoot, body))
      }
      if (request.method === 'GET' && parts.length === 0) {
        if (!workspaceId) return sendJson(response, 400, { error: 'workspaceId is required' })
        return sendJson(response, 200, listManagedJobs(dataRoot, workspaceId))
      }
      if (!workspaceId || !parts[0]) return sendJson(response, 400, { error: 'workspaceId and job id are required' })
      if (request.method === 'GET' && parts[1] === 'events') {
        const record = getManagedJob(dataRoot, workspaceId, parts[0])
        response.writeHead(200, {
          'Content-Type': 'text/event-stream; charset=utf-8', 'Cache-Control': 'no-cache, no-transform',
          Connection: 'keep-alive', 'X-Accel-Buffering': 'no',
        })
        response.write(`data: ${JSON.stringify(record)}\n\n`)
        if (['succeeded', 'failed', 'cancelled'].includes(record.status)) return response.end()
        const key = `${workspaceId}:${record.id}`
        const listeners = subscribers.get(key) || new Set<ServerResponse>()
        listeners.add(response)
        subscribers.set(key, listeners)
        request.on('close', () => {
          listeners.delete(response)
          if (!listeners.size) subscribers.delete(key)
        })
        return
      }
      if (request.method === 'GET' && parts.length === 1) {
        return sendJson(response, 200, getManagedJob(dataRoot, workspaceId, parts[0]))
      }
      if (request.method === 'DELETE' && parts.length === 1) {
        return sendJson(response, 202, cancelManagedJob(dataRoot, workspaceId, parts[0]))
      }
      return next()
    } catch (error: any) {
      return sendJson(response, errorStatus(error), { error: error?.message || String(error) })
    }
  })
}
