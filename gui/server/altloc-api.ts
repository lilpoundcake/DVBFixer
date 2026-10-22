/** Workspace-scoped CLI adapter. No selection science is implemented in Node. */
import crypto from 'node:crypto'
import fs from 'node:fs'
import path from 'node:path'
import { Type } from '@sinclair/typebox'
import { Value } from '@sinclair/typebox/value'
import type { IncomingMessage } from 'node:http'
import type { ApiRouteHost } from './http-types'
import { runDvbfixerArgs } from './dvbfixer-runner'
import { readRequestBody, errorStatus } from './request-body'
import { acquireWorkspaceRunLock } from './workspace-lock'
import { assertWorkspaceQuota } from './storage-quota'
import { assertWorkspaceAccess, loadWorkspace, resolveWorkspaceFile, updateWorkspace, workspaceRoot, type WorkspaceArtifact } from './workspace-api'
import type { AltlocRequest, AltlocReport, AltlocResponse } from '../src/lib/altloc-types'

const label = Type.String({ minLength: 1, maxLength: 1, pattern: '^\\S$' })
const digest = Type.String({ pattern: '^[a-f0-9]{64}$' })
const selection = Type.Record(Type.String({ maxLength: 200 }), label, { maxProperties: 10000 })
const policy = Type.Union(['manual', 'occupancy', 'first', 'prefer'].map(value => Type.Literal(value)))
const fallback = Type.Union(['error', 'occupancy', 'first'].map(value => Type.Literal(value)))
const requestSchema = Type.Object({
  inputArtifactId: Type.String({ minLength: 1, maxLength: 200 }),
  action: Type.Union([Type.Literal('analyze'), Type.Literal('apply')]), policy,
  prefer: Type.Optional(label), fallback,
  plan: Type.Optional(Type.Object({ schemaVersion: Type.Literal(1), sourceSha256: digest, selections: selection }, { additionalProperties: false })),
}, { additionalProperties: false })
const strings = Type.Array(Type.String(), { maxItems: 100000 })
const serials = Type.Array(Type.Integer({ minimum: 1, maximum: 99999 }), { maxItems: 100000 })
const rowSchema = Type.Object({
  key: Type.String(), model: Type.Integer({ minimum: 1 }), chain: Type.String({ minLength: 1, maxLength: 1 }),
  residueNumber: Type.String(), insertionCode: Type.String({ maxLength: 1 }), occurrence: Type.Integer({ minimum: 1 }),
  selected: Type.Union([label, Type.Null()]), method: Type.String(), sharedSerials: serials, warnings: strings,
  states: Type.Array(Type.Object({ label, residueName: Type.String(), occupancy: Type.Union([Type.Number({ minimum: 0, maximum: 1 }), Type.Null()]),
    atomNames: strings, serials, atomCount: Type.Integer({ minimum: 1 }) }), { minItems: 1, maxItems: 100 }),
})
const reportSchema = Type.Object({
  schemaVersion: Type.Literal(1), sourceSha256: digest, normalizedSha256: digest, policy,
  prefer: Type.Union([label, Type.Null()]), fallback, overrides: selection,
  residues: Type.Array(rowSchema, { maxItems: 10000 }), warnings: strings, unresolved: strings,
  inputAtoms: Type.Integer({ minimum: 1 }), removedAtoms: Type.Union([Type.Integer({ minimum: 0 }), Type.Null()]),
  outputSha256: Type.Union([digest, Type.Null()]), outputWritten: Type.Boolean(),
  normalizedPdb: Type.Optional(Type.String()), toolVersion: Type.String(), chainMap: Type.Record(Type.String(), Type.String()),
})

class AltlocApiError extends Error {
  statusCode: number
  constructor(statusCode: number, message: string) { super(message); this.statusCode = statusCode }
}
const sha256 = (bytes: Buffer | string) => crypto.createHash('sha256').update(bytes).digest('hex')
const MAX_SOURCE = 10 * 1024 * 1024
function readFile(file: string, max: number): Buffer {
  const fd = fs.openSync(file, fs.constants.O_RDONLY | fs.constants.O_NOFOLLOW)
  try {
    const stat = fs.fstatSync(fd)
    if (!stat.isFile() || stat.size > max) throw new AltlocApiError(413, 'File exceeds the altloc operation limit')
    const bytes = fs.readFileSync(fd)
    if (bytes.length > max) throw new AltlocApiError(413, 'File exceeds the altloc operation limit')
    return bytes
  } finally { fs.closeSync(fd) }
}

function replacePrivateFile(file: string, bytes: Buffer): void {
  const temporary = path.join(path.dirname(file), `.validated-${crypto.randomUUID()}`)
  try {
    fs.writeFileSync(temporary, bytes, { flag: 'wx', mode: 0o600 })
    fs.renameSync(temporary, file)
  } finally { fs.rmSync(temporary, { force: true }) }
}

export async function executeAltlocSelection(
  dataRoot: string, workspaceId: string, request: AltlocRequest,
  runner = runDvbfixerArgs,
  authorize: (workspace: ReturnType<typeof loadWorkspace>) => void = () => {},
  signal?: AbortSignal,
): Promise<AltlocResponse> {
  if (!Value.Check(requestSchema, request)) throw new AltlocApiError(400, 'Invalid altloc request')
  if (request.action === 'apply' && !request.plan) throw new AltlocApiError(400, 'Apply requires a reviewed selection plan')
  const workspace = loadWorkspace(dataRoot, workspaceId)
  authorize(workspace)
  const artifact = workspace.artifacts.find(item => item.id === request.inputArtifactId && !item.hidden && item.kind === 'structure')
  if (!artifact) throw new AltlocApiError(404, 'Input structure artifact not found')
  const ext = path.extname(artifact.file).toLowerCase()
  if (!['.pdb', '.ent', '.cif', '.mmcif'].includes(ext)) throw new AltlocApiError(415, 'Altlocs requires PDB or CIF input')
  const source = resolveWorkspaceFile(dataRoot, workspaceId, artifact.file)
  const bytes = readFile(source, MAX_SOURCE)
  const sourceDigest = sha256(bytes)
  if (request.plan && request.plan.sourceSha256 !== sourceDigest) throw new AltlocApiError(409, 'Input changed; analyze again')
  const root = workspaceRoot(dataRoot, workspaceId)
  assertWorkspaceQuota(root, bytes.length * 4 + 1024 * 1024)
  const release = acquireWorkspaceRunLock(dataRoot, workspaceId)
  let operation: string | undefined
  let published = false
  try {
    const runs = path.join(root, 'runs')
    fs.mkdirSync(runs, { recursive: true })
    operation = fs.mkdtempSync(path.join(runs, 'altlocs_'))
    fs.chmodSync(operation, 0o700)
    const snapshot = path.join(operation, `source${ext}`)
    const output = path.join(operation, 'selected.pdb')
    const reportFile = path.join(operation, 'selection-report.json')
    fs.writeFileSync(snapshot, bytes, { mode: 0o600, flag: 'wx' })
    const args = [snapshot, '--policy', request.policy, '--fallback', request.fallback, '--report-json', reportFile]
    if (request.prefer) args.push('--prefer', request.prefer)
    if (request.action === 'analyze') args.push('--include-structure')
    else args.push('-o', output)
    if (request.plan) {
      const plan = path.join(operation, 'selection-plan.json')
      fs.writeFileSync(plan, JSON.stringify(request.plan), { mode: 0o600, flag: 'wx' })
      args.push('--selection-plan', plan)
    }
    const run = await runner('altlocs', args, operation, { timeoutMs: 60000, maxOutputBytes: 1024 * 1024, signal })
    if (run.code !== 0 || run.failure) throw new AltlocApiError(run.failure === 'overloaded' ? 503 : 422,
      run.failure === 'overloaded' ? 'Server busy' : `Altloc selection failed: ${run.stderr.slice(-3000) || run.stdout.slice(-3000)}`)
    const reportBytes = readFile(reportFile, 25 * 1024 * 1024)
    const reportValue: unknown = JSON.parse(reportBytes.toString('utf8'))
    if (!Value.Check(reportSchema, reportValue)) throw new AltlocApiError(500, 'Invalid altloc CLI report')
    const report = reportValue as AltlocReport
    if (report.sourceSha256 !== sourceDigest || report.policy !== request.policy || report.fallback !== request.fallback ||
        report.prefer !== (request.prefer || null) || JSON.stringify(Object.entries(report.overrides).sort()) !== JSON.stringify(Object.entries(request.plan?.selections || {}).sort())) {
      throw new AltlocApiError(500, 'Altloc CLI report/request mismatch')
    }
    if (report.normalizedPdb && sha256(report.normalizedPdb) !== report.normalizedSha256) throw new AltlocApiError(500, 'Altloc preview digest mismatch')
    const response: AltlocResponse = { report, outputArtifact: null }
    // Revalidate source identity and permissions even for an analysis response.
    const checkSource = (latest: ReturnType<typeof loadWorkspace>) => {
      authorize(latest)
      const current = latest.artifacts.find(item => item.id === artifact.id && !item.hidden)
      if (!current || current.file !== artifact.file || sha256(readFile(resolveWorkspaceFile(dataRoot, workspaceId, current.file), MAX_SOURCE)) !== sourceDigest) {
        throw new AltlocApiError(409, 'Input artifact changed; analyze again')
      }
    }
    checkSource(loadWorkspace(dataRoot, workspaceId))
    if (request.action === 'analyze') {
      if (report.outputWritten || fs.existsSync(output) || !report.normalizedPdb) throw new AltlocApiError(500, 'Analysis unexpectedly wrote output or omitted preview')
      return response
    }
    if (signal?.aborted) throw new AltlocApiError(499, 'Request cancelled')
    if (!report.outputWritten || report.unresolved.length || !report.outputSha256) throw new AltlocApiError(500, 'Altloc output is unresolved')
    if (Object.keys(request.plan!.selections).length !== report.residues.length ||
        report.residues.some(row => request.plan!.selections[row.key] !== row.selected)) {
      throw new AltlocApiError(422, 'Apply requires exactly one reviewed state for every alternate residue')
    }
    const outputBytes = readFile(output, MAX_SOURCE)
    if (!outputBytes.length || sha256(outputBytes) !== report.outputSha256) throw new AltlocApiError(500, 'Altloc output digest mismatch')
    // Detach validated outputs from any child-created hardlink before publication.
    replacePrivateFile(output, outputBytes)
    replacePrivateFile(reportFile, reportBytes)
    replacePrivateFile(path.join(operation, 'stdout.log'), Buffer.from(run.stdout))
    replacePrivateFile(path.join(operation, 'stderr.log'), Buffer.from(run.stderr))
    fs.unlinkSync(snapshot)
    const relative = path.relative(root, output).replace(/\\/g, '/')
    const result: WorkspaceArtifact = {
      id: crypto.randomUUID(), file: relative, name: `${path.basename(artifact.file, ext)}_selected.pdb`,
      kind: 'structure', parent: artifact.file, command: 'altlocs', artifactType: 'altloc-selection',
      altlocProvenance: { sourceArtifactId: artifact.id, sourceSha256: sourceDigest,
        outputSha256: report.outputSha256, report: path.relative(root, reportFile).replace(/\\/g, '/'),
        selections: request.plan!.selections, toolVersion: report.toolVersion },
    }
    updateWorkspace(dataRoot, workspaceId, latest => {
      checkSource(latest)
      assertWorkspaceQuota(root)
      latest.artifacts.push(result)
      return latest
    })
    published = true
    response.outputArtifact = { id: result.id, file: result.file, name: result.name, kind: 'structure' }
    response.reportFile = result.altlocProvenance!.report
    return response
  } finally {
    if (operation && !published) fs.rmSync(operation, { recursive: true, force: true })
    release()
  }
}

export function registerAltlocApi(server: ApiRouteHost, dataRoot: string,
  principal: (req: IncomingMessage) => string, legacyOwner: string): void {
  server.middlewares.use('/api/workspaces', async (req, res, next) => {
    const match = (req.url || '').match(/^\/([a-zA-Z0-9_-]+)\/altlocs$/)
    if (!match) return next()
    const send = (status: number, body: unknown) => {
      res.statusCode = status
      res.setHeader('Content-Type', 'application/json')
      res.setHeader('Cache-Control', 'no-store')
      res.end(JSON.stringify(body))
    }
    if (req.method !== 'POST') return send(405, { error: 'POST required' })
    if (!String(req.headers['content-type']).startsWith('application/json')) return send(415, { error: 'JSON required' })
    const controller = new AbortController()
    const abort = () => { if (!res.writableFinished) controller.abort() }
    res.once('close', abort)
    try {
      const authorize = (workspace: ReturnType<typeof loadWorkspace>) => assertWorkspaceAccess(workspace, principal(req), 'writer')
      authorize(loadWorkspace(dataRoot, match[1], legacyOwner))
      let value: unknown
      try { value = JSON.parse((await readRequestBody(req)).toString('utf8')) } catch (error) {
        if (errorStatus(error) === 413) throw error
        throw new AltlocApiError(400, 'Invalid JSON')
      }
      if (!Value.Check(requestSchema, value)) return send(400, { error: 'Invalid altloc request' })
      const result = await executeAltlocSelection(dataRoot, match[1], value as AltlocRequest, runDvbfixerArgs, authorize, controller.signal)
      return send(result.outputArtifact ? 201 : 200, result)
    } catch (error) {
      return send(error instanceof AltlocApiError ? error.statusCode : errorStatus(error),
        { error: error instanceof AltlocApiError ? error.message : 'Altloc operation failed' })
    } finally { res.removeListener('close', abort) }
  })
}
