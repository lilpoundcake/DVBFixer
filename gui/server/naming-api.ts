import type { IncomingMessage, ServerResponse } from 'node:http'
import type { ViteDevServer } from 'vite'
import crypto from 'node:crypto'
import fs from 'node:fs'
import path from 'node:path'
import { Value } from '@sinclair/typebox/value'
import { acquireWorkspaceRun } from './managed-jobs'
import {
  ApiErrorSchema,
  NAMING_OPENAPI_DOCUMENT,
  NamingCliReportSchema,
  NamingConversionRequestSchema,
  NamingConversionResponseSchema,
  type NamingCliReport,
  type NamingConversionRequest,
  type NamingConversionResponse,
} from './naming-api-schema'
import { runDvbfixerArgs } from './dvbfixer-runner'
import { errorStatus, readRequestBody } from './request-body'
import {
  loadWorkspace,
  resolveWorkspaceFile,
  saveWorkspace,
  workspaceRoot,
  writeJsonAtomic,
  type NamingConversionProvenance,
  type WorkspaceArtifact,
} from './workspace-api'

const MAX_SOURCE_BYTES = Number(process.env.DVBFIXER_NAMING_MAX_SOURCE_BYTES || 50 * 1024 * 1024)
const TIMEOUT_MS = Number(process.env.DVBFIXER_NAMING_TIMEOUT_MS || 60_000)
const MAX_OUTPUT_BYTES = 1024 * 1024

function maxReportBytes(): number {
  return Number(process.env.DVBFIXER_NAMING_MAX_REPORT_BYTES || 20 * 1024 * 1024)
}

type NamingRunner = typeof runDvbfixerArgs

class NamingApiError extends Error {
  readonly statusCode: number
  readonly code: string
  readonly details?: unknown

  constructor(
    statusCode: number,
    code: string,
    message: string,
    details?: unknown,
  ) {
    super(message)
    this.name = 'NamingApiError'
    this.statusCode = statusCode
    this.code = code
    this.details = details
  }
}

function sha256(content: Buffer): string {
  return crypto.createHash('sha256').update(content).digest('hex')
}

function schemaError(schema: Parameters<typeof Value.Errors>[0], value: unknown): string {
  const first = Value.Errors(schema, value).First()
  return first ? `${first.path || '/'} ${first.message}` : 'value does not match schema'
}

function parseRequest(value: unknown): NamingConversionRequest {
  if (!Value.Check(NamingConversionRequestSchema, value)) {
    throw new NamingApiError(400, 'INVALID_REQUEST', schemaError(NamingConversionRequestSchema, value))
  }
  const identities = new Set<string>()
  for (const override of value.variantOverrides || []) {
    const identity = JSON.stringify([override.chainId, override.residueNumber, override.insertionCode])
    if (identities.has(identity)) {
      throw new NamingApiError(400, 'INVALID_REQUEST', 'variantOverrides contains a duplicate residue identity')
    }
    identities.add(identity)
  }
  return value
}

function parseReport(file: string): { report: NamingCliReport; bytes: Buffer } {
  let bytes: Buffer
  let value: unknown
  try {
    const stat = fs.statSync(file)
    if (!stat.isFile() || stat.size > maxReportBytes()) {
      throw new NamingApiError(500, 'INVALID_CLI_REPORT', 'atom-names report is missing or exceeds the configured limit')
    }
    bytes = fs.readFileSync(file)
    value = JSON.parse(bytes.toString('utf8'))
  } catch (error) {
    if (error instanceof NamingApiError) throw error
    throw new NamingApiError(500, 'INVALID_CLI_REPORT', 'atom-names did not produce a readable JSON report')
  }
  if (!Value.Check(NamingCliReportSchema, value)) {
    throw new NamingApiError(500, 'INVALID_CLI_REPORT', schemaError(NamingCliReportSchema, value))
  }
  return { report: value, bytes }
}

function validateReportRequest(report: NamingCliReport, request: NamingConversionRequest): void {
  const overrides = request.variantOverrides || []
  if (report.request.targetForceField !== request.target.forceField ||
      report.request.profile !== request.target.profile ||
      report.request.dryRun !== (request.dryRun || false) ||
      JSON.stringify(report.request.variantOverrides) !== JSON.stringify(overrides)) {
    throw new NamingApiError(500, 'CLI_REPORT_MISMATCH', 'atom-names report does not match the API request')
  }
  if (report.status === 'success' && (!report.result || report.error !== null ||
      report.output.bytes === null || report.output.sha256 === null)) {
    throw new NamingApiError(500, 'INVALID_CLI_REPORT', 'atom-names success report is internally inconsistent')
  }
  if (report.status === 'error' && (report.result !== null || report.error === null || report.output.written)) {
    throw new NamingApiError(500, 'INVALID_CLI_REPORT', 'atom-names error report is internally inconsistent')
  }
}

function outputName(request: NamingConversionRequest, sourceFile: string): string {
  if (request.outputName) return request.outputName
  const stem = path.basename(sourceFile, path.extname(sourceFile))
  return `${stem}_${request.target.forceField}_gromacs.pdb`
}

function moveToFailureArea(operationDirectory: string, workspaceDirectory: string): void {
  if (!fs.existsSync(operationDirectory)) return
  const failed = path.join(workspaceDirectory, 'runs', '_failed')
  fs.mkdirSync(failed, { recursive: true })
  fs.renameSync(operationDirectory, path.join(failed, path.basename(operationDirectory)))
}

function removeSuccessHelpers(operationDirectory: string, output: string): void {
  for (const entry of fs.readdirSync(operationDirectory, { withFileTypes: true })) {
    const candidate = path.join(operationDirectory, entry.name)
    if (candidate === output) continue
    fs.rmSync(candidate, { recursive: entry.isDirectory(), force: true })
  }
}

export async function executeNamingConversion(
  dataRoot: string,
  workspaceId: string,
  request: NamingConversionRequest,
  runner: NamingRunner = runDvbfixerArgs,
): Promise<{ statusCode: 200 | 201; response: NamingConversionResponse }> {
  let initial
  try { initial = loadWorkspace(dataRoot, workspaceId) } catch {
    throw new NamingApiError(404, 'WORKSPACE_NOT_FOUND', 'workspace does not exist')
  }
  const sourceArtifact = initial.artifacts.find(artifact => artifact.id === request.inputArtifactId && !artifact.hidden)
  if (!sourceArtifact) throw new NamingApiError(404, 'ARTIFACT_NOT_FOUND', 'input artifact does not exist')
  if (sourceArtifact.kind !== 'structure') {
    throw new NamingApiError(415, 'UNSUPPORTED_SOURCE_FORMAT', 'input artifact is not a structure')
  }
  if (!['.pdb', '.ent'].includes(path.extname(sourceArtifact.file).toLowerCase())) {
    throw new NamingApiError(415, 'UNSUPPORTED_SOURCE_FORMAT', 'atom naming V1 accepts PDB artifacts only')
  }
  let source: string
  try { source = resolveWorkspaceFile(dataRoot, workspaceId, sourceArtifact.file) } catch {
    throw new NamingApiError(404, 'ARTIFACT_NOT_FOUND', 'input artifact file does not exist')
  }
  const stat = fs.statSync(source)
  if (!stat.isFile()) throw new NamingApiError(404, 'ARTIFACT_NOT_FOUND', 'input artifact is not a file')
  if (stat.size > MAX_SOURCE_BYTES) {
    throw new NamingApiError(413, 'SOURCE_TOO_LARGE', `input artifact exceeds ${MAX_SOURCE_BYTES} bytes`)
  }
  const sourceSha256 = sha256(fs.readFileSync(source))

  const operationId = crypto.randomUUID()
  const release = acquireWorkspaceRun(workspaceId, `naming-${operationId}`)
  if (!release) throw new NamingApiError(409, 'WORKSPACE_BUSY', 'another DVBFixer operation is running in this workspace')
  const root = workspaceRoot(dataRoot, workspaceId)
  const operationDirectory = path.join(root, 'runs', `naming_${operationId}`)
  const output = path.join(operationDirectory, outputName(request, sourceArtifact.file))
  const reportFile = path.join(operationDirectory, 'report.json')
  try {
    fs.mkdirSync(operationDirectory, { recursive: false })
    const args = [
      source, '-o', output,
      '--target-ff', request.target.forceField,
      '--profile', request.target.profile,
      '--report-json', reportFile,
    ]
    if (request.variantOverrides?.length) {
      const overridesFile = path.join(operationDirectory, 'variant-overrides.json')
      writeJsonAtomic(overridesFile, request.variantOverrides)
      args.push('--variant-overrides', overridesFile)
    }
    if (request.dryRun) args.push('--dry-run')

    const run = await runner('atom-names', args, operationDirectory, {
      timeoutMs: TIMEOUT_MS, maxOutputBytes: MAX_OUTPUT_BYTES, killGraceMs: 2_000,
    })
    fs.writeFileSync(path.join(operationDirectory, 'stdout.log'), run.stdout)
    fs.writeFileSync(path.join(operationDirectory, 'stderr.log'), run.stderr)
    const parsed = parseReport(reportFile)
    validateReportRequest(parsed.report, request)
    if (run.code !== 0 || parsed.report.status !== 'success' || !parsed.report.result) {
      const cliError = parsed.report.error
      const isConversionError = cliError?.category === 'conversion'
      throw new NamingApiError(
        isConversionError ? 422 : 500,
        isConversionError ? cliError.code : 'SUBPROCESS_FAILED',
        isConversionError ? cliError.message : 'atom-names failed',
        isConversionError ? cliError.details : undefined,
      )
    }
    const report = parsed.report
    const conversionResult = report.result
    if (!conversionResult) {
      throw new NamingApiError(500, 'INVALID_CLI_REPORT', 'successful atom-names report has no result')
    }
    if (report.output.path !== output) {
      throw new NamingApiError(500, 'CLI_REPORT_MISMATCH', 'atom-names reported an unexpected output path')
    }

    let outputArtifact: WorkspaceArtifact | null = null
    let manifestToPublish: ReturnType<typeof loadWorkspace> | null = null
    if (request.dryRun) {
      if (report.output.written || fs.existsSync(output)) {
        throw new NamingApiError(500, 'CLI_REPORT_MISMATCH', 'dry run unexpectedly published an output')
      }
    } else {
      if (!report.output.written || !fs.existsSync(output) || !fs.statSync(output).isFile()) {
        throw new NamingApiError(500, 'MISSING_OUTPUT', 'atom-names did not publish the expected output')
      }
      const outputBytes = fs.readFileSync(output)
      if (!outputBytes.length || report.output.bytes !== outputBytes.length || report.output.sha256 !== sha256(outputBytes)) {
        throw new NamingApiError(500, 'OUTPUT_DIGEST_MISMATCH', 'atom-names output does not match its report')
      }
      const relativeFolder = path.relative(root, operationDirectory).replace(/\\/g, '/')
      const relativeOutput = path.relative(root, output).replace(/\\/g, '/')
      const provenance: NamingConversionProvenance = {
        operation: 'pdb-force-field-naming', operationId, apiVersion: 1,
        reportSchemaVersion: 1, sourceArtifactId: sourceArtifact.id,
        sourceFile: sourceArtifact.file, sourceSha256, command: 'atom-names',
        targetForceField: request.target.forceField, profile: request.target.profile,
        variantOverrides: request.variantOverrides || [], reportSha256: sha256(parsed.bytes),
        outputSha256: report.output.sha256, dvbfixerVersion: report.tool.version,
      }
      outputArtifact = {
        id: crypto.randomUUID(), file: relativeOutput, name: path.basename(output),
        kind: 'structure', folder: relativeFolder, parent: sourceArtifact.file,
        command: 'atom-names', artifactType: 'naming-conversion', namingProvenance: provenance,
      }
      const latest = loadWorkspace(dataRoot, workspaceId)
      const currentSource = latest.artifacts.find(artifact => artifact.id === sourceArtifact.id)
      if (!currentSource || currentSource.hidden || currentSource.kind !== 'structure' ||
          currentSource.file !== sourceArtifact.file || sha256(fs.readFileSync(source)) !== sourceSha256) {
        throw new NamingApiError(409, 'SOURCE_ARTIFACT_CHANGED', 'input artifact changed during conversion')
      }
      if (latest.artifacts.some(artifact => artifact.file === relativeOutput)) {
        throw new NamingApiError(409, 'OUTPUT_CONFLICT', 'output artifact path already exists')
      }
      latest.artifacts.push(outputArtifact)
      manifestToPublish = latest
    }

    const response: NamingConversionResponse = {
      operationId, status: 'succeeded', sourceArtifactId: sourceArtifact.id,
      outputArtifact: outputArtifact && {
        id: outputArtifact.id, file: outputArtifact.file, name: outputArtifact.name, kind: 'structure',
      },
      output: { bytes: report.output.bytes!, sha256: report.output.sha256! },
      result: conversionResult,
    }
    if (!Value.Check(NamingConversionResponseSchema, response)) {
      throw new NamingApiError(500, 'INVALID_RESPONSE', schemaError(NamingConversionResponseSchema, response))
    }
    if (manifestToPublish) {
      removeSuccessHelpers(operationDirectory, output)
      saveWorkspace(dataRoot, manifestToPublish)
    }
    if (request.dryRun) fs.rmSync(operationDirectory, { recursive: true, force: true })
    return { statusCode: request.dryRun ? 200 : 201, response }
  } catch (error) {
    try { moveToFailureArea(operationDirectory, root) } catch { /* preserve the original failure */ }
    throw error
  } finally {
    release()
  }
}

function sendJson(res: ServerResponse, status: number, body: unknown): void {
  res.statusCode = status
  res.setHeader('Content-Type', 'application/json')
  res.setHeader('Cache-Control', 'no-store')
  res.end(JSON.stringify(body))
}

function sendError(res: ServerResponse, error: unknown, requestId: string): void {
  const apiError = error instanceof NamingApiError ? error : null
  const status = apiError?.statusCode || errorStatus(error)
  const body = {
    error: {
      code: apiError?.code || (status === 413 ? 'PAYLOAD_TOO_LARGE' : 'INTERNAL_ERROR'),
      message: apiError?.message || 'unexpected server error',
      ...(apiError?.details === undefined ? {} : { details: apiError.details }),
      requestId,
    },
  }
  if (!Value.Check(ApiErrorSchema, body)) {
    return sendJson(res, 500, { error: { code: 'INTERNAL_ERROR', message: 'invalid error response', requestId } })
  }
  sendJson(res, status, body)
}

export function registerNamingApi(
  server: ViteDevServer,
  dataRoot: string,
  runner: NamingRunner = runDvbfixerArgs,
): void {
  server.middlewares.use('/api/v1', async (req: IncomingMessage, res: ServerResponse, next: () => void) => {
    const route = (req.url || '').split('?')[0]
    if (req.method === 'GET' && route === '/openapi.json') {
      return sendJson(res, 200, NAMING_OPENAPI_DOCUMENT)
    }
    const match = route.match(/^\/workspaces\/([^/]+)\/naming-conversions$/)
    if (!match) return next()
    const requestId = `req_${crypto.randomUUID()}`
    res.setHeader('X-Request-Id', requestId)
    if (req.method !== 'POST') return sendError(res, new NamingApiError(405, 'METHOD_NOT_ALLOWED', 'method not allowed'), requestId)
    if (!String(req.headers['content-type'] || '').toLowerCase().startsWith('application/json')) {
      return sendError(res, new NamingApiError(415, 'UNSUPPORTED_MEDIA_TYPE', 'Content-Type must be application/json'), requestId)
    }
    try {
      let workspaceId: string
      try { workspaceId = decodeURIComponent(match[1]) } catch {
        throw new NamingApiError(400, 'INVALID_WORKSPACE_ID', 'workspace ID is not valid URL encoding')
      }
      if (!/^[a-zA-Z0-9_-]+$/.test(workspaceId)) {
        throw new NamingApiError(400, 'INVALID_WORKSPACE_ID', 'workspace ID contains unsupported characters')
      }
      let value: unknown
      try { value = JSON.parse((await readRequestBody(req)).toString('utf8')) } catch (error) {
        if (errorStatus(error) === 413) throw error
        throw new NamingApiError(400, 'INVALID_JSON', 'request body is not valid JSON')
      }
      const request = parseRequest(value)
      const result = await executeNamingConversion(dataRoot, workspaceId, request, runner)
      return sendJson(res, result.statusCode, result.response)
    } catch (error) {
      return sendError(res, error, requestId)
    }
  })
}
