import { afterEach, describe, expect, it } from 'vitest'
import fs from 'node:fs'
import os from 'node:os'
import path from 'node:path'
import { Readable } from 'node:stream'
import { Value } from '@sinclair/typebox/value'
import { executeNamingConversion, registerNamingApi } from './naming-api'
import {
  NAMING_OPENAPI_DOCUMENT,
  NAMING_REQUEST_EXAMPLE,
  NamingConversionRequestSchema,
  type NamingConversionRequest,
} from './naming-api-schema'
import { loadWorkspace, saveWorkspace, workspaceRoot, type WorkspaceManifest } from './workspace-api'
import {
  DEFAULT_MAX_UPLOAD_BYTES, DEFAULT_WORKSPACE_QUOTA_BYTES, initializeStorageQuota,
  workspaceLogicalBytes,
} from './storage-quota'

const temporaryDirectories: string[] = []

function temp(): string {
  const directory = fs.mkdtempSync(path.join(os.tmpdir(), 'dvbfixer-naming-api-test-'))
  temporaryDirectories.push(directory)
  return directory
}

function workspace(dataRoot: string, id = 'workspace-a'): WorkspaceManifest {
  const root = workspaceRoot(dataRoot, id)
  fs.mkdirSync(path.join(root, 'files'), { recursive: true })
  fs.writeFileSync(path.join(root, 'files', 'input.pdb'), 'ATOM      1  CA  ALA A   1      10.000  10.000  10.000  1.00 20.00           C  \nEND\n')
  return saveWorkspace(dataRoot, {
    version: 1, revision: 0, id, name: id,
    createdAt: new Date(0).toISOString(), updatedAt: new Date(0).toISOString(),
    primaryFile: null, secondaryFile: null,
    artifacts: [{ id: 'source', file: 'files/input.pdb', name: 'input.pdb', kind: 'structure' }],
    toolState: {},
  })
}

const request: NamingConversionRequest = {
  inputArtifactId: 'source',
  target: { forceField: 'amber', profile: 'gromacs' },
  variantOverrides: [],
  outputName: 'named.pdb',
  dryRun: false,
}

function successfulRunner(options: {
  dryRun?: boolean
  beforeWrite?: () => void
  outputSha256?: string
} = {}) {
  return async (_command: string, args: string[]) => {
    options.beforeWrite?.()
    const output = args[args.indexOf('-o') + 1]
    const reportFile = args[args.indexOf('--report-json') + 1]
    const outputBytes = Buffer.from('ATOM      1  CA  ALA A   1      10.000  10.000  10.000  1.00 20.00           C  \nEND\n')
    const digest = options.outputSha256 || await import('node:crypto').then(({ default: crypto }) =>
      crypto.createHash('sha256').update(outputBytes).digest('hex'))
    if (!options.dryRun) fs.writeFileSync(output, outputBytes)
    fs.writeFileSync(reportFile, JSON.stringify({
      schemaVersion: 1,
      operation: 'pdb-force-field-naming',
      status: 'success',
      tool: { name: 'dvbfixer', version: '0.8.5' },
      request: {
        targetForceField: 'amber', profile: 'gromacs', dryRun: Boolean(options.dryRun), variantOverrides: [],
      },
      output: { path: output, written: !options.dryRun, bytes: outputBytes.length, sha256: digest },
      result: {
        model: 1,
        summary: { changedAtoms: 0, changedResidues: 0, variantChanges: 0, atomNameChanges: 0 },
        changes: [], diagnostics: [],
      },
      error: null,
    }))
    return { code: 0, stdout: '', stderr: '' }
  }
}

async function apiRequest(
  dataRoot: string,
  method: string,
  url: string,
  body = '',
  contentType = 'application/json',
  runner?: Parameters<typeof executeNamingConversion>[3],
) {
  let middleware: ((req: any, res: any, next: () => void) => Promise<void>) | undefined
  const server = { middlewares: { use: (_route: string, handler: typeof middleware) => { middleware = handler } } }
  registerNamingApi(server as any, dataRoot, runner)
  const incoming = Readable.from([body]) as Readable & {
    method: string; url: string; headers: Record<string, string>
  }
  incoming.method = method
  incoming.url = url
  incoming.headers = { 'content-type': contentType }
  let responseBody = ''
  const headers: Record<string, string> = {}
  const response = {
    statusCode: 0,
    setHeader: (name: string, value: unknown) => { headers[name.toLowerCase()] = String(value) },
    end: (value?: unknown) => { responseBody = value === undefined ? '' : String(value) },
  }
  await middleware!(incoming as any, response as any, () => {})
  return {
    status: response.statusCode,
    headers,
    body: responseBody ? JSON.parse(responseBody) : null,
  }
}

afterEach(() => {
  initializeStorageQuota({
    maxUploadBytes: DEFAULT_MAX_UPLOAD_BYTES,
    workspaceQuotaBytes: DEFAULT_WORKSPACE_QUOTA_BYTES,
  })
  temporaryDirectories.splice(0).forEach(directory =>
    fs.rmSync(directory, { recursive: true, force: true }))
})

describe('naming conversion application boundary', () => {
  it('registers exactly one output with complete provenance', async () => {
    const dataRoot = temp()
    workspace(dataRoot)

    const result = await executeNamingConversion(dataRoot, 'workspace-a', request, successfulRunner())

    expect(result.statusCode).toBe(201)
    expect(result.response.outputArtifact?.name).toBe('named.pdb')
    const saved = loadWorkspace(dataRoot, 'workspace-a')
    expect(saved.artifacts).toHaveLength(2)
    const output = saved.artifacts[1]
    expect(output.command).toBe('atom-names')
    expect(output.namingProvenance).toMatchObject({
      operation: 'pdb-force-field-naming', sourceArtifactId: 'source',
      targetForceField: 'amber', profile: 'gromacs', dvbfixerVersion: '0.8.5',
    })
    expect(output.namingProvenance?.reportSha256).toMatch(/^[a-f0-9]{64}$/)
    expect(output.namingProvenance?.sourceSha256).toMatch(/^[a-f0-9]{64}$/)
    expect(output.namingProvenance?.outputSha256).toBe(result.response.output.sha256)
    expect(fs.readdirSync(path.join(workspaceRoot(dataRoot, 'workspace-a'), output.folder!)))
      .toEqual(['named.pdb'])
  })

  it('returns a dry-run report without changing the manifest or retaining output', async () => {
    const dataRoot = temp()
    const before = workspace(dataRoot)
    const dryRequest = { ...request, dryRun: true }

    const result = await executeNamingConversion(
      dataRoot, 'workspace-a', dryRequest, successfulRunner({ dryRun: true }),
    )

    expect(result.statusCode).toBe(200)
    expect(result.response.outputArtifact).toBeNull()
    const saved = loadWorkspace(dataRoot, 'workspace-a')
    expect(saved.revision).toBe(before.revision)
    expect(saved.artifacts).toHaveLength(1)
    expect(fs.readdirSync(path.join(workspaceRoot(dataRoot, 'workspace-a'), 'runs'))).toEqual([])
  })

  it('preserves a concurrent manifest update made while the command runs', async () => {
    const dataRoot = temp()
    workspace(dataRoot)
    const runner = successfulRunner({
      beforeWrite: () => {
        const latest = loadWorkspace(dataRoot, 'workspace-a')
        latest.toolState.concurrent = 'preserved'
        latest.artifacts.push({ id: 'concurrent', file: 'files/input.pdb', name: 'other', kind: 'structure' })
        saveWorkspace(dataRoot, latest)
      },
    })

    await executeNamingConversion(dataRoot, 'workspace-a', request, runner)

    const saved = loadWorkspace(dataRoot, 'workspace-a')
    expect(saved.toolState.concurrent).toBe('preserved')
    expect(saved.artifacts.map(artifact => artifact.id)).toContain('concurrent')
    expect(saved.artifacts.filter(artifact => artifact.command === 'atom-names')).toHaveLength(1)
  })

  it('does not register an output when its digest disagrees with the report', async () => {
    const dataRoot = temp()
    const before = workspace(dataRoot)

    await expect(executeNamingConversion(
      dataRoot, 'workspace-a', request, successfulRunner({ outputSha256: '0'.repeat(64) }),
    )).rejects.toMatchObject({ code: 'OUTPUT_DIGEST_MISMATCH' })

    const saved = loadWorkspace(dataRoot, 'workspace-a')
    expect(saved.revision).toBe(before.revision)
    expect(saved.artifacts).toHaveLength(1)
    expect(fs.readdirSync(path.join(workspaceRoot(dataRoot, 'workspace-a'), 'runs', '_failed'))).toHaveLength(1)
  })

  it('deletes operation data instead of retaining it when workspace quota is exhausted', async () => {
    const dataRoot = temp()
    const before = workspace(dataRoot)
    const root = workspaceRoot(dataRoot, 'workspace-a')
    const usedBefore = workspaceLogicalBytes(root)
    initializeStorageQuota({
      maxUploadBytes: DEFAULT_MAX_UPLOAD_BYTES,
      workspaceQuotaBytes: usedBefore + 256,
    })

    await expect(executeNamingConversion(dataRoot, 'workspace-a', request, successfulRunner()))
      .rejects.toThrow(/workspace quota exceeded/)

    expect(loadWorkspace(dataRoot, 'workspace-a')).toMatchObject({ revision: before.revision })
    expect(fs.existsSync(path.join(root, 'runs', '_failed'))).toBe(false)
    expect(workspaceLogicalBytes(root)).toBe(usedBefore)
  })

  it('does not publish when the source artifact changes during conversion', async () => {
    const dataRoot = temp()
    workspace(dataRoot)
    const runner = successfulRunner({
      beforeWrite: () => {
        const latest = loadWorkspace(dataRoot, 'workspace-a')
        latest.artifacts[0].file = 'files/replaced.pdb'
        fs.writeFileSync(path.join(workspaceRoot(dataRoot, 'workspace-a'), 'files', 'replaced.pdb'), 'END\n')
        saveWorkspace(dataRoot, latest)
      },
    })

    await expect(executeNamingConversion(dataRoot, 'workspace-a', request, runner))
      .rejects.toMatchObject({ statusCode: 409, code: 'SOURCE_ARTIFACT_CHANGED' })

    expect(loadWorkspace(dataRoot, 'workspace-a').artifacts).toHaveLength(1)
  })

  it('maps a stable CLI conversion error without registering an artifact', async () => {
    const dataRoot = temp()
    const before = workspace(dataRoot)
    const runner = async (_command: string, args: string[]) => {
      const reportFile = args[args.indexOf('--report-json') + 1]
      const output = args[args.indexOf('-o') + 1]
      fs.writeFileSync(reportFile, JSON.stringify({
        schemaVersion: 1, operation: 'pdb-force-field-naming', status: 'error',
        tool: { name: 'dvbfixer', version: '0.8.5' },
        request: { targetForceField: 'amber', profile: 'gromacs', dryRun: false, variantOverrides: [] },
        output: { path: output, written: false, bytes: null, sha256: null },
        result: null,
        error: { code: 'NAMING_COLLISION', category: 'conversion', message: 'collision', details: {} },
      }))
      return { code: 2, stdout: '', stderr: 'collision' }
    }

    await expect(executeNamingConversion(dataRoot, 'workspace-a', request, runner))
      .rejects.toMatchObject({ statusCode: 422, code: 'NAMING_COLLISION' })
    expect(loadWorkspace(dataRoot, 'workspace-a')).toMatchObject({ revision: before.revision })
  })

  it('rejects oversized and contradictory CLI reports', async () => {
    const dataRoot = temp()
    workspace(dataRoot)
    const previousLimit = process.env.DVBFIXER_NAMING_MAX_REPORT_BYTES
    process.env.DVBFIXER_NAMING_MAX_REPORT_BYTES = '64'
    try {
      const oversized = async (_command: string, args: string[]) => {
        fs.writeFileSync(args[args.indexOf('--report-json') + 1], 'x'.repeat(65))
        return { code: 0, stdout: '', stderr: '' }
      }
      await expect(executeNamingConversion(dataRoot, 'workspace-a', request, oversized))
        .rejects.toMatchObject({ code: 'INVALID_CLI_REPORT' })
    } finally {
      if (previousLimit === undefined) delete process.env.DVBFIXER_NAMING_MAX_REPORT_BYTES
      else process.env.DVBFIXER_NAMING_MAX_REPORT_BYTES = previousLimit
    }

    const contradictory = async (_command: string, args: string[]) => {
      const output = args[args.indexOf('-o') + 1]
      const outputBytes = Buffer.from('END\n')
      fs.writeFileSync(output, outputBytes)
      const digest = await import('node:crypto').then(({ default: crypto }) =>
        crypto.createHash('sha256').update(outputBytes).digest('hex'))
      fs.writeFileSync(args[args.indexOf('--report-json') + 1], JSON.stringify({
        schemaVersion: 1, operation: 'pdb-force-field-naming', status: 'success',
        tool: { name: 'dvbfixer', version: '0.8.5' },
        request: { targetForceField: 'amber', profile: 'gromacs', dryRun: false, variantOverrides: [] },
        output: { path: output, written: true, bytes: outputBytes.length, sha256: digest },
        result: {
          model: 1, summary: { changedAtoms: 0, changedResidues: 0, variantChanges: 0, atomNameChanges: 0 },
          changes: [], diagnostics: [],
        },
        error: { code: 'IMPOSSIBLE', category: 'conversion', message: 'contradiction', details: {} },
      }))
      return { code: 0, stdout: '', stderr: '' }
    }
    await expect(executeNamingConversion(dataRoot, 'workspace-a', request, contradictory))
      .rejects.toMatchObject({ code: 'INVALID_CLI_REPORT' })
  })

  it('reports queue saturation before attempting to parse a CLI report', async () => {
    const dataRoot = temp()
    workspace(dataRoot)
    await expect(executeNamingConversion(dataRoot, 'workspace-a', request, async () => ({
      code: -1, stdout: '', stderr: 'DVBfixer process queue is full',
      started: false, failure: 'overloaded' as const,
    }))).rejects.toMatchObject({ statusCode: 503, code: 'SERVER_BUSY' })
  })

  it('rejects source symlinks that escape the workspace', async () => {
    const dataRoot = temp()
    workspace(dataRoot)
    const source = path.join(workspaceRoot(dataRoot, 'workspace-a'), 'files', 'input.pdb')
    const outside = path.join(temp(), 'outside.pdb')
    fs.writeFileSync(outside, 'END\n')
    fs.unlinkSync(source)
    fs.symlinkSync(outside, source)

    await expect(executeNamingConversion(dataRoot, 'workspace-a', request, successfulRunner()))
      .rejects.toMatchObject({ statusCode: 404, code: 'ARTIFACT_NOT_FOUND' })
  })

  it('rejects an artifact ID belonging to another workspace', async () => {
    const dataRoot = temp()
    workspace(dataRoot, 'workspace-a')
    const other = workspace(dataRoot, 'workspace-b')
    other.artifacts[0].id = 'other-source'
    saveWorkspace(dataRoot, other)

    await expect(executeNamingConversion(
      dataRoot, 'workspace-a', { ...request, inputArtifactId: 'other-source' }, successfulRunner(),
    )).rejects.toMatchObject({ statusCode: 404, code: 'ARTIFACT_NOT_FOUND' })
  })
})

describe('naming API transport contract', () => {
  it('returns 507 and deletes generated data when publication exhausts quota', async () => {
    const dataRoot = temp()
    workspace(dataRoot)
    const root = workspaceRoot(dataRoot, 'workspace-a')
    const usedBefore = workspaceLogicalBytes(root)
    initializeStorageQuota({
      maxUploadBytes: DEFAULT_MAX_UPLOAD_BYTES,
      workspaceQuotaBytes: usedBefore + 256,
    })

    const response = await apiRequest(
      dataRoot,
      'POST',
      '/workspaces/workspace-a/naming-conversions',
      JSON.stringify(request),
      'application/json',
      successfulRunner(),
    )

    expect(response).toMatchObject({
      status: 507,
      body: { error: { code: 'WORKSPACE_QUOTA_EXCEEDED' } },
    })
    expect(workspaceLogicalBytes(root)).toBe(usedBefore)
    expect(fs.existsSync(path.join(root, 'runs', '_failed'))).toBe(false)
  })

  it('executes a valid request and returns the registered artifact', async () => {
    const dataRoot = temp()
    workspace(dataRoot)
    const response = await apiRequest(
      dataRoot,
      'POST',
      '/workspaces/workspace-a/naming-conversions',
      JSON.stringify(request),
      'application/json',
      successfulRunner(),
    )
    expect(response.status).toBe(201)
    expect(response.body).toMatchObject({ status: 'succeeded', sourceArtifactId: 'source' })
    expect(response.body.outputArtifact.file).toMatch(/^runs\/naming_/)
    expect(loadWorkspace(dataRoot, 'workspace-a').artifacts).toHaveLength(2)
  })

  it('publishes OpenAPI from the runtime request schema', async () => {
    const response = await apiRequest(temp(), 'GET', '/openapi.json')
    expect(response.status).toBe(200)
    expect(response.body).toEqual(JSON.parse(JSON.stringify(NAMING_OPENAPI_DOCUMENT)))
    expect(Value.Check(NamingConversionRequestSchema, NAMING_REQUEST_EXAMPLE)).toBe(true)
    expect(NAMING_OPENAPI_DOCUMENT.components.securitySchemes.bearerAuth).toEqual({
      type: 'http', scheme: 'bearer',
    })
    expect(NAMING_OPENAPI_DOCUMENT.paths['/api/v1/workspaces/{workspaceId}/naming-conversions']
      .post.security).toEqual([{ bearerAuth: [] }])
  })

  it('rejects unsupported media types and malformed or unknown request fields', async () => {
    const dataRoot = temp()
    const media = await apiRequest(dataRoot, 'POST', '/workspaces/workspace-a/naming-conversions', '{}', 'text/plain')
    expect(media.status).toBe(415)
    expect(media.body.error.code).toBe('UNSUPPORTED_MEDIA_TYPE')
    expect(media.headers['x-request-id']).toMatch(/^req_/)

    const malformed = await apiRequest(dataRoot, 'POST', '/workspaces/workspace-a/naming-conversions', '{')
    expect(malformed.status).toBe(400)
    expect(malformed.body.error.code).toBe('INVALID_JSON')

    const unknown = await apiRequest(dataRoot, 'POST', '/workspaces/workspace-a/naming-conversions', JSON.stringify({
      ...request, unexpected: true,
    }))
    expect(unknown.status).toBe(400)
    expect(unknown.body.error.code).toBe('INVALID_REQUEST')
  })
})
