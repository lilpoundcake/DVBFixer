/**
 * Host-neutral registration for the GUI and versioned HTTP APIs.
 *
 *   POST   /api/v1/workspaces/:id/jobs — start a managed DVBFixer workflow.
 *   GET    /api/mutations           — list mutations
 *   POST   /api/mutations           — create
 *   PUT    /api/mutations/:id       — update
 *   DELETE /api/mutations/:id       — delete
 *
 * Env vars:
 *   DVBFIXER_EXECUTABLE  default 'dvbfixer' — CLI executable
 *   DVBFIXER_ARGS        optional JSON string array prepended to CLI arguments
 *   DATABASE_URL   postgres connection string
 */

import type { IncomingMessage, ServerResponse } from 'node:http'
import path from 'node:path'
import fs from 'node:fs'
import crypto from 'node:crypto'
import { COMMANDS } from './dvbfixer-spec'
import {
  registerManagedJobApi, resetManagedJobAdmission, shutdownManagedJobs,
} from './managed-jobs'
import { registerHomologyApi } from './homology-api'
import { registerNamingApi } from './naming-api'
import {
  initializeDvbfixerProcessAdmission, resetDvbfixerProcessAdmission, runDvbfixerArgs,
  shutdownDvbfixerProcesses, getDvbfixerProcessAdmissionSnapshot,
} from './dvbfixer-runner'
import {
  assertWorkspaceAccess, loadWorkspace, registerWorkspaceApi, resolveWorkspaceFile, updateWorkspace,
  migrateWorkspaceOwnership, workspaceRoot, writeJsonAtomic,
} from './workspace-api'
import { errorStatus, readRequestBody } from './request-body'
import type { ApiRouteHost } from './http-types'
import {
  createAuthMiddleware, getApiPrincipal, parseAuthConfig, resolveLegacyWorkspaceOwner,
  type AuthConfig,
} from './auth'
import { createCorsMiddleware, parseCorsAllowedOrigins } from './cors'
import { createRateLimitMiddleware, parseRateLimitConfig, type RateLimitConfig } from './rate-limit'
import { assertWorkspaceQuota, configureStorageQuota, type StorageQuotaOptions } from './storage-quota'
import {
  createRequestObservabilityMiddleware, ensureRequestId, parseAccessLogConfig, type AccessLogConfig,
} from './request-observability'
import { ApiMetrics, parseMetricsConfig, type MetricsConfig } from './metrics'
import { createAuditSink, parseAuditLogConfig, type AuditLogConfig } from './audit-log'
export { runDvbfixer } from './dvbfixer-runner'
export { buildArgs } from './command-args'

// Defer the pg import so the plugin loads even if pg is missing or DB is unset.
type PgClient = {
  query: (sql: string, params?: any[]) => Promise<{ rows: any[] }>
  end?: () => Promise<void>
}
let pgPool: PgClient | null = null
let pgInitialization: Promise<PgClient | null> | null = null

export async function closeApiResources(): Promise<void> {
  const pool = pgPool
  pgPool = null
  if (pool?.end) await pool.end()
}

export async function shutdownApiWork(): Promise<void> {
  shutdownManagedJobs()
  await shutdownDvbfixerProcesses()
}

export function resetApiShutdown(): void {
  resetManagedJobAdmission()
  resetDvbfixerProcessAdmission()
}

// Configurable path for the mutations backup. Development defaults to the
// git-tracked project-root file; hardened deployments keep it on bounded state.
let mutationsBackupFile: string | null = null
let mutationsBackupWrites: Promise<void> = Promise.resolve()

function mutationsBackupPath(): string {
  return mutationsBackupFile ?? path.join(process.cwd(), 'mutations.json')
}

/**
 * Snapshot the full mutations table to the configured JSON backup.
 * Called after every successful CRUD write so the git-tracked backup
 * stays in sync with the live DB. Failures are warned but not fatal —
 * a dump miss is better than crashing the API.
 */
async function dumpMutationsToBackup(pg: PgClient) {
  const write = mutationsBackupWrites.then(async () => {
    const { rows } = await pg.query(
      'SELECT id, chain, mutation_name, mutations, igg_subclass, properties, display_order FROM mutations ORDER BY display_order ASC, id ASC'
    )
    writeJsonAtomic(mutationsBackupPath(), rows)
  }).catch(err => {
    console.warn('[api] failed to dump mutations backup:', err)
  })
  mutationsBackupWrites = write
  await write
}

/**
 * If the mutations table is empty and a configured backup file exists,
 * seed the table from it. Runs once on first connection — subsequent
 * runs find the table non-empty and skip. This means a fresh clone with
 * the committed `mutations.json` gets the team's mutation library
 * automatically.
 */
async function seedMutationsFromBackup(pg: PgClient) {
  const filePath = mutationsBackupPath()
  if (!fs.existsSync(filePath)) return
  try {
    const { rows: countRows } = await pg.query('SELECT COUNT(*)::int AS n FROM mutations')
    if ((countRows[0]?.n ?? 0) > 0) return
    const raw = fs.readFileSync(filePath, 'utf-8').trim()
    if (!raw) return
    const data = JSON.parse(raw) as Array<{
      id?: number
      chain?: string
      mutation_name?: string
      mutations?: string
      igg_subclass?: string
      properties?: string
      display_order?: number
    }>
    if (!Array.isArray(data) || data.length === 0) return
    for (let i = 0; i < data.length; i++) {
      const row = data[i]
      // Preserve the order in the JSON file: row at index i gets
      // display_order i+1 if the JSON doesn't specify one explicitly.
      const order = typeof row.display_order === 'number' ? row.display_order : i + 1
      if (row.id !== undefined) {
        await pg.query(
          'INSERT INTO mutations (id, chain, mutation_name, mutations, igg_subclass, properties, display_order) VALUES ($1, $2, $3, $4, $5, $6, $7)',
          [row.id, row.chain ?? '', row.mutation_name ?? '', row.mutations ?? '', row.igg_subclass ?? '', row.properties ?? '', order]
        )
      } else {
        await pg.query(
          'INSERT INTO mutations (chain, mutation_name, mutations, igg_subclass, properties, display_order) VALUES ($1, $2, $3, $4, $5, $6)',
          [row.chain ?? '', row.mutation_name ?? '', row.mutations ?? '', row.igg_subclass ?? '', row.properties ?? '', order]
        )
      }
    }
    // Bump the auto-id sequence above any explicitly-inserted ids so new
    // rows don't collide.
    await pg.query("SELECT setval('mutations_id_seq', COALESCE((SELECT MAX(id) FROM mutations), 1))")
    console.log(`[api] seeded ${data.length} mutations from mutations.json`)
  } catch (err) {
    console.warn('[api] failed to seed mutations from backup:', err)
  }
}

export async function getPg(): Promise<PgClient | null> {
  if (pgPool) return pgPool
  const url = process.env.DATABASE_URL
  if (!url) return null
  if (!pgInitialization) pgInitialization = initializePg(url)
  try { return await pgInitialization } finally { pgInitialization = null }
}

async function initializePg(url: string): Promise<PgClient | null> {
  let candidate: PgClient | null = null
  try {
    const pg = await import('pg')
    const { Pool } = pg.default ?? pg
    const pool = new Pool({ connectionString: url })
    candidate = pool as unknown as PgClient
    // Auto-create schema. `igg_subclass` is the per-row antibody
    // subclass tag (e.g. IgG1 / IgG2 / IgG4) — empty string by default.
    await pool.query(`
      CREATE TABLE IF NOT EXISTS mutations (
        id SERIAL PRIMARY KEY,
        chain TEXT NOT NULL,
        mutation_name TEXT NOT NULL,
        mutations TEXT NOT NULL,
        igg_subclass TEXT NOT NULL DEFAULT '',
        properties TEXT NOT NULL DEFAULT '',
        display_order INTEGER NOT NULL DEFAULT 0,
        created_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
      )
    `)
    // Migrate older deployments whose table predates igg_subclass / display_order.
    await pool.query(
      `ALTER TABLE mutations ADD COLUMN IF NOT EXISTS igg_subclass TEXT NOT NULL DEFAULT ''`
    )
    await pool.query(
      `ALTER TABLE mutations ADD COLUMN IF NOT EXISTS display_order INTEGER NOT NULL DEFAULT 0`
    )
    await pool.query(
      `ALTER TABLE mutations ADD COLUMN IF NOT EXISTS properties TEXT NOT NULL DEFAULT ''`
    )
    // Seed display_order from id where it's still 0 (e.g. just-added column).
    // Preserves the existing visible order.
    await pool.query(
      `UPDATE mutations SET display_order = id WHERE display_order = 0`
    )
    // Seed from the git-tracked backup file if the table is empty.
    await seedMutationsFromBackup(candidate)
    pgPool = candidate
    return pgPool
  } catch (err) {
    if (candidate?.end) await candidate.end().catch(() => {})
    console.error('[api] postgres init failed:', err)
    return null
  }
}

/* ────────────────────────────────────────────────────────────────────────
 * SSE helpers (used by /api/antibody-engineer/run and any future
 * long-running orchestrator that wants real-time progress).
 * ──────────────────────────────────────────────────────────────────────── */

/** Set SSE headers and flush so the browser begins reading. Call once
 *  before the first sseSend(). */
export function writeSSEHeaders(res: ServerResponse): void {
  res.writeHead(200, {
    'Content-Type': 'text/event-stream; charset=utf-8',
    'Cache-Control': 'no-cache, no-transform',
    Connection: 'keep-alive',
    // Defeat upstream proxy buffering (nginx, etc.) so events arrive promptly.
    'X-Accel-Buffering': 'no',
  })
  // res.flushHeaders exists on Node's ServerResponse — try-cast for compat.
  if (typeof (res as any).flushHeaders === 'function') (res as any).flushHeaders()
}

/** Emit one SSE message. Single-channel: client switches on payload.status. */
export function sseSend(res: ServerResponse, payload: unknown): void {
  res.write(`data: ${JSON.stringify(payload)}\n\n`)
}

/** Read the entire request body as a string. */
function readBody(req: IncomingMessage): Promise<string> {
  return readRequestBody(req).then(body => body.toString('utf8'))
}

function sendJson(res: ServerResponse, status: number, body: any) {
  res.statusCode = status
  res.setHeader('Content-Type', 'application/json')
  res.setHeader('Cache-Control', 'no-store')
  res.end(JSON.stringify(body))
}

export function readServiceVersion(projectRoot: string): string {
  try {
    const value = JSON.parse(fs.readFileSync(path.join(projectRoot, 'package.json'), 'utf8')).version
    if (typeof value === 'string' && value.length > 0) return value
  } catch (error) {
    throw new Error('GUI package metadata must contain a service version', { cause: error })
  }
  throw new Error('GUI package metadata must contain a service version')
}

export interface ApiRouteOptions {
  projectRoot: string
  dataRoot?: string
  mutationsBackupFile?: string
  authConfig?: AuthConfig
  legacyWorkspaceOwner?: string
  corsAllowedOrigins?: readonly string[]
  storageQuota?: StorageQuotaOptions
  maxConcurrentProcesses?: number
  maxQueuedProcesses?: number
  rateLimit?: RateLimitConfig
  accessLog?: AccessLogConfig
  metrics?: MetricsConfig
  auditLog?: AuditLogConfig
}

export function registerApiRoutes(server: ApiRouteHost, options: ApiRouteOptions): void {
      const structuresDir = path.resolve(options.dataRoot || process.env.DVBFIXER_GUI_DATA_DIR || path.join(options.projectRoot, 'structures'))
      const serviceVersion = readServiceVersion(options.projectRoot)
      const authConfig = options.authConfig || parseAuthConfig(process.env)
      const legacyWorkspaceOwner = options.legacyWorkspaceOwner || resolveLegacyWorkspaceOwner(authConfig)
      const corsAllowedOrigins = options.corsAllowedOrigins || parseCorsAllowedOrigins(process.env)
      const rateLimit = options.rateLimit || parseRateLimitConfig(process.env)
      const accessLog = options.accessLog || parseAccessLogConfig(process.env)
      const metricsConfig = options.metrics || parseMetricsConfig(process.env)
      const auditLog = options.auditLog || parseAuditLogConfig(process.env)
      const metrics = new ApiMetrics()
      configureStorageQuota(options.storageQuota)
      initializeDvbfixerProcessAdmission(
        options.maxConcurrentProcesses === undefined
          ? process.env.DVBFIXER_MAX_CONCURRENT_PROCESSES
          : String(options.maxConcurrentProcesses),
        options.maxQueuedProcesses === undefined
          ? process.env.DVBFIXER_MAX_QUEUED_PROCESSES
          : String(options.maxQueuedProcesses),
      )
      fs.mkdirSync(structuresDir, { recursive: true })
      mutationsBackupFile = path.resolve(options.projectRoot, options.mutationsBackupFile || 'mutations.json')
      const rateLimitMiddleware = createRateLimitMiddleware(rateLimit)
      server.middlewares.use('/api', createRequestObservabilityMiddleware(
        accessLog, undefined, undefined, metricsConfig.enabled ? metrics : undefined,
        createAuditSink(structuresDir, auditLog),
      ))
      server.middlewares.use('/api', createCorsMiddleware(corsAllowedOrigins, rateLimitMiddleware))
      server.middlewares.use('/api', rateLimitMiddleware)
      server.middlewares.use('/api', createAuthMiddleware(authConfig))
      if (metricsConfig.enabled) server.middlewares.use('/api/metrics', (req, res) => {
        if (req.method !== 'GET' && req.method !== 'HEAD') {
          res.setHeader('Allow', 'GET, HEAD')
          return sendJson(res, 405, { error: 'method not allowed' })
        }
        const body = Buffer.from(metrics.render(getDvbfixerProcessAdmissionSnapshot()))
        res.statusCode = 200
        res.setHeader('Content-Type', 'text/plain; version=0.0.4; charset=utf-8')
        res.setHeader('Cache-Control', 'no-store')
        res.setHeader('X-Content-Type-Options', 'nosniff')
        res.setHeader('Content-Length', String(body.length))
        res.end(req.method === 'HEAD' ? undefined : body)
      })
      migrateWorkspaceOwnership(
        structuresDir,
        legacyWorkspaceOwner,
        authConfig.enabled,
      )
      server.middlewares.use('/api/health', (req, res, next) => {
        if (req.method !== 'GET') return next()
        let storageWritable = false
        try { fs.accessSync(structuresDir, fs.constants.R_OK | fs.constants.W_OK); storageWritable = true } catch { /* reported below */ }
        return sendJson(res, storageWritable ? 200 : 503, {
          status: storageWritable ? 'ready' : 'degraded', version: serviceVersion, storageWritable,
          databaseConfigured: Boolean(process.env.DATABASE_URL), authenticationRequired: authConfig.enabled,
        })
      })
      server.middlewares.use('/api/session', (req, res, next) => {
        if (req.method !== 'GET') return next()
        return sendJson(res, 200, { principal: getApiPrincipal(req) })
      })
      const principalId = (request: IncomingMessage) => getApiPrincipal(request).id
      registerWorkspaceApi(server, structuresDir, principalId, legacyWorkspaceOwner, !authConfig.enabled)
      registerHomologyApi(server, structuresDir, principalId, legacyWorkspaceOwner)
      registerManagedJobApi(server, structuresDir, principalId, legacyWorkspaceOwner, serviceVersion)
      registerNamingApi(
        server, structuresDir, runDvbfixerArgs, principalId, legacyWorkspaceOwner, serviceVersion,
      )

      // ── Mutations CRUD ────────────────────────────────────────────────
      server.middlewares.use('/api/mutations', async (req, res, next) => {
        const pg = await getPg()
        if (!pg) {
          return sendJson(res, 503, { error: 'DATABASE_URL not configured' })
        }
        try {
          const url = req.url || ''
          // Strip query string and split. /api/mutations/123 → ['', '123']
          const pathOnly = url.split('?')[0]
          const idMatch = pathOnly.match(/^\/(\d+)$/)
          const id = idMatch ? parseInt(idMatch[1], 10) : null

          // PATCH /api/mutations/reorder — atomic reorder via full id list.
          // Frontend computes the new sequence via @dnd-kit's arrayMove and
          // POSTs the entire ordered id array; backend rewrites every row's
          // display_order in one transaction.
          if (pathOnly === '/reorder' && (req.method === 'PATCH' || req.method === 'POST')) {
            const body = JSON.parse(await readBody(req) || '{}') as { ids?: number[] }
            if (!Array.isArray(body.ids)) return sendJson(res, 400, { error: 'ids[] required' })
            // Single UPDATE ... FROM (VALUES …) is the simplest transactional
            // form. Postgres-specific; pg client allows multi-statement.
            // Using a CTE keeps it readable.
            const values: any[] = []
            const tuples: string[] = body.ids.map((rowId, i) => {
              values.push(rowId, i + 1)
              return `($${values.length - 1}::int, $${values.length}::int)`
            })
            if (tuples.length === 0) return sendJson(res, 200, { ok: true })
            await pg.query(
              `UPDATE mutations SET display_order = v.new_order
               FROM (VALUES ${tuples.join(', ')}) AS v(id, new_order)
               WHERE mutations.id = v.id`,
              values
            )
            await dumpMutationsToBackup(pg)
            return sendJson(res, 200, { ok: true })
          }

          if (req.method === 'GET' && !id) {
            const { rows } = await pg.query(
              'SELECT id, chain, mutation_name, mutations, igg_subclass, properties, display_order FROM mutations ORDER BY display_order ASC, id ASC'
            )
            return sendJson(res, 200, rows)
          }

          if (req.method === 'POST' && !id) {
            const body = JSON.parse(await readBody(req) || '{}')
            const { chain = '', mutation_name = '', mutations = '', igg_subclass = '', properties = '' } = body
            // New row goes to the bottom of the list: display_order = max+1.
            const { rows } = await pg.query(
              `INSERT INTO mutations (chain, mutation_name, mutations, igg_subclass, properties, display_order)
               VALUES ($1, $2, $3, $4, $5, COALESCE((SELECT MAX(display_order) + 1 FROM mutations), 1))
               RETURNING id, chain, mutation_name, mutations, igg_subclass, properties, display_order`,
              [chain, mutation_name, mutations, igg_subclass, properties]
            )
            await dumpMutationsToBackup(pg)
            return sendJson(res, 201, rows[0])
          }

          if (req.method === 'PUT' && id) {
            const body = JSON.parse(await readBody(req) || '{}')
            const { chain, mutation_name, mutations, igg_subclass, properties } = body
            const { rows } = await pg.query(
              `UPDATE mutations SET
                 chain = COALESCE($2, chain),
                 mutation_name = COALESCE($3, mutation_name),
                 mutations = COALESCE($4, mutations),
                 igg_subclass = COALESCE($5, igg_subclass),
                 properties = COALESCE($6, properties)
               WHERE id = $1
               RETURNING id, chain, mutation_name, mutations, igg_subclass, properties, display_order`,
              [id, chain, mutation_name, mutations, igg_subclass, properties]
            )
            await dumpMutationsToBackup(pg)
            return sendJson(res, 200, rows[0])
          }

          if (req.method === 'DELETE' && id) {
            await pg.query('DELETE FROM mutations WHERE id = $1', [id])
            await dumpMutationsToBackup(pg)
            return sendJson(res, 204, {})
          }

          return next()
        } catch (err: any) {
          sendJson(res, errorStatus(err), { error: err.message ?? String(err) })
        }
      })

      // ── Spec exposure (so frontend doesn't import server/) ────────────
      server.middlewares.use('/api/dvbfixer-spec', async (req, res, next) => {
        if (req.method !== 'GET') return next()
        sendJson(res, 200, COMMANDS)
      })

      // ── Antibody Engineer pipeline (SSE) ─────────────────────────────
      // POST /api/antibody-engineer/run streams progress for the
      // multi-step DVBFixer pipeline that applies one or more selected
      // Mutations DB rows (with equivalent-chain fan-out) to a structure.
      server.middlewares.use('/api/antibody-engineer/run', async (req, res, next) => {
        if (req.method !== 'POST') return next()
        try {
          const bodyText = await readBody(req)
          const body = JSON.parse(bodyText || '{}') as {
            workspaceId?: string
            inputFile?: string
            mutationIds?: number[]
            equivalentChainsMap?: Record<string, string[]>
            /** Per-mutation-id override of target chains; bypasses
             *  equivalent-chains expansion for that row. Used by the AE
             *  panel when a Mutations DB row has an empty `chain`
             *  field. */
            manualChainsByMutationId?: Record<number, string[]>
            hasGlycan?: boolean
            scheme?: 'EU' | 'Kabat'
          }
          if (!body.workspaceId) return sendJson(res, 400, { error: 'workspaceId required' })
          if (!body.inputFile) return sendJson(res, 400, { error: 'inputFile required' })
          if (!Array.isArray(body.mutationIds) || body.mutationIds.length === 0) {
            return sendJson(res, 400, { error: 'mutationIds required (non-empty array)' })
          }
          if (typeof body.hasGlycan !== 'boolean') return sendJson(res, 400, { error: 'hasGlycan required' })
          if (body.scheme !== 'EU' && body.scheme !== 'Kabat') return sendJson(res, 400, { error: 'scheme must be EU or Kabat' })

          assertWorkspaceAccess(
            loadWorkspace(structuresDir, body.workspaceId, legacyWorkspaceOwner),
            principalId(req),
            'writer',
          )

          const activeRoot = workspaceRoot(structuresDir, body.workspaceId)
          resolveWorkspaceFile(structuresDir, body.workspaceId, body.inputFile)

          let aborted = false
          const controller = new AbortController()
          const abortPipeline = () => {
            if (!res.writableFinished) {
              aborted = true
              controller.abort()
            }
          }
          res.on('close', abortPipeline)
          if (res.destroyed || res.closed) abortPipeline()

          // Dynamic imports so SSE-specific deps stay out of cold-path code.
          const { runEngineerPipeline, engineerChecksum, findCachedEntry } =
            await import('./antibody-pipeline')
          if (aborted) return

          // Compute checksum + lookup cache. Cache hit = no pipeline run.
          const checksum = engineerChecksum({
            inputFile: body.inputFile,
            mutationIds: body.mutationIds,
            hasGlycan: body.hasGlycan,
            scheme: body.scheme,
          })

          assertWorkspaceQuota(activeRoot)
          writeSSEHeaders(res)

          const workspaceBeforeRun = loadWorkspace(structuresDir, body.workspaceId, legacyWorkspaceOwner)
          const inputArtifact = workspaceBeforeRun.artifacts.find(artifact => artifact.file === body.inputFile)
          const cached = findCachedEntry(activeRoot, workspaceBeforeRun.artifacts, body.inputFile, checksum)
          if (cached) {
            sseSend(res, { step: 0, total: 0, name: 'cached', status: 'done', outputFile: cached.file })
            sseSend(res, { step: 0, total: 0, status: 'complete', outputFile: cached.file, entry: cached })
            res.end()
            return
          }

          // Look up the selected mutation rows from postgres. The Mutations
          // DB is the source of truth for what to mutate; the engineer tool
          // only references them by id.
          const pg = await getPg()
          if (!pg) {
            sseSend(res, { step: 0, total: 0, name: 'validate', status: 'error', stderr: 'DATABASE_URL not configured — mutations table unavailable.' })
            res.end()
            return
          }
          const { rows } = await pg.query(
            'SELECT id, chain, mutation_name, mutations FROM mutations WHERE id = ANY($1::int[]) ORDER BY id ASC',
            [body.mutationIds]
          )
          if (rows.length !== body.mutationIds.length) {
            const found = new Set(rows.map((r: any) => r.id))
            const missing = body.mutationIds.filter(id => !found.has(id))
            sseSend(res, { step: 0, total: 0, name: 'validate', status: 'error', stderr: `Mutation row(s) not found: ${missing.join(', ')}` })
            res.end()
            return
          }

          const generated = await runEngineerPipeline({
            structuresDir: activeRoot,
            inputFile: body.inputFile,
            mutationRows: rows as any,
            mutationIds: body.mutationIds,
            equivalentChainsMap: body.equivalentChainsMap ?? {},
            manualChainsByMutationId: body.manualChainsByMutationId ?? {},
            hasGlycan: body.hasGlycan,
            scheme: body.scheme,
            checksum,
            inputMetadata: inputArtifact,
            onEvent: (e) => sseSend(res, e),
            isAborted: () => aborted,
            signal: controller.signal,
            assertStorage: () => assertWorkspaceQuota(activeRoot),
          })

          if (generated.length) {
            updateWorkspace(structuresDir, body.workspaceId, workspace => {
              assertWorkspaceAccess(workspace, principalId(req), 'writer')
              for (const entry of generated.filter(item => item.file && !workspace.artifacts.some(existing => existing.file === item.file))) {
                workspace.artifacts.push({ id: crypto.randomUUID(), file: entry.file, name: entry.name || path.basename(entry.file),
                  kind: /\.(pdb|cif|mmcif)$/i.test(entry.file) ? 'structure' : 'artifact', command: entry.command,
                  parent: entry.parent, description: entry.description, allotype: entry.allotype,
                  iggSubtype: entry.iggSubtype,
                  ...('engineerChecksum' in entry ? {
                    engineerChecksum: entry.engineerChecksum, mutationIds: entry.mutationIds,
                    mutationsResolved: entry.mutationsResolved, hasGlycan: entry.hasGlycan, scheme: entry.scheme,
                  } : {}),
                })
              }
              return workspace
            }, legacyWorkspaceOwner)
          }

          res.end()
        } catch (err: any) {
          if (!res.headersSent) return sendJson(res, errorStatus(err), { error: err?.message ?? String(err) })
          sseSend(res, { step: 0, total: 0, name: 'fatal', status: 'error', stderr: err?.message ?? String(err) })
          res.end()
        }
      })

      // ── Health / config ───────────────────────────────────────────────
      server.middlewares.use('/api/status', async (req, res, next) => {
        if (req.method !== 'GET') return next()
        const pg = await getPg()
        const dvbfixer = process.env.DVBFIXER_EXECUTABLE || 'dvbfixer'
        sendJson(res, 200, {
          dvbfixer,
          databaseConfigured: !!process.env.DATABASE_URL,
          databaseConnected: !!pg,
        })
      })

      server.middlewares.use('/api/v1', (req, res) => sendJson(res, 404, {
        error: {
          code: 'ROUTE_NOT_FOUND', message: 'V1 route not found',
          requestId: ensureRequestId(req, res),
        },
      }))
}
