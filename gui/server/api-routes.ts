/**
 * Host-neutral registration for the GUI and versioned HTTP APIs.
 *
 *   POST   /api/v1/workspaces/:id/jobs — start a managed DVBFixer workflow.
 *
 * Env vars:
 *   DVBFIXER_EXECUTABLE  default 'dvbfixer' — CLI executable
 *   DVBFIXER_ARGS        optional JSON string array prepended to CLI arguments
 */

import type { IncomingMessage, ServerResponse } from 'node:http'
import path from 'node:path'
import fs from 'node:fs'
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
import { registerWorkspaceApi, migrateWorkspaceOwnership } from './workspace-api'
import type { ApiRouteHost } from './http-types'
import {
  createAuthMiddleware, getApiPrincipal, parseAuthConfig, resolveLegacyWorkspaceOwner,
  type AuthConfig,
} from './auth'
import { createCorsMiddleware, parseCorsAllowedOrigins } from './cors'
import { createRateLimitMiddleware, parseRateLimitConfig, type RateLimitConfig } from './rate-limit'
import { configureStorageQuota, type StorageQuotaOptions } from './storage-quota'
import {
  createRequestObservabilityMiddleware, ensureRequestId, parseAccessLogConfig, type AccessLogConfig,
} from './request-observability'
import { ApiMetrics, parseMetricsConfig, type MetricsConfig } from './metrics'
import { createAuditSink, parseAuditLogConfig, type AuditLogConfig } from './audit-log'
export { runDvbfixer } from './dvbfixer-runner'
export { buildArgs } from './command-args'

export async function shutdownApiWork(): Promise<void> {
  shutdownManagedJobs()
  await shutdownDvbfixerProcesses()
}

export function resetApiShutdown(): void {
  resetManagedJobAdmission()
  resetDvbfixerProcessAdmission()
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
          authenticationRequired: authConfig.enabled,
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

      // ── Spec exposure (so frontend doesn't import server/) ────────────
      server.middlewares.use('/api/dvbfixer-spec', async (req, res, next) => {
        if (req.method !== 'GET') return next()
        sendJson(res, 200, COMMANDS)
      })

      // ── Health / config ───────────────────────────────────────────────
      server.middlewares.use('/api/v1', (req, res) => sendJson(res, 404, {
        error: {
          code: 'ROUTE_NOT_FOUND', message: 'V1 route not found',
          requestId: ensureRequestId(req, res),
        },
      }))
}
