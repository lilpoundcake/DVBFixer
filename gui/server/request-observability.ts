import { randomUUID } from 'node:crypto'
import type { IncomingMessage, ServerResponse } from 'node:http'
import type { ApiMiddleware } from './http-types'
import type { AuditEvent } from './audit-log'

interface RequestContext {
  requestId: string
  principalId: string | null
}

export interface AccessLogConfig {
  enabled: boolean
}

export type AccessLogSink = (line: string) => void

export interface HttpObservationStart { method: string; route: string }
export interface HttpObservation extends HttpObservationStart {
  status: number | null
  durationSeconds: number
  outcome: 'completed' | 'aborted'
}
export interface HttpObserver {
  start(observation: HttpObservationStart): void
  finish(observation: HttpObservation): void
}

const contexts = new WeakMap<IncomingMessage, RequestContext>()

export function parseAccessLogConfig(environment: NodeJS.ProcessEnv = process.env): AccessLogConfig {
  const value = environment.DVBFIXER_ACCESS_LOG ?? 'off'
  if (value !== 'off' && value !== 'json') throw new Error('DVBFIXER_ACCESS_LOG must be off or json')
  return { enabled: value === 'json' }
}

export function ensureRequestId(request: IncomingMessage, response: ServerResponse): string {
  let context = contexts.get(request)
  if (!context) {
    context = { requestId: `req_${randomUUID()}`, principalId: null }
    contexts.set(request, context)
  }
  response.setHeader('X-Request-Id', context.requestId)
  return context.requestId
}

export function getRequestId(request: IncomingMessage): string {
  const context = contexts.get(request)
  if (!context) throw new Error('No request context is attached to this request')
  return context.requestId
}

export function setRequestPrincipal(request: IncomingMessage, principalId: string): void {
  const context = contexts.get(request)
  if (context) context.principalId = principalId
}

export function routeLabel(url: string | undefined): string {
  const pathname = (url || '').split('?', 1)[0]
  if (pathname === '/health' || pathname === '/session' || pathname === '/status' ||
      pathname === '/dvbfixer-spec' || pathname === '/metrics' || pathname === '/v1/openapi.json') return `/api${pathname}`
  if (/^\/v1\/workspaces\/[^/]+\/naming-conversions$/.test(pathname)) {
    return '/api/v1/workspaces/:workspaceId/naming-conversions'
  }
  if (/^\/v1\/workspaces\/[^/]+\/jobs(?:\/[^/]+(?:\/events)?)?$/.test(pathname)) {
    return '/api/v1/workspaces/:workspaceId/jobs'
  }
  const first = pathname.split('/').filter(Boolean)[0]
  const known = new Set(['workspaces', 'homology', 'jobs', 'dvbfixer', 'mutations', 'antibody-engineer'])
  if (first && known.has(first)) return `/api/${first}`
  return pathname.startsWith('/v1/') ? '/api/v1/unmatched' : '/api/unmatched'
}

export function createRequestObservabilityMiddleware(
  config: AccessLogConfig,
  sink: AccessLogSink = line => console.log(line),
  now: () => number = () => performance.now(),
  observer?: HttpObserver,
  auditSink?: (event: AuditEvent) => void,
): ApiMiddleware {
  return (request, response, next) => {
    const requestId = ensureRequestId(request, response)
    const started = now()
    const route = routeLabel(request.url)
    const method = /^[A-Z]{1,16}$/.test(request.method || '') ? request.method! : 'OTHER'
    observer?.start({ method, route })
    if (config.enabled || observer) {
      let emitted = false
      const emit = (outcome: 'completed' | 'aborted') => {
        if (emitted) return
        emitted = true
        const context = contexts.get(request)
        const status = response.headersSent || outcome === 'completed' ? response.statusCode : null
        const durationSeconds = Math.max(0, (now() - started) / 1000)
        observer?.finish({ method, route, status, durationSeconds, outcome })
        const record = {
          schema: 'dvbfixer.http_access.v1',
          event: 'http.server.request',
          timestamp: new Date().toISOString(),
          requestId,
          method,
          route,
          status,
          durationMs: Math.round(durationSeconds * 1_000_000) / 1000,
          outcome,
          principalId: context?.principalId ?? null,
        }
        if (config.enabled) try { sink(JSON.stringify(record)) } catch { /* logging must not affect requests */ }
        if (auditSink && (route.startsWith('/api/v1/') || method !== 'GET' || status === 401 || status === 403)) {
          try {
            auditSink({
              schema: 'dvbfixer.audit.v1', timestamp: record.timestamp, requestId,
              method, route, status, outcome, principalId: record.principalId,
            })
          } catch (error) {
            console.error('[audit] failed to persist event:', error)
          }
        }
      }
      response.once('finish', () => emit('completed'))
      response.once('close', () => emit(response.writableFinished ? 'completed' : 'aborted'))
    }
    next()
  }
}
