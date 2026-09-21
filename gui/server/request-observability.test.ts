import { EventEmitter } from 'node:events'
import type { IncomingMessage, ServerResponse } from 'node:http'
import { describe, expect, it, vi } from 'vitest'
import {
  createRequestObservabilityMiddleware, ensureRequestId, getRequestId, parseAccessLogConfig,
  setRequestPrincipal,
} from './request-observability'

function exchange(url = '/v1/workspaces/secret-id/naming-conversions?token=secret') {
  const request = Object.assign(new EventEmitter(), {
    method: 'POST', url, headers: { 'x-request-id': 'attacker-value' },
  }) as unknown as IncomingMessage
  const headers = new Map<string, string>()
  const response = Object.assign(new EventEmitter(), {
    statusCode: 201, headersSent: true, writableFinished: true,
    setHeader: (name: string, value: string) => headers.set(name.toLowerCase(), value),
  }) as unknown as ServerResponse
  return { request, response, headers }
}

describe('request observability', () => {
  it('validates configuration', () => {
    expect(parseAccessLogConfig({})).toEqual({ enabled: false })
    expect(parseAccessLogConfig({ DVBFIXER_ACCESS_LOG: 'json' })).toEqual({ enabled: true })
    expect(() => parseAccessLogConfig({ DVBFIXER_ACCESS_LOG: 'text' })).toThrow(/off or json/)
  })

  it('generates one server-owned ID and ignores a supplied header', () => {
    const { request, response, headers } = exchange()
    const first = ensureRequestId(request, response)
    expect(first).toMatch(/^req_[0-9a-f-]{36}$/)
    expect(first).not.toBe('attacker-value')
    expect(ensureRequestId(request, response)).toBe(first)
    expect(getRequestId(request)).toBe(first)
    expect(headers.get('x-request-id')).toBe(first)
  })

  it('logs one redacted completion record with the authenticated principal', () => {
    const { request, response } = exchange()
    const lines: string[] = []
    let clock = 10
    createRequestObservabilityMiddleware({ enabled: true }, line => lines.push(line), () => clock)(
      request, response, () => setRequestPrincipal(request, 'CaseSensitive'),
    )
    clock = 12.5
    response.emit('finish')
    response.emit('close')

    expect(lines).toHaveLength(1)
    expect(JSON.parse(lines[0])).toMatchObject({
      schema: 'dvbfixer.http_access.v1', requestId: getRequestId(request), method: 'POST',
      route: '/api/v1/workspaces/:workspaceId/naming-conversions', status: 201,
      durationMs: 2.5, outcome: 'completed', principalId: 'CaseSensitive',
    })
    expect(lines[0]).not.toContain('secret-id')
    expect(lines[0]).not.toContain('token=secret')
    expect(lines[0]).not.toContain('attacker-value')
  })

  it('logs an abort without inventing a successful status', () => {
    const { request, response } = exchange('/unknown?secret=yes')
    Object.assign(response, { headersSent: false, writableFinished: false, statusCode: 200 })
    const lines: string[] = []
    createRequestObservabilityMiddleware({ enabled: true }, line => lines.push(line))(
      request, response, () => {},
    )
    response.emit('close')
    expect(JSON.parse(lines[0])).toMatchObject({
      route: '/api/unmatched', status: null, outcome: 'aborted', principalId: null,
    })
  })

  it('does not let sink failures escape into request handling', () => {
    const { request, response } = exchange()
    const next = vi.fn()
    createRequestObservabilityMiddleware({ enabled: true }, () => { throw new Error('sink failed') })(
      request, response, next,
    )
    expect(() => response.emit('finish')).not.toThrow()
    expect(next).toHaveBeenCalledOnce()
  })

  it('observes completion even when access logging is disabled', () => {
    const { request, response } = exchange('/health')
    const starts: unknown[] = []
    const finishes: unknown[] = []
    createRequestObservabilityMiddleware(
      { enabled: false }, undefined, () => 1000,
      { start: value => starts.push(value), finish: value => finishes.push(value) },
    )(request, response, () => {})
    response.emit('finish')
    expect(starts).toEqual([{ method: 'POST', route: '/api/health' }])
    expect(finishes).toEqual([{
      method: 'POST', route: '/api/health', status: 201, durationSeconds: 0, outcome: 'completed',
    }])
  })
})
