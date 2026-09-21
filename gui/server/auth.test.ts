import { createHash } from 'node:crypto'
import type { IncomingMessage, ServerResponse } from 'node:http'
import { Readable } from 'node:stream'
import { describe, expect, it, vi } from 'vitest'
import {
  AuthConfigError,
  createAuthMiddleware,
  getApiPrincipal,
  parseAuthConfig,
  resolveLegacyWorkspaceOwner,
} from './auth'

const TOKEN_A = Buffer.alloc(32, 0x11).toString('base64url')
const TOKEN_B = Buffer.alloc(32, 0x22).toString('base64url')

function tokenHash(token: string): string {
  return createHash('sha256').update(token, 'ascii').digest('hex')
}

function configured(principals = [
  { id: 'alice', tokenSha256: tokenHash(TOKEN_A) },
  { id: 'Alice', tokenSha256: tokenHash(TOKEN_B) },
]) {
  return JSON.stringify({ version: 1, principals })
}

function request(options: {
  method?: string
  url?: string
  authorization?: string
  rawAuthorizations?: string[]
} = {}): IncomingMessage {
  const stream = Readable.from(['request body']) as unknown as IncomingMessage
  stream.method = options.method ?? 'POST'
  stream.url = options.url ?? '/v1/workspaces'
  stream.headers = options.authorization === undefined
    ? {}
    : { authorization: options.authorization }
  stream.rawHeaders = (options.rawAuthorizations ?? []).flatMap(value => ['Authorization', value])
  return stream
}

interface ResponseCapture {
  readonly response: ServerResponse
  readonly headers: Map<string, string>
  statusCode: number
  body?: Buffer
}

function response(): ResponseCapture {
  const headers = new Map<string, string>()
  const capture: ResponseCapture = {
    headers,
    statusCode: 200,
    response: undefined as unknown as ServerResponse,
  }
  capture.response = {
    get statusCode() { return capture.statusCode },
    set statusCode(value: number) { capture.statusCode = value },
    setHeader(name: string, value: string | number | readonly string[]) {
      headers.set(name.toLowerCase(), String(value))
      return this
    },
    end(chunk?: Buffer) {
      capture.body = chunk
      return this
    },
  } as ServerResponse
  return capture
}

function runMiddleware(
  config: ReturnType<typeof parseAuthConfig>,
  req: IncomingMessage,
): { next: ReturnType<typeof vi.fn>; result: ResponseCapture } {
  const next = vi.fn()
  const result = response()
  createAuthMiddleware(config)(req, result.response, next)
  return { next, result }
}

describe('authentication configuration', () => {
  it('uses an immutable local principal only when configuration is absent', () => {
    const direct = parseAuthConfig(undefined)
    const fromEnvironment = parseAuthConfig({})

    expect(direct).toEqual({ mode: 'disabled', enabled: false, principal: { id: 'local' } })
    expect(fromEnvironment).toEqual(direct)
    expect(Object.isFrozen(direct)).toBe(true)
    if (direct.mode === 'disabled') expect(Object.isFrozen(direct.principal)).toBe(true)
  })

  it('parses, preserves case-sensitive IDs, and freezes configured principals', () => {
    const config = parseAuthConfig({ DVBFIXER_AUTH_PRINCIPALS: configured() })

    expect(config.mode).toBe('enabled')
    if (config.mode !== 'enabled') throw new Error('Expected enabled authentication')
    expect(config.principals.map(principal => principal.id)).toEqual(['alice', 'Alice'])
    expect(Object.isFrozen(config)).toBe(true)
    expect(Object.isFrozen(config.principals)).toBe(true)
    expect(config.principals.every(Object.isFrozen)).toBe(true)
  })

  it('rejects blank, malformed, oversized, empty, and unknown-field configuration', () => {
    const invalid = [
      '',
      '   ',
      '{secret',
      JSON.stringify({ version: 1, principals: [] }),
      JSON.stringify({ version: 1, principals: [], extra: true }),
      JSON.stringify({ version: 1, principals: [{
        id: 'alice', tokenSha256: tokenHash(TOKEN_A), token: TOKEN_A,
      }] }),
      ' '.repeat(64 * 1024 + 1),
    ]

    for (const value of invalid) expect(() => parseAuthConfig(value)).toThrow(AuthConfigError)
  })

  it('enforces principal identifiers, lowercase hashes, uniqueness, and count limits', () => {
    const validHash = tokenHash(TOKEN_A)
    const invalidPrincipals = [
      [{ id: '-alice', tokenSha256: validHash }],
      [{ id: 'a'.repeat(65), tokenSha256: validHash }],
      [{ id: 'alice', tokenSha256: validHash.toUpperCase() }],
      [{ id: 'alice', tokenSha256: validHash.slice(1) }],
      [{ id: 'alice', tokenSha256: validHash }, { id: 'alice', tokenSha256: tokenHash(TOKEN_B) }],
      [{ id: 'alice', tokenSha256: validHash }, { id: 'bob', tokenSha256: validHash }],
      Array.from({ length: 257 }, (_, index) => ({
        id: `p${index}`,
        tokenSha256: index.toString(16).padStart(64, '0'),
      })),
    ]

    for (const principals of invalidPrincipals) {
      expect(() => parseAuthConfig(configured(principals))).toThrow(AuthConfigError)
    }
  })

  it('never includes configuration contents in parse errors', () => {
    const secret = 'plaintext-token-that-must-not-leak'
    let error: unknown
    try {
      parseAuthConfig(JSON.stringify({
        version: 1,
        principals: [{ id: 'alice', tokenSha256: tokenHash(TOKEN_A), token: secret }],
      }))
    } catch (caught) {
      error = caught
    }
    expect(String(error)).not.toContain(secret)
  })

  it('resolves legacy ownership deterministically', () => {
    expect(resolveLegacyWorkspaceOwner(parseAuthConfig(undefined), {})).toBe('local')
    const single = parseAuthConfig(configured([{ id: 'alice', tokenSha256: tokenHash(TOKEN_A) }]))
    expect(resolveLegacyWorkspaceOwner(single, {})).toBe('alice')
    const multiple = parseAuthConfig(configured())
    expect(resolveLegacyWorkspaceOwner(multiple, { DVBFIXER_LEGACY_WORKSPACE_OWNER: 'Alice' })).toBe('Alice')
    expect(() => resolveLegacyWorkspaceOwner(multiple, {})).toThrow(/LEGACY_WORKSPACE_OWNER/)
    expect(() => resolveLegacyWorkspaceOwner(multiple, {
      DVBFIXER_LEGACY_WORKSPACE_OWNER: 'missing',
    })).toThrow(/configured principal/)
  })
})

describe('authentication middleware', () => {
  it('attaches the local principal to every request when disabled', () => {
    const req = request({ method: 'GET', url: '/health' })
    const { next } = runMiddleware(parseAuthConfig(undefined), req)

    expect(next).toHaveBeenCalledOnce()
    expect(getApiPrincipal(req)).toEqual({ id: 'local' })
  })

  it('authenticates one canonical bearer token and attaches its principal', () => {
    const req = request({ authorization: `Bearer ${TOKEN_B}` })
    const { next } = runMiddleware(parseAuthConfig(configured()), req)

    expect(next).toHaveBeenCalledOnce()
    expect(getApiPrincipal(req)).toEqual({ id: 'Alice' })
    expect(Object.isFrozen(getApiPrincipal(req))).toBe(true)
  })

  it('allows only the two exact public GET requests without a principal', () => {
    const config = parseAuthConfig(configured())
    for (const url of ['/health', '/v1/openapi.json']) {
      const req = request({ method: 'GET', url })
      const { next } = runMiddleware(config, req)
      expect(next).toHaveBeenCalledOnce()
      expect(() => getApiPrincipal(req)).toThrow(/No API principal/)
    }

    for (const [method, url] of [
      ['POST', '/health'],
      ['GET', '/health?full=1'],
      ['GET', '/v1/openapi.json/'],
    ]) {
      const { next, result } = runMiddleware(config, request({ method, url }))
      expect(next).not.toHaveBeenCalled()
      expect(result.statusCode).toBe(401)
    }
  })

  it('rejects malformed, noncanonical, duplicate, and incorrect credentials', () => {
    const config = parseAuthConfig(configured())
    const invalidRequests = [
      request(),
      request({ authorization: `Basic ${TOKEN_A}` }),
      request({ authorization: `Bearer ${TOKEN_A}=` }),
      request({ authorization: 'Bearer short' }),
      request({ authorization: `Bearer ${Buffer.alloc(32, 0x33).toString('base64url')}` }),
      request({
        authorization: `Bearer ${TOKEN_A}`,
        rawAuthorizations: [`Bearer ${TOKEN_A}`, `Bearer ${TOKEN_B}`],
      }),
    ]

    for (const req of invalidRequests) {
      const resume = vi.spyOn(req, 'resume')
      const { next, result } = runMiddleware(config, req)
      expect(next).not.toHaveBeenCalled()
      expect(resume).toHaveBeenCalledOnce()
      expect(result.statusCode).toBe(401)
      expect(result.headers.get('content-type')).toBe('application/json; charset=utf-8')
      expect(result.headers.get('cache-control')).toBe('no-store')
      expect(result.headers.get('www-authenticate')).toBe('Bearer realm="dvbfixer"')
      const body = JSON.parse(result.body?.toString('utf8') ?? '')
      expect(body).toMatchObject({
        error: { code: 'UNAUTHORIZED', message: 'Missing or invalid bearer credential' },
      })
      expect(body.error.requestId).toMatch(/^req_/)
      expect(() => getApiPrincipal(req)).toThrow(/No API principal/)
    }
  })
})
