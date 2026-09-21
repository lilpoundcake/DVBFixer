import type { IncomingMessage, ServerResponse } from 'node:http'
import { Readable } from 'node:stream'
import { describe, expect, it, vi } from 'vitest'
import {
  CorsConfigError, createCorsMiddleware, parseCorsAllowedOrigins,
} from './cors'
import type { ApiNext } from './http-types'

function request(options: {
  method?: string
  url?: string
  encrypted?: boolean
  headers?: Array<[string, string]>
} = {}): IncomingMessage {
  const stream = Readable.from(['body']) as unknown as IncomingMessage
  stream.method = options.method ?? 'GET'
  stream.url = options.url ?? '/workspaces'
  stream.rawHeaders = (options.headers ?? []).flat()
  stream.headers = {}
  for (const [name, value] of options.headers ?? []) {
    stream.headers[name.toLowerCase()] = value
  }
  Object.defineProperty(stream, 'socket', { value: { encrypted: options.encrypted === true } })
  return stream
}

interface ResponseCapture {
  readonly response: ServerResponse
  readonly headers: Map<string, string>
  statusCode: number
  ended: boolean
  body?: Buffer
}

function response(initialVary?: string): ResponseCapture {
  const headers = new Map<string, string>()
  if (initialVary !== undefined) headers.set('vary', initialVary)
  const capture: ResponseCapture = {
    headers,
    statusCode: 200,
    ended: false,
    response: undefined as unknown as ServerResponse,
  }
  capture.response = {
    get statusCode() { return capture.statusCode },
    set statusCode(value: number) { capture.statusCode = value },
    getHeader(name: string) { return headers.get(name.toLowerCase()) },
    setHeader(name: string, value: string | number | readonly string[]) {
      headers.set(name.toLowerCase(), Array.isArray(value) ? value.join(', ') : String(value))
      return this
    },
    end(chunk?: Buffer) {
      capture.ended = true
      capture.body = chunk
      return this
    },
  } as ServerResponse
  return capture
}

function run(
  req: IncomingMessage,
  allowedOrigins: readonly string[] = [],
  initialVary?: string,
) {
  const result = response(initialVary)
  const next = vi.fn()
  createCorsMiddleware(allowedOrigins)(req, result.response, next)
  return { next, result }
}

describe('CORS configuration', () => {
  it('treats an absent setting as an immutable empty allowlist', () => {
    const direct = parseCorsAllowedOrigins(undefined)
    expect(direct).toEqual([])
    expect(parseCorsAllowedOrigins({})).toEqual([])
    expect(Object.isFrozen(direct)).toBe(true)
  })

  it('accepts only unique canonical HTTP(S) origins', () => {
    const origins = parseCorsAllowedOrigins(JSON.stringify([
      'https://api.example.test',
      'http://localhost:5173',
      'https://[2001:db8::1]:8443',
    ]))
    expect(origins).toEqual([
      'https://api.example.test',
      'http://localhost:5173',
      'https://[2001:db8::1]:8443',
    ])
    expect(Object.isFrozen(origins)).toBe(true)
  })

  it('rejects blank, malformed, oversized, overlong, duplicate, and noncanonical settings', () => {
    const invalid = [
      '',
      '  ',
      '{bad',
      '{}',
      JSON.stringify(['https://example.test', 'https://example.test']),
      JSON.stringify(['null']),
      JSON.stringify(['ftp://example.test']),
      JSON.stringify(['https://EXAMPLE.test']),
      JSON.stringify(['https://example.test/']),
      JSON.stringify(['https://example.test:443']),
      JSON.stringify(['https://user@example.test']),
      JSON.stringify(Array.from({ length: 65 }, (_, index) => `https://x${index}.test`)),
      JSON.stringify(['https://example.test', 'x'.repeat(16 * 1024)]),
    ]
    for (const value of invalid) {
      expect(() => parseCorsAllowedOrigins(value)).toThrow(CorsConfigError)
    }
  })
})

describe('CORS middleware', () => {
  it('allows requests without Origin and adds a deduplicated Vary field', () => {
    const { next, result } = run(request(), [], 'Accept-Encoding, origin, Accept-Encoding')
    expect(next).toHaveBeenCalledOnce()
    expect(result.headers.get('vary')).toBe('Accept-Encoding, origin')
    expect(result.headers.has('access-control-allow-origin')).toBe(false)
  })

  it('accepts exact same-origin requests using socket encryption and strict Host parsing', () => {
    const httpResult = run(request({ headers: [
      ['Host', 'example.test:8080'],
      ['Origin', 'http://example.test:8080'],
    ] }))
    const httpsResult = run(request({ encrypted: true, headers: [
      ['Host', 'secure.example.test'],
      ['Origin', 'https://secure.example.test'],
    ] }))

    expect(httpResult.next).toHaveBeenCalledOnce()
    expect(httpResult.result.headers.get('access-control-allow-origin')).toBe('http://example.test:8080')
    expect(httpsResult.next).toHaveBeenCalledOnce()
    expect(httpsResult.result.headers.get('access-control-allow-origin')).toBe('https://secure.example.test')
  })

  it('accepts exact allowlisted origins independently of the request Host', () => {
    const { next, result } = run(request({ headers: [
      ['Host', 'service.example.test'],
      ['Origin', 'https://client.example.test'],
    ] }), ['https://client.example.test'])

    expect(next).toHaveBeenCalledOnce()
    expect(result.headers.get('access-control-allow-origin')).toBe('https://client.example.test')
    expect(result.headers.get('access-control-expose-headers')).toContain('RateLimit-Remaining')
  })

  it('rejects null, malformed, untrusted, and duplicate Origin values', () => {
    const cases = [
      [['Host', 'example.test'], ['Origin', 'null']],
      [['Host', 'example.test'], ['Origin', 'https://example.test/']],
      [['Host', 'example.test'], ['Origin', 'https://other.test']],
      [['Host', 'example.test'], ['Origin', 'http://example.test'], ['Origin', 'http://example.test']],
      [['Host', 'example.test'], ['Host', 'example.test'], ['Origin', 'http://example.test']],
      [['Host', 'bad host'], ['Origin', 'http://bad']],
    ]

    for (const headers of cases) {
      const req = request({ headers: headers as Array<[string, string]> })
      const resume = vi.spyOn(req, 'resume')
      const { next, result } = run(req)
      expect(next).not.toHaveBeenCalled()
      expect(resume).toHaveBeenCalledOnce()
      expect(result.statusCode).toBe(403)
      expect(result.headers.get('vary')).toBe('Origin')
      expect(result.headers.get('cache-control')).toBe('no-store')
    }
  })

  it('uses the common error envelope for V1 CORS rejection', () => {
    const { result } = run(request({
      url: '/v1/workspaces/example/jobs',
      headers: [['Host', 'api.example.test'], ['Origin', 'https://other.test']],
    }))
    expect(JSON.parse(result.body!.toString())).toMatchObject({
      error: { code: 'CORS_FORBIDDEN', message: 'Origin is not allowed', requestId: expect.stringMatching(/^req_/) },
    })
  })

  it('answers valid preflight before downstream middleware with exact reflected permissions', () => {
    const { next, result } = run(request({ method: 'OPTIONS', headers: [
      ['Host', 'api.example.test'],
      ['Origin', 'https://client.example.test'],
      ['Access-Control-Request-Method', 'PATCH'],
      ['Access-Control-Request-Headers', 'Authorization, X-File-Name'],
    ] }), ['https://client.example.test'])

    expect(next).not.toHaveBeenCalled()
    expect(result.statusCode).toBe(204)
    expect(result.ended).toBe(true)
    expect(result.body).toBeUndefined()
    expect(result.headers.get('access-control-allow-origin')).toBe('https://client.example.test')
    expect(result.headers.get('access-control-allow-methods')).toBe('PATCH')
    expect(result.headers.get('access-control-allow-headers')).toBe('Authorization, X-File-Name')
    expect(result.headers.has('access-control-allow-credentials')).toBe(false)
    expect(result.headers.has('access-control-allow-private-network')).toBe(false)
  })

  it('supports every allowed method and an omitted requested-header list', () => {
    for (const method of ['GET', 'POST', 'PUT', 'PATCH', 'DELETE']) {
      const { result } = run(request({ method: 'OPTIONS', headers: [
        ['Host', 'api.example.test'],
        ['Origin', 'http://api.example.test'],
        ['Access-Control-Request-Method', method],
      ] }))
      expect(result.statusCode).toBe(204)
      expect(result.headers.get('access-control-allow-methods')).toBe(method)
      expect(result.headers.has('access-control-allow-headers')).toBe(false)
    }
  })

  it('rejects invalid preflight requests without permissive response headers', () => {
    const cases: Array<Array<[string, string]>> = [
      [['Host', 'api.example.test']],
      [['Host', 'api.example.test'], ['Origin', 'http://api.example.test']],
      [['Host', 'api.example.test'], ['Origin', 'http://api.example.test'], ['Access-Control-Request-Method', 'OPTIONS']],
      [['Host', 'api.example.test'], ['Origin', 'http://api.example.test'], ['Access-Control-Request-Method', 'get']],
      [['Host', 'api.example.test'], ['Origin', 'http://api.example.test'], ['Access-Control-Request-Method', 'GET'], ['Access-Control-Request-Headers', 'Cookie']],
      [['Host', 'api.example.test'], ['Origin', 'http://api.example.test'], ['Access-Control-Request-Method', 'GET'], ['Access-Control-Request-Headers', 'Accept, accept']],
      [['Host', 'api.example.test'], ['Origin', 'http://api.example.test'], ['Access-Control-Request-Method', 'GET'], ['Access-Control-Request-Headers', 'Accept'], ['Access-Control-Request-Headers', 'Content-Type']],
      [['Host', 'api.example.test'], ['Origin', 'http://api.example.test'], ['Access-Control-Request-Method', 'GET'], ['Access-Control-Request-Private-Network', 'true']],
    ]

    for (const headers of cases) {
      const { next, result } = run(request({ method: 'OPTIONS', headers }))
      expect(next).not.toHaveBeenCalled()
      expect(result.statusCode).toBe(403)
      expect(result.headers.has('access-control-allow-origin')).toBe(false)
      expect(result.headers.has('access-control-allow-methods')).toBe(false)
      expect(result.headers.has('access-control-allow-headers')).toBe(false)
    }
  })

  it('admits rejected requests to the limiter but exempts valid preflight', () => {
    const admission = vi.fn((_req: IncomingMessage, _res: ServerResponse, next: ApiNext) => next())
    const validResponse = response()
    createCorsMiddleware(['https://client.example.test'], admission)(request({
      method: 'OPTIONS',
      headers: [
        ['Host', 'api.example.test'],
        ['Origin', 'https://client.example.test'],
        ['Access-Control-Request-Method', 'GET'],
      ],
    }), validResponse.response, vi.fn())
    expect(validResponse.statusCode).toBe(204)
    expect(admission).not.toHaveBeenCalled()

    const invalidResponse = response()
    createCorsMiddleware([], admission)(request({
      method: 'OPTIONS', headers: [['Host', 'api.example.test']],
    }), invalidResponse.response, vi.fn())
    expect(admission).toHaveBeenCalledOnce()
    expect(invalidResponse.statusCode).toBe(403)
  })
})
