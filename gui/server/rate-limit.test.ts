import { describe, expect, it } from 'vitest'
import { Readable } from 'node:stream'
import {
  createRateLimitMiddleware, DEFAULT_RATE_LIMIT_MAX_KEYS, DEFAULT_RATE_LIMIT_REQUESTS,
  DEFAULT_RATE_LIMIT_WINDOW_MS, parseRateLimitConfig, type RateLimitConfig,
} from './rate-limit'

function invoke(
  middleware: ReturnType<typeof createRateLimitMiddleware>,
  address: string | undefined,
  method = 'GET',
) {
  const request = Readable.from([]) as any
  request.method = method
  request.socket = { remoteAddress: address }
  const headers = new Map<string, string>()
  let body = ''
  let nextCalled = false
  const response = {
    statusCode: 200,
    setHeader(name: string, value: unknown) { headers.set(name.toLowerCase(), String(value)) },
    end(value?: unknown) { body = value === undefined ? '' : Buffer.from(value as any).toString() },
  }
  middleware(request, response as any, () => { nextCalled = true })
  return { status: response.statusCode, headers, body, nextCalled }
}

const config = (overrides: Partial<RateLimitConfig> = {}): RateLimitConfig => ({
  requests: 2,
  windowMs: 10_000,
  maxKeys: 10,
  ...overrides,
})

describe('rate limit configuration', () => {
  it('uses bounded defaults and accepts zero only for disabling requests', () => {
    expect(parseRateLimitConfig({})).toEqual({
      requests: DEFAULT_RATE_LIMIT_REQUESTS,
      windowMs: DEFAULT_RATE_LIMIT_WINDOW_MS,
      maxKeys: DEFAULT_RATE_LIMIT_MAX_KEYS,
    })
    expect(parseRateLimitConfig({
      DVBFIXER_RATE_LIMIT_REQUESTS: '0',
      DVBFIXER_RATE_LIMIT_WINDOW_MS: '1000',
      DVBFIXER_RATE_LIMIT_MAX_KEYS: '1',
    })).toEqual({ requests: 0, windowMs: 1000, maxKeys: 1 })
  })

  it.each(['', '-1', '+1', '01', '1.5', ' 1', '1000001'])(
    'rejects invalid request limits: %j',
    value => expect(() => parseRateLimitConfig({ DVBFIXER_RATE_LIMIT_REQUESTS: value })).toThrow(),
  )
})

describe('rate limit middleware', () => {
  it('allows the configured count and returns retry metadata when exhausted', () => {
    let time = 1_000
    const middleware = createRateLimitMiddleware(config(), () => time)
    expect(invoke(middleware, '192.0.2.1')).toMatchObject({ status: 200, nextCalled: true })
    const second = invoke(middleware, '192.0.2.1')
    expect(second.headers.get('ratelimit-remaining')).toBe('0')
    const rejected = invoke(middleware, '192.0.2.1')
    expect(rejected).toMatchObject({ status: 429, nextCalled: false })
    expect(rejected.headers.get('retry-after')).toBe('10')
    expect(JSON.parse(rejected.body)).toEqual({ error: 'rate limit exceeded' })
    time = 11_000
    expect(invoke(middleware, '192.0.2.1')).toMatchObject({ status: 200, nextCalled: true })
  })

  it('normalizes IPv4-mapped IPv6 and never trusts forwarding headers', () => {
    const middleware = createRateLimitMiddleware(config({ requests: 1 }))
    expect(invoke(middleware, '::ffff:192.0.2.8').status).toBe(200)
    expect(invoke(middleware, '192.0.2.8').status).toBe(429)
  })

  it('counts OPTIONS presented directly and supports a disabled limiter', () => {
    const limited = createRateLimitMiddleware(config({ requests: 1 }))
    expect(invoke(limited, '192.0.2.9', 'OPTIONS').nextCalled).toBe(true)
    expect(invoke(limited, '192.0.2.9').status).toBe(429)
    const disabled = createRateLimitMiddleware(config({ requests: 0 }))
    expect(invoke(disabled, '192.0.2.9')).toMatchObject({ status: 200, nextCalled: true })
    expect(invoke(disabled, '192.0.2.9')).toMatchObject({ status: 200, nextCalled: true })
  })

  it('fails closed for new clients when active tracker capacity is exhausted', () => {
    let time = 0
    const middleware = createRateLimitMiddleware(config({ maxKeys: 1 }), () => time)
    expect(invoke(middleware, '192.0.2.10').status).toBe(200)
    expect(invoke(middleware, '192.0.2.11')).toMatchObject({ status: 503, nextCalled: false })
    time = 10_000
    expect(invoke(middleware, '192.0.2.11').status).toBe(200)
  })

  it('does not rescan active entries for each rejected overflow key', () => {
    let clockReads = 0
    const middleware = createRateLimitMiddleware(config({ maxKeys: 1 }), () => {
      clockReads += 1
      return 0
    })
    expect(invoke(middleware, '192.0.2.20').status).toBe(200)
    expect(invoke(middleware, '192.0.2.21').status).toBe(503)
    expect(invoke(middleware, '192.0.2.22').status).toBe(503)
    expect(clockReads).toBe(3)
  })

  it('shares one bounded key for requests without a remote address', () => {
    const middleware = createRateLimitMiddleware(config({ requests: 1 }))
    expect(invoke(middleware, undefined).status).toBe(200)
    expect(invoke(middleware, 'not-an-ip').status).toBe(429)
  })
})
