import type { IncomingMessage, ServerResponse } from 'node:http'
import { isIP } from 'node:net'
import type { ApiMiddleware } from './http-types'

export const DEFAULT_RATE_LIMIT_REQUESTS = 120
export const DEFAULT_RATE_LIMIT_WINDOW_MS = 60_000
export const DEFAULT_RATE_LIMIT_MAX_KEYS = 10_000

export interface RateLimitConfig {
  requests: number
  windowMs: number
  maxKeys: number
}

type RateLimitEnvironment = Readonly<Record<string, string | undefined>>

interface WindowEntry {
  count: number
  resetAt: number
}

function decimalSetting(
  name: string,
  value: string | undefined,
  fallback: number,
  minimum: number,
  maximum: number,
): number {
  if (value === undefined) return fallback
  if (!/^(0|[1-9]\d*)$/.test(value)) {
    throw new Error(`${name} must be a canonical non-negative integer`)
  }
  const parsed = Number(value)
  if (!Number.isSafeInteger(parsed) || parsed < minimum || parsed > maximum) {
    throw new Error(`${name} must be between ${minimum} and ${maximum}`)
  }
  return parsed
}

export function parseRateLimitConfig(
  environment: RateLimitEnvironment = process.env,
): RateLimitConfig {
  return Object.freeze({
    requests: decimalSetting(
      'DVBFIXER_RATE_LIMIT_REQUESTS',
      environment.DVBFIXER_RATE_LIMIT_REQUESTS,
      DEFAULT_RATE_LIMIT_REQUESTS,
      0,
      1_000_000,
    ),
    windowMs: decimalSetting(
      'DVBFIXER_RATE_LIMIT_WINDOW_MS',
      environment.DVBFIXER_RATE_LIMIT_WINDOW_MS,
      DEFAULT_RATE_LIMIT_WINDOW_MS,
      1_000,
      3_600_000,
    ),
    maxKeys: decimalSetting(
      'DVBFIXER_RATE_LIMIT_MAX_KEYS',
      environment.DVBFIXER_RATE_LIMIT_MAX_KEYS,
      DEFAULT_RATE_LIMIT_MAX_KEYS,
      1,
      100_000,
    ),
  })
}

function clientKey(request: IncomingMessage): string {
  const address = request.socket.remoteAddress
  if (!address) return '<unknown>'
  const mapped = /^::ffff:(\d{1,3}(?:\.\d{1,3}){3})$/i.exec(address)
  const candidate = mapped?.[1] || address
  if (isIP(candidate) === 0) return '<unknown>'
  return candidate.toLowerCase()
}

function secondsUntil(resetAt: number, now: number): number {
  return Math.max(1, Math.ceil((resetAt - now) / 1000))
}

function setRateHeaders(
  response: ServerResponse,
  config: RateLimitConfig,
  remaining: number,
  resetAt: number,
  now: number,
): void {
  response.setHeader('RateLimit-Limit', String(config.requests))
  response.setHeader('RateLimit-Remaining', String(Math.max(0, remaining)))
  response.setHeader('RateLimit-Reset', String(secondsUntil(resetAt, now)))
}

function sendJson(
  request: IncomingMessage,
  response: ServerResponse,
  status: number,
  body: { error: string },
): void {
  const encoded = Buffer.from(JSON.stringify(body))
  request.resume()
  response.statusCode = status
  response.setHeader('Content-Type', 'application/json; charset=utf-8')
  response.setHeader('Content-Length', String(encoded.length))
  response.setHeader('Cache-Control', 'no-store')
  response.end(encoded)
}

export function createRateLimitMiddleware(
  config: RateLimitConfig,
  now: () => number = Date.now,
): ApiMiddleware {
  const windows = new Map<string, WindowEntry>()
  const expirations: Array<{ key: string; entry: WindowEntry }> = []
  let expirationIndex = 0
  let lastTime = 0
  return (request, response, next) => {
    if (config.requests === 0) return next()

    const currentTime = Math.max(lastTime, now())
    lastTime = currentTime
    while (expirationIndex < expirations.length) {
      const expired = expirations[expirationIndex]
      if (expired.entry.resetAt > currentTime) break
      if (windows.get(expired.key) === expired.entry) windows.delete(expired.key)
      expirationIndex += 1
    }
    if (expirationIndex >= 1024 && expirationIndex * 2 >= expirations.length) {
      expirations.splice(0, expirationIndex)
      expirationIndex = 0
    }
    const key = clientKey(request)
    let entry = windows.get(key)
    if (entry && entry.resetAt <= currentTime) {
      windows.delete(key)
      entry = undefined
    }

    if (!entry) {
      if (windows.size >= config.maxKeys) {
        return sendJson(request, response, 503, { error: 'rate limit tracker capacity exhausted' })
      }
      entry = { count: 0, resetAt: currentTime + config.windowMs }
      windows.set(key, entry)
      expirations.push({ key, entry })
    }

    if (entry.count >= config.requests) {
      setRateHeaders(response, config, 0, entry.resetAt, currentTime)
      response.setHeader('Retry-After', String(secondsUntil(entry.resetAt, currentTime)))
      return sendJson(request, response, 429, { error: 'rate limit exceeded' })
    }

    entry.count += 1
    setRateHeaders(response, config, config.requests - entry.count, entry.resetAt, currentTime)
    return next()
  }
}
