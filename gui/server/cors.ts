import type { IncomingMessage, ServerResponse } from 'node:http'
import type { ApiMiddleware } from './http-types'

const CORS_ENV_NAME = 'DVBFIXER_CORS_ALLOWED_ORIGINS'
const MAX_CONFIG_BYTES = 16 * 1024
const MAX_ORIGINS = 64
const ALLOWED_METHODS = new Set(['GET', 'POST', 'PUT', 'PATCH', 'DELETE'])
const ALLOWED_HEADERS = new Set(['authorization', 'content-type', 'accept', 'x-file-name'])
const HEADER_NAME_PATTERN = /^[!#$%&'*+.^_`|~0-9A-Za-z-]+$/

type CorsEnvironment = Readonly<Record<string, string | undefined>>

export class CorsConfigError extends Error {
  constructor(message: string) {
    super(message)
    this.name = 'CorsConfigError'
  }
}

function configText(source: CorsEnvironment | string | undefined): string | undefined {
  if (typeof source === 'string' || source === undefined) return source
  if (!Object.hasOwn(source, CORS_ENV_NAME)) return undefined
  return source[CORS_ENV_NAME] ?? ''
}

function isCanonicalHttpOrigin(value: string): boolean {
  try {
    const parsed = new URL(value)
    return (parsed.protocol === 'http:' || parsed.protocol === 'https:')
      && parsed.origin === value
      && parsed.username === ''
      && parsed.password === ''
      && parsed.pathname === '/'
      && parsed.search === ''
      && parsed.hash === ''
  } catch {
    return false
  }
}

export function parseCorsAllowedOrigins(
  source: CorsEnvironment | string | undefined,
): readonly string[] {
  const text = configText(source)
  if (text === undefined) return Object.freeze([])
  if (text.trim() === '') throw new CorsConfigError(`${CORS_ENV_NAME} must not be blank`)
  if (Buffer.byteLength(text, 'utf8') > MAX_CONFIG_BYTES) {
    throw new CorsConfigError(`${CORS_ENV_NAME} exceeds the 16 KiB limit`)
  }

  let parsed: unknown
  try {
    parsed = JSON.parse(text)
  } catch {
    throw new CorsConfigError(`${CORS_ENV_NAME} must be valid JSON`)
  }
  if (!Array.isArray(parsed)) {
    throw new CorsConfigError(`${CORS_ENV_NAME} must be a JSON array`)
  }
  if (parsed.length > MAX_ORIGINS) {
    throw new CorsConfigError(`${CORS_ENV_NAME} must contain at most 64 origins`)
  }

  const seen = new Set<string>()
  const origins = parsed.map((candidate, index) => {
    if (typeof candidate !== 'string' || !isCanonicalHttpOrigin(candidate)) {
      throw new CorsConfigError(`Origin ${index} must be a canonical HTTP(S) origin`)
    }
    if (seen.has(candidate)) throw new CorsConfigError(`Origin ${index} is duplicated`)
    seen.add(candidate)
    return candidate
  })
  return Object.freeze(origins)
}

function headerValues(request: IncomingMessage, name: string): string[] {
  const values: string[] = []
  for (let index = 0; index < request.rawHeaders.length; index += 2) {
    if (request.rawHeaders[index]?.toLowerCase() === name) {
      values.push(request.rawHeaders[index + 1] ?? '')
    }
  }
  if (values.length > 0) return values

  const value = request.headers[name]
  if (Array.isArray(value)) return value
  return value === undefined ? [] : [value]
}

function appendVaryOrigin(response: ServerResponse): void {
  const existing = response.getHeader('Vary')
  const candidates = (Array.isArray(existing) ? existing : existing === undefined ? [] : [existing])
    .flatMap(value => String(value).split(','))
    .map(value => value.trim())
    .filter(Boolean)
  const seen = new Set<string>()
  const fields = candidates.filter(value => {
    const normalized = value.toLowerCase()
    if (seen.has(normalized)) return false
    seen.add(normalized)
    return true
  })
  if (!seen.has('origin')) fields.push('Origin')
  response.setHeader('Vary', fields.join(', '))
}

function requestOrigin(request: IncomingMessage): string | null | undefined {
  const values = headerValues(request, 'origin')
  if (values.length === 0) return undefined
  if (values.length !== 1 || !isCanonicalHttpOrigin(values[0])) return null
  return values[0]
}

function requestAuthority(request: IncomingMessage): string | null {
  const values = headerValues(request, 'host')
  if (values.length !== 1) return null
  const authority = values[0]
  if (authority === '' || authority !== authority.trim()
    || /[\s,@/\\?#]/.test(authority) || authority.endsWith(':')) return null

  const encrypted = (request.socket as typeof request.socket & { encrypted?: boolean }).encrypted === true
  const protocol = encrypted ? 'https:' : 'http:'
  try {
    const parsed = new URL(`${protocol}//${authority}`)
    if (parsed.username || parsed.password || parsed.pathname !== '/' || parsed.search || parsed.hash) {
      return null
    }
    return parsed.origin
  } catch {
    return null
  }
}

function singleHeader(request: IncomingMessage, name: string): string | null {
  const values = headerValues(request, name)
  return values.length === 1 ? values[0] : null
}

function reject(request: IncomingMessage, response: ServerResponse): void {
  const body = Buffer.from(JSON.stringify({ error: 'Forbidden' }))
  request.resume()
  response.statusCode = 403
  response.setHeader('Content-Type', 'application/json; charset=utf-8')
  response.setHeader('Content-Length', String(body.length))
  response.setHeader('Cache-Control', 'no-store')
  response.end(body)
}

function handlePreflight(request: IncomingMessage, response: ServerResponse, origin: string): boolean {
  const requestedMethod = singleHeader(request, 'access-control-request-method')
  if (requestedMethod === null || !ALLOWED_METHODS.has(requestedMethod)) return false

  const requestedHeaders = headerValues(request, 'access-control-request-headers')
  if (requestedHeaders.length > 1) return false
  if (headerValues(request, 'access-control-request-private-network').length > 0) return false

  let reflectedHeaders: string | undefined
  if (requestedHeaders.length === 1) {
    reflectedHeaders = requestedHeaders[0]
    const fields = reflectedHeaders.split(',').map(value => value.trim())
    const seen = new Set<string>()
    if (fields.length === 0 || fields.some(field => {
      const normalized = field.toLowerCase()
      if (!HEADER_NAME_PATTERN.test(field) || !ALLOWED_HEADERS.has(normalized) || seen.has(normalized)) {
        return true
      }
      seen.add(normalized)
      return false
    })) return false
  }

  response.statusCode = 204
  response.setHeader('Access-Control-Allow-Origin', origin)
  response.setHeader('Access-Control-Allow-Methods', requestedMethod)
  if (reflectedHeaders !== undefined) {
    response.setHeader('Access-Control-Allow-Headers', reflectedHeaders)
  }
  response.end()
  return true
}

export function createCorsMiddleware(
  allowedOrigins: readonly string[],
  rejectedRequestAdmission?: ApiMiddleware,
): ApiMiddleware {
  const allowlist = new Set(allowedOrigins)
  return (request, response, next) => {
    const rejectAfterAdmission = () => {
      if (rejectedRequestAdmission) {
        return rejectedRequestAdmission(request, response, () => reject(request, response))
      }
      return reject(request, response)
    }
    appendVaryOrigin(response)
    const origin = requestOrigin(request)

    if (origin === undefined) {
      if (request.method === 'OPTIONS') rejectAfterAdmission()
      else next()
      return
    }
    if (origin === null) {
      rejectAfterAdmission()
      return
    }

    const sameOrigin = requestAuthority(request) === origin
    if (!sameOrigin && !allowlist.has(origin)) {
      rejectAfterAdmission()
      return
    }

    if (request.method === 'OPTIONS') {
      if (!handlePreflight(request, response, origin)) rejectAfterAdmission()
      return
    }

    response.setHeader('Access-Control-Allow-Origin', origin)
    response.setHeader(
      'Access-Control-Expose-Headers',
      'RateLimit-Limit, RateLimit-Remaining, RateLimit-Reset, Retry-After',
    )
    next()
  }
}
