import { createHash, timingSafeEqual } from 'node:crypto'
import type { IncomingMessage, ServerResponse } from 'node:http'
import type { ApiMiddleware } from './http-types'
import { ensureRequestId, setRequestPrincipal } from './request-observability'

const AUTH_ENV_NAME = 'DVBFIXER_AUTH_PRINCIPALS'
const LEGACY_OWNER_ENV_NAME = 'DVBFIXER_LEGACY_WORKSPACE_OWNER'
const MAX_CONFIG_BYTES = 64 * 1024
const MAX_PRINCIPALS = 256
const PRINCIPAL_ID_PATTERN = /^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$/
const SHA256_PATTERN = /^[0-9a-f]{64}$/
const BEARER_TOKEN_PATTERN = /^[A-Za-z0-9_-]{43}$/

export interface ApiPrincipal {
  readonly id: string
}

export interface ConfiguredAuthPrincipal extends ApiPrincipal {
  readonly tokenSha256: string
}

export type AuthConfig =
  | Readonly<{
    mode: 'disabled'
    enabled: false
    principal: ApiPrincipal
  }>
  | Readonly<{
    mode: 'enabled'
    enabled: true
    principals: readonly ConfiguredAuthPrincipal[]
  }>

export class AuthConfigError extends Error {
  constructor(message: string) {
    super(message)
    this.name = 'AuthConfigError'
  }
}

type AuthEnvironment = Readonly<Record<string, string | undefined>>

const requestPrincipals = new WeakMap<IncomingMessage, ApiPrincipal>()

function frozenPrincipal(id: string): ApiPrincipal {
  return Object.freeze({ id })
}

function configText(source: AuthEnvironment | string | undefined): string | undefined {
  if (typeof source === 'string' || source === undefined) return source
  if (!Object.hasOwn(source, AUTH_ENV_NAME)) return undefined
  return source[AUTH_ENV_NAME] ?? ''
}

function isRecord(value: unknown): value is Record<string, unknown> {
  return typeof value === 'object' && value !== null && !Array.isArray(value)
}

function hasExactFields(value: Record<string, unknown>, fields: readonly string[]): boolean {
  const keys = Object.keys(value)
  return keys.length === fields.length && fields.every(field => Object.hasOwn(value, field))
}

export function parseAuthConfig(source: AuthEnvironment | string | undefined): AuthConfig {
  const text = configText(source)
  if (text === undefined) {
    return Object.freeze({
      mode: 'disabled',
      enabled: false,
      principal: frozenPrincipal('local'),
    })
  }
  if (text.trim() === '') throw new AuthConfigError(`${AUTH_ENV_NAME} must not be blank`)
  if (Buffer.byteLength(text, 'utf8') > MAX_CONFIG_BYTES) {
    throw new AuthConfigError(`${AUTH_ENV_NAME} exceeds the 64 KiB limit`)
  }

  let parsed: unknown
  try {
    parsed = JSON.parse(text)
  } catch {
    throw new AuthConfigError(`${AUTH_ENV_NAME} must be valid JSON`)
  }
  if (!isRecord(parsed) || !hasExactFields(parsed, ['version', 'principals'])) {
    throw new AuthConfigError(`${AUTH_ENV_NAME} has an invalid top-level shape`)
  }
  if (parsed.version !== 1 || !Array.isArray(parsed.principals)) {
    throw new AuthConfigError(`${AUTH_ENV_NAME} must contain version 1 and a principals array`)
  }
  if (parsed.principals.length === 0 || parsed.principals.length > MAX_PRINCIPALS) {
    throw new AuthConfigError(`${AUTH_ENV_NAME} must contain between 1 and 256 principals`)
  }

  const ids = new Set<string>()
  const hashes = new Set<string>()
  const principals = parsed.principals.map((candidate, index): ConfiguredAuthPrincipal => {
    if (!isRecord(candidate) || !hasExactFields(candidate, ['id', 'tokenSha256'])) {
      throw new AuthConfigError(`Principal ${index} has an invalid shape`)
    }
    const { id, tokenSha256 } = candidate
    if (typeof id !== 'string' || !PRINCIPAL_ID_PATTERN.test(id)) {
      throw new AuthConfigError(`Principal ${index} has an invalid id`)
    }
    if (ids.has(id)) throw new AuthConfigError(`Principal ${index} has a duplicate id`)
    if (typeof tokenSha256 !== 'string' || !SHA256_PATTERN.test(tokenSha256)) {
      throw new AuthConfigError(`Principal ${index} has an invalid tokenSha256`)
    }
    if (hashes.has(tokenSha256)) {
      throw new AuthConfigError(`Principal ${index} has a duplicate tokenSha256`)
    }
    ids.add(id)
    hashes.add(tokenSha256)
    return Object.freeze({ id, tokenSha256 })
  })

  return Object.freeze({
    mode: 'enabled',
    enabled: true,
    principals: Object.freeze(principals),
  })
}

export function resolveLegacyWorkspaceOwner(
  config: AuthConfig,
  environment: AuthEnvironment = process.env,
): string {
  if (config.mode === 'disabled') return config.principal.id
  const configured = environment[LEGACY_OWNER_ENV_NAME]?.trim()
  if (configured) {
    if (!config.principals.some(principal => principal.id === configured)) {
      throw new AuthConfigError(`${LEGACY_OWNER_ENV_NAME} must name a configured principal`)
    }
    return configured
  }
  if (config.principals.length === 1) return config.principals[0].id
  throw new AuthConfigError(
    `${LEGACY_OWNER_ENV_NAME} is required when multiple principals are configured`,
  )
}

function authorizationValues(request: IncomingMessage): string[] {
  const values: string[] = []
  for (let index = 0; index < request.rawHeaders.length; index += 2) {
    if (request.rawHeaders[index]?.toLowerCase() === 'authorization') {
      values.push(request.rawHeaders[index + 1] ?? '')
    }
  }
  if (values.length > 0) return values

  const header = request.headers.authorization
  if (Array.isArray(header)) return header
  return header === undefined ? [] : [header]
}

function bearerToken(request: IncomingMessage): string | undefined {
  const values = authorizationValues(request)
  if (values.length !== 1) return undefined
  const match = /^Bearer ([A-Za-z0-9_-]{43})$/i.exec(values[0])
  if (!match || !BEARER_TOKEN_PATTERN.test(match[1])) return undefined
  const decoded = Buffer.from(match[1], 'base64url')
  if (decoded.length !== 32 || decoded.toString('base64url') !== match[1]) return undefined
  return match[1]
}

function isPublicRequest(request: IncomingMessage): boolean {
  return request.method === 'GET'
    && (request.url === '/health' || request.url === '/v1/openapi.json')
}

function rejectUnauthorized(request: IncomingMessage, response: ServerResponse): void {
  const versioned = (request.url || '').startsWith('/v1/')
  const requestId = versioned ? ensureRequestId(request, response) : null
  const body = Buffer.from(JSON.stringify(versioned ? {
    error: { code: 'UNAUTHORIZED', message: 'Missing or invalid bearer credential', requestId },
  } : { error: 'Unauthorized' }))
  request.resume()
  response.statusCode = 401
  response.setHeader('Content-Type', 'application/json; charset=utf-8')
  response.setHeader('Content-Length', String(body.length))
  response.setHeader('Cache-Control', 'no-store')
  response.setHeader('WWW-Authenticate', 'Bearer realm="dvbfixer"')
  if (requestId) response.setHeader('X-Request-Id', requestId)
  response.end(body)
}

export function createAuthMiddleware(config: AuthConfig): ApiMiddleware {
  if (config.mode === 'disabled') {
    return (request, _response, next) => {
      requestPrincipals.set(request, config.principal)
      setRequestPrincipal(request, config.principal.id)
      next()
    }
  }

  const configured = config.principals.map(({ id, tokenSha256 }) => ({
    principal: frozenPrincipal(id),
    tokenHash: Buffer.from(tokenSha256, 'hex'),
  }))

  return (request, response, next) => {
    if (isPublicRequest(request)) {
      next()
      return
    }

    const token = bearerToken(request)
    if (token === undefined) {
      rejectUnauthorized(request, response)
      return
    }

    const presentedHash = createHash('sha256').update(token, 'ascii').digest()
    const matches = configured.map(entry => timingSafeEqual(presentedHash, entry.tokenHash))
    const matchIndex = matches.findIndex(Boolean)
    if (matchIndex < 0) {
      rejectUnauthorized(request, response)
      return
    }

    requestPrincipals.set(request, configured[matchIndex].principal)
    setRequestPrincipal(request, configured[matchIndex].principal.id)
    next()
  }
}

export function getApiPrincipal(request: IncomingMessage): ApiPrincipal {
  const principal = requestPrincipals.get(request)
  if (principal === undefined) throw new Error('No API principal is attached to this request')
  return principal
}
