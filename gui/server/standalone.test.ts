import { afterEach, describe, expect, it, vi } from 'vitest'
import fs from 'node:fs'
import http from 'node:http'
import os from 'node:os'
import path from 'node:path'
import crypto from 'node:crypto'
import {
  createStandaloneApplication,
  createStandaloneServer,
  loadStandaloneConfig,
  type StandaloneConfig,
} from './standalone'
import { parseAuthConfig } from './auth'

const temporaryDirectories: string[] = []

function temp(): string {
  const directory = fs.mkdtempSync(path.join(os.tmpdir(), 'dvbfixer-standalone-test-'))
  temporaryDirectories.push(directory)
  return directory
}

function fixture(): StandaloneConfig {
  const projectRoot = temp()
  const staticRoot = path.join(projectRoot, 'dist')
  fs.mkdirSync(path.join(staticRoot, 'assets'), { recursive: true })
  fs.writeFileSync(path.join(projectRoot, 'package.json'), JSON.stringify({ version: '9.9.9' }))
  fs.writeFileSync(path.join(staticRoot, 'index.html'), '<!doctype html><title>DVBFixer</title>')
  fs.writeFileSync(path.join(staticRoot, 'assets', 'app.js'), 'console.log("ok")')
  return {
    ...loadStandaloneConfig({}, projectRoot),
    port: 0,
    shutdownGraceMs: 100,
  }
}

afterEach(() => temporaryDirectories.splice(0).forEach(directory =>
  fs.rmSync(directory, { recursive: true, force: true })))

describe('standalone configuration', () => {
  it('loads loopback defaults and resolves paths from the project root', () => {
    const root = temp()
    const config = loadStandaloneConfig({}, root)
    expect(config).toMatchObject({
      host: '127.0.0.1', port: 5173,
      projectRoot: root,
      dataRoot: path.join(root, 'structures'),
      staticRoot: path.join(root, 'dist'),
      mutationsBackupFile: path.join(root, 'mutations.json'),
      deploymentResources: { required: false },
      accessLog: { enabled: false },
      metrics: { enabled: false },
    })
  })

  it('rejects invalid ports and remote binding without explicit acknowledgement', () => {
    expect(() => loadStandaloneConfig({ DVBFIXER_PORT: 'nope' }, temp())).toThrow(/DVBFIXER_PORT/)
    expect(() => loadStandaloneConfig({ DVBFIXER_PORT: '65536' }, temp())).toThrow(/65535/)
    expect(() => loadStandaloneConfig({ DVBFIXER_HOST: '0.0.0.0' }, temp())).toThrow(/INSECURE_REMOTE/)
    expect(() => loadStandaloneConfig({
      DVBFIXER_HOST: '0.0.0.0', DVBFIXER_ALLOW_INSECURE_REMOTE: '1',
    }, temp())).toThrow(/AUTH_PRINCIPALS/)
    expect(loadStandaloneConfig({
      DVBFIXER_HOST: '0.0.0.0', DVBFIXER_ALLOW_INSECURE_REMOTE: '1',
      DVBFIXER_AUTH_PRINCIPALS: JSON.stringify({
        version: 1, principals: [{ id: 'owner', tokenSha256: 'a'.repeat(64) }],
      }),
    }, temp()).host).toBe('0.0.0.0')
  })

  it('validates CORS and resource limits before startup', () => {
    expect(() => loadStandaloneConfig({
      DVBFIXER_CORS_ALLOWED_ORIGINS: '',
    }, temp())).toThrow(/must not be blank/)
    expect(() => loadStandaloneConfig({
      DVBFIXER_GUI_MAX_UPLOAD_BYTES: '0',
    }, temp())).toThrow(/positive safe integer/)
    expect(() => loadStandaloneConfig({
      DVBFIXER_MAX_CONCURRENT_PROCESSES: '65',
    }, temp())).toThrow(/at most 64/)
    expect(() => loadStandaloneConfig({
      DVBFIXER_MAX_QUEUED_PROCESSES: '257',
    }, temp())).toThrow(/at most 256/)
    expect(() => loadStandaloneConfig({
      DVBFIXER_RATE_LIMIT_WINDOW_MS: '999',
    }, temp())).toThrow(/between 1000 and 3600000/)
    expect(() => loadStandaloneConfig({
      DVBFIXER_OS_RESOURCE_LIMITS_REQUIRED: '1',
    }, temp())).toThrow(/TMPDIR/)
    expect(() => loadStandaloneConfig({ DVBFIXER_ACCESS_LOG: 'text' }, temp())).toThrow(/off or json/)
    expect(() => loadStandaloneConfig({ DVBFIXER_METRICS: 'json' }, temp())).toThrow(/off or prometheus/)
  })

  it('fails before listening when the static client is absent', () => {
    const config = fixture()
    fs.unlinkSync(path.join(config.staticRoot, 'index.html'))
    expect(() => createStandaloneApplication(config)).toThrow(/build:client/)
  })

  it('rejects static roots that expose workspace data', () => {
    const config = fixture()
    expect(() => createStandaloneApplication({
      ...config, staticRoot: config.projectRoot,
    })).toThrow(/must not overlap/)
    expect(() => createStandaloneApplication({
      ...config, dataRoot: path.join(config.staticRoot, 'workspaces'),
    })).toThrow(/must not overlap/)
    const alias = path.join(config.projectRoot, 'static-alias')
    fs.symlinkSync(config.staticRoot, alias)
    expect(() => createStandaloneApplication({
      ...config, dataRoot: path.join(alias, 'future-workspaces'),
    })).toThrow(/symlinks/)
  })
})

describe('standalone HTTP server', () => {
  it('serves protected process-local Prometheus metrics when enabled', async () => {
    const token = Buffer.alloc(32, 7).toString('base64url')
    const config = fixture()
    config.metrics = { enabled: true }
    config.authConfig = parseAuthConfig(JSON.stringify({
      version: 1,
      principals: [{ id: 'operator', tokenSha256: crypto.createHash('sha256').update(token).digest('hex') }],
    }))
    config.legacyWorkspaceOwner = 'operator'
    const instance = createStandaloneServer(config)
    const address = await instance.start()
    const base = `http://127.0.0.1:${address.port}`
    try {
      expect((await fetch(`${base}/api/metrics`)).status).toBe(401)
      await fetch(`${base}/api/health`)
      const response = await fetch(`${base}/api/metrics`, {
        headers: { Authorization: `Bearer ${token}` },
      })
      expect(response.status).toBe(200)
      expect(response.headers.get('content-type')).toBe('text/plain; version=0.0.4; charset=utf-8')
      expect(response.headers.get('cache-control')).toBe('no-store')
      const body = await response.text()
      expect(body).toContain('dvbfixer_http_requests_total')
      expect(body).toContain('route="/api/health"')
      expect(body).not.toContain(token)
      expect(body).not.toContain('operator')
    } finally {
      await instance.close()
    }
  })

  it('emits a redacted access record using the ingress route label', async () => {
    const config = fixture()
    config.accessLog = { enabled: true }
    const lines: string[] = []
    const log = vi.spyOn(console, 'log').mockImplementation(line => lines.push(String(line)))
    const instance = createStandaloneServer(config)
    const address = await instance.start()
    try {
      const response = await fetch(`http://127.0.0.1:${address.port}/api/health?secret=value`)
      expect(response.status).toBe(200)
      await new Promise(resolve => setImmediate(resolve))
      const records = lines.map(line => { try { return JSON.parse(line) } catch { return null } })
        .filter(record => record?.schema === 'dvbfixer.http_access.v1')
      expect(records).toHaveLength(1)
      expect(records[0]).toMatchObject({ route: '/api/health', status: 200, outcome: 'completed' })
      expect(lines.join('\n')).not.toContain('secret=value')
    } finally {
      await instance.close()
      log.mockRestore()
    }
  })

  it('rate-limits by direct client address before authentication', async () => {
    const config = fixture()
    config.rateLimit = { requests: 2, windowMs: 60_000, maxKeys: 10 }
    const instance = createStandaloneServer(config)
    const address = await instance.start()
    const base = `http://127.0.0.1:${address.port}`
    try {
      const denied = await fetch(`${base}/api/health`, {
        headers: { Origin: 'https://evil.example' },
      })
      expect(denied.status).toBe(403)
      expect(denied.headers.get('x-request-id')).toMatch(/^req_/)
      const healthy = await fetch(`${base}/api/health`, { headers: { 'X-Request-Id': 'client-value' } })
      expect(healthy.status).toBe(200)
      expect(healthy.headers.get('x-request-id')).toMatch(/^req_/)
      expect(healthy.headers.get('x-request-id')).not.toBe('client-value')
      const limited = await fetch(`${base}/api/health`)
      expect(limited.status).toBe(429)
      expect(limited.headers.get('x-request-id')).toMatch(/^req_/)
      expect(limited.headers.get('retry-after')).toBe('60')
      expect(await limited.json()).toEqual({ error: 'rate limit exceeded' })
    } finally {
      await instance.close()
    }
  })

  it('enforces browser origins before authentication and handles preflight', async () => {
    const config = fixture()
    config.corsAllowedOrigins = ['https://client.example']
    const instance = createStandaloneServer(config)
    const address = await instance.start()
    const base = `http://127.0.0.1:${address.port}`
    try {
      const denied = await fetch(`${base}/api/health`, {
        headers: { Origin: 'https://evil.example' },
      })
      expect(denied.status).toBe(403)
      const preflight = await fetch(`${base}/api/workspaces`, {
        method: 'OPTIONS',
        headers: {
          Origin: 'https://client.example',
          'Access-Control-Request-Method': 'POST',
          'Access-Control-Request-Headers': 'Authorization, Content-Type',
        },
      })
      expect(preflight.status).toBe(204)
      expect(preflight.headers.get('access-control-allow-origin')).toBe('https://client.example')
      expect(preflight.headers.get('access-control-allow-methods')).toBe('POST')
      const allowed = await fetch(`${base}/api/health`, { headers: { Origin: 'https://client.example' } })
      expect(allowed.headers.get('access-control-expose-headers')).toContain('X-Request-Id')
    } finally {
      await instance.close()
    }
  })

  it('enforces workspace ownership and ACLs per authenticated principal', async () => {
    const aliceToken = Buffer.alloc(32, 8).toString('base64url')
    const bobToken = Buffer.alloc(32, 9).toString('base64url')
    const config = fixture()
    config.authConfig = parseAuthConfig(JSON.stringify({
      version: 1,
      principals: [
        { id: 'alice', tokenSha256: crypto.createHash('sha256').update(aliceToken).digest('hex') },
        { id: 'bob', tokenSha256: crypto.createHash('sha256').update(bobToken).digest('hex') },
      ],
    }))
    config.legacyWorkspaceOwner = 'alice'
    const instance = createStandaloneServer(config)
    const address = await instance.start()
    const base = `http://127.0.0.1:${address.port}`
    const headers = (token: string) => ({
      Authorization: `Bearer ${token}`, 'Content-Type': 'application/json',
    })

    const createdResponse = await fetch(`${base}/api/workspaces`, {
      method: 'POST', headers: headers(aliceToken), body: JSON.stringify({ name: 'Alice workspace' }),
    })
    expect(createdResponse.status).toBe(201)
    const created = await createdResponse.json() as { id: string; revision: number; effectiveRole: string }
    expect(created.effectiveRole).toBe('owner')

    const bobList = await fetch(`${base}/api/workspaces`, { headers: headers(bobToken) })
    expect(await bobList.json()).toEqual([])
    expect((await fetch(`${base}/api/workspaces/${created.id}`, { headers: headers(bobToken) })).status).toBe(404)

    const acl = await fetch(`${base}/api/workspaces/${created.id}/acl`, {
      method: 'PATCH', headers: headers(aliceToken),
      body: JSON.stringify({ revision: created.revision, acl: [{ principalId: 'bob', role: 'reader' }] }),
    })
    expect(acl.status).toBe(200)
    const shared = await acl.json() as { revision: number }
    const bobRead = await fetch(`${base}/api/workspaces/${created.id}`, { headers: headers(bobToken) })
    expect(bobRead.status).toBe(200)
    expect(await bobRead.json()).toMatchObject({ effectiveRole: 'reader' })
    const bobWrite = await fetch(`${base}/api/workspaces/${created.id}`, {
      method: 'PATCH', headers: headers(bobToken),
      body: JSON.stringify({ revision: shared.revision, toolState: {} }),
    })
    expect(bobWrite.status).toBe(403)

    const bobJob = await fetch(`${base}/api/v1/workspaces/${created.id}/jobs`, {
      method: 'POST', headers: headers(bobToken), body: JSON.stringify({ command: 'diagnose' }),
    })
    expect(bobJob.status).toBe(403)
    expect(await bobJob.json()).toMatchObject({
      error: { code: expect.any(String), requestId: expect.stringMatching(/^req_/) },
    })

    await instance.close()
  })

  it('requires a configured bearer token for protected APIs', async () => {
    const token = Buffer.alloc(32, 7).toString('base64url')
    const config = fixture()
    config.authConfig = parseAuthConfig(JSON.stringify({
      version: 1,
      principals: [{ id: 'alice', tokenSha256: crypto.createHash('sha256').update(token).digest('hex') }],
    }))
    config.legacyWorkspaceOwner = 'alice'
    const instance = createStandaloneServer(config)
    const address = await instance.start()
    const base = `http://127.0.0.1:${address.port}`

    const health = await fetch(`${base}/api/health`)
    expect(await health.json()).toMatchObject({ authenticationRequired: true })
    const missing = await fetch(`${base}/api/session`)
    expect(missing.status).toBe(401)
    expect(missing.headers.get('www-authenticate')).toContain('Bearer')
    const missingV1 = await fetch(`${base}/api/v1/workspaces/private/jobs`)
    expect(missingV1.status).toBe(401)
    expect(await missingV1.json()).toMatchObject({
      error: { code: 'UNAUTHORIZED', requestId: expect.stringMatching(/^req_/) },
    })
    const session = await fetch(`${base}/api/session`, {
      headers: { Authorization: `Bearer ${token}` },
    })
    expect(session.status).toBe(200)
    expect(await session.json()).toEqual({ principal: { id: 'alice' } })

    await instance.close()
  })

  it('serves the client and composed APIs and closes idempotently', async () => {
    const instance = createStandaloneServer(fixture())
    const address = await instance.start()
    const base = `http://127.0.0.1:${address.port}`

    const health = await fetch(`${base}/api/health`)
    expect(health.status).toBe(200)
    expect(await health.json()).toMatchObject({ status: 'ready', version: '9.9.9' })

    const openApi = await fetch(`${base}/api/v1/openapi.json`)
    expect(openApi.status).toBe(200)
    expect(await openApi.json()).toMatchObject({ openapi: '3.1.0' })

    const unknownV1 = await fetch(`${base}/api/v1/unknown`)
    expect(unknownV1.status).toBe(404)
    expect(await unknownV1.json()).toMatchObject({
      error: { code: 'ROUTE_NOT_FOUND', requestId: expect.stringMatching(/^req_/) },
    })

    const reboundStatus = await new Promise<number | undefined>((resolve, reject) => {
      const request = http.request({
        host: '127.0.0.1', port: address.port, path: '/api/health',
        headers: { Host: 'attacker.example' },
      }, response => {
        response.resume()
        response.on('end', () => resolve(response.statusCode))
      })
      request.on('error', reject)
      request.end()
    })
    expect(reboundStatus).toBe(403)

    const asset = await fetch(`${base}/assets/app.js`)
    expect(asset.status).toBe(200)
    expect(asset.headers.get('cache-control')).toContain('immutable')

    const shell = await fetch(`${base}/workspace/example`, { headers: { Accept: 'text/html' } })
    expect(shell.status).toBe(200)
    expect(await shell.text()).toContain('<title>DVBFixer</title>')

    const missingApi = await fetch(`${base}/api/unknown`)
    expect(missingApi.status).toBe(404)
    expect(missingApi.headers.get('content-type')).toContain('application/json')

    await instance.close()
    await instance.close()
    expect(instance.server.listening).toBe(false)
  })
})
