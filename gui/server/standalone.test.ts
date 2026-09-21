import { afterEach, describe, expect, it } from 'vitest'
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
    host: '127.0.0.1', port: 0, projectRoot,
    dataRoot: path.join(projectRoot, 'structures'), staticRoot, shutdownGraceMs: 100,
    authConfig: parseAuthConfig(undefined), legacyWorkspaceOwner: 'local',
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
