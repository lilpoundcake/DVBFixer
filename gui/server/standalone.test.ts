import { afterEach, describe, expect, it } from 'vitest'
import fs from 'node:fs'
import http from 'node:http'
import os from 'node:os'
import path from 'node:path'
import {
  createStandaloneApplication,
  createStandaloneServer,
  loadStandaloneConfig,
  type StandaloneConfig,
} from './standalone'

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
    expect(loadStandaloneConfig({
      DVBFIXER_HOST: '0.0.0.0', DVBFIXER_ALLOW_INSECURE_REMOTE: '1',
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
