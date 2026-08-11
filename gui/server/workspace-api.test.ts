import { afterEach, describe, expect, it } from 'vitest'
import fs from 'node:fs'
import os from 'node:os'
import path from 'node:path'
import { Readable } from 'node:stream'
import {
  applyArtifactMetadataPatch, applyWorkspacePatch, artifactResponseHeaders, assertWorkspaceRevision,
  ensureRetiredWorkspaceIndexMigration, ensureWorkspaceMetadataMigration, ensureWorkspaceMigration, listWorkspaces, loadWorkspace,
  mergeClientWorkspaceUpdate, resolveWorkspaceFile, saveWorkspace,
  registerWorkspaceApi, WorkspaceRevisionConflictError, writeJsonAtomic, type WorkspaceManifest,
} from './workspace-api'

const temporaryDirectories: string[] = []
function temp(): string {
  const directory = fs.mkdtempSync(path.join(os.tmpdir(), 'dvbfixer-workspace-test-'))
  temporaryDirectories.push(directory)
  return directory
}

async function apiRequest(root: string, method: string, url: string, body: unknown) {
  let middleware: ((req: any, res: any, next: () => void) => Promise<void>) | undefined
  const server = { middlewares: { use: (_route: string, handler: typeof middleware) => { middleware = handler } } }
  registerWorkspaceApi(server as any, root)
  const request = Readable.from([JSON.stringify(body)]) as Readable & {
    method: string; url: string; headers: Record<string, string>
  }
  request.method = method
  request.url = url
  request.headers = { 'content-type': 'application/json' }
  let responseBody = ''
  let finish: (() => void) | undefined
  const ended = new Promise<void>(resolve => { finish = resolve })
  const response = {
    statusCode: 0,
    setHeader: () => {},
    end: (value?: unknown) => { responseBody = value === undefined ? '' : String(value); finish?.() },
  }
  await middleware!(request as any, response as any, () => {})
  await ended
  return { status: response.statusCode, body: responseBody ? JSON.parse(responseBody) : null }
}
afterEach(() => temporaryDirectories.splice(0).forEach(directory => fs.rmSync(directory, { recursive: true, force: true })))

describe('workspace storage', () => {
  it('migrates top-level folders to projects and ungrouped files to Unsorted', () => {
    const root = temp()
    fs.mkdirSync(path.join(root, 'legacy'))
    fs.writeFileSync(path.join(root, 'legacy', 'template.pdb'), 'END\n')
    fs.writeFileSync(path.join(root, 'target.fasta'), '>A\nAGS\n')
    fs.writeFileSync(path.join(root, 'index.json'), JSON.stringify([
      { id: '__root__', kind: 'folder', name: '__root__', children: ['folder-one', 'target.fasta'] },
      { id: 'folder-one', kind: 'folder', name: 'Project One', children: ['legacy/template.pdb'] },
      { id: 'legacy/template', file: 'legacy/template.pdb', name: 'Template', kind: 'structure' },
      { id: 'target', file: 'target.fasta', name: 'Target', kind: 'artifact' },
    ]))
    ensureWorkspaceMigration(root)
    const projects = listWorkspaces(root)
    expect(projects.map(project => project.name).sort()).toEqual(['Project One', 'Unsorted'])
    const project = projects.find(item => item.name === 'Project One')!
    expect(project.artifacts[0].file).toMatch(/^files\//)
    expect(fs.existsSync(resolveWorkspaceFile(root, project.id, project.artifacts[0].file))).toBe(true)
    expect(() => resolveWorkspaceFile(root, project.id, '../index.json')).toThrow(/workspace artifact/)
  })

  it('persists versioned workflow state', () => {
    const root = temp()
    ensureWorkspaceMigration(root)
    const workspace = listWorkspaces(root)[0]
    const saved = saveWorkspace(root, { ...workspace, toolState: { dvbfixer: { inputFile: 'a.pdb' } } })
    expect(loadWorkspace(root, saved.id).toolState).toEqual({ dvbfixer: { inputFile: 'a.pdb' } })
    expect(saved.revision).toBe(workspace.revision + 1)
  })

  it('migrates legacy index metadata onto a uniquely matching artifact once', () => {
    const root = temp()
    fs.mkdirSync(path.join(root, 'projects', 'existing', 'files'), { recursive: true })
    fs.writeFileSync(path.join(root, 'projects', '.migrated-v1'), 'done\n')
    fs.writeFileSync(path.join(root, 'index.json'), JSON.stringify([
      { file: 'legacy/target.pdb', name: 'Legacy target', organism: 'Human', method: 'X-ray',
        resolution: 2.1, description: 'Imported', iggSubtype: 'IgG1', allotype: 'G1m3',
        equivalentChains: [['H', 'A']], starred: true },
    ]))
    const now = new Date().toISOString()
    saveWorkspace(root, {
      version: 1, revision: 0, id: 'existing', name: 'Existing', createdAt: now, updatedAt: now,
      primaryFile: null, secondaryFile: null, toolState: {}, artifacts: [
        { id: 'target', file: 'files/target.pdb', name: 'Target', kind: 'structure' },
      ],
    })
    ensureWorkspaceMetadataMigration(root)
    const artifact = loadWorkspace(root, 'existing').artifacts[0]
    expect(artifact).toMatchObject({
      organism: 'Human', method: 'X-ray', resolution: '2.1', description: 'Imported',
      iggSubtype: 'IgG1', allotype: 'G1m3', equivalentChains: [['H', 'A']],
      legacySourceFile: 'legacy/target.pdb',
    })
    expect(artifact).not.toHaveProperty('starred')
    const revision = loadWorkspace(root, 'existing').revision
    ensureWorkspaceMetadataMigration(root)
    expect(loadWorkspace(root, 'existing').revision).toBe(revision)
  })

  it('imports retired per-workspace index artifacts without deleting the recovery source', () => {
    const root = temp()
    fs.mkdirSync(path.join(root, 'projects', 'existing', 'dvb_prepare'), { recursive: true })
    fs.writeFileSync(path.join(root, 'projects', '.migrated-v1'), 'done\n')
    fs.writeFileSync(path.join(root, 'projects', '.metadata-migrated-v1'), 'done\n')
    fs.writeFileSync(path.join(root, 'projects', 'existing', 'dvb_prepare', 'out.pdb'), 'END\n')
    const now = new Date().toISOString()
    saveWorkspace(root, {
      version: 1, revision: 0, id: 'existing', name: 'Existing', createdAt: now, updatedAt: now,
      primaryFile: null, secondaryFile: null, toolState: {}, artifacts: [],
    })
    const legacyIndex = path.join(root, 'projects', 'existing', 'index.json')
    fs.writeFileSync(legacyIndex, JSON.stringify([{
      file: 'dvb_prepare/out.pdb', name: 'Engineered', parent: 'files/input.pdb',
      command: 'antibody-engineer', _engineerChecksum: 'abc', mutationIds: [1, 2],
      mutationsResolved: 'H:MET1ALA', hasGlycan: false, scheme: 'EU', starred: true,
    }]))

    ensureRetiredWorkspaceIndexMigration(root)
    expect(loadWorkspace(root, 'existing').artifacts[0]).toMatchObject({
      file: 'dvb_prepare/out.pdb', parent: 'files/input.pdb', command: 'antibody-engineer',
      engineerChecksum: 'abc', mutationIds: [1, 2], mutationsResolved: 'H:MET1ALA',
      hasGlycan: false, scheme: 'EU',
    })
    expect(loadWorkspace(root, 'existing').artifacts[0]).not.toHaveProperty('starred')
    expect(fs.existsSync(legacyIndex)).toBe(true)
  })

  it('loads pre-revision manifests at revision zero and advances them on save', () => {
    const root = temp()
    const project = path.join(root, 'projects', 'legacy')
    fs.mkdirSync(project, { recursive: true })
    const now = new Date().toISOString()
    fs.writeFileSync(path.join(project, 'workspace.json'), JSON.stringify({
      version: 1, id: 'legacy', name: 'Legacy', createdAt: now, updatedAt: now,
      primaryFile: null, secondaryFile: null, artifacts: [], toolState: {},
    }))
    const loaded = loadWorkspace(root, 'legacy')
    expect(loaded.revision).toBe(0)
    expect(saveWorkspace(root, loaded).revision).toBe(1)
  })

  it('writes JSON atomically without leaving temporary files', () => {
    const root = temp()
    const file = path.join(root, 'state.json')
    writeJsonAtomic(file, { revision: 1 })
    writeJsonAtomic(file, { revision: 2 })
    expect(JSON.parse(fs.readFileSync(file, 'utf8'))).toEqual({ revision: 2 })
    expect(fs.readdirSync(root)).toEqual(['state.json'])
  })

  it('rejects workspace symlinks that resolve outside the workspace', () => {
    const root = temp()
    ensureWorkspaceMigration(root)
    const workspace = listWorkspaces(root)[0]
    const outside = path.join(root, 'outside.pdb')
    fs.writeFileSync(outside, 'END\n')
    const link = path.join(root, 'projects', workspace.id, 'files', 'outside.pdb')
    try { fs.symlinkSync(outside, link) } catch { return }
    expect(() => resolveWorkspaceFile(root, workspace.id, 'files/outside.pdb')).toThrow(/workspace artifact/)
  })

  it('preserves server-owned artifacts when applying client workspace updates', () => {
    const now = new Date().toISOString()
    const current: WorkspaceManifest = {
      version: 1, revision: 4, id: 'safe', name: 'Safe', createdAt: now, updatedAt: now,
      primaryFile: 'files/a.pdb', secondaryFile: null, toolState: {}, artifacts: [
        { id: 'a', file: 'files/a.pdb', name: 'A', kind: 'structure' },
        { id: 'b', file: 'files/b.txt', name: 'B', kind: 'artifact' },
      ],
    }
    const incoming: WorkspaceManifest = {
      ...current, name: 'Injected', primaryFile: 'workspace.json', toolState: { dvbfixer: { input: 'a' } },
      artifacts: [
        { id: 'evil', file: 'workspace.json', name: 'Manifest', kind: 'artifact' },
        { id: 'b', file: 'files/b.txt', name: 'Changed', kind: 'artifact' },
      ],
    }
    const merged = mergeClientWorkspaceUpdate(current, incoming)
    expect(merged.name).toBe('Safe')
    expect(merged.primaryFile).toBeNull()
    expect(merged.artifacts.map(artifact => [artifact.id, artifact.name])).toEqual([['b', 'B'], ['a', 'A']])
    expect(merged.toolState).toEqual(incoming.toolState)
  })

  it('applies only supported narrow patches and enforces revisions', () => {
    const now = new Date().toISOString()
    const current: WorkspaceManifest = {
      version: 1, revision: 7, id: 'safe', name: 'Safe', createdAt: now, updatedAt: now,
      primaryFile: 'files/a.pdb', secondaryFile: null, toolState: {}, artifacts: [
        { id: 'a', file: 'files/a.pdb', name: 'A', kind: 'structure' },
        { id: 'b', file: 'files/b.pdb', name: 'B', kind: 'structure' },
        { id: 'notes', file: 'files/notes.txt', name: 'Notes', kind: 'artifact' },
      ],
    }
    const patched = applyWorkspacePatch(current, {
      revision: 7, primaryFile: 'files/b.pdb', toolState: { dvbfixer: { input: 'files/b.pdb' } },
      artifactOrder: ['notes', 'b', 'a'],
    })
    expect(patched.primaryFile).toBe('files/b.pdb')
    expect(patched.artifacts.map(artifact => artifact.id)).toEqual(['notes', 'b', 'a'])
    expect(patched.toolState).toEqual({ dvbfixer: { input: 'files/b.pdb' } })
    expect(() => assertWorkspaceRevision(current, 6)).toThrow(WorkspaceRevisionConflictError)
    try { assertWorkspaceRevision(current, 6) } catch (error) {
      expect((error as WorkspaceRevisionConflictError).statusCode).toBe(409)
    }
    expect(() => applyWorkspacePatch(current, { revision: 7, name: 'Nope' } as any)).toThrow(/unsupported workspace patch field/)
    expect(() => applyWorkspacePatch(current, { revision: 7, primaryFile: 'files/notes.txt' })).toThrow(/structure artifact/)
  })

  it('validates and applies revision-safe artifact metadata patches', async () => {
    const root = temp()
    ensureWorkspaceMigration(root)
    const workspace = listWorkspaces(root)[0]
    const seeded = saveWorkspace(root, {
      ...workspace,
      artifacts: [{ id: 'structure', file: 'files/a.pdb', name: 'A', kind: 'structure' }],
    })
    const patched = applyArtifactMetadataPatch(seeded, 'structure', {
      revision: seeded.revision, name: 'Renamed', organism: 'Human', resolution: 1.8,
      equivalentChains: [['H', 'A']],
    })
    expect(patched.artifacts[0]).toMatchObject({
      name: 'Renamed', organism: 'Human', resolution: '1.8', equivalentChains: [['H', 'A']],
    })
    expect(() => applyArtifactMetadataPatch(seeded, 'structure', {
      revision: seeded.revision - 1, organism: 'Mouse',
    })).toThrow(WorkspaceRevisionConflictError)
    expect(() => applyArtifactMetadataPatch(seeded, 'structure', {
      revision: seeded.revision, equivalentChains: ['bad'] as any,
    })).toThrow(/equivalentChains/)

    const response = await apiRequest(root, 'PATCH', `/${seeded.id}/artifacts/structure/metadata`, {
      revision: seeded.revision, description: 'Updated through API', allotype: null,
    })
    expect(response.status).toBe(200)
    expect(response.body.revision).toBe(seeded.revision + 1)
    expect(response.body.artifacts[0].description).toBe('Updated through API')
  })

  it('returns 409 for stale PUT/PATCH and advances a successful PATCH revision', async () => {
    const root = temp()
    ensureWorkspaceMigration(root)
    const workspace = listWorkspaces(root)[0]
    const stale = await apiRequest(root, 'PATCH', `/${workspace.id}`, {
      revision: workspace.revision - 1, toolState: { stale: true },
    })
    expect(stale.status).toBe(409)
    expect(stale.body.error).toMatch(/revision conflict/)

    const updated = await apiRequest(root, 'PATCH', `/${workspace.id}`, {
      revision: workspace.revision, toolState: { dvbfixer: { input: 'a.pdb' } },
    })
    expect(updated.status).toBe(200)
    expect(updated.body.revision).toBe(workspace.revision + 1)

    const stalePut = await apiRequest(root, 'PUT', `/${workspace.id}`, workspace)
    expect(stalePut.status).toBe(409)
  })

  it('drops protected control files from malformed persisted artifact lists', () => {
    const root = temp()
    ensureWorkspaceMigration(root)
    const workspace = listWorkspaces(root)[0]
    saveWorkspace(root, {
      ...workspace,
      artifacts: [{ id: 'evil', file: 'workspace.json', name: 'Manifest', kind: 'artifact' }],
    })
    expect(loadWorkspace(root, workspace.id).artifacts).toEqual([])
  })

  it('forces active content to download with restrictive headers', () => {
    expect(artifactResponseHeaders('report.html')).toMatchObject({
      'Content-Type': 'text/html',
      'Content-Disposition': 'attachment; filename="report.html"',
      'Content-Security-Policy': "sandbox; default-src 'none'",
      'X-Content-Type-Options': 'nosniff',
    })
    expect(artifactResponseHeaders('model.pdb')['Content-Disposition']).toBe('inline; filename="model.pdb"')
  })
})
