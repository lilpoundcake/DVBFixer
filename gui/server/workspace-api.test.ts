import { afterEach, describe, expect, it, vi } from 'vitest'
import fs from 'node:fs'
import os from 'node:os'
import path from 'node:path'
import { Readable } from 'node:stream'
import {
  DEFAULT_MAX_UPLOAD_BYTES, DEFAULT_WORKSPACE_QUOTA_BYTES, initializeStorageQuota,
  workspaceLogicalBytes,
} from './storage-quota'
import {
  applyArtifactMetadataPatch, applyWorkspaceAclPatch, applyWorkspacePatch, artifactResponseHeaders,
  assertWorkspaceAccess, assertWorkspaceRevision, canReadWorkspace, canWriteWorkspace, createWorkspace,
  deleteWorkspace, deleteWorkspaceArtifact,
  ensureRetiredWorkspaceIndexMigration, ensureWorkspaceMetadataMigration, ensureWorkspaceMigration, isWorkspaceOwner,
  listWorkspaces, loadWorkspace, mergeClientWorkspaceUpdate, migrateWorkspaceOwnership, reorderWorkspaces,
  resolveWorkspaceFile, saveWorkspace,
  registerWorkspaceApi, WorkspaceAuthorizationError, WorkspaceRevisionConflictError, workspaceRole,
  workspaceOrderPath, workspaceRoot, writeJsonAtomic, type WorkspaceManifest,
} from './workspace-api'

const temporaryDirectories: string[] = []
function temp(): string {
  const directory = fs.mkdtempSync(path.join(os.tmpdir(), 'dvbfixer-workspace-test-'))
  temporaryDirectories.push(directory)
  return directory
}

async function apiRequest(root: string, method: string, url: string, body: unknown, principalId = 'local') {
  let middleware: ((req: any, res: any, next: () => void) => Promise<void>) | undefined
  const server = { middlewares: { use: (_route: string, handler: typeof middleware) => { middleware = handler } } }
  registerWorkspaceApi(server as any, root, principalId)
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

async function importRequest(
  root: string,
  workspaceId: string,
  content: Buffer,
  quotaOptions: { maxUploadBytes: number; workspaceQuotaBytes: number },
) {
  initializeStorageQuota(quotaOptions)
  let middleware: ((req: any, res: any, next: () => void) => Promise<void>) | undefined
  const server = { middlewares: { use: (_route: string, handler: typeof middleware) => { middleware = handler } } }
  registerWorkspaceApi(server as any, root, 'local', 'local', false)
  const request = Readable.from([content]) as Readable & {
    method: string; url: string; headers: Record<string, string>
  }
  request.method = 'POST'
  request.url = `/${workspaceId}/import`
  request.headers = {
    'content-length': String(content.length),
    'content-type': 'application/octet-stream',
    'x-file-name': encodeURIComponent('uploaded.pdb'),
  }
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

afterEach(() => {
  vi.restoreAllMocks()
  initializeStorageQuota({
    maxUploadBytes: DEFAULT_MAX_UPLOAD_BYTES,
    workspaceQuotaBytes: DEFAULT_WORKSPACE_QUOTA_BYTES,
  })
  temporaryDirectories.splice(0).forEach(directory => fs.rmSync(directory, { recursive: true, force: true }))
})

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

  it('transfers legacy-index workspaces created in local mode when authentication is enabled', () => {
    const root = temp()
    fs.writeFileSync(path.join(root, 'input.pdb'), 'END\n')
    fs.writeFileSync(path.join(root, 'index.json'), JSON.stringify([
      { id: '__root__', kind: 'folder', name: '__root__', children: ['input.pdb'] },
      { id: 'input', file: 'input.pdb', name: 'Input', kind: 'structure' },
    ]))
    ensureWorkspaceMigration(root, 'local', true)
    const local = listWorkspaces(root)[0]
    expect(local).toMatchObject({ ownerPrincipalId: 'local', provisionalOwner: true })

    migrateWorkspaceOwnership(root, 'alice', true)
    expect(loadWorkspace(root, local.id)).toMatchObject({ ownerPrincipalId: 'alice' })
    expect(loadWorkspace(root, local.id).provisionalOwner).toBeUndefined()
  })

  it('persists versioned workflow state', () => {
    const root = temp()
    ensureWorkspaceMigration(root)
    const workspace = listWorkspaces(root)[0]
    const saved = saveWorkspace(root, { ...workspace, toolState: { dvbfixer: { inputFile: 'a.pdb' } } })
    expect(loadWorkspace(root, saved.id).toolState).toEqual({ dvbfixer: { inputFile: 'a.pdb' } })
    expect(saved.revision).toBe(workspace.revision + 1)
  })

  it('compares the on-disk revision while holding the manifest lock', () => {
    const root = temp()
    ensureWorkspaceMigration(root)
    const first = listWorkspaces(root)[0]
    const stale = loadWorkspace(root, first.id)
    const saved = saveWorkspace(root, { ...first, name: 'first writer' })

    expect(() => saveWorkspace(root, { ...stale, name: 'stale writer' }))
      .toThrow(WorkspaceRevisionConflictError)
    expect(loadWorkspace(root, first.id)).toMatchObject({
      revision: saved.revision,
      name: 'first writer',
    })
  })

  it('checks replacement-aware manifest growth before writing', () => {
    const root = temp()
    ensureWorkspaceMigration(root)
    const workspace = listWorkspaces(root)[0]
    const project = path.join(root, 'projects', workspace.id)
    const before = workspaceLogicalBytes(project)
    initializeStorageQuota({ maxUploadBytes: 1024, workspaceQuotaBytes: before + 16 })
    expect(() => saveWorkspace(root, { ...workspace, toolState: { oversized: 'x'.repeat(1000) } }))
      .toThrow(/workspace quota exceeded/)
    expect(loadWorkspace(root, workspace.id).revision).toBe(workspace.revision)
  })

  it('distinguishes upload body limits from workspace quota failures without publishing files', async () => {
    const root = temp()
    ensureWorkspaceMigration(root)
    const workspace = listWorkspaces(root)[0]
    const project = path.join(root, 'projects', workspace.id)
    const originalRevision = workspace.revision

    const tooLarge = await importRequest(root, workspace.id, Buffer.from('12345'), {
      maxUploadBytes: 4,
      workspaceQuotaBytes: DEFAULT_WORKSPACE_QUOTA_BYTES,
    })
    expect(tooLarge.status).toBe(413)

    const usedBytes = workspaceLogicalBytes(project)
    const overQuota = await importRequest(root, workspace.id, Buffer.from('END\n'), {
      maxUploadBytes: 1024,
      workspaceQuotaBytes: usedBytes + 4,
    })
    expect(overQuota).toMatchObject({
      status: 507,
      body: { code: 'WORKSPACE_QUOTA_EXCEEDED' },
    })
    expect(loadWorkspace(root, workspace.id)).toMatchObject({ revision: originalRevision, artifacts: [] })
    const files = fs.readdirSync(project, { recursive: true, encoding: 'utf8' })
    expect(files.some(file => file.endsWith('uploaded.pdb') || file.endsWith('.tmp'))).toBe(false)
  })

  it('removes the published artifact and temporary files when manifest publication fails', async () => {
    const root = temp()
    ensureWorkspaceMigration(root)
    const workspace = listWorkspaces(root)[0]
    const project = path.join(root, 'projects', workspace.id)
    const renameSync = fs.renameSync
    vi.spyOn(fs, 'renameSync').mockImplementation((oldPath, newPath) => {
      if (path.resolve(String(newPath)) === path.join(project, 'workspace.json')) {
        throw new Error('simulated manifest publication failure')
      }
      return renameSync(oldPath, newPath)
    })

    const response = await importRequest(root, workspace.id, Buffer.from('END\n'), {
      maxUploadBytes: 1024,
      workspaceQuotaBytes: DEFAULT_WORKSPACE_QUOTA_BYTES,
    })
    expect(response).toMatchObject({ status: 500, body: { error: 'simulated manifest publication failure' } })
    expect(loadWorkspace(root, workspace.id)).toMatchObject({ revision: workspace.revision, artifacts: [] })
    const files = fs.readdirSync(project, { recursive: true, encoding: 'utf8' })
    expect(files.some(file => file.endsWith('uploaded.pdb') || file.endsWith('.tmp'))).toBe(false)
  })

  it('restores an artifact when locked deletion cannot publish its manifest', () => {
    const root = temp()
    const created = createWorkspace(root, 'Delete rollback')
    const project = workspaceRoot(root, created.id)
    const source = path.join(project, 'files', 'input.pdb')
    fs.writeFileSync(source, 'END\n')
    const workspace = saveWorkspace(root, {
      ...created,
      artifacts: [{ id: 'input', file: 'files/input.pdb', name: 'input.pdb', kind: 'structure' }],
    })
    const renameSync = fs.renameSync
    vi.spyOn(fs, 'renameSync').mockImplementation((oldPath, newPath) => {
      if (path.resolve(String(newPath)) === path.join(project, 'workspace.json')) {
        throw new Error('simulated delete publication failure')
      }
      return renameSync(oldPath, newPath)
    })

    expect(() => deleteWorkspaceArtifact(root, created.id, 'input', 'local'))
      .toThrow('simulated delete publication failure')
    expect(fs.readFileSync(source, 'utf8')).toBe('END\n')
    expect(loadWorkspace(root, created.id)).toMatchObject({
      revision: workspace.revision,
      artifacts: [{ id: 'input' }],
    })
  })

  it('does not let a stale revision-zero manifest resurrect a deleted workspace', () => {
    const root = temp()
    const created = createWorkspace(root, 'Deleted')
    const stale = { ...created, revision: 0 }

    deleteWorkspace(root, created.id, 'local')

    expect(() => saveWorkspace(root, stale)).toThrow(WorkspaceRevisionConflictError)
    expect(fs.existsSync(path.join(root, 'projects', created.id))).toBe(false)
  })

  it('finishes a committed workspace deletion interrupted before its move', () => {
    const root = temp()
    const created = createWorkspace(root, 'Interrupted delete')
    const trashName = `${created.id}-pending-delete`
    writeJsonAtomic(path.join(root, '.workspace-deletions', `${created.id}.json`), {
      workspaceId: created.id,
      trashName,
      deletedAt: new Date().toISOString(),
    })

    registerWorkspaceApi({ middlewares: { use: () => {} } } as any, root)

    expect(fs.existsSync(path.join(root, 'projects', created.id))).toBe(false)
    expect(fs.existsSync(path.join(root, '_workspace_trash', trashName, 'workspace.json'))).toBe(true)
  })

  it('quarantines directories recreated by stale work after workspace deletion', () => {
    const root = temp()
    const created = createWorkspace(root, 'Late recreation')
    deleteWorkspace(root, created.id, 'local')
    const lateFile = path.join(root, 'projects', created.id, 'files', 'imports', 'late.pdb')
    fs.mkdirSync(path.dirname(lateFile), { recursive: true })
    fs.writeFileSync(lateFile, 'END\n')

    expect(() => registerWorkspaceApi({ middlewares: { use: () => {} } } as any, root)).not.toThrow()

    expect(fs.existsSync(path.join(root, 'projects', created.id))).toBe(false)
    const trashEntries = fs.readdirSync(path.join(root, '_workspace_trash'))
      .filter(entry => entry.startsWith(`${created.id}-`))
    expect(trashEntries).toHaveLength(2)
    expect(trashEntries.some(entry => fs.existsSync(
      path.join(root, '_workspace_trash', entry, 'files', 'imports', 'late.pdb'),
    ))).toBe(true)
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

  it('atomically migrates ownerless v1 manifests to v2 without changing revision or timestamps', () => {
    const root = temp()
    const project = path.join(root, 'projects', 'legacy')
    fs.mkdirSync(project, { recursive: true })
    const createdAt = '2024-01-02T03:04:05.000Z'
    const updatedAt = '2024-06-07T08:09:10.000Z'
    fs.writeFileSync(path.join(project, 'workspace.json'), JSON.stringify({
      version: 1, revision: 17, id: 'legacy', name: 'Legacy', createdAt, updatedAt,
      primaryFile: null, secondaryFile: null, artifacts: [], toolState: {},
    }))

    const loaded = loadWorkspace(root, 'legacy', 'legacy-owner')
    expect(loaded).toMatchObject({
      version: 2, revision: 17, createdAt, updatedAt, ownerPrincipalId: 'legacy-owner', acl: [],
    })
    expect(JSON.parse(fs.readFileSync(path.join(project, 'workspace.json'), 'utf8'))).toMatchObject({
      version: 2, revision: 17, createdAt, updatedAt, ownerPrincipalId: 'legacy-owner', acl: [],
    })
    expect(fs.readdirSync(project)).toEqual(['workspace.json'])
  })

  it('transfers only the provisional local owner when configured authentication is enabled', () => {
    const root = temp()
    const local = createWorkspace(root, 'Local', 'local', true)
    const intentionalLocal = createWorkspace(root, 'Intentional local', 'local')
    const explicit = createWorkspace(root, 'Explicit', 'bob')
    migrateWorkspaceOwnership(root, 'alice', true)
    expect(loadWorkspace(root, local.id).ownerPrincipalId).toBe('alice')
    expect(loadWorkspace(root, local.id).provisionalOwner).toBeUndefined()
    expect(loadWorkspace(root, intentionalLocal.id).ownerPrincipalId).toBe('local')
    expect(loadWorkspace(root, explicit.id).ownerPrincipalId).toBe('bob')
  })

  it('rejects malformed v2 security metadata without repairing or broadening access', () => {
    const root = temp()
    const project = path.join(root, 'projects', 'invalid')
    fs.mkdirSync(project, { recursive: true })
    const now = new Date().toISOString()
    const base = {
      version: 2, revision: 1, id: 'invalid', name: 'Invalid', createdAt: now, updatedAt: now,
      primaryFile: null, secondaryFile: null, artifacts: [], toolState: {}, ownerPrincipalId: 'owner', acl: [],
    }
    const write = (security: Record<string, unknown>) => {
      fs.writeFileSync(path.join(project, 'workspace.json'), JSON.stringify({ ...base, ...security }))
    }

    write({ ownerPrincipalId: undefined })
    expect(() => loadWorkspace(root, 'invalid')).toThrow(/ownerPrincipalId/)
    write({ acl: [{ principalId: 'reader', role: 'admin' }] })
    expect(() => loadWorkspace(root, 'invalid')).toThrow(/acl entry/)
    write({ acl: [{ principalId: 'reader', role: 'reader' }, { principalId: 'reader', role: 'writer' }] })
    expect(() => loadWorkspace(root, 'invalid')).toThrow(/duplicate/)
    write({ acl: [{ principalId: 'owner', role: 'writer' }] })
    expect(() => loadWorkspace(root, 'invalid')).toThrow(/owner must not appear/)
    write({ provisionalOwner: true })
    expect(() => loadWorkspace(root, 'invalid')).toThrow(/provisional workspace owner/)
    expect(() => migrateWorkspaceOwnership(root, 'alice', true)).toThrow(/provisional workspace owner/)
    fs.writeFileSync(path.join(project, 'workspace.json'), JSON.stringify({
      ...base, version: 1, ownerPrincipalId: undefined, acl: undefined, provisionalOwner: true,
    }))
    expect(() => loadWorkspace(root, 'invalid')).toThrow(/security fields/)
  })

  it('uses case-sensitive principals and distinguishes hidden from insufficient access', () => {
    const root = temp()
    const workspace = saveWorkspace(root, {
      version: 2, revision: 0, id: 'secured', name: 'Secured', createdAt: new Date().toISOString(),
      updatedAt: new Date().toISOString(), primaryFile: null, secondaryFile: null, artifacts: [], toolState: {},
      ownerPrincipalId: 'Owner', acl: [
        { principalId: 'Reader', role: 'reader' },
        { principalId: 'reader', role: 'writer' },
      ],
    })

    expect(workspaceRole(workspace, 'Owner')).toBe('owner')
    expect(workspaceRole(workspace, 'owner')).toBeNull()
    expect(canReadWorkspace(workspace, 'Reader')).toBe(true)
    expect(canWriteWorkspace(workspace, 'Reader')).toBe(false)
    expect(canWriteWorkspace(workspace, 'reader')).toBe(true)
    expect(isWorkspaceOwner(workspace, 'Owner')).toBe(true)
    expect(() => assertWorkspaceAccess(workspace, 'Reader', 'writer')).toThrow(WorkspaceAuthorizationError)
    try { assertWorkspaceAccess(workspace, 'Reader', 'writer') } catch (error) {
      expect((error as WorkspaceAuthorizationError).statusCode).toBe(403)
    }
    try { assertWorkspaceAccess(workspace, 'missing') } catch (error) {
      expect((error as WorkspaceAuthorizationError).statusCode).toBe(404)
    }
  })

  it('filters listings and preserves inaccessible workspace positions during visible reordering', () => {
    const root = temp()
    fs.mkdirSync(path.join(root, 'projects'), { recursive: true })
    fs.writeFileSync(path.join(root, 'projects', '.migrated-v1'), 'done\n')
    fs.writeFileSync(path.join(root, 'projects', '.metadata-migrated-v1'), 'done\n')
    fs.writeFileSync(path.join(root, 'projects', '.artifact-index-migrated-v1'), 'done\n')
    const alice = createWorkspace(root, 'Alice', 'alice')
    const shared = saveWorkspace(root, {
      ...createWorkspace(root, 'Shared', 'bob'), acl: [{ principalId: 'alice', role: 'reader' }],
    })
    const hidden = createWorkspace(root, 'Hidden', 'carol')
    fs.writeFileSync(workspaceOrderPath(root, 'alice'), JSON.stringify([hidden.id, alice.id, shared.id]))

    expect(listWorkspaces(root, 'alice').map(workspace => workspace.id)).toEqual([alice.id, shared.id])
    expect(listWorkspaces(root, 'carol').map(workspace => workspace.id)).toEqual([hidden.id])
    expect(reorderWorkspaces(root, [shared.id, alice.id], 'alice')).toEqual([shared.id, alice.id])
    expect(JSON.parse(fs.readFileSync(workspaceOrderPath(root, 'alice'), 'utf8'))).toEqual([
      hidden.id, shared.id, alice.id,
    ])
    expect(listWorkspaces(root, 'alice').map(workspace => workspace.id)).toEqual([shared.id, alice.id])
  })

  it('applies revisioned ACL replacement only for the owner and enforces it through the route', async () => {
    const root = temp()
    const workspace = saveWorkspace(root, {
      ...createWorkspace(root, 'Secured', 'owner'), acl: [{ principalId: 'reader', role: 'reader' }],
    })
    expect(() => applyWorkspaceAclPatch(workspace, 'reader', {
      revision: workspace.revision, acl: [],
    })).toThrow(WorkspaceAuthorizationError)

    const forbidden = await apiRequest(root, 'PATCH', `/${workspace.id}/acl`, {
      revision: workspace.revision, acl: [],
    }, 'reader')
    expect(forbidden.status).toBe(403)
    const hidden = await apiRequest(root, 'GET', `/${workspace.id}`, {}, 'outsider')
    expect(hidden.status).toBe(404)

    const updated = await apiRequest(root, 'PATCH', `/${workspace.id}/acl`, {
      revision: workspace.revision,
      acl: [{ principalId: 'writer', role: 'writer' }],
    }, 'owner')
    expect(updated.status).toBe(200)
    expect(updated.body.revision).toBe(workspace.revision + 1)
    expect(updated.body.acl).toEqual([{ principalId: 'writer', role: 'writer' }])
    expect((await apiRequest(root, 'PATCH', `/${workspace.id}/acl`, {
      revision: workspace.revision, acl: [],
    }, 'owner')).status).toBe(409)
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
