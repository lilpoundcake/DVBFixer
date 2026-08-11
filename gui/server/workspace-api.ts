import type { IncomingMessage, ServerResponse } from 'node:http'
import type { ViteDevServer } from 'vite'
import crypto from 'node:crypto'
import fs from 'node:fs'
import path from 'node:path'

export interface WorkspaceArtifact {
  id: string
  file: string
  name: string
  kind: 'structure' | 'artifact'
  folder?: string
  parent?: string
  command?: string
  artifactType?: string
  /** Internal run bookkeeping retained on disk but hidden from user pickers. */
  hidden?: boolean
}

export interface WorkspaceManifest {
  version: 1
  id: string
  name: string
  createdAt: string
  updatedAt: string
  primaryFile: string | null
  secondaryFile: string | null
  artifacts: WorkspaceArtifact[]
  toolState: Record<string, unknown>
}

const ALLOWED = new Set([
  '.pdb', '.cif', '.mmcif', '.gro', '.fasta', '.fa', '.faa', '.pir', '.aln',
  '.xtc', '.trr', '.dcd', '.tpr', '.mol2', '.sdf', '.top', '.itp', '.dat',
  '.json', '.txt', '.html', '.csv', '.log', '.md', '.yaml', '.yml', '.toml',
  '.xml', '.pml', '.py', '.sh', '.mdp',
])

function sendJson(res: ServerResponse, status: number, body: unknown): void {
  res.statusCode = status
  res.setHeader('Content-Type', 'application/json')
  res.setHeader('Cache-Control', 'no-store')
  res.end(JSON.stringify(body))
}

function readBody(req: IncomingMessage): Promise<Buffer> {
  return new Promise((resolve, reject) => {
    const chunks: Buffer[] = []
    req.on('data', chunk => chunks.push(Buffer.from(chunk)))
    req.on('end', () => resolve(Buffer.concat(chunks)))
    req.on('error', reject)
  })
}

function safeId(id: string): string {
  if (!/^[a-zA-Z0-9_-]+$/.test(id)) throw new Error('invalid workspace id')
  return id
}

export function workspacesRoot(dataRoot: string): string {
  const root = path.join(dataRoot, 'projects')
  fs.mkdirSync(root, { recursive: true })
  return root
}

export function workspaceRoot(dataRoot: string, id: string): string {
  return path.join(workspacesRoot(dataRoot), safeId(id))
}

function manifestPath(dataRoot: string, id: string): string {
  return path.join(workspaceRoot(dataRoot, id), 'workspace.json')
}

export function loadWorkspace(dataRoot: string, id: string): WorkspaceManifest {
  const file = manifestPath(dataRoot, id)
  if (!fs.existsSync(file)) throw new Error(`workspace not found: ${id}`)
  return JSON.parse(fs.readFileSync(file, 'utf8')) as WorkspaceManifest
}

export function saveWorkspace(dataRoot: string, manifest: WorkspaceManifest): WorkspaceManifest {
  const root = workspaceRoot(dataRoot, manifest.id)
  fs.mkdirSync(path.join(root, 'files'), { recursive: true })
  fs.mkdirSync(path.join(root, 'runs'), { recursive: true })
  fs.mkdirSync(path.join(root, 'homology'), { recursive: true })
  const next = { ...manifest, version: 1 as const, updatedAt: new Date().toISOString() }
  fs.writeFileSync(manifestPath(dataRoot, manifest.id), `${JSON.stringify(next, null, 2)}\n`)
  return next
}

export function resolveWorkspaceFile(dataRoot: string, workspaceId: string, relative: string): string {
  if (!relative || path.isAbsolute(relative)) throw new Error(`invalid workspace path: ${relative}`)
  const root = workspaceRoot(dataRoot, workspaceId)
  const resolved = path.resolve(root, relative)
  const rel = path.relative(root, resolved)
  if (rel.startsWith('..') || path.isAbsolute(rel) || !fs.existsSync(resolved)) {
    throw new Error(`workspace artifact not found: ${relative}`)
  }
  return resolved
}

function uniqueId(name: string): string {
  const stem = name.toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/^-|-$/g, '') || 'workspace'
  return `${stem}-${crypto.randomUUID().slice(0, 8)}`
}

function artifactFor(file: string, folder?: string): WorkspaceArtifact {
  const extension = path.extname(file).toLowerCase()
  const structure = ['.pdb', '.cif', '.mmcif', '.gro'].includes(extension)
  const basename = path.basename(file)
  return {
    id: crypto.randomUUID(), file: file.replace(/\\/g, '/'),
    name: path.basename(file, extension), kind: structure ? 'structure' : 'artifact', folder,
    hidden: ['workspace.json', 'index.json', 'run.json', 'stdout.log', 'stderr.log'].includes(basename) || basename.startsWith('_'),
  }
}

function createWorkspace(dataRoot: string, name: string): WorkspaceManifest {
  const now = new Date().toISOString()
  return saveWorkspace(dataRoot, {
    version: 1, id: uniqueId(name), name, createdAt: now, updatedAt: now,
    primaryFile: null, secondaryFile: null, artifacts: [], toolState: {},
  })
}

/**
 * One-time compatibility migration. Files are copied before manifests are
 * published, so a failed migration never damages the legacy library.
 */
export function ensureWorkspaceMigration(dataRoot: string): void {
  const marker = path.join(workspacesRoot(dataRoot), '.migrated-v1')
  if (fs.existsSync(marker)) return
  const indexFile = path.join(dataRoot, 'index.json')
  let entries: any[] = []
  if (fs.existsSync(indexFile)) {
    try { entries = JSON.parse(fs.readFileSync(indexFile, 'utf8')) } catch { entries = [] }
  }
  const folders = entries.filter(entry => entry.kind === 'folder')
  const root = folders.find(entry => entry.id === '__root__')
  const byId = new Map(entries.map(entry => [entry.kind === 'folder' ? entry.id : entry.file, entry]))
  const assigned = new Set<string>()

  const collect = (id: string, prefix = ''): Array<{ entry: any; folder: string }> => {
    const entry = byId.get(id)
    if (!entry) return []
    if (entry.kind !== 'folder') return [{ entry, folder: prefix }]
    const nextPrefix = entry.id === '__root__' ? prefix : [prefix, entry.name].filter(Boolean).join('/')
    return (entry.children || []).flatMap((child: string) => collect(child, nextPrefix))
  }

  const groups: Array<{ name: string; files: Array<{ entry: any; folder: string }> }> = []
  for (const child of root?.children || []) {
    const entry = byId.get(child)
    if (entry?.kind === 'folder') groups.push({ name: entry.name, files: (entry.children || []).flatMap((nested: string) => collect(nested, '')) })
  }
  const explicitlyGrouped = new Set(groups.flatMap(group => group.files.map(item => item.entry.file)))
  const remaining = entries.filter(entry => entry.kind !== 'folder' && !explicitlyGrouped.has(entry.file))
    .map(entry => ({ entry, folder: '' }))
  if (remaining.length || groups.length === 0) groups.push({ name: 'Unsorted', files: remaining })

  const ownership = new Map<string, { workspaceId: string; relative: string }>()
  const createdWorkspaceIds: string[] = []
  for (const group of groups) {
    const workspace = createWorkspace(dataRoot, group.name)
    createdWorkspaceIds.push(workspace.id)
    const artifacts: WorkspaceArtifact[] = []
    for (const { entry, folder } of group.files) {
      if (!entry.file || assigned.has(entry.file)) continue
      const source = path.resolve(dataRoot, entry.file)
      if (!source.startsWith(path.resolve(dataRoot)) || !fs.existsSync(source) || !fs.statSync(source).isFile()) continue
      assigned.add(entry.file)
      const relative = path.join('files', folder, path.basename(entry.file)).replace(/\\/g, '/')
      const destination = path.join(workspaceRoot(dataRoot, workspace.id), relative)
      fs.mkdirSync(path.dirname(destination), { recursive: true })
      fs.copyFileSync(source, destination)
      ownership.set(entry.file, { workspaceId: workspace.id, relative })
      artifacts.push({ ...artifactFor(relative, folder || undefined), name: entry.name || artifactFor(relative).name,
        parent: entry.parent, command: entry.command, artifactType: entry.artifactType })
    }
    saveWorkspace(dataRoot, { ...workspace, artifacts, primaryFile: artifacts.find(a => a.kind === 'structure')?.file || null })
  }

  const legacyHomology = path.join(dataRoot, 'homology_projects')
  if (fs.existsSync(legacyHomology)) {
    for (const directory of fs.readdirSync(legacyHomology, { withFileTypes: true }).filter(entry => entry.isDirectory())) {
      const sourceProject = path.join(legacyHomology, directory.name, 'homology-project.json')
      if (!fs.existsSync(sourceProject)) continue
      try {
        const project = JSON.parse(fs.readFileSync(sourceProject, 'utf8')) as { id: string; templates?: Array<{ file: string }> }
        const counts = new Map<string, number>()
        for (const template of project.templates || []) {
          const owner = ownership.get(template.file)
          if (owner) counts.set(owner.workspaceId, (counts.get(owner.workspaceId) || 0) + 1)
        }
        const targetId = [...counts].sort((a, b) => b[1] - a[1])[0]?.[0] || createdWorkspaceIds[0]
        if (!targetId) continue
        const targetWorkspace = loadWorkspace(dataRoot, targetId)
        for (const template of project.templates || []) {
          const owner = ownership.get(template.file)
          if (owner?.workspaceId === targetId) { template.file = owner.relative; continue }
          const legacySource = path.resolve(dataRoot, template.file)
          if (!legacySource.startsWith(path.resolve(dataRoot)) || !fs.existsSync(legacySource)) continue
          const relative = path.join('files', 'legacy_homology', path.basename(template.file)).replace(/\\/g, '/')
          const destination = path.join(workspaceRoot(dataRoot, targetId), relative)
          fs.mkdirSync(path.dirname(destination), { recursive: true })
          fs.copyFileSync(legacySource, destination)
          template.file = relative
          if (!targetWorkspace.artifacts.some(artifact => artifact.file === relative)) targetWorkspace.artifacts.push(artifactFor(relative, 'legacy_homology'))
        }
        const destination = path.join(workspaceRoot(dataRoot, targetId), 'homology', safeId(project.id), 'homology-project.json')
        fs.mkdirSync(path.dirname(destination), { recursive: true })
        fs.writeFileSync(destination, `${JSON.stringify(project, null, 2)}\n`)
        saveWorkspace(dataRoot, targetWorkspace)
      } catch { /* keep malformed legacy workflows untouched */ }
    }
  }
  fs.writeFileSync(marker, `${new Date().toISOString()}\n`)
}

export function listWorkspaces(dataRoot: string): WorkspaceManifest[] {
  ensureWorkspaceMigration(dataRoot)
  return fs.readdirSync(workspacesRoot(dataRoot), { withFileTypes: true })
    .filter(entry => entry.isDirectory())
    .flatMap(entry => { try { return [loadWorkspace(dataRoot, entry.name)] } catch { return [] } })
    .sort((a, b) => b.updatedAt.localeCompare(a.updatedAt))
}

export function registerWorkspaceApi(server: ViteDevServer, dataRoot: string): void {
  ensureWorkspaceMigration(dataRoot)
  server.middlewares.use('/api/workspaces', async (req, res, next) => {
    try {
      const parts = (req.url || '').split('?')[0].split('/').filter(Boolean).map(decodeURIComponent)
      if (!parts.length && req.method === 'GET') {
        return sendJson(res, 200, listWorkspaces(dataRoot).map(({ id, name, updatedAt, artifacts }) => ({ id, name, updatedAt, artifactCount: artifacts.filter(artifact => !artifact.hidden).length })))
      }
      if (!parts.length && req.method === 'POST') {
        const body = JSON.parse((await readBody(req)).toString('utf8') || '{}')
        return sendJson(res, 201, createWorkspace(dataRoot, String(body.name || 'Untitled workspace')))
      }
      const id = parts[0]
      if (parts.length === 1 && req.method === 'GET') return sendJson(res, 200, loadWorkspace(dataRoot, id))
      if (parts.length === 1 && req.method === 'DELETE') {
        const source = workspaceRoot(dataRoot, id)
        if (!fs.existsSync(source)) return sendJson(res, 404, { error: `workspace not found: ${id}` })
        const trash = path.join(dataRoot, '_workspace_trash')
        fs.mkdirSync(trash, { recursive: true })
        const destination = path.join(trash, `${safeId(id)}-${new Date().toISOString().replace(/[:.]/g, '-')}`)
        fs.renameSync(source, destination)
        res.statusCode = 204
        return res.end()
      }
      if (parts.length === 1 && req.method === 'PUT') {
        const incoming = JSON.parse((await readBody(req)).toString('utf8')) as WorkspaceManifest
        const current = loadWorkspace(dataRoot, id)
        return sendJson(res, 200, saveWorkspace(dataRoot, { ...incoming, id: current.id, createdAt: current.createdAt, version: 1 }))
      }
      if (parts[1] === 'rename' && req.method === 'PATCH') {
        const workspace = loadWorkspace(dataRoot, id)
        const body = JSON.parse((await readBody(req)).toString('utf8') || '{}')
        const name = String(body.name || '').trim()
        if (!name) return sendJson(res, 400, { error: 'workspace name cannot be empty' })
        return sendJson(res, 200, saveWorkspace(dataRoot, { ...workspace, name }))
      }
      if (parts[1] === 'import' && req.method === 'POST') {
        const workspace = loadWorkspace(dataRoot, id)
        const encoded = Array.isArray(req.headers['x-file-name']) ? req.headers['x-file-name'][0] : req.headers['x-file-name']
        const original = decodeURIComponent(encoded || '')
        const safe = path.basename(original).replace(/[^a-zA-Z0-9._-]+/g, '_')
        if (!safe || !ALLOWED.has(path.extname(safe).toLowerCase())) return sendJson(res, 400, { error: `unsupported artifact type: ${original}` })
        const content = await readBody(req)
        if (!content.length) return sendJson(res, 400, { error: 'uploaded artifact is empty' })
        const directory = path.join(workspaceRoot(dataRoot, id), 'files', 'imports')
        fs.mkdirSync(directory, { recursive: true })
        let destination = path.join(directory, safe)
        if (fs.existsSync(destination)) destination = path.join(directory, `${path.basename(safe, path.extname(safe))}_${Date.now()}${path.extname(safe)}`)
        fs.writeFileSync(destination, content)
        const relative = path.relative(workspaceRoot(dataRoot, id), destination).replace(/\\/g, '/')
        const artifact = artifactFor(relative, 'imports')
        saveWorkspace(dataRoot, { ...workspace, artifacts: [...workspace.artifacts, artifact] })
        return sendJson(res, 201, artifact)
      }
      if (parts[1] === 'artifacts' && parts[2] && req.method === 'DELETE') {
        const workspace = loadWorkspace(dataRoot, id)
        const artifactIndex = workspace.artifacts.findIndex(artifact => artifact.id === parts[2])
        if (artifactIndex < 0) return sendJson(res, 404, { error: `artifact not found: ${parts[2]}` })
        const artifact = workspace.artifacts[artifactIndex]
        const root = workspaceRoot(dataRoot, id)
        const source = path.resolve(root, artifact.file)
        const relativeSource = path.relative(root, source)
        if (relativeSource.startsWith('..') || path.isAbsolute(relativeSource)) throw new Error(`invalid workspace path: ${artifact.file}`)
        if (fs.existsSync(source)) {
          const trash = path.join(root, '.trash')
          fs.mkdirSync(trash, { recursive: true })
          let destination = path.join(trash, `${Date.now()}-${path.basename(source)}`)
          while (fs.existsSync(destination)) destination = path.join(trash, `${Date.now()}-${crypto.randomUUID().slice(0, 5)}-${path.basename(source)}`)
          fs.renameSync(source, destination)
        }
        workspace.artifacts.splice(artifactIndex, 1)
        if (workspace.primaryFile === artifact.file) workspace.primaryFile = null
        if (workspace.secondaryFile === artifact.file) workspace.secondaryFile = null
        saveWorkspace(dataRoot, workspace)
        res.statusCode = 204
        return res.end()
      }
      if (parts[1] === 'artifacts' && parts[2] && parts[3] === 'rename' && req.method === 'PATCH') {
        const workspace = loadWorkspace(dataRoot, id)
        const artifact = workspace.artifacts.find(item => item.id === parts[2])
        if (!artifact) return sendJson(res, 404, { error: `artifact not found: ${parts[2]}` })
        const body = JSON.parse((await readBody(req)).toString('utf8') || '{}')
        const name = String(body.name || '').trim()
        if (!name) return sendJson(res, 400, { error: 'file name cannot be empty' })
        artifact.name = name
        return sendJson(res, 200, saveWorkspace(dataRoot, workspace))
      }
      if (parts[1] === 'files' && parts.length >= 3 && req.method === 'GET') {
        const relative = parts.slice(2).join('/')
        const file = resolveWorkspaceFile(dataRoot, id, relative)
        return fs.createReadStream(file).pipe(res)
      }
      return next()
    } catch (error: any) {
      return sendJson(res, /not found/.test(error?.message || '') ? 404 : 500, { error: error?.message || String(error) })
    }
  })
}
