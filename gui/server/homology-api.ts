import type { IncomingMessage, ServerResponse } from 'node:http'
import type { ViteDevServer } from 'vite'
import crypto from 'node:crypto'
import fs from 'node:fs'
import path from 'node:path'
import { runDvbfixerArgs } from './dvbfixer-runner'

export interface TemplateSelection {
  id: string
  file: string
  chain: string
  targetChain: string
  label?: string
}

export interface AlignmentSpan { start: number; end: number }
export interface AlignmentRow {
  id: string
  kind: 'target' | 'template'
  sequence: string
  templateId?: string
}
export interface AlignmentGroup {
  chainId: string
  rows: AlignmentRow[]
  masks: Record<string, AlignmentSpan[]>
}
export interface HomologyProject {
  version: 1
  id: string
  name: string
  targetFasta: string
  templates: TemplateSelection[]
  engine: 'mafft' | 'muscle' | 'clustalo'
  alignmentGroups: AlignmentGroup[]
  structuralAlignment?: string
  modelOptions: Record<string, string | number | boolean>
  createdAt: string
  updatedAt: string
}

const AA3: Record<string, string> = {
  ALA: 'A', ARG: 'R', ASN: 'N', ASP: 'D', CYS: 'C', GLN: 'Q', GLU: 'E', GLY: 'G',
  HIS: 'H', ILE: 'I', LEU: 'L', LYS: 'K', MET: 'M', PHE: 'F', PRO: 'P', SER: 'S',
  THR: 'T', TRP: 'W', TYR: 'Y', VAL: 'V', MSE: 'M', HID: 'H', HIE: 'H', HIP: 'H',
  ASH: 'D', GLH: 'E', LYN: 'K', CYX: 'C', CYM: 'C',
}

function sendJson(res: ServerResponse, status: number, body: unknown): void {
  res.statusCode = status
  res.setHeader('Content-Type', 'application/json')
  res.setHeader('Cache-Control', 'no-store')
  res.end(JSON.stringify(body))
}

function readBody(req: IncomingMessage): Promise<string> {
  return new Promise((resolve, reject) => {
    const chunks: Buffer[] = []
    req.on('data', chunk => chunks.push(chunk))
    req.on('end', () => resolve(Buffer.concat(chunks).toString('utf8')))
    req.on('error', reject)
  })
}

function safeId(id: string): string {
  if (!/^[a-zA-Z0-9_-]+$/.test(id)) throw new Error('invalid project id')
  return id
}

function resolveArtifact(root: string, relative: string): string {
  if (!relative || path.isAbsolute(relative)) throw new Error(`invalid artifact path: ${relative}`)
  const resolved = path.resolve(root, relative)
  const rel = path.relative(root, resolved)
  if (rel.startsWith('..') || path.isAbsolute(rel) || !fs.existsSync(resolved)) {
    throw new Error(`artifact not found: ${relative}`)
  }
  return resolved
}

function projectRoot(root: string): string {
  const directory = path.join(root, 'homology_projects')
  fs.mkdirSync(directory, { recursive: true })
  return directory
}

function projectPath(root: string, id: string): string {
  return path.join(projectRoot(root), safeId(id), 'homology-project.json')
}

function saveProject(root: string, project: HomologyProject): HomologyProject {
  const file = projectPath(root, project.id)
  fs.mkdirSync(path.dirname(file), { recursive: true })
  const next = { ...project, version: 1 as const, updatedAt: new Date().toISOString() }
  fs.writeFileSync(file, `${JSON.stringify(next, null, 2)}\n`)
  return next
}

function loadProject(root: string, id: string): HomologyProject {
  const file = projectPath(root, id)
  if (!fs.existsSync(file)) throw new Error(`project not found: ${id}`)
  return JSON.parse(fs.readFileSync(file, 'utf8')) as HomologyProject
}

export function parseFasta(text: string): Array<{ id: string; sequence: string }> {
  const records: Array<{ id: string; sequence: string }> = []
  let current: { id: string; sequence: string } | null = null
  for (const raw of text.split(/\r?\n/)) {
    const line = raw.trim()
    if (!line) continue
    if (line.startsWith('>')) {
      if (current) records.push(current)
      const id = line.slice(1).trim().split(/\s+/)[0]
      if (!id) throw new Error('target FASTA header is empty')
      current = { id, sequence: '' }
    } else {
      if (!current) throw new Error('target sequence appears before a FASTA header')
      current.sequence += line.replace(/\s+/g, '').toUpperCase()
    }
  }
  if (current) records.push(current)
  if (!records.length) throw new Error('target FASTA is empty')
  if (new Set(records.map(record => record.id)).size !== records.length) throw new Error('target chain IDs must be unique')
  for (const record of records) {
    if (!record.sequence || /[^A-Z*?-]/.test(record.sequence)) throw new Error(`invalid target sequence: ${record.id}`)
  }
  return records
}

function parseAlignedFasta(file: string): Array<{ id: string; sequence: string }> {
  return parseFasta(fs.readFileSync(file, 'utf8'))
}

function listFiles(directory: string): string[] {
  const output: string[] = []
  const walk = (current: string) => {
    for (const entry of fs.readdirSync(current, { withFileTypes: true })) {
      const full = path.join(current, entry.name)
      if (entry.isDirectory()) walk(full)
      else if (entry.isFile()) output.push(full)
    }
  }
  walk(directory)
  return output
}

async function alignProject(root: string, project: HomologyProject): Promise<HomologyProject> {
  const records = parseFasta(project.targetFasta)
  const directory = path.dirname(projectPath(root, project.id))
  const groups: AlignmentGroup[] = []
  for (const record of records) {
    const templates = project.templates.filter(template => template.targetChain === record.id)
    if (!templates.length) throw new Error(`target chain ${record.id} has no assigned template`)
    const input = path.join(directory, `_target_${record.id}.fasta`)
    const output = path.join(directory, `_alignment_${record.id}.fasta`)
    fs.writeFileSync(input, `>${record.id}\n${record.sequence}\n`)
    const args = [input, '-o', output, '--engine', project.engine]
    for (const template of templates) {
      const source = resolveArtifact(root, template.file)
      args.push('--template', `${source}:${template.chain}`)
    }
    const result = await runDvbfixerArgs('msa', args, directory)
    if (result.code !== 0) throw new Error(result.stderr || result.stdout || `MSA failed for ${record.id}`)
    const aligned = parseAlignedFasta(output)
    const rows: AlignmentRow[] = aligned.map(row => {
      if (row.id === record.id) return { ...row, kind: 'target' as const }
      const match = /^template_(\d+)_/.exec(row.id)
      const index = match ? Number(match[1]) - 1 : -1
      const template = templates[index]
      if (!template) throw new Error(`cannot map aligned template row ${row.id}`)
      return { ...row, id: template.id, kind: 'template' as const, templateId: template.id }
    })
    groups.push({ chainId: record.id, rows, masks: {} })
  }
  return saveProject(root, { ...project, alignmentGroups: groups })
}

interface PdbResidue { key: string; chain: string; number: string; aa: string; lines: string[] }

function pdbResidues(file: string, chain: string): PdbResidue[] {
  const residues: PdbResidue[] = []
  const byKey = new Map<string, PdbResidue>()
  for (const line of fs.readFileSync(file, 'utf8').split(/\r?\n/)) {
    if (!line.startsWith('ATOM  ') && !line.startsWith('HETATM')) continue
    if ((line[21] || ' ').trim() !== chain) continue
    const residueName = line.slice(17, 20).trim().toUpperCase()
    const aa = AA3[residueName]
    if (!aa) continue
    const number = line.slice(22, 27).trim()
    const key = `${chain}:${number}`
    let residue = byKey.get(key)
    if (!residue) {
      residue = { key, chain, number, aa, lines: [] }
      byKey.set(key, residue)
      residues.push(residue)
    }
    residue.lines.push(`${line}\n`)
  }
  if (!residues.length) throw new Error(`no protein residues found for chain ${chain} in ${file}`)
  return residues
}

function sanitize(value: string): string {
  return value.replace(/[^a-zA-Z0-9_-]+/g, '_')
}

export function materializeModelInputs(root: string, project: HomologyProject, runDir: string): {
  fasta: string; pir: string; fragments: string[]
} {
  if (!project.alignmentGroups.length) throw new Error('generate an alignment before modeling')
  for (const group of project.alignmentGroups) {
    if (new Set(group.rows.map(row => row.sequence.length)).size !== 1) {
      throw new Error(`alignment rows for target chain ${group.chainId} have different lengths`)
    }
  }
  const targetRecords = parseFasta(project.targetFasta)
  const groupLengths = project.alignmentGroups.map(group => group.rows[0]?.sequence.length || 0)
  if (groupLengths.some(length => length === 0)) throw new Error('alignment contains an empty chain group')
  const fasta = path.join(runDir, 'target.fasta')
  fs.writeFileSync(fasta, targetRecords.map(record => `>${record.id}\n${record.sequence.replace(/[-*]/g, '')}\n`).join(''))

  const fragments: string[] = []
  const pirEntries: string[] = []
  for (let groupIndex = 0; groupIndex < project.alignmentGroups.length; groupIndex++) {
    const group = project.alignmentGroups[groupIndex]
    for (const row of group.rows.filter(item => item.kind === 'template')) {
      const template = project.templates.find(item => item.id === row.templateId)
      if (!template) throw new Error(`template metadata missing for ${row.id}`)
      const source = resolveArtifact(root, template.file)
      const residues = pdbResidues(source, template.chain)
      const rowResidueByColumn: Array<PdbResidue | null> = []
      let residueIndex = 0
      for (const character of row.sequence) {
        rowResidueByColumn.push(character === '-' ? null : (residues[residueIndex++] || null))
      }
      if (residueIndex !== residues.length) {
        throw new Error(`alignment row ${row.id} does not match ${template.file}:${template.chain}`)
      }
      const spans = group.masks[row.id]?.length ? group.masks[row.id] : [{ start: 0, end: row.sequence.length }]
      for (let spanIndex = 0; spanIndex < spans.length; spanIndex++) {
        const span = spans[spanIndex]
        const selected = rowResidueByColumn.slice(span.start, span.end).filter((item): item is PdbResidue => !!item)
        if (!selected.length) continue
        const code = `${sanitize(template.id)}_${groupIndex + 1}_${spanIndex + 1}`
        const fragment = path.join(runDir, `${code}.pdb`)
        fs.writeFileSync(fragment, `${selected.flatMap(residue => residue.lines).join('')}TER\nEND\n`)
        fragments.push(fragment)
        const blocks = groupLengths.map(length => '-'.repeat(length))
        blocks[groupIndex] = row.sequence.split('').map((character, column) => (
          column >= span.start && column < span.end ? character : '-'
        )).join('')
        pirEntries.push(
          `>P1;${code}\nstructureX:${code}:${selected[0].number}:${template.chain}:` +
          `${selected[selected.length - 1].number}:${template.chain}::::\n${blocks.join('/')}*\n`,
        )
      }
    }
  }
  if (!fragments.length) throw new Error('template masks select no residues')
  const targetBlocks = project.alignmentGroups.map(group => {
    const row = group.rows.find(item => item.kind === 'target')
    if (!row) throw new Error(`target row missing for chain ${group.chainId}`)
    return row.sequence
  })
  const pir = path.join(runDir, 'alignment.pir')
  fs.writeFileSync(pir, `${pirEntries.join('')}>P1;target\nsequence:target::::::::\n${targetBlocks.join('/')}*\n`)
  return { fasta, pir, fragments }
}

function registerRun(root: string, project: HomologyProject, files: string[]): string {
  const relativeFiles = files.map(file => path.relative(root, file).replace(/\\/g, '/'))
  const primary = relativeFiles.find(file => file.endsWith('.pdb')) || relativeFiles[0]
  if (!primary) return ''
  const indexPath = path.join(root, 'index.json')
  let entries: any[] = []
  if (fs.existsSync(indexPath)) try { entries = JSON.parse(fs.readFileSync(indexPath, 'utf8')) } catch {}
  entries.push({
    id: primary, file: primary, name: `${project.name} → homology`,
    kind: primary.endsWith('.pdb') ? 'structure' : 'artifact', artifactType: 'homology-model',
    artifacts: relativeFiles, command: 'homology', organism: '', chains: 0, residues: 0,
    description: `DVBfixer GUI homology project ${project.id}`,
  })
  fs.writeFileSync(indexPath, JSON.stringify(entries, null, 2))
  return primary
}

async function modelProject(root: string, project: HomologyProject): Promise<unknown> {
  const runDir = path.join(root, `dvb_homology_${new Date().toISOString().replace(/[:.]/g, '-').slice(0, 19)}`)
  fs.mkdirSync(runDir, { recursive: true })
  const materialized = materializeModelInputs(root, project, runDir)
  const prefix = path.join(runDir, sanitize(project.name || 'target'))
  const args = [materialized.fasta]
  for (const fragment of materialized.fragments) args.push('--template', fragment)
  args.push('--alignment', materialized.pir, '-o', prefix)
  for (const [flag, value] of Object.entries(project.modelOptions || {})) {
    if (value === '' || value === false || value === null || value === undefined) continue
    if (value === true) args.push(flag)
    else args.push(flag, String(value))
  }
  const result = await runDvbfixerArgs('homology', args, runDir)
  fs.writeFileSync(path.join(runDir, 'run.json'), `${JSON.stringify({ projectId: project.id, args, exitCode: result.code }, null, 2)}\n`)
  fs.writeFileSync(path.join(runDir, 'stdout.log'), result.stdout)
  fs.writeFileSync(path.join(runDir, 'stderr.log'), result.stderr)
  const files = listFiles(runDir)
  const outputFile = result.code === 0 ? registerRun(root, project, files) : ''
  return { ok: result.code === 0, exitCode: result.code, stdout: result.stdout, stderr: result.stderr,
    outputFile, outputDir: path.relative(root, runDir).replace(/\\/g, '/'),
    artifacts: files.map(file => path.relative(root, file).replace(/\\/g, '/')) }
}

export function registerHomologyApi(server: ViteDevServer, root: string): void {
  server.middlewares.use('/api/homology', async (req, res, next) => {
    try {
      const parts = (req.url || '').split('?')[0].split('/').filter(Boolean)
      if (req.method === 'GET' && parts[0] === 'projects' && !parts[1]) {
        const projects = fs.readdirSync(projectRoot(root), { withFileTypes: true })
          .filter(entry => entry.isDirectory())
          .flatMap(entry => {
            try { return [loadProject(root, entry.name)] } catch { return [] }
          })
          .sort((a, b) => b.updatedAt.localeCompare(a.updatedAt))
          .map(project => ({ id: project.id, name: project.name, updatedAt: project.updatedAt }))
        return sendJson(res, 200, projects)
      }
      if (req.method === 'GET' && parts[0] === 'engines') {
        const result = await runDvbfixerArgs('msa', ['--list-engines'])
        const engines = Object.fromEntries(result.stdout.trim().split(/\r?\n/).filter(Boolean).map(line => {
          const [name, ...rest] = line.split(':')
          const executable = rest.join(':').trim()
          return [name, { available: executable !== 'MISSING', path: executable === 'MISSING' ? null : executable }]
        }))
        return sendJson(res, result.code === 0 ? 200 : 500, engines)
      }
      if (req.method === 'POST' && parts[0] === 'projects' && !parts[1]) {
        const now = new Date().toISOString()
        const project: HomologyProject = {
          version: 1, id: crypto.randomUUID(), name: 'Untitled homology project', targetFasta: '',
          templates: [], engine: 'mafft', alignmentGroups: [], modelOptions: {
            '--num-models': 5, '--md-level': 'fast',
          }, createdAt: now, updatedAt: now,
        }
        return sendJson(res, 201, saveProject(root, project))
      }
      if (parts[0] === 'projects' && parts[1] && req.method === 'GET') {
        return sendJson(res, 200, loadProject(root, parts[1]))
      }
      if (parts[0] === 'projects' && parts[1] && req.method === 'PUT') {
        const incoming = JSON.parse(await readBody(req)) as HomologyProject
        const current = loadProject(root, parts[1])
        return sendJson(res, 200, saveProject(root, { ...incoming, id: current.id, createdAt: current.createdAt, version: 1 }))
      }
      if (req.method === 'POST' && parts[0] === 'align') {
        const project = JSON.parse(await readBody(req)) as HomologyProject
        return sendJson(res, 200, await alignProject(root, project))
      }
      if (req.method === 'POST' && parts[0] === 'salign') {
        const project = JSON.parse(await readBody(req)) as HomologyProject
        if (project.templates.length < 2) return sendJson(res, 400, { error: 'select at least two templates' })
        const directory = path.dirname(projectPath(root, project.id))
        const output = path.join(directory, 'structural_alignment.pir')
        const fitDir = path.join(directory, 'fitted')
        const args = project.templates.map(template => `${resolveArtifact(root, template.file)}:${template.chain}`)
        args.push('-o', output, '--fit-dir', fitDir)
        const result = await runDvbfixerArgs('salign', args, directory)
        if (result.code !== 0) return sendJson(res, 500, { error: result.stderr || result.stdout })
        const relative = path.relative(root, output).replace(/\\/g, '/')
        return sendJson(res, 200, saveProject(root, { ...project, structuralAlignment: relative }))
      }
      if (req.method === 'POST' && parts[0] === 'model') {
        const project = JSON.parse(await readBody(req)) as HomologyProject
        const result = await modelProject(root, saveProject(root, project)) as { ok: boolean }
        return sendJson(res, result.ok ? 200 : 500, result)
      }
      return next()
    } catch (error: any) {
      return sendJson(res, 500, { error: error?.message || String(error) })
    }
  })
}
