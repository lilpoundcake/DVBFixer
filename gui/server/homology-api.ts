import type { IncomingMessage, ServerResponse } from 'node:http'
import type { ViteDevServer } from 'vite'
import crypto from 'node:crypto'
import fs from 'node:fs'
import path from 'node:path'
import { runDvbfixerArgs } from './dvbfixer-runner'
import { workspaceRoot } from './workspace-api'

export interface TemplateSelection {
  id: string
  file: string
  chain: string
  targetChain: string
  label?: string
  fittedFile?: string
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
  maskModes?: Record<string, 'all' | 'ranges' | 'none'>
}
export interface HomologyProject {
  version: 1
  id: string
  name: string
  targetFasta: string
  templates: TemplateSelection[]
  engine: 'mafft' | 'muscle' | 'clustalo'
  alignmentGroups: AlignmentGroup[]
  structuralAlignment?: string | Record<string, string>
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

export function parseSequenceArtifact(file: string): Array<{ id: string; sequence: string; length: number }> {
  const extension = path.extname(file).toLowerCase()
  const text = fs.readFileSync(file, 'utf8')
  if (['.fasta', '.fa', '.faa', '.aln', '.txt'].includes(extension) || text.trimStart().startsWith('>')) {
    if (/^>P1;/m.test(text)) {
      const records: Array<{ id: string; sequence: string; length: number }> = []
      const chunks = text.split(/^>P1;/m).slice(1)
      for (const chunk of chunks) {
        const lines = chunk.split(/\r?\n/)
        const id = (lines.shift() || '').trim()
        lines.shift()
        const sequence = lines.join('').replace(/[\s*/]/g, '').toUpperCase()
        if (id && sequence) records.push({ id, sequence, length: sequence.length })
      }
      if (records.length) return records
    }
    return parseFasta(text).map(record => ({ ...record, length: record.sequence.replace(/[-*]/g, '').length }))
  }
  if (extension === '.pdb') {
    const seqres = new Map<string, string[]>()
    const seqresOrder: string[] = []
    for (const line of text.split(/\r?\n/)) {
      if (!line.startsWith('SEQRES')) continue
      const chain = (line[11] || ' ').trim() || '_'
      if (!seqres.has(chain)) { seqres.set(chain, []); seqresOrder.push(chain) }
      for (const residue of line.slice(19).trim().split(/\s+/)) {
        const aa = AA3[residue.toUpperCase()]
        if (aa) seqres.get(chain)!.push(aa)
      }
    }
    const seqresRecords = seqresOrder.map(id => ({ id, sequence: seqres.get(id)!.join(''), length: seqres.get(id)!.length })).filter(record => record.length)
    if (seqresRecords.length) return seqresRecords
    const chainOrder: string[] = []
    const residues = new Map<string, Map<string, string>>()
    for (const line of text.split(/\r?\n/)) {
      if (!line.startsWith('ATOM  ') && !line.startsWith('HETATM')) continue
      const aa = AA3[line.slice(17, 20).trim().toUpperCase()]
      if (!aa) continue
      const chain = (line[21] || ' ').trim() || '_'
      if (!residues.has(chain)) { residues.set(chain, new Map()); chainOrder.push(chain) }
      residues.get(chain)!.set(line.slice(22, 27).trim(), aa)
    }
    const records = chainOrder.map(id => ({ id, sequence: [...residues.get(id)!.values()].join(''), length: residues.get(id)!.size }))
    if (records.length) return records
  }
  if (extension === '.cif' || extension === '.mmcif') {
    const tokens = text.match(/'(?:[^']|'')*'|"(?:[^"]|"")*"|[^\s]+/g)?.map(token => token.replace(/^(['"])(.*)\1$/, '$2')) || []
    for (let index = 0; index < tokens.length; index++) {
      if (tokens[index] !== 'loop_') continue
      const headers: string[] = []
      while (tokens[index + 1]?.startsWith('_')) headers.push(tokens[++index])
      if (!headers.includes('_atom_site.group_PDB')) continue
      const groupAt = headers.indexOf('_atom_site.group_PDB')
      const chainAt = headers.indexOf('_atom_site.auth_asym_id') >= 0 ? headers.indexOf('_atom_site.auth_asym_id') : headers.indexOf('_atom_site.label_asym_id')
      const seqAt = headers.indexOf('_atom_site.auth_seq_id') >= 0 ? headers.indexOf('_atom_site.auth_seq_id') : headers.indexOf('_atom_site.label_seq_id')
      const compAt = headers.indexOf('_atom_site.auth_comp_id') >= 0 ? headers.indexOf('_atom_site.auth_comp_id') : headers.indexOf('_atom_site.label_comp_id')
      if (chainAt < 0 || seqAt < 0 || compAt < 0) break
      const rows = new Map<string, Map<string, string>>()
      const order: string[] = []
      while (index + headers.length < tokens.length) {
        const row = tokens.slice(index + 1, index + 1 + headers.length)
        if (row[0] === 'loop_' || row[0] === 'stop_' || row[0]?.startsWith('_') || row[0]?.startsWith('data_')) break
        index += headers.length
        if (row[groupAt] !== 'ATOM' && row[groupAt] !== 'HETATM') continue
        const aa = AA3[row[compAt]?.toUpperCase()]
        if (!aa) continue
        const chain = row[chainAt] === '.' || row[chainAt] === '?' ? '_' : row[chainAt]
        if (!rows.has(chain)) { rows.set(chain, new Map()); order.push(chain) }
        rows.get(chain)!.set(row[seqAt], aa)
      }
      const records = order.map(id => ({ id, sequence: [...rows.get(id)!.values()].join(''), length: rows.get(id)!.size }))
      if (records.length) return records
    }
  }
  throw new Error(`sequence parsing is not supported for ${path.basename(file)}`)
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
      // MSA correspondence must always come from the user's source template.
      // Per-chain fitted artifacts can live in different reference frames and
      // are for inspection/export only, never for rebuilding model inputs.
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
    groups.push({ chainId: record.id, rows, masks: {}, maskModes: Object.fromEntries(rows.filter(row => row.kind === 'template').map(row => [row.id, 'all'])) })
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

  const compositeBlocks: string[] = []
  const compositeLines: string[] = []
  let atomSerial = 1
  let firstPosition: { number: number; chain: string } | null = null
  let lastPosition: { number: number; chain: string } | null = null
  for (let groupIndex = 0; groupIndex < project.alignmentGroups.length; groupIndex++) {
    const group = project.alignmentGroups[groupIndex]
    const targetRow = group.rows.find(item => item.kind === 'target')
    if (!targetRow) throw new Error(`target row missing for chain ${group.chainId}`)
    const chosen: Array<PdbResidue | null> = Array(groupLengths[groupIndex]).fill(null)
    for (const row of group.rows.filter(item => item.kind === 'template')) {
      const template = project.templates.find(item => item.id === row.templateId)
      if (!template) throw new Error(`template metadata missing for ${row.id}`)
      const source = resolveArtifact(root, template.fittedFile || template.file)
      const residues = pdbResidues(source, template.chain)
      const alignedSequence = row.sequence.replace(/-/g, '')
      const coordinateSequence = residues.map(residue => residue.aa).join('')
      if (alignedSequence !== coordinateSequence) {
        let mismatch = 0
        while (mismatch < Math.min(alignedSequence.length, coordinateSequence.length) &&
               alignedSequence[mismatch] === coordinateSequence[mismatch]) mismatch++
        throw new Error(
          `alignment row ${row.id} is not mapped to ${template.file}:${template.chain} ` +
          `(first difference at residue ${mismatch + 1}); regenerate the alignment`,
        )
      }
      const rowResidueByColumn: Array<PdbResidue | null> = []
      let residueIndex = 0
      for (const character of row.sequence) {
        rowResidueByColumn.push(character === '-' ? null : (residues[residueIndex++] || null))
      }
      if (residueIndex !== residues.length) {
        throw new Error(`alignment row ${row.id} does not match ${template.file}:${template.chain}`)
      }
      const mode = group.maskModes?.[row.id] || (group.masks[row.id]?.length ? 'ranges' : 'all')
      const spans = mode === 'none' ? [] : mode === 'ranges' ? (group.masks[row.id] || []) : [{ start: 0, end: row.sequence.length }]
      for (const span of spans) {
        for (let column = Math.max(0, span.start); column < Math.min(span.end, chosen.length); column++) {
          // Template insertions have no target residue to receive them. When
          // masks overlap, earlier template rows have explicit precedence.
          if (targetRow.sequence[column] !== '-' && !chosen[column] && rowResidueByColumn[column]) {
            chosen[column] = rowResidueByColumn[column]
          }
        }
      }
    }
    let targetOrdinal = 0
    const block: string[] = []
    for (let column = 0; column < targetRow.sequence.length; column++) {
      if (targetRow.sequence[column] !== '-') targetOrdinal++
      const residue = chosen[column]
      block.push(residue ? residue.aa : '-')
      if (!residue) continue
      firstPosition ||= { number: targetOrdinal, chain: group.chainId }
      lastPosition = { number: targetOrdinal, chain: group.chainId }
      for (const raw of residue.lines) {
        const line = raw.replace(/\r?\n$/, '').padEnd(80, ' ')
        compositeLines.push(
          `${line.slice(0, 6)}${String(atomSerial++).padStart(5)}${line.slice(11, 21)}` +
          `${group.chainId.slice(0, 1) || ' '}${String(targetOrdinal).padStart(4)} ${line.slice(27)}\n`,
        )
      }
    }
    compositeBlocks.push(block.join(''))
    compositeLines.push('TER\n')
  }
  if (!firstPosition || !lastPosition) throw new Error('template masks select no residues')
  const compositeCode = 'selected_template_mosaic'
  const composite = path.join(runDir, `${compositeCode}.pdb`)
  fs.writeFileSync(composite, `${compositeLines.join('')}END\n`)
  const fragments = [composite]
  const targetBlocks = project.alignmentGroups.map(group => {
    const row = group.rows.find(item => item.kind === 'target')
    if (!row) throw new Error(`target row missing for chain ${group.chainId}`)
    return row.sequence
  })
  const pir = path.join(runDir, 'alignment.pir')
  fs.writeFileSync(pir,
    `>P1;${compositeCode}\nstructureX:${compositeCode}:${firstPosition.number}:${firstPosition.chain}:` +
    `${lastPosition.number}:${lastPosition.chain}::::\n${compositeBlocks.join('/')}*\n` +
    `>P1;target\nsequence:target::::::::\n${targetBlocks.join('/')}*\n`)
  return { fasta, pir, fragments }
}

function registerRun(root: string, project: HomologyProject, files: string[]): string {
  const relativeFiles = files.map(file => path.relative(root, file).replace(/\\/g, '/'))
  const primary = relativeFiles.find(file => file.endsWith('.pdb')) || relativeFiles[0]
  if (!primary) return ''
  const manifest = JSON.parse(fs.readFileSync(path.join(root, 'workspace.json'), 'utf8'))
  for (const relative of relativeFiles) manifest.artifacts.push({
    id: crypto.randomUUID(), file: relative,
    name: relative === primary ? `${project.name} → homology` : path.basename(relative),
    kind: relative.endsWith('.pdb') ? 'structure' : 'artifact', artifactType: 'homology-model',
    command: 'homology', folder: path.dirname(relative),
    hidden: ['run.json', 'stdout.log', 'stderr.log', 'target.fasta', 'alignment.pir'].includes(path.basename(relative)) ||
      path.basename(relative).startsWith('_') || (relative.endsWith('.pdb') && relative !== primary && !path.basename(relative).startsWith(sanitize(project.name || 'target'))),
  })
  manifest.toolState = { ...manifest.toolState, lastHomologyRun: { projectId: project.id, files: relativeFiles } }
  fs.writeFileSync(path.join(root, 'workspace.json'), `${JSON.stringify(manifest, null, 2)}\n`)
  return primary
}

async function modelProject(root: string, project: HomologyProject): Promise<unknown> {
  const runDir = path.join(root, 'runs', `dvb_homology_${new Date().toISOString().replace(/[:.]/g, '-').slice(0, 19)}`)
  fs.mkdirSync(runDir, { recursive: true })
  const fasta = path.join(runDir, 'target.fasta')
  fs.writeFileSync(fasta, parseFasta(project.targetFasta).map(record => `>${record.id}\n${record.sequence.replace(/[-*]/g, '')}\n`).join(''))
  const plan = path.join(runDir, 'template-plan.json')
  fs.writeFileSync(plan, `${JSON.stringify({
    templates: project.templates.map(template => ({
      id: template.id, path: resolveArtifact(root, template.file), chain: template.chain,
      targetChain: template.targetChain,
    })),
    alignmentGroups: project.alignmentGroups,
  }, null, 2)}\n`)
  const prefix = path.join(runDir, sanitize(project.name || 'target'))
  const args = [fasta, '--template-plan', plan, '-o', prefix]
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

async function structurallyAlignGroups(root: string, project: HomologyProject): Promise<HomologyProject> {
  const directory = path.dirname(projectPath(root, project.id))
  const templates: TemplateSelection[] = project.templates.map(template => ({ ...template, fittedFile: undefined }))
  const outputs: Record<string, string> = {}
  const targetChains = [...new Set(templates.map(template => template.targetChain))]
  const globalReferenceFile = templates[0]?.file
  for (const targetChain of targetChains) {
    const members = templates.filter(template => template.targetChain === targetChain)
    const sharedReference = members.find(template => template.file === globalReferenceFile)
    if (targetChains.length > 1 && !sharedReference) {
      throw new Error(
        `target chain ${targetChain} has no template from the common reference structure ${globalReferenceFile}; ` +
        'add that structure/chain to every target-chain group so fitted coordinates share one frame',
      )
    }
    const group = sharedReference
      ? [sharedReference, ...members.filter(template => template.id !== sharedReference.id)]
      : members
    if (group.length < 2) continue
    const safeChain = sanitize(targetChain)
    const output = path.join(directory, `structural_alignment_${safeChain}.pir`)
    const fitDir = path.join(directory, 'fitted', safeChain)
    const args = group.map(template => `${resolveArtifact(root, template.file)}:${template.chain}`)
    args.push('-o', output, '--fit-dir', fitDir, '--engine', 'biopython')
    const result = await runDvbfixerArgs('salign', args, directory)
    if (result.code !== 0) throw new Error(result.stderr || result.stdout || `structural alignment failed for target chain ${targetChain}`)
    outputs[targetChain] = path.relative(root, output).replace(/\\/g, '/')
    group.forEach((template, index) => {
      const candidate = path.join(fitDir, `template_${index + 1}_${template.chain || 'all'}_fit.pdb`)
      if (!fs.existsSync(candidate)) throw new Error(`structural fitting did not produce ${path.basename(candidate)}`)
      template.fittedFile = path.relative(root, candidate).replace(/\\/g, '/')
    })
  }
  return saveProject(root, { ...project, templates, structuralAlignment: outputs })
}

export function registerHomologyApi(server: ViteDevServer, root: string): void {
  server.middlewares.use('/api/homology', async (req, res, next) => {
    try {
      const requestUrl = new URL(req.url || '/', 'http://localhost')
      const parts = requestUrl.pathname.split('/').filter(Boolean)
      const queryWorkspace = requestUrl.searchParams.get('workspaceId') || ''
      if (req.method === 'GET' && parts[0] === 'projects' && !parts[1]) {
        if (!queryWorkspace) return sendJson(res, 400, { error: 'workspaceId is required' })
        const storageRoot = workspaceRoot(root, queryWorkspace)
        const projects = fs.readdirSync(projectRoot(storageRoot), { withFileTypes: true })
          .filter(entry => entry.isDirectory())
          .flatMap(entry => {
            try { return [loadProject(storageRoot, entry.name)] } catch { return [] }
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
      if (req.method === 'POST' && (parts[0] === 'parse-sequence' || parts[0] === 'chains')) {
        const body = JSON.parse(await readBody(req) || '{}') as { workspaceId?: string; file?: string }
        if (!body.workspaceId || !body.file) return sendJson(res, 400, { error: 'workspaceId and file are required' })
        const file = resolveArtifact(workspaceRoot(root, body.workspaceId), body.file)
        const records = parseSequenceArtifact(file)
        return sendJson(res, 200, parts[0] === 'chains' ? records.map(record => ({ id: record.id, length: record.length })) : { records })
      }
      if (req.method === 'POST' && parts[0] === 'projects' && !parts[1]) {
        const request = JSON.parse(await readBody(req) || '{}') as { workspaceId?: string }
        if (!request.workspaceId) return sendJson(res, 400, { error: 'workspaceId is required' })
        const storageRoot = workspaceRoot(root, request.workspaceId)
        const now = new Date().toISOString()
        const project: HomologyProject = {
          version: 1, id: crypto.randomUUID(), name: 'Untitled homology project', targetFasta: '',
          templates: [], engine: 'mafft', alignmentGroups: [], modelOptions: {
            '--num-models': 5, '--md-level': 'fast',
          }, createdAt: now, updatedAt: now,
        }
        return sendJson(res, 201, saveProject(storageRoot, project))
      }
      if (parts[0] === 'projects' && parts[1] && req.method === 'GET') {
        if (!queryWorkspace) return sendJson(res, 400, { error: 'workspaceId is required' })
        return sendJson(res, 200, loadProject(workspaceRoot(root, queryWorkspace), parts[1]))
      }
      if (parts[0] === 'projects' && parts[1] && req.method === 'PUT') {
        const incoming = JSON.parse(await readBody(req)) as HomologyProject
        const workspaceId = (incoming as any).workspaceId as string
        if (!workspaceId) return sendJson(res, 400, { error: 'workspaceId is required' })
        const storageRoot = workspaceRoot(root, workspaceId)
        const current = loadProject(storageRoot, parts[1])
        const { workspaceId: _workspaceId, ...clean } = incoming as any
        void _workspaceId
        return sendJson(res, 200, saveProject(storageRoot, { ...clean, id: current.id, createdAt: current.createdAt, version: 1 }))
      }
      if (req.method === 'POST' && parts[0] === 'align') {
        const incoming = JSON.parse(await readBody(req)) as HomologyProject & { workspaceId?: string }
        if (!incoming.workspaceId) return sendJson(res, 400, { error: 'workspaceId is required' })
        const { workspaceId, ...project } = incoming
        return sendJson(res, 200, await alignProject(workspaceRoot(root, workspaceId), project))
      }
      if (req.method === 'POST' && parts[0] === 'salign') {
        const incoming = JSON.parse(await readBody(req)) as HomologyProject & { workspaceId?: string }
        if (!incoming.workspaceId) return sendJson(res, 400, { error: 'workspaceId is required' })
        const { workspaceId, ...project } = incoming
        const storageRoot = workspaceRoot(root, workspaceId)
        if (project.templates.length < 2) return sendJson(res, 400, { error: 'select at least two templates' })
        return sendJson(res, 200, await structurallyAlignGroups(storageRoot, project))
      }
      if (req.method === 'POST' && parts[0] === 'model') {
        const incoming = JSON.parse(await readBody(req)) as HomologyProject & { workspaceId?: string }
        if (!incoming.workspaceId) return sendJson(res, 400, { error: 'workspaceId is required' })
        const { workspaceId, ...project } = incoming
        const storageRoot = workspaceRoot(root, workspaceId)
        const result = await modelProject(storageRoot, project) as { ok: boolean }
        return sendJson(res, result.ok ? 200 : 500, result)
      }
      return next()
    } catch (error: any) {
      return sendJson(res, 500, { error: error?.message || String(error) })
    }
  })
}
