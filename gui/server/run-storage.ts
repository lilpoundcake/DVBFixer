/** Common transport/storage policy; never infer scientific acceptance here. */
import crypto from 'node:crypto'
import fs from 'node:fs'
import path from 'node:path'
import { Readable } from 'node:stream'
import { pipeline } from 'node:stream/promises'
import { createGzip, createGunzip } from 'node:zlib'
import { tarHeader, writeJsonAtomic } from './workspace-api'

export const RUN_ARCHIVE = 'diagnostics.tar.gz'
export const RUN_SUMMARY = 'run-summary.json'
const CONTROL = new Set(['job.json', 'run.json', 'stdout.log', 'stderr.log', 'cancel.requested', RUN_ARCHIVE, RUN_SUMMARY])
const CACHE = new Set(['.cache', '.home', '.tmp'])

export type RunStatus = 'succeeded' | 'failed' | 'cancelled' | 'validation_failed'
export interface RunStorageSummary {
  version: 1
  command: string
  status: RunStatus
  outputs: string[]
  diagnostics: { archive: string; files: Array<{ file: string; bytes: number; sha256: string }> } | null
  warning?: string
}

/** No traversal through links or special files, including hidden directories. */
export function runFiles(root: string): string[] {
  if (fs.lstatSync(root).isSymbolicLink()) throw new Error('run root must not be a symlink')
  const files: string[] = []
  const walk = (directory: string) => {
    for (const entry of fs.readdirSync(directory, { withFileTypes: true })) {
      const full = path.join(directory, entry.name)
      if (entry.isDirectory()) walk(full)
      else if (entry.isFile() && fs.lstatSync(full).nlink === 1) files.push(path.relative(root, full).split(path.sep).join('/'))
      else throw new Error('run contains a link or special file')
    }
  }
  walk(root)
  return files.sort()
}

/** Copy native runtime evidence, not runtime caches; retain failures as well. */
export function copyNativeRun(source: string, destination: string): void {
  if (fs.existsSync(destination)) runFiles(destination)
  fs.cpSync(source, destination, {
    recursive: true, force: false, errorOnExist: true,
    filter: file => {
      const relative = path.relative(source, file)
      if (relative.split(path.sep).some(part => CACHE.has(part))) return false
      const stat = fs.lstatSync(file)
      if (stat.isSymbolicLink() || (!stat.isDirectory() && (!stat.isFile() || stat.nlink !== 1))) {
        throw new Error('native run contains a link or special file')
      }
      return true
    },
  })
}

export function outputFiles(files: string[], declared: string | null, bundleDirectories: string[] = []): string[] {
  const isInternal = (file: string) => file.split('/').some(part => part.startsWith('.') || part.startsWith('_'))
  return files.filter(file => {
    if (CONTROL.has(file)) return false
    // A published bundle is indivisible: retain all its files and relative paths.
    if (bundleDirectories.some(dir => file.startsWith(`${dir}/`))) return true
    if (!declared) return false
    if (file === declared || file.startsWith(`${declared}/`)) return true
    if (isInternal(file)) return false
    // TOP packages include per-molecule/restraint files whose names need not
    // share the main topology's stem. Never archive those live dependencies.
    if (/\.top$/i.test(declared) && (/\.itp$/i.test(file) || file.split('/').some(part => part.endsWith('.ff')))) return true
    // Prefix outputs and numbered candidates, including their companion files.
    const stem = declared.replace(/\.[^/.]+$/, '')
    if (file.startsWith(`${stem}.`)) return true
    if (!file.startsWith(`${stem}_`)) return false
    // File-mode PDB tools may publish numbered candidates. Arbitrary suffixes
    // such as _renum/_prepare/_minimize are pipeline intermediates, not outputs.
    return !/\.pdb$/i.test(declared) || /^\d+(?:\.|$)/.test(file.slice(stem.length + 1))
  })
}

async function digestFile(file: string): Promise<string> {
  const hash = crypto.createHash('sha256')
  for await (const chunk of fs.createReadStream(file)) hash.update(chunk)
  return hash.digest('hex')
}

async function archiveDiagnostics(root: string, files: string[]): Promise<NonNullable<RunStorageSummary['diagnostics']>> {
  const destination = path.join(root, RUN_ARCHIVE)
  const temporary = path.join(root, `.diagnostics-${crypto.randomUUID()}.tmp`)
  if (fs.existsSync(destination)) throw new Error('diagnostics archive already exists')
  if (files.includes('diagnostics-manifest.json')) throw new Error('reserved diagnostics manifest name already exists')
  const entries: NonNullable<RunStorageSummary['diagnostics']>['files'] = []
  const originals = new Map<string, fs.Stats>()
  const rawHash = crypto.createHash('sha256')
  async function* tar(): AsyncGenerator<Buffer> {
    for (const file of files) {
      if (file.includes('\\') || file.split('/').some(part => part === '..' || part === '.' || /^[A-Za-z]:/.test(part))) {
        throw new Error('unsafe diagnostic archive path')
      }
      const full = path.join(root, file)
      const stat = fs.lstatSync(full)
      if (!stat.isFile() || stat.nlink !== 1) throw new Error('diagnostic file is not private and regular')
      originals.set(file, stat)
      yield tarHeader(file, stat.size, Math.max(0, stat.mtimeMs))
      const hash = crypto.createHash('sha256')
      let bytes = 0
      for await (const chunk of fs.createReadStream(full)) {
        const buffer = Buffer.from(chunk)
        bytes += buffer.length
        hash.update(buffer)
        yield buffer
      }
      if (bytes !== stat.size) throw new Error('diagnostic file changed during archival')
      entries.push({ file, bytes, sha256: hash.digest('hex') })
      if (stat.size % 512) yield Buffer.alloc(512 - stat.size % 512)
    }
    const manifest = Buffer.from(JSON.stringify({ version: 1, files: entries }, null, 2))
    yield tarHeader('diagnostics-manifest.json', manifest.length, 0)
    yield manifest
    if (manifest.length % 512) yield Buffer.alloc(512 - manifest.length % 512)
    yield Buffer.alloc(1024)
  }
  async function* measured(): AsyncGenerator<Buffer> {
    for await (const chunk of tar()) { rawHash.update(chunk); yield chunk }
  }
  try {
    await pipeline(Readable.from(measured()), createGzip(), fs.createWriteStream(temporary, { flags: 'wx', mode: 0o600 }))
    // Verify the entire decompressed stream before publishing or removing sources.
    const verified = crypto.createHash('sha256')
    await pipeline(fs.createReadStream(temporary), createGunzip(), async stream => {
      for await (const chunk of stream) verified.update(chunk)
    })
    if (verified.digest('hex') !== rawHash.digest('hex')) throw new Error('diagnostics archive verification failed')
    for (const entry of entries) {
      const full = path.join(root, entry.file)
      const before = originals.get(entry.file)!
      const after = fs.lstatSync(full)
      if (!after.isFile() || after.nlink !== 1 || before.ino !== after.ino || before.size !== after.size ||
          before.mtimeMs !== after.mtimeMs || await digestFile(full) !== entry.sha256) {
        throw new Error('diagnostic file changed before archival commit')
      }
    }
    const descriptor = fs.openSync(temporary, 'r')
    try { fs.fsyncSync(descriptor) } finally { fs.closeSync(descriptor) }
    fs.renameSync(temporary, destination)
    for (const entry of entries) {
      try { fs.unlinkSync(path.join(root, entry.file)) } catch { /* The verified archive exists; keep any undeletable original. */ }
    }
    const prune = (directory: string) => {
      for (const entry of fs.readdirSync(directory, { withFileTypes: true })) {
        if (!entry.isDirectory()) continue
        const full = path.join(directory, entry.name)
        prune(full)
        if (fs.readdirSync(full).length === 0) {
          try { fs.rmdirSync(full) } catch { /* Only empty, run-owned directories are candidates. */ }
        }
      }
    }
    prune(root)
    return { archive: RUN_ARCHIVE, files: entries }
  } finally {
    fs.rmSync(temporary, { force: true })
  }
}

/** Archive only unpublished extras. Unknown declared outputs remain unpacked. */
export async function compactRun(root: string, command: string, status: RunStatus, outputs: string[]): Promise<RunStorageSummary> {
  const summary: RunStorageSummary = { version: 1, command, status, outputs, diagnostics: null }
  try {
    const files = runFiles(root).filter(file => !CONTROL.has(file) && !outputs.includes(file))
    if (files.length) summary.diagnostics = await archiveDiagnostics(root, files)
    if (files.some(file => fs.existsSync(path.join(root, file)))) summary.warning = 'Some archived diagnostics remain unpacked.'
  } catch (error) {
    // Packaging is not scientific failure. Never discard originals on failure.
    summary.warning = `Diagnostics kept unpacked: ${error instanceof Error ? error.message : String(error)}`
  }
  writeJsonAtomic(path.join(root, RUN_SUMMARY), summary)
  return summary
}
