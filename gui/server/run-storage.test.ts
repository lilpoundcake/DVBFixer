import { afterEach, describe, expect, it, vi } from 'vitest'
import fs from 'node:fs'
import os from 'node:os'
import path from 'node:path'
import crypto from 'node:crypto'
import { gunzipSync } from 'node:zlib'
import { compactRun, copyNativeRun, outputFiles, RUN_ARCHIVE, RUN_SUMMARY, runFiles } from './run-storage'

const roots: string[] = []
function root() { const directory = fs.mkdtempSync(path.join(os.tmpdir(), 'run-storage-')); roots.push(directory); return directory }
function write(directory: string, file: string, text = file) {
  fs.mkdirSync(path.dirname(path.join(directory, file)), { recursive: true })
  fs.writeFileSync(path.join(directory, file), text)
}
function unpack(directory: string) {
  const tar = gunzipSync(fs.readFileSync(path.join(directory, RUN_ARCHIVE)))
  const result: Record<string, Buffer> = {}
  let offset = 0
  while (tar.subarray(offset, offset + 512).some(byte => byte !== 0)) {
    const header = tar.subarray(offset, offset + 512)
    const field = (start: number, length: number) => header.subarray(start, start + length).toString().replace(/\0.*$/, '')
    const name = [field(345, 155), field(0, 100)].filter(Boolean).join('/')
    const size = parseInt(field(124, 12), 8)
    result[name] = tar.subarray(offset + 512, offset + 512 + size)
    offset += 512 + Math.ceil(size / 512) * 512
  }
  return result
}
afterEach(() => { vi.restoreAllMocks(); for (const directory of roots.splice(0)) fs.rmSync(directory, { recursive: true, force: true }) })

describe('shared GUI run storage', () => {
  it('preserves output sets, numbered candidates, sidecars and complete bundles', () => {
    const files = ['result/model.pdb', 'result/model.dat', 'result/bundle.json', 'result/_required.json',
      'private/motif.pdb', 'stdout.log', 'target_model_1.pdb', 'target_model_1.dat', 'target_model_prepare.pdb', '_input.pdb']
    expect(outputFiles(files, 'target_model.pdb', ['result'])).toEqual([
      'result/model.pdb', 'result/model.dat', 'result/bundle.json', 'result/_required.json', 'target_model_1.pdb', 'target_model_1.dat',
    ])
    expect(outputFiles(files, 'result')).toContain('result/_required.json')
    expect(outputFiles(['final.top', 'posre_Chain_A.itp', 'LIG.itp', 'interchain_ss.itp', 'amber.ff/ff.dat', 'scratch.txt'], 'final.top'))
      .toEqual(['final.top', 'posre_Chain_A.itp', 'LIG.itp', 'interchain_ss.itp', 'amber.ff/ff.dat'])
  })

  it.each(['succeeded', 'failed', 'cancelled', 'validation_failed'] as const)('archives and verifies diagnostics for %s without touching declared outputs or logs', async status => {
    const directory = root()
    write(directory, 'result/final.pdb', 'ATOM preserved\n')
    write(directory, 'result/final.dat', '{"schema":1}')
    write(directory, 'private/gap-1/motif.pdb', 'motif\n'.repeat(2000))
    write(directory, 'private/request.json', '{"mask":[true,false]}')
    write(directory, 'stdout.log', 'progress')
    write(directory, 'job.json', '{"id":"job"}')
    const outputs = ['result/final.pdb', 'result/final.dat']
    const summary = await compactRun(directory, 'model', status, outputs)
    expect(summary).toMatchObject({ status, outputs, diagnostics: { archive: RUN_ARCHIVE } })
    expect(summary.warning).toBeUndefined()
    expect(fs.readFileSync(path.join(directory, outputs[0]), 'utf8')).toBe('ATOM preserved\n')
    expect(fs.readFileSync(path.join(directory, 'stdout.log'), 'utf8')).toBe('progress')
    expect(fs.existsSync(path.join(directory, 'private'))).toBe(false)
    const contents = unpack(directory)
    for (const entry of summary.diagnostics!.files) {
      expect(crypto.createHash('sha256').update(contents[entry.file]).digest('hex')).toBe(entry.sha256)
      expect(contents[entry.file].length).toBe(entry.bytes)
    }
    expect(JSON.parse(contents['diagnostics-manifest.json'].toString()).files).toEqual(summary.diagnostics!.files)
    expect(JSON.parse(fs.readFileSync(path.join(directory, RUN_SUMMARY), 'utf8'))).toEqual(summary)
    expect(fs.statSync(path.join(directory, RUN_ARCHIVE)).size).toBeLessThan(2000)
  })

  it('keeps originals and records a warning if archive publication fails', async () => {
    const directory = root()
    write(directory, 'evidence.txt', 'do not lose this')
    vi.spyOn(fs, 'renameSync').mockImplementationOnce(() => { throw new Error('injected archive failure') })
    const summary = await compactRun(directory, 'prepare', 'failed', [])
    expect(summary.warning).toContain('Diagnostics kept unpacked')
    expect(fs.readFileSync(path.join(directory, 'evidence.txt'), 'utf8')).toBe('do not lose this')
    expect(fs.existsSync(path.join(directory, RUN_ARCHIVE))).toBe(false)
    expect(fs.readdirSync(directory).some(file => file.endsWith('.tmp'))).toBe(false)
  })

  it('does not follow symlinks or archive hardlinks', async () => {
    const directory = root(), outside = root()
    write(outside, 'secret', 'not run data')
    fs.symlinkSync(path.join(outside, 'secret'), path.join(directory, 'link'))
    expect(() => runFiles(directory)).toThrow(/link/)
    expect((await compactRun(directory, 'model', 'failed', [])).warning).toContain('link')
    expect(() => copyNativeRun(directory, root())).toThrow(/link/)
    expect(fs.readFileSync(path.join(outside, 'secret'), 'utf8')).toBe('not run data')
    const hardlinked = root()
    fs.linkSync(path.join(outside, 'secret'), path.join(hardlinked, 'linked-file'))
    expect((await compactRun(hardlinked, 'prepare', 'failed', [])).warning).toContain('link')
    expect(fs.existsSync(path.join(hardlinked, 'linked-file'))).toBe(true)
  })

  it('leaves unsafe cross-platform archive names unpacked', async () => {
    const directory = root()
    write(directory, '..\\outside.txt', 'keep safely on this filesystem')
    const summary = await compactRun(directory, 'model', 'failed', [])
    expect(summary.warning).toContain('unsafe diagnostic archive path')
    expect(fs.existsSync(path.join(directory, '..\\outside.txt'))).toBe(true)
    expect(fs.existsSync(path.join(directory, RUN_ARCHIVE))).toBe(false)
  })

  it('copies native failures but omits environment/cache trees', () => {
    const source = root(), destination = root()
    write(source, '.dvbfixer-run/raw/candidate.pdb')
    write(source, '.dvbfixer-run/.cache/huge.cache')
    write(source, '.dvbfixer-run/.home/settings')
    write(source, '.dvbfixer-run/.tmp/scratch')
    copyNativeRun(source, destination)
    expect(runFiles(destination)).toEqual(['.dvbfixer-run/raw/candidate.pdb'])
  })

  it('merges existing directories without overwriting any existing file', () => {
    const source = root(), destination = root(), conflict = root()
    write(source, 'nested/raw.pdb', 'new evidence')
    write(destination, 'nested/control.json', 'keep this control')
    copyNativeRun(source, destination)
    expect(fs.readFileSync(path.join(destination, 'nested/raw.pdb'), 'utf8')).toBe('new evidence')
    expect(fs.readFileSync(path.join(destination, 'nested/control.json'), 'utf8')).toBe('keep this control')
    write(conflict, 'nested/control.json', 'replacement')
    expect(() => copyNativeRun(conflict, destination)).toThrow()
    expect(fs.readFileSync(path.join(destination, 'nested/control.json'), 'utf8')).toBe('keep this control')
  })

  it('rejects copying into the source tree before creating a destination', () => {
    const source = root()
    write(source, 'raw.pdb')
    expect(() => copyNativeRun(source, path.join(source, 'nested'))).toThrow(/inside its source/)
    expect(fs.existsSync(path.join(source, 'nested'))).toBe(false)
  })

  it('does not overwrite an existing archive or create ambiguous manifest entries', async () => {
    const directory = root()
    write(directory, RUN_ARCHIVE, 'existing archive')
    write(directory, 'raw.pdb')
    const summary = await compactRun(directory, 'model', 'failed', [])
    expect(summary.warning).toContain('already exists')
    expect(fs.readFileSync(path.join(directory, RUN_ARCHIVE), 'utf8')).toBe('existing archive')
    expect(fs.existsSync(path.join(directory, 'raw.pdb'))).toBe(true)
  })
})
