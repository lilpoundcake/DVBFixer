import { afterEach, expect, it } from 'vitest'
import fs from 'node:fs'
import os from 'node:os'
import path from 'node:path'
import { compactRun } from './run-storage'
import { inspectionInventory, inspectionStructures, readInspectionStructure } from './run-inspection'

const directories: string[] = []
afterEach(() => { for (const dir of directories.splice(0)) fs.rmSync(dir, { recursive: true, force: true }) })
function root() { const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'inspection-')); directories.push(dir); return dir }

it('previews exact archived bytes without extraction, including long paths and preceding large entries', async () => {
  const dir = root()
  const file = `private/${'nested/'.repeat(15)}seed-7.pdb`
  fs.mkdirSync(path.dirname(path.join(dir, file)), { recursive: true })
  fs.writeFileSync(path.join(dir, 'aaa.log'), 'log\n'.repeat(100000))
  fs.writeFileSync(path.join(dir, file), 'ATOM rejected coordinates\n')
  await compactRun(dir, 'model', 'failed', [])
  const inventory = inspectionInventory(path.join(dir, 'run-summary.json'))
  expect(inspectionStructures(inventory)).toEqual([file])
  const before = fs.readdirSync(dir)
  expect((await readInspectionStructure(path.join(dir, 'diagnostics.tar.gz'), inventory, file)).toString()).toBe('ATOM rejected coordinates\n')
  expect(fs.readdirSync(dir)).toEqual(before)
  inventory.diagnostics!.files.find(entry => entry.file === file)!.sha256 = 'a'.repeat(64)
  await expect(readInspectionStructure(path.join(dir, 'diagnostics.tar.gz'), inventory, file)).rejects.toThrow('checksum')
  await expect(readInspectionStructure(path.join(dir, 'diagnostics.tar.gz'), inventory, '../secret.pdb')).rejects.toThrow('not found')
})

it('lists rejected bundle coordinates but not accepted outputs or unrelated files', () => {
  expect(inspectionStructures({ status: 'succeeded', outputs: ['result/accepted.pdb', 'result/rejected/seed-19.pdb', 'report.json'] }))
    .toEqual(['result/rejected/seed-19.pdb'])
  expect(inspectionStructures({ status: 'validation_failed', outputs: ['result/seed-7.pdb'] })).toEqual(['result/seed-7.pdb'])
  expect(inspectionStructures({ status: 'failed', outputs: [] })).toEqual([])
})

it.each(['../escape.pdb', '/absolute.pdb', 'nested/../../escape.pdb', 'C:/escape.pdb'])('refuses unsafe inventory paths: %s', file => {
  const dir = root(), summary = path.join(dir, 'run-summary.json')
  fs.writeFileSync(summary, JSON.stringify({ version: 1, outputs: [], diagnostics: {
    archive: 'diagnostics.tar.gz', files: [{ file, bytes: 10, sha256: 'a'.repeat(64) }],
  } }))
  expect(() => inspectionInventory(summary)).toThrow('invalid run summary')
})
