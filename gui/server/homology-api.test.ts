import { afterEach, describe, expect, it } from 'vitest'
import fs from 'node:fs'
import os from 'node:os'
import path from 'node:path'
import { materializeModelInputs, parseFasta, type HomologyProject } from './homology-api'

const temporaryDirectories: string[] = []

function temporaryDirectory(): string {
  const directory = fs.mkdtempSync(path.join(os.tmpdir(), 'dvbfixer-gui-test-'))
  temporaryDirectories.push(directory)
  return directory
}

afterEach(() => {
  for (const directory of temporaryDirectories.splice(0)) {
    fs.rmSync(directory, { recursive: true, force: true })
  }
})

describe('homology project materialization', () => {
  it('parses multiple target chains and rejects duplicate ids', () => {
    expect(parseFasta('>A\nAGS\n>B\nTT\n')).toEqual([
      { id: 'A', sequence: 'AGS' }, { id: 'B', sequence: 'TT' },
    ])
    expect(() => parseFasta('>A\nAA\n>A\nGG\n')).toThrow(/unique/)
  })

  it('turns a selected template span into a fragment and padded PIR row', () => {
    const root = temporaryDirectory()
    const run = path.join(root, 'run')
    fs.mkdirSync(run)
    const pdb = path.join(root, 'template.pdb')
    const atom = (serial: number, residue: string, number: number) =>
      `ATOM  ${String(serial).padStart(5)}  CA  ${residue} A${String(number).padStart(4)}       0.000   0.000   0.000  1.00 20.00           C  \n`
    fs.writeFileSync(pdb, atom(1, 'ALA', 1) + atom(2, 'GLY', 2) + atom(3, 'SER', 3))
    const now = new Date().toISOString()
    const project: HomologyProject = {
      version: 1, id: 'project', name: 'target', targetFasta: '>A\nAGS\n', engine: 'mafft',
      templates: [{ id: 'tpl', file: 'template.pdb', chain: 'A', targetChain: 'A' }],
      alignmentGroups: [{
        chainId: 'A',
        rows: [
          { id: 'A', kind: 'target', sequence: 'AGS' },
          { id: 'tpl', kind: 'template', templateId: 'tpl', sequence: 'AGS' },
        ],
        masks: { tpl: [{ start: 1, end: 3 }] },
      }],
      modelOptions: {}, createdAt: now, updatedAt: now,
    }
    const result = materializeModelInputs(root, project, run)
    expect(result.fragments).toHaveLength(1)
    const fragment = fs.readFileSync(result.fragments[0], 'utf8')
    expect(fragment).not.toContain('ALA A   1')
    expect(fragment).toContain('GLY A   2')
    expect(fragment).toContain('SER A   3')
    const pir = fs.readFileSync(result.pir, 'utf8')
    expect(pir).toContain('-GS*')
    expect(pir).toContain('>P1;target')
    expect(pir).toContain('AGS*')
  })
})
