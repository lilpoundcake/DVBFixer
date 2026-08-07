import { afterEach, describe, expect, it } from 'vitest'
import fs from 'node:fs'
import os from 'node:os'
import path from 'node:path'
import { materializeModelInputs, parseFasta, parseSequenceArtifact, type HomologyProject } from './homology-api'

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

  it('extracts target sequences from PDB SEQRES and mmCIF atom_site records', () => {
    const root = temporaryDirectory()
    const pdb = path.join(root, 'target.pdb')
    fs.writeFileSync(pdb, 'SEQRES   1 A    3  ALA GLY SER\nEND\n')
    expect(parseSequenceArtifact(pdb)).toEqual([{ id: 'A', sequence: 'AGS', length: 3 }])

    const cif = path.join(root, 'target.cif')
    fs.writeFileSync(cif, `data_target
loop_
_atom_site.group_PDB
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.auth_comp_id
ATOM A 1 ALA
ATOM A 2 GLY
ATOM B 1 SER
#
`)
    expect(parseSequenceArtifact(cif)).toEqual([
      { id: 'A', sequence: 'AG', length: 2 },
      { id: 'B', sequence: 'S', length: 1 },
    ])
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
        masks: { tpl: [{ start: 1, end: 3 }] }, maskModes: { tpl: 'ranges' },
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

  it('distinguishes a cleared template selection from whole-template default', () => {
    const root = temporaryDirectory()
    const run = path.join(root, 'run')
    fs.mkdirSync(run)
    const pdb = path.join(root, 'template.pdb')
    fs.writeFileSync(pdb, 'ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n')
    const now = new Date().toISOString()
    const project: HomologyProject = {
      version: 1, id: 'cleared', name: 'cleared', targetFasta: '>A\nA\n', engine: 'mafft',
      templates: [{ id: 'tpl', file: 'template.pdb', chain: 'A', targetChain: 'A' }],
      alignmentGroups: [{ chainId: 'A', rows: [
        { id: 'A', kind: 'target', sequence: 'A' },
        { id: 'tpl', kind: 'template', templateId: 'tpl', sequence: 'A' },
      ], masks: { tpl: [] }, maskModes: { tpl: 'none' } }],
      modelOptions: {}, createdAt: now, updatedAt: now,
    }
    expect(() => materializeModelInputs(root, project, run)).toThrow(/select no residues/)
  })

  it('rejects a template row mapped to the wrong coordinate sequence', () => {
    const root = temporaryDirectory()
    const run = path.join(root, 'run')
    fs.mkdirSync(run)
    fs.writeFileSync(path.join(root, 'template.pdb'),
      'ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n')
    const now = new Date().toISOString()
    const project: HomologyProject = {
      version: 1, id: 'mismatch', name: 'mismatch', targetFasta: '>A\nG\n', engine: 'mafft',
      templates: [{ id: 'tpl', file: 'template.pdb', chain: 'A', targetChain: 'A' }],
      alignmentGroups: [{ chainId: 'A', rows: [
        { id: 'A', kind: 'target', sequence: 'G' },
        { id: 'tpl', kind: 'template', templateId: 'tpl', sequence: 'G' },
      ], masks: {}, maskModes: { tpl: 'all' } }],
      modelOptions: {}, createdAt: now, updatedAt: now,
    }
    expect(() => materializeModelInputs(root, project, run)).toThrow(/first difference at residue 1/)
  })

  it('cuts Modeller fragments from structurally fitted coordinates', () => {
    const root = temporaryDirectory()
    const run = path.join(root, 'run')
    fs.mkdirSync(run)
    const atom = (x: string) => `ATOM      1  CA  ALA A   1      ${x}   0.000   0.000  1.00 20.00           C  \n`
    fs.writeFileSync(path.join(root, 'original.pdb'), atom(' 0.000'))
    fs.writeFileSync(path.join(root, 'fitted.pdb'), atom('42.000'))
    const now = new Date().toISOString()
    const project: HomologyProject = {
      version: 1, id: 'fitted', name: 'fitted', targetFasta: '>A\nA\n', engine: 'mafft',
      templates: [{ id: 'tpl', file: 'original.pdb', fittedFile: 'fitted.pdb', chain: 'A', targetChain: 'A' }],
      alignmentGroups: [{ chainId: 'A', rows: [
        { id: 'A', kind: 'target', sequence: 'A' },
        { id: 'tpl', kind: 'template', templateId: 'tpl', sequence: 'A' },
      ], masks: {}, maskModes: { tpl: 'all' } }],
      modelOptions: {}, createdAt: now, updatedAt: now,
    }
    const result = materializeModelInputs(root, project, run)
    expect(fs.readFileSync(result.fragments[0], 'utf8')).toContain('42.000')
  })

  it('combines selected parts into one coordinate-preserving mosaic known', () => {
    const root = temporaryDirectory()
    const run = path.join(root, 'run')
    fs.mkdirSync(run)
    const atom = (serial: number, residue: string, number: number, x: number) =>
      `ATOM  ${String(serial).padStart(5)}  CA  ${residue} A${String(number).padStart(4)}    ${x.toFixed(3).padStart(8)}   0.000   0.000  1.00 20.00           C  \n`
    fs.writeFileSync(path.join(root, 'left.pdb'), atom(1, 'ALA', 1, 1) + atom(2, 'GLY', 2, 2))
    fs.writeFileSync(path.join(root, 'right.pdb'), atom(1, 'ALA', 1, 31) + atom(2, 'GLY', 2, 32))
    const now = new Date().toISOString()
    const project: HomologyProject = {
      version: 1, id: 'mosaic', name: 'mosaic', targetFasta: '>A\nAG\n', engine: 'mafft',
      templates: [
        { id: 'left', file: 'left.pdb', chain: 'A', targetChain: 'A' },
        { id: 'right', file: 'right.pdb', chain: 'A', targetChain: 'A' },
      ],
      alignmentGroups: [{ chainId: 'A', rows: [
        { id: 'A', kind: 'target', sequence: 'AG' },
        { id: 'left', kind: 'template', templateId: 'left', sequence: 'AG' },
        { id: 'right', kind: 'template', templateId: 'right', sequence: 'AG' },
      ], masks: { left: [{ start: 0, end: 1 }], right: [{ start: 1, end: 2 }] },
      maskModes: { left: 'ranges', right: 'ranges' } }],
      modelOptions: {}, createdAt: now, updatedAt: now,
    }
    const result = materializeModelInputs(root, project, run)
    expect(result.fragments).toHaveLength(1)
    const mosaic = fs.readFileSync(result.fragments[0], 'utf8')
    expect(mosaic).toContain('   1.000')
    expect(mosaic).toContain('  32.000')
    expect(fs.readFileSync(result.pir, 'utf8')).toContain('>P1;selected_template_mosaic')
  })
})
