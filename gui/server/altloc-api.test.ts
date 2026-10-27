import { afterEach, describe, expect, it, vi } from 'vitest'
import fs from 'node:fs'
import os from 'node:os'
import path from 'node:path'
import crypto from 'node:crypto'
import { executeAltlocSelection } from './altloc-api'
import { loadWorkspace, saveWorkspace, workspaceRoot, updateWorkspace } from './workspace-api'
import type { AltlocReport, AltlocRequest } from '../src/lib/altloc-types'

const temps: string[] = []
const sha = (bytes: Buffer | string) => crypto.createHash('sha256').update(bytes).digest('hex')
const pdb = 'ATOM      1  CA AALA C   1       1.000   2.000   3.000  0.60 20.00           C  \n'
const selected = pdb.slice(0, 16) + ' ' + pdb.slice(17)
const key = '[1,"C","1","",1]'
const request: AltlocRequest = { inputArtifactId: 'input', action: 'analyze', policy: 'occupancy', fallback: 'error' }
function workspace() {
  const root = fs.mkdtempSync(path.join(os.tmpdir(), 'altloc-api-test-')); temps.push(root)
  const directory = workspaceRoot(root, 'w')
  fs.mkdirSync(path.join(directory, 'files'), { recursive: true })
  fs.writeFileSync(path.join(directory, 'files/input.pdb'), pdb)
  saveWorkspace(root, { version: 1, revision: 0, id: 'w', name: 'Workspace', createdAt: '', updatedAt: '',
    primaryFile: null, secondaryFile: null, toolState: {}, artifacts: [{ id: 'input', file: 'files/input.pdb', name: 'Input', kind: 'structure' }] })
  return root
}
function runner(hook?: (report: AltlocReport, output?: string) => void) {
  return vi.fn(async (_command: string, args: string[]) => {
    const source = args[0]
    const output = args.includes('-o') ? args[args.indexOf('-o') + 1] : undefined
    const plan = args.includes('--selection-plan') ? JSON.parse(fs.readFileSync(args[args.indexOf('--selection-plan') + 1], 'utf8')) : null
    const report: AltlocReport = {
      schemaVersion: 1, sourceSha256: sha(fs.readFileSync(source)), normalizedSha256: sha(pdb),
      policy: args[args.indexOf('--policy') + 1] as AltlocReport['policy'], prefer: null, fallback: 'error',
      overrides: plan?.selections || {}, residues: [{ key, model: 1, chain: 'C', residueNumber: '1', insertionCode: '', occurrence: 1,
        selected: 'A', method: 'occupancy', sharedSerials: [], warnings: [], states: [{ label: 'A', residueName: 'ALA', occupancy: .6, atomNames: ['CA'], serials: [1], atomCount: 1 }] }],
      warnings: [], unresolved: [], inputAtoms: 1, removedAtoms: 0, outputSha256: sha(selected), outputWritten: !!output,
      toolVersion: 'test', chainMap: {}, ...(output ? {} : { normalizedPdb: pdb }),
    }
    if (output) fs.writeFileSync(output, selected)
    hook?.(report, output)
    fs.writeFileSync(args[args.indexOf('--report-json') + 1], JSON.stringify(report))
    return { code: 0, stdout: 'Selected', stderr: '', started: true }
  })
}
afterEach(() => temps.splice(0).forEach(root => fs.rmSync(root, { recursive: true, force: true })))

describe('altloc CLI workspace boundary', () => {
  it('analysis registers nothing; apply registers only the digest-verified new structure', async () => {
    const root = workspace()
    const run = runner()
    const analysis = await executeAltlocSelection(root, 'w', request, run)
    expect(analysis.report.residues).toHaveLength(1)
    expect(loadWorkspace(root, 'w').artifacts).toHaveLength(1)
    const result = await executeAltlocSelection(root, 'w', { ...request, action: 'apply',
      plan: { schemaVersion: 1, sourceSha256: sha(pdb), selections: { [key]: 'A' } } }, run)
    const latest = loadWorkspace(root, 'w')
    expect(latest.artifacts).toHaveLength(2)
    expect(latest.artifacts[1].altlocProvenance?.sourceSha256).toBe(sha(pdb))
    expect(fs.readFileSync(path.join(workspaceRoot(root, 'w'), result.outputArtifact!.file), 'utf8')).toBe(selected)
    expect(fs.readFileSync(path.join(workspaceRoot(root, 'w'), 'files/input.pdb'), 'utf8')).toBe(pdb)
  })
  it('rejects stale source plans before invoking Python', async () => {
    const root = workspace(); const run = runner()
    await expect(executeAltlocSelection(root, 'w', { ...request, action: 'apply', plan: { schemaVersion: 1, sourceSha256: '0'.repeat(64), selections: {} } }, run)).rejects.toThrow('changed')
    expect(run).not.toHaveBeenCalled()
  })
  it.each(['digest', 'source', 'artifact', 'permission', 'symlink', 'missing-plan'])('publishes nothing after %s failure', async failure => {
    const root = workspace()
    const run = runner((report, output) => {
      if (failure === 'digest') report.outputSha256 = '0'.repeat(64)
      if (failure === 'source') fs.appendFileSync(path.join(workspaceRoot(root, 'w'), 'files/input.pdb'), 'END\n')
      if (failure === 'artifact') updateWorkspace(root, 'w', latest => { latest.artifacts = []; return latest })
      if (failure === 'symlink' && output) { fs.unlinkSync(output); fs.symlinkSync(path.join(workspaceRoot(root, 'w'), 'files/input.pdb'), output) }
    })
    let authorizations = 0
    await expect(executeAltlocSelection(root, 'w', { ...request, action: 'apply', ...(failure === 'missing-plan' ? {} : {
      plan: { schemaVersion: 1 as const, sourceSha256: sha(pdb), selections: { [key]: 'A' } },
    }) }, run, () => { if (++authorizations > 1 && failure === 'permission') throw new Error('permission revoked') })).rejects.toThrow()
    expect(loadWorkspace(root, 'w').artifacts.filter(item => item.command === 'altlocs')).toHaveLength(0)
    expect(fs.readdirSync(path.join(workspaceRoot(root, 'w'), 'runs'))).toEqual([])
  })
  it('rebases publication onto the current manifest without dropping concurrent artifacts', async () => {
    const root = workspace()
    await executeAltlocSelection(root, 'w', { ...request, action: 'apply', plan: { schemaVersion: 1, sourceSha256: sha(pdb), selections: { [key]: 'A' } } }, runner(() => {
      updateWorkspace(root, 'w', latest => { latest.artifacts.push({ id: 'other', file: 'files/other.json', name: 'Other', kind: 'artifact' }); return latest })
    }))
    expect(loadWorkspace(root, 'w').artifacts.map(item => item.id)).toContain('other')
    expect(loadWorkspace(root, 'w').artifacts).toHaveLength(3)
  })
  it('rejects incomplete reviewed plans even when the automatic CLI policy resolves the input', async () => {
    const root = workspace()
    await expect(executeAltlocSelection(root, 'w', { ...request, action: 'apply',
      plan: { schemaVersion: 1, sourceSha256: sha(pdb), selections: {} } }, runner())).rejects.toThrow('every alternate residue')
    expect(loadWorkspace(root, 'w').artifacts).toHaveLength(1)
  })
  it('detaches child-created hardlinks and never follows a child-created log symlink', async () => {
    const root = workspace()
    const external = path.join(root, 'external.pdb')
    fs.writeFileSync(external, selected)
    const result = await executeAltlocSelection(root, 'w', { ...request, action: 'apply',
      plan: { schemaVersion: 1, sourceSha256: sha(pdb), selections: { [key]: 'A' } } }, runner((_report, output) => {
      fs.unlinkSync(output!)
      fs.linkSync(external, output!)
      fs.symlinkSync(external, path.join(path.dirname(output!), 'stdout.log'))
    }))
    fs.writeFileSync(external, 'changed externally')
    expect(fs.readFileSync(path.join(workspaceRoot(root, 'w'), result.outputArtifact!.file), 'utf8')).toBe(selected)
    expect(fs.lstatSync(path.join(workspaceRoot(root, 'w'), path.dirname(result.outputArtifact!.file), 'stdout.log')).isSymbolicLink()).toBe(false)
  })
})
