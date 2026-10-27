import { afterEach, describe, expect, it, vi } from 'vitest'
import fs from 'node:fs'
import os from 'node:os'
import path from 'node:path'
import crypto from 'node:crypto'
import { EventEmitter } from 'node:events'
import { Readable } from 'node:stream'
import { registerHomologyApi, type HomologyProject } from './homology-api'
import { ensureWorkspaceMigration, listWorkspaces, loadWorkspace, workspaceRoot } from './workspace-api'

const { runner } = vi.hoisted(() => ({ runner: vi.fn() }))
vi.mock('./dvbfixer-runner', () => ({ runDvbfixerArgs: runner }))
const roots: string[] = []
afterEach(() => { runner.mockReset(); for (const root of roots.splice(0)) fs.rmSync(root, { recursive: true, force: true }) })

async function request(dataRoot: string, route: string, body: unknown) {
  let middleware: any
  registerHomologyApi({ middlewares: { use: (_route: string, handler: any) => { middleware = handler } } } as any, dataRoot)
  const incoming = Readable.from([JSON.stringify(body)]) as any
  incoming.method = 'POST'; incoming.url = route; incoming.headers = { 'content-type': 'application/json' }
  let text = ''
  const response = Object.assign(new EventEmitter(), { statusCode: 200, writableFinished: false,
    setHeader: () => {}, end(value: string) { text = value; this.writableFinished = true } })
  await middleware(incoming, response, () => { throw new Error('unhandled route') })
  return { status: response.statusCode, body: JSON.parse(text) }
}

describe('Homology shared run storage', () => {
  it.each(['model', 'diffusion', 'failed-diffusion', 'align', 'salign'])('packages %s without registering intermediate structures', async mode => {
    const dataRoot = fs.mkdtempSync(path.join(os.tmpdir(), 'homology-runs-')); roots.push(dataRoot)
    ensureWorkspaceMigration(dataRoot)
    const workspace = listWorkspaces(dataRoot)[0]
    const root = workspaceRoot(dataRoot, workspace.id)
    fs.writeFileSync(path.join(root, 'template.pdb'), 'END\n')
    fs.writeFileSync(path.join(root, 'other.pdb'), 'END\n')
    const project: HomologyProject & { workspaceId: string } = {
      workspaceId: workspace.id, version: 1, id: 'project', name: 'target', targetFasta: '>A\nAGS\n', engine: 'mafft',
      templates: [
        { id: 'one', file: 'template.pdb', chain: 'A', targetChain: 'A' },
        { id: 'two', file: 'other.pdb', chain: 'A', targetChain: 'A' },
      ], alignmentGroups: [], modelOptions: mode.includes('diffusion') ? { '--backend': 'diffusion' } : {},
      createdAt: '2026-01-01T00:00:00Z', updatedAt: '2026-01-01T00:00:00Z',
    }
    runner.mockImplementation(async (command: string, args: string[], cwd: string) => {
      const output = args[args.indexOf('-o') + 1]
      if (mode.includes('diffusion')) {
        const native = args[args.indexOf('--diffusion-work-parent') + 1]
        fs.mkdirSync(path.join(native, '.private', 'raw'), { recursive: true })
        fs.writeFileSync(path.join(native, '.private', 'raw', 'candidate.pdb'), 'raw')
        if (mode === 'diffusion') {
          const bundle = path.join(native, 'result'); fs.mkdirSync(bundle)
          const references: Record<string, unknown> = {}
          for (const label of ['pdb', 'dat', 'provenance', 'sampler_trace']) {
            const data = label === 'provenance' ? '{}' : label
            fs.writeFileSync(path.join(bundle, `final.${label}`), data)
            references[label] = { path: `final.${label}`, sha256: crypto.createHash('sha256').update(data).digest('hex') }
          }
          fs.writeFileSync(path.join(bundle, 'bundle.json'), JSON.stringify({ schema_version: 2, status: 'success', candidates: [
            { candidate_id: 'one', validation_passed: true, ...references },
          ] }))
        }
      } else if (command === 'homology') {
        fs.writeFileSync(`${output}.pdb`, 'END\n'); fs.writeFileSync(`${output}.dat`, '{}')
      } else if (command === 'msa') {
        fs.writeFileSync(output, '>A\nAGS\n>template_1_source\nAGS\n>template_2_source\nAGS\n')
      } else if (command === 'salign') {
        fs.writeFileSync(output, 'alignment')
        const fit = args[args.indexOf('--fit-dir') + 1]; fs.mkdirSync(fit, { recursive: true })
        for (const i of [1, 2]) fs.writeFileSync(path.join(fit, `template_${i}_A_fit.pdb`), 'END\n')
      }
      fs.writeFileSync(path.join(cwd, 'scientific-intermediate.tmp'), 'evidence')
      return { code: mode === 'failed-diffusion' ? 2 : 0, stdout: 'progress', stderr: '' }
    })
    const response = await request(dataRoot, mode === 'align' || mode === 'salign' ? `/${mode}` : '/model', project)
    expect(response.status).toBe(mode === 'failed-diffusion' ? 500 : 200)
    const artifacts = loadWorkspace(dataRoot, workspace.id).artifacts
    expect(artifacts.some(item => item.file.endsWith('/run-summary.json'))).toBe(true)
    expect(artifacts.some(item => item.file.endsWith('/diagnostics.tar.gz'))).toBe(true)
    expect(artifacts.some(item => /scientific-intermediate|\/raw\//.test(item.file))).toBe(false)
    if (mode === 'failed-diffusion') expect(artifacts.every(item => item.kind !== 'structure')).toBe(true)
    else if (mode === 'align') expect(artifacts.some(item => item.file.endsWith('alignment_A.fasta'))).toBe(true)
    else expect(artifacts.some(item => item.kind === 'structure')).toBe(true)
    for (const artifact of artifacts) expect(fs.existsSync(path.join(root, artifact.file))).toBe(true)
  })
})
