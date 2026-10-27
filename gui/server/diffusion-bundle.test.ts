import crypto from 'node:crypto'
import fs from 'node:fs'
import os from 'node:os'
import path from 'node:path'
import { afterEach, expect, it } from 'vitest'
import { readDiffusionPublication } from './diffusion-bundle'

const roots: string[] = []
afterEach(() => { for (const root of roots.splice(0)) fs.rmSync(root, { recursive: true, force: true }) })

function bundle(passed = true) {
  const root = fs.mkdtempSync(path.join(os.tmpdir(), 'diffusion-publication-'))
  roots.push(root)
  const candidates = ['ranked-z', 'ranked-a'].map(candidate_id => {
    const references = Object.fromEntries(['pdb', 'dat', 'provenance', 'sampler_trace'].map(label => {
      const name = `${candidate_id}.${label}`
      const data = label === 'provenance' ? '{}' : label
      fs.writeFileSync(path.join(root, name), data)
      return [label, { path: name, sha256: crypto.createHash('sha256').update(data).digest('hex') }]
    }))
    return { candidate_id, validation_passed: passed, ...references }
  })
  const index = path.join(root, 'bundle.json')
  fs.writeFileSync(index, JSON.stringify({ schema_version: 2, status: passed ? 'success' : 'validation_failed', candidates }))
  return { root, index }
}

it('selects publication rank rather than filename order', () => {
  const { root, index } = bundle()
  expect(readDiffusionPublication(index).primary).toBe(path.join(root, 'ranked-z.pdb'))
})

it('retains rejected candidates for inspection without automatic primary output', () => {
  const { index } = bundle(false)
  expect(readDiffusionPublication(index)).toMatchObject({ status: 'validation_failed', primary: null })
})

it('rejects artifact corruption before workspace registration', () => {
  const { root, index } = bundle()
  fs.writeFileSync(path.join(root, 'ranked-a.dat'), 'corrupt')
  expect(() => readDiffusionPublication(index)).toThrow(/digest mismatch/)
})

it('verifies declared staged artifacts without interpreting geometry', () => {
  const { root, index } = bundle()
  const data = JSON.parse(fs.readFileSync(index, 'utf8'))
  fs.writeFileSync(path.join(root, 'raw.pdb'), 'raw coordinates')
  const provenance = JSON.stringify({ artifacts: [{ role: 'raw', path: 'raw.pdb',
    sha256: crypto.createHash('sha256').update('raw coordinates').digest('hex') }] })
  fs.writeFileSync(path.join(root, data.candidates[0].provenance.path), provenance)
  data.candidates[0].provenance.sha256 = crypto.createHash('sha256').update(provenance).digest('hex')
  fs.writeFileSync(index, JSON.stringify(data))
  expect(readDiffusionPublication(index).primary).toBeTruthy()
  fs.writeFileSync(path.join(root, 'raw.pdb'), 'tampered')
  expect(() => readDiffusionPublication(index)).toThrow(/digest mismatch/)
})
