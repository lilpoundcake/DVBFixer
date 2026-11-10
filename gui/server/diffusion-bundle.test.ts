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

function withOutcomes() {
  const { root, index } = bundle()
  const document = JSON.parse(fs.readFileSync(index, 'utf8'))
  const write = (name: string, data: string) => {
    fs.mkdirSync(path.dirname(path.join(root, name)), { recursive: true })
    fs.writeFileSync(path.join(root, name), data)
    return { path: name, sha256: crypto.createHash('sha256').update(data).digest('hex') }
  }
  for (const [i, candidate] of document.candidates.entries()) {
    candidate.seed = [7, 19][i]
    candidate.provenance = write(candidate.provenance.path, JSON.stringify({ request_sha256: 'a'.repeat(64) }))
  }
  const report = { version: 1, request_sha256: 'a'.repeat(64), requested_seeds: [42, 7, 19, 101], outcomes: [
    { seed: 42, candidate_id: 'failed-42', status: 'rejected', validation: { passed: false, hard_gate_failures: ['angle'],
      warnings: ['Resolved raw issue', 'Final validation: d/ASP28A: angle = 79.1°; allowed 95–125°.'] },
      coordinate: write('rejected/seed-42.pdb', 'REJECTED COORDINATES\n') },
    ...document.candidates.map((candidate: any) => ({ seed: candidate.seed, candidate_id: candidate.candidate_id,
      status: 'accepted', validation: { passed: true, hard_gate_failures: [], warnings: [] }, coordinate: null })),
    { seed: 101, candidate_id: null, status: 'not_returned', validation: null, coordinate: null },
  ] }
  const save = () => {
    document.seed_outcomes = write('seed-outcomes.json', JSON.stringify(report))
    fs.writeFileSync(index, JSON.stringify(document))
  }
  save()
  return { root, index, document, report, save }
}

it('transports all seed outcomes while selecting only the accepted ranked structure', () => {
  const { root, index } = withOutcomes()
  const result = readDiffusionPublication(index)
  expect(result.primary).toBe(path.join(root, 'ranked-z.pdb'))
  expect(result.candidates).toHaveLength(2)
  expect(result.rejectedCoordinates).toEqual([path.join(root, 'rejected/seed-42.pdb')])
  expect(result.seedOutcomes?.report).toBe(path.join(root, 'seed-outcomes.json'))
  expect(result.seedOutcomes?.outcomes.map(row => [row.seed, row.status])).toEqual([
    [42, 'rejected'], [7, 'accepted'], [19, 'accepted'], [101, 'not_returned'],
  ])
  expect(result.seedOutcomes?.outcomes[0].findings).toEqual(['d/ASP28A: angle = 79.1°; allowed 95–125°.'])
})

it.each(['binding', 'omit', 'promote', 'path', 'repeat', 'corrupt'])('rejects invalid seed outcome custody (%s)', mutation => {
  const { root, index, report, save } = withOutcomes()
  if (mutation === 'binding') report.request_sha256 = 'b'.repeat(64)
  if (mutation === 'omit') report.outcomes.pop()
  if (mutation === 'promote') report.outcomes[0].status = 'accepted'
  if (mutation === 'path') report.outcomes[0].coordinate.path = '../escape.pdb'
  if (mutation === 'repeat') report.requested_seeds[1] = report.requested_seeds[0]
  save()
  if (mutation === 'corrupt') fs.writeFileSync(path.join(root, 'rejected/seed-42.pdb'), 'tampered')
  expect(() => readDiffusionPublication(index)).toThrow()
})
