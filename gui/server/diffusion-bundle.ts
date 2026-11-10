import crypto from 'node:crypto'
import fs from 'node:fs'
import path from 'node:path'
import type { DiffusionOutcomes } from './diffusion-outcomes'

function verifyArtifact(root: string, reference: { path: string; sha256: string }, maxBytes = 1_000_000_000): { artifact: string; data: Buffer } {
  if (!reference || typeof reference.path !== 'string' || typeof reference.sha256 !== 'string' ||
      path.isAbsolute(reference.path) || reference.path.split(/[\\/]/).includes('..')) {
    throw new Error('invalid diffusion artifact reference')
  }
  const artifact = path.resolve(root, reference.path)
  if (!artifact.startsWith(`${path.resolve(root)}${path.sep}`)) throw new Error('invalid diffusion artifact path')
  let component = root
  for (const part of reference.path.split(/[\\/]/)) {
    component = path.join(component, part)
    if (fs.lstatSync(component).isSymbolicLink()) throw new Error('diffusion artifact must not be a symlink')
  }
  const info = fs.lstatSync(artifact)
  if (!info.isFile() || info.nlink !== 1 || info.size > maxBytes) throw new Error('invalid diffusion artifact file')
  const data = fs.readFileSync(artifact)
  if (crypto.createHash('sha256').update(data).digest('hex') !== reference.sha256) throw new Error('diffusion artifact digest mismatch')
  return { artifact, data }
}

/** Transport verification only: Python owns scientific validation and ranking. */
export function readDiffusionPublication(indexPath: string): {
  status: 'success' | 'validation_failed'
  primary: string | null
  candidates: string[]
  sidecars: Record<string, string>
  seedOutcomes?: DiffusionOutcomes
  rejectedCoordinates?: string[]
} {
  const root = path.dirname(indexPath)
  const stat = fs.lstatSync(indexPath)
  if (!stat.isFile() || stat.isSymbolicLink() || stat.size > 10_000_000) throw new Error('invalid diffusion bundle index')
  const index = JSON.parse(fs.readFileSync(indexPath, 'utf8'))
  if (index.schema_version !== 2 || !['success', 'validation_failed'].includes(index.status) ||
      !Array.isArray(index.candidates) || index.candidates.length < 1 || index.candidates.length > 1000) {
    throw new Error('invalid diffusion publication schema')
  }
  const candidates: string[] = []
  const sidecars: Record<string, string> = {}
  const seen = new Set<string>()
  const requestDigests = new Set<string>()
  for (const candidate of index.candidates) {
    if (candidate.validation_passed !== (index.status === 'success')) throw new Error('inconsistent diffusion validation status')
    for (const label of ['pdb', 'dat', 'provenance', 'sampler_trace']) {
      const { artifact, data } = verifyArtifact(root, candidate[label])
      if (seen.has(artifact)) throw new Error('invalid diffusion artifact path')
      seen.add(artifact)
      if (label === 'provenance') {
        if (data.length > 16_000_000) throw new Error('diffusion provenance exceeds size limit')
        const provenance = JSON.parse(data.toString('utf8'))
        if (!provenance || typeof provenance !== 'object' || Array.isArray(provenance)) throw new Error('invalid diffusion provenance')
        requestDigests.add(provenance.request_sha256)
        if (provenance.artifacts !== undefined) {
          if (!Array.isArray(provenance.artifacts) || provenance.artifacts.length > 1000) throw new Error('invalid diffusion artifact manifest')
          // Custody only: Python owns scientific meaning and acceptance. This
          // includes staged raw/stage1/evidence/report files when declared.
          for (const reference of provenance.artifacts) verifyArtifact(root, reference)
        }
      }
      if (label === 'pdb') candidates.push(artifact)
    }
    sidecars[path.resolve(root, candidate.pdb.path)] = path.resolve(root, candidate.dat.path)
  }
  const result = { status: index.status, primary: index.status === 'success' ? candidates[0] : null, candidates, sidecars }
  if (index.seed_outcomes === undefined) return result
  const { artifact: reportPath, data } = verifyArtifact(root, index.seed_outcomes, 10_000_000)
  if (data.length > 10_000_000 || seen.has(reportPath)) throw new Error('invalid seed outcome report')
  const report = JSON.parse(data.toString('utf8'))
  if (report.version !== 1 || requestDigests.size !== 1 || !/^[a-f0-9]{64}$/.test(report.request_sha256) ||
    !requestDigests.has(report.request_sha256) || !Array.isArray(report.requested_seeds) ||
    !Array.isArray(report.outcomes) || report.outcomes.length < 1 || report.outcomes.length > 1000 ||
    report.requested_seeds.length !== report.outcomes.length ||
    report.requested_seeds.some((seed: unknown) => !Number.isSafeInteger(seed) || (seed as number) < 0) ||
    new Set(report.requested_seeds).size !== report.requested_seeds.length) throw new Error('invalid seed outcome inventory or binding')
  const published = new Map(index.candidates.map((candidate: any) => [candidate.seed, candidate])) as Map<number, any>
  if (published.size !== index.candidates.length || [...published.keys()].some(seed => !report.requested_seeds.includes(seed))) {
    throw new Error('seed report omits or repeats published candidates')
  }
  const rejectedCoordinates: string[] = []
  const outcomes: DiffusionOutcomes['outcomes'] = report.outcomes.map((row: any, i: number) => {
    if (!row || row.seed !== report.requested_seeds[i]) throw new Error('seed outcome order or identity mismatch')
    const candidate = published.get(row.seed)
    if (row.status === 'not_returned') {
      if (candidate || row.validation !== null || row.coordinate !== null || row.candidate_id !== null) throw new Error('not-returned seed claims a candidate')
      return { seed: row.seed, candidateId: null, status: row.status, findings: [], hardGateFailures: [], omittedFindings: 0 }
    }
    const summary = row.validation
    if (!['accepted', 'rejected'].includes(row.status) || typeof row.candidate_id !== 'string' || !row.candidate_id ||
      !summary || summary.passed !== (row.status === 'accepted') || !Array.isArray(summary.hard_gate_failures) ||
      !Array.isArray(summary.warnings) || [...summary.hard_gate_failures, ...summary.warnings].some(text => typeof text !== 'string') ||
      Boolean(summary.hard_gate_failures.length) !== (row.status === 'rejected')) throw new Error('inconsistent seed validation outcome')
    if (candidate) {
      if (candidate.candidate_id !== row.candidate_id || candidate.validation_passed !== summary.passed || row.coordinate !== null) {
        throw new Error('seed outcome differs from published candidate')
      }
    } else {
      if (row.status !== 'rejected' || !row.coordinate) throw new Error('unpublished seed must retain rejected coordinates')
      if (row.coordinate.path !== `rejected/seed-${row.seed}.pdb`) throw new Error('rejected coordinate has an unexpected path')
      const { artifact } = verifyArtifact(root, row.coordinate)
      if (seen.has(artifact)) throw new Error('rejected coordinate repeats an existing artifact')
      seen.add(artifact)
      rejectedCoordinates.push(artifact)
    }
    const findings = summary.warnings.filter((text: string) => text.startsWith('Final validation: '))
      .map((text: string) => text.slice('Final validation: '.length)) as string[]
    return { seed: row.seed, candidateId: row.candidate_id, status: row.status,
      findings: findings.slice(0, 20).map(text => text.slice(0, 2048)),
      hardGateFailures: summary.hard_gate_failures.slice(0, 100), omittedFindings: Math.max(0, findings.length - 20) }
  })
  return { ...result, seedOutcomes: { report: reportPath, outcomes }, rejectedCoordinates }
}
