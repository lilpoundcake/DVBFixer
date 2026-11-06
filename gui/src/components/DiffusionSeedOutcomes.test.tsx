import { renderToStaticMarkup } from 'react-dom/server'
import { expect, it } from 'vitest'
import { DiffusionSeedOutcomes } from './DiffusionSeedOutcomes'

it('shows every seed, final findings and a retained report without treating omissions as rejections', () => {
  const html = renderToStaticMarkup(<DiffusionSeedOutcomes diffusion={{ report: 'runs/result/seed-outcomes.json', outcomes: [
    { seed: 7, candidateId: 'pass-7', status: 'accepted', findings: [], hardGateFailures: [], omittedFindings: 0 },
    { seed: 19, candidateId: 'fail-19', status: 'rejected', findings: ['d/ASP28A: angle = 79.1°; limit 95–125°.'], hardGateFailures: ['angle'], omittedFindings: 0 },
    { seed: 101, candidateId: null, status: 'not_returned', findings: [], hardGateFailures: [], omittedFindings: 0 },
  ] }} onOpenReport={() => {}} />)
  expect(html).toContain('Seed 7: Accepted')
  expect(html).toContain('Seed 19: Rejected')
  expect(html).toContain('Seed 101: Not returned for final validation')
  expect(html).toContain('d/ASP28A')
  expect(html).toContain('79.1°')
  expect(html).toContain('inspection only')
  expect(html).toContain('Open per-seed validation report')
  expect(html).not.toContain('Seed 101: Rejected')
})

it('does not add empty candidate banners for legacy or non-diffusion runs', () => {
  expect(renderToStaticMarkup(<DiffusionSeedOutcomes />)).toBe('')
})
