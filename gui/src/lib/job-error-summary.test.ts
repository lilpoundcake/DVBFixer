import { expect, it } from 'vitest'
import { jobErrorSummary } from './tool-output'

it('extracts a concise timeout cause instead of a full diagnostic summary', () => {
  const error = '=== dvbfixer 0.10.5 ===\n!!! ERROR !!! ERROR: external diffusion runner timed out after 900 seconds\n' +
    '!!!!\nDIAGNOSTIC SUMMARY: 1 error(s), 18 warning(s)\n[WARNING] input clash\n[ERROR] ERROR: external diffusion runner timed out after 900 seconds'
  expect(jobErrorSummary(error)).toBe('external diffusion runner timed out after 900 seconds')
})
it('bounds unknown multiline errors and does not invent validation outcomes', () => {
  expect(jobErrorSummary('Traceback\n' + 'x'.repeat(10000))).toBe('Traceback')
  expect(jobErrorSummary('x'.repeat(10000)).length).toBeLessThanOrEqual(240)
})
