import { expect, it } from 'vitest'
import { managedJobProgress, type ManagedJobRecord } from './managed-jobs'

it('reports stage-local counts only and falls back to indeterminate for legacy records', () => {
  const job = {} as ManagedJobRecord
  expect(managedJobProgress(job)).toEqual({ label: 'Waiting for stage events' })
  job.progress = { version: 1, stage: 'runner' }
  expect(managedJobProgress(job)).toEqual({ label: 'Runner (sampling / refinement)' })
  job.progress = { version: 1, stage: 'sampling', completed: 1, total: 4 }
  expect(managedJobProgress(job)).toEqual({ label: 'Sampling · 1/4', percent: 25 })
  job.progress.completed = 5
  expect(managedJobProgress(job)).toEqual({ label: 'Sampling' })
  job.progress = { version: 1, stage: 'refinement', seed: 19 }
  expect(managedJobProgress(job)).toEqual({ label: 'Refinement · Seed 19' })
  job.progress = { version: 1, stage: 'publication', completed: 1, total: 1 }
  expect(managedJobProgress(job)).toEqual({ label: 'Publication · 1/1', percent: 100 })
  expect(job.status).toBeUndefined() // No lifecycle or scientific outcome inferred.
})
