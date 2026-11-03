import { expect, it } from 'vitest'
import { formatJobTimestamp } from './managed-jobs'

it('formats job times without ISO separators or milliseconds, in the requested local zone', () => {
  expect(formatJobTimestamp('2026-11-03T17:05:05.823Z', 'en-GB', 'UTC')).toBe('03/11/2026, 17:05:05')
  expect(formatJobTimestamp('2026-11-03T17:05:05.823Z', 'en-GB', 'Europe/Moscow')).toBe('03/11/2026, 20:05:05')
  expect(formatJobTimestamp(null)).toBe('Not started')
  expect(formatJobTimestamp('bad')).toBe('Unknown time')
})
