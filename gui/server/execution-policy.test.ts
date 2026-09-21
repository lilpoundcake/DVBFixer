import { describe, expect, it } from 'vitest'
import { COMMANDS } from './dvbfixer-spec'
import { API_EXECUTION_POLICY, executionPolicy } from './execution-policy'

describe('API execution policy', () => {
  it('classifies every public command exactly once', () => {
    expect(Object.keys(API_EXECUTION_POLICY).sort()).toEqual(COMMANDS.map(command => command.name).sort())
    expect(new Set(Object.values(API_EXECUTION_POLICY)))
      .toEqual(new Set(['synchronous-transform', 'managed-workflow']))
  })

  it('keeps long-running scientific operations asynchronous', () => {
    for (const command of ['prepare', 'minimize', 'model', 'zbs', 'parametrize']) {
      expect(executionPolicy(command)).toBe('managed-workflow')
    }
  })

  it('routes naming only through its dedicated synchronous contract', () => {
    expect(executionPolicy('atom-names')).toBe('synchronous-transform')
    expect(executionPolicy('diagnose')).toBe('managed-workflow')
    expect(executionPolicy('conect')).toBe('managed-workflow')
    expect(executionPolicy('renumber')).toBe('managed-workflow')
    expect(executionPolicy('missing')).toBeNull()
  })
})
