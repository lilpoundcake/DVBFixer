import { describe, expect, it } from 'vitest'
import { diffusionSeeds } from './diffusion-seeds'
import { buildArgs } from '../../server/command-args'
import { workflowValues } from './workflow-backend-options'

describe('diffusion seeds', () => {
  it.each(['1, 17, 148', '1 17 148', '1,17\n148', [1, 17, 148]])('normalizes %j in order', value => {
    expect(diffusionSeeds(value)).toEqual([1, 17, 148])
    expect(buildArgs('model', { '--diffusion-seed': value })).toEqual([
      '--diffusion-seed', '1', '--diffusion-seed', '17', '--diffusion-seed', '148',
    ])
  })
  it.each([undefined, null, '', '  ', []])('defaults %j to 7', value => expect(diffusionSeeds(value)).toEqual([7]))
  it.each(['-1', '1.5', '1, x', '1 01', ',', '9007199254740992', true, [null]])('rejects %j', value => {
    expect(() => diffusionSeeds(value)).toThrow()
    expect(() => buildArgs('model', { '--diffusion-seed': value })).toThrow()
  })
  it('normalizes active seeds and ignores inactive stale input', () => {
    expect(workflowValues('model', { '--backend': 'diffusion', '--diffusion-seed': '1 17' })['--diffusion-seed']).toEqual([1, 17])
    expect(workflowValues('model', { '--backend': 'modeller', '--diffusion-seed': 'bad' })).toEqual({ '--backend': 'modeller' })
  })
})
