import { expect, it } from 'vitest'
import { workflowValues } from './workflow-backend-options'

it('only exposes context overrides for explicitly selected Protenix', () => {
  for (const model of ['', 'protpardelle']) {
    expect(workflowValues('model', { '--backend': 'diffusion', '--diffusion-model': model, '--diffusion-context': 'full' }))
      .not.toHaveProperty('--diffusion-context')
  }
  expect(workflowValues('model', { '--backend': 'diffusion', '--diffusion-model': 'protenix', '--diffusion-context': 'full' }))
    .toHaveProperty('--diffusion-context', 'full')
  expect(workflowValues('model', { '--backend': 'diffusion', '--diffusion-model': 'protenix', '--diffusion-context': 'none' }))
    .not.toHaveProperty('--diffusion-context')
})

it('removes inactive Modeller defaults when selecting diffusion', () => {
  expect(workflowValues('model', { '--backend': 'diffusion', '--md-level': 'fast', '--num-loops': 2, '--diffusion-seed': '7,11' }))
    .toEqual({ '--backend': 'diffusion', '--diffusion-seed': '7,11' })
})

it('keeps ZBS prepare backend separate from its model backend', () => {
  expect(workflowValues('zbs', { '--model-backend': 'diffusion', '--backend': 'legacy', '--num-output': 1, '--diffusion-model': 'protenix' }))
    .toEqual({ '--model-backend': 'diffusion', '--backend': 'legacy', '--diffusion-model': 'protenix' })
})
