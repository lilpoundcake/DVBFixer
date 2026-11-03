import { describe, expect, it } from 'vitest'
import { orderToolTabs } from './tool-tab-order'

describe('tool tab order', () => {
  it('places Diagnose before Cluster without modifying source metadata', () => {
    const commands = [{ name: 'cluster' }, { name: 'diagnose' }]
    expect(orderToolTabs(commands).map(command => command.name)).toEqual(['diagnose', 'cluster'])
    expect(commands.map(command => command.name)).toEqual(['cluster', 'diagnose'])
  })

  it('preserves unrelated tabs and an already-correct order', () => {
    const commands = [{ name: 'diagnose' }, { name: 'other' }, { name: 'cluster' }]
    expect(orderToolTabs(commands)).toEqual(commands)
    expect(orderToolTabs([{ name: 'prepare' }, { name: 'minimize' }])).toEqual([{ name: 'prepare' }, { name: 'minimize' }])
  })
  it('puts Diagnose first and Altloc immediately before Model', () => {
    expect(orderToolTabs([{ name: 'pull' }, { name: 'diagnose' }, { name: 'rename' }]).map(item => item.name))
      .toEqual(['diagnose', 'pull', 'rename'])
    expect(orderToolTabs([{ name: 'split' }, { name: 'model' }, { name: 'prepare' }, { name: 'altlocs' }]).map(item => item.name))
      .toEqual(['split', 'altlocs', 'model', 'prepare'])
  })
})
