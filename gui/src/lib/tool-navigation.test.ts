import { describe, expect, it } from 'vitest'
import { COMMANDS } from '../../server/dvbfixer-spec'

describe('shared tool navigation registry', () => {
  const tools = (category: string) => COMMANDS.filter(command => command.category === category).map(command => command.name)
  it('keeps the requested utility order and removes Analysis', () => {
    const utilities = tools('Utilities')
    expect(utilities[0]).toBe('diagnose')
    expect(utilities.slice(utilities.indexOf('conect'), utilities.indexOf('pull') + 1)).toEqual(['conect', 'puppet', 'pull'])
    expect(tools('Analysis')).toEqual([])
  })
  it('keeps Altloc before Model and Convert before Transplant/Cluster', () => {
    const preparation = tools('Structure preparation')
    expect(preparation[preparation.indexOf('model') - 1]).toBe('altlocs')
    expect(tools('Glycoprotein preparation')).toEqual(['convert', 'transplant', 'cluster'])
  })
})
