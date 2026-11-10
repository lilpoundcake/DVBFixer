import { describe, expect, it } from 'vitest'
import { COMMANDS } from '../../server/dvbfixer-spec'
import { resolveToolCategory } from './tool-navigation'

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
    expect(tools('Glycans')).toEqual(['convert', 'transplant', 'cluster'])
    expect(tools('Glycoprotein preparation')).toEqual([])
  })
  it('restores renamed categories without discarding the selected tool group', () => {
    const categories = [...new Set(COMMANDS.map(command => command.category))]
    expect(resolveToolCategory('Glycoprotein preparation', categories)).toBe('Glycans')
    expect(resolveToolCategory('Analysis', categories)).toBe('Utilities')
    expect(resolveToolCategory('Glycans', categories)).toBe('Glycans')
    expect(resolveToolCategory('Unknown historical group', categories)).toBe(categories[0])
    expect(resolveToolCategory('Glycoprotein preparation', [])).toBe('')
  })
})
