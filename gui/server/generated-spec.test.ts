import { describe, expect, it } from 'vitest'
import { GENERATED_COMMANDS } from './generated-dvbfixer-spec'
import { buildArgs } from './api-plugin'

describe('generated DVBfixer command schema', () => {
  it('contains every current command and semantic groups', () => {
    expect(GENERATED_COMMANDS).toHaveLength(21)
    expect(GENERATED_COMMANDS.map(command => command.name)).toContain('salign')
    expect(GENERATED_COMMANDS.map(command => command.name)).toContain('msa')
    const prepare = GENERATED_COMMANDS.find(command => command.name === 'prepare')!
    expect(prepare.groups.map(group => group.name)).toContain('Mutations')
    expect(prepare.flags.map(field => field.flag)).toContain('--smiles')
    expect(prepare.flags.find(field => field.flag === '--strip-heterogens')?.default).toBe(false)
    const model = GENERATED_COMMANDS.find(command => command.name === 'model')!
    expect(model.flags.find(field => field.flag === '--pin-input')).toMatchObject({
      default: true, falseFlag: '--no-pin-input',
    })
  })

  it('expands repeatable fixed-cardinality arguments correctly', () => {
    expect(buildArgs('pull', {
      '--bond': 'A:1:SG B:2:SG,C:3:SG D:4:SG',
    })).toEqual([
      '--bond', 'A:1:SG', 'B:2:SG', '--bond', 'C:3:SG', 'D:4:SG',
    ])
  })
})
