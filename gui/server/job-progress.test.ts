import { describe, expect, it } from 'vitest'
import { createJobProgressDecoder, type JobProgress } from './job-progress'

describe('versioned job progress pipe', () => {
  it('decodes fragmented frames and multiple events independently of logs', () => {
    const events: JobProgress[] = []
    const decode = createJobProgressDecoder(event => events.push(event))
    decode(Buffer.from('{"version":1,"sta'))
    decode(Buffer.from('ge":"sampling"}\n{"version":1,"stage":"refinement","completed":1,"total":2,"seed":19}\n'))
    expect(events).toEqual([{ version: 1, stage: 'sampling' }, { version: 1, stage: 'refinement', completed: 1, total: 2, seed: 19 }])
  })

  it('ignores malformed, oversized, unknown-version and inconsistent events, then recovers', () => {
    const events: JobProgress[] = []
    const decode = createJobProgressDecoder(event => events.push(event))
    for (const event of [
      { version: 2, stage: 'sampling' }, { version: 1, stage: 'accepted' },
      { version: 1, stage: 'sampling', completed: 1 },
      { version: 1, stage: 'sampling', completed: 4, total: 3 },
      { version: 1, stage: 'sampling', completed: 0, total: 0 },
      { version: 1, stage: 'sampling', completed: 0.5, total: 2 },
      { version: 1, stage: 'sampling', completed: 0, total: 2 ** 53 },
      { version: 1, stage: 'sampling', scientificAcceptance: true },
    ]) decode(Buffer.from(`${JSON.stringify(event)}\n`))
    decode(Buffer.from('not json\n' + 'x'.repeat(65536)))
    decode(Buffer.from('\n{"version":1,"stage":"validation"}\n'))
    expect(events).toEqual([{ version: 1, stage: 'validation' }])
  })

  it('bounds the number of delivered events even for an overactive emitter', () => {
    let count = 0
    const decode = createJobProgressDecoder(() => count++)
    decode(Buffer.from('{"version":1,"stage":"sampling"}\n'.repeat(5000)))
    decode(Buffer.from('{"version":1,"stage":"publication"}\n'))
    expect(count).toBe(4096)
  })
})
