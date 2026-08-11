import { afterEach, describe, expect, it } from 'vitest'
import { runDvbfixerArgs } from './dvbfixer-runner'

const originalEnvironment = {
  executable: process.env.DVBFIXER_EXECUTABLE,
  args: process.env.DVBFIXER_ARGS,
}

afterEach(() => {
  for (const [key, value] of Object.entries({
    DVBFIXER_EXECUTABLE: originalEnvironment.executable,
    DVBFIXER_ARGS: originalEnvironment.args,
  })) {
    if (value === undefined) delete process.env[key]
    else process.env[key] = value
  }
})

function useNode(script: string): void {
  process.env.DVBFIXER_EXECUTABLE = process.execPath
  process.env.DVBFIXER_ARGS = JSON.stringify(['-e', script])
}

describe('DVBfixer runner', () => {
  it('keeps the executable separate from its JSON argument prefix', async () => {
    useNode('console.log(JSON.stringify(process.argv.slice(1)))')
    const result = await runDvbfixerArgs('doctor', ['--format', 'json'])
    expect(result.code).toBe(0)
    expect(JSON.parse(result.stdout.trim())).toEqual(['doctor', '--format', 'json'])
  })

  it('caps captured process output', async () => {
    useNode("process.stdout.write('x'.repeat(4096))")
    const result = await runDvbfixerArgs('doctor', [], process.cwd(), { maxOutputBytes: 128 })
    expect(Buffer.byteLength(result.stdout)).toBeLessThan(180)
    expect(result.stdout).toContain('[output truncated]')
  })

  it('terminates a timed-out process', async () => {
    useNode('setInterval(() => {}, 1000)')
    const result = await runDvbfixerArgs('doctor', [], process.cwd(), { timeoutMs: 30 })
    expect(result.code).not.toBe(0)
    expect(result.stderr).toContain('timed out')
  })
})
