import { spawn } from 'node:child_process'

export interface DvbfixerRunOptions {
  timeoutMs?: number
  maxOutputBytes?: number
  signal?: AbortSignal
}

function commandConfiguration(): { executable: string; prefix: string[] } {
  const executable = process.env.DVBFIXER_EXECUTABLE || 'dvbfixer'
  let prefix: unknown = []
  if (process.env.DVBFIXER_ARGS) {
    try { prefix = JSON.parse(process.env.DVBFIXER_ARGS) } catch {
      throw new Error('DVBFIXER_ARGS must be a JSON array of strings')
    }
  }
  if (!Array.isArray(prefix) || !prefix.every(item => typeof item === 'string')) {
    throw new Error('DVBFIXER_ARGS must be a JSON array of strings')
  }
  return { executable, prefix }
}

function appendBounded(current: string, chunk: Buffer, limit: number): string {
  if (Buffer.byteLength(current) >= limit) return current
  const remaining = limit - Buffer.byteLength(current)
  const text = chunk.subarray(0, remaining).toString()
  return current + text + (chunk.length > remaining ? '\n[output truncated]\n' : '')
}

export function runDvbfixerArgs(
  command: string,
  args: string[],
  cwd = process.cwd(),
  options: DvbfixerRunOptions = {},
): Promise<{ code: number; stdout: string; stderr: string }> {
  const { executable, prefix } = commandConfiguration()
  const commandArgs = [...prefix, command, ...args]
  const timeoutMs = options.timeoutMs ?? Number(process.env.DVBFIXER_TIMEOUT_MS || 30 * 60_000)
  const maxOutputBytes = options.maxOutputBytes ?? Number(process.env.DVBFIXER_MAX_OUTPUT_BYTES || 10 * 1024 * 1024)
  return new Promise((resolve) => {
    const child = spawn(executable, commandArgs, { cwd, detached: process.platform !== 'win32' })
    let stdout = ''
    let stderr = ''
    let settled = false
    let terminationReason = ''
    const finish = (code: number) => {
      if (settled) return
      settled = true
      clearTimeout(timer)
      options.signal?.removeEventListener('abort', abort)
      resolve({ code, stdout, stderr: terminationReason ? `${stderr}\n${terminationReason}`.trim() : stderr })
    }
    const terminate = (reason: string) => {
      if (settled || child.killed) return
      terminationReason = reason
      if (process.platform !== 'win32' && child.pid) {
        try { process.kill(-child.pid, 'SIGTERM') } catch { child.kill('SIGTERM') }
      } else child.kill('SIGTERM')
    }
    const abort = () => terminate('DVBfixer run cancelled')
    const timer = setTimeout(() => terminate(`DVBfixer run timed out after ${timeoutMs} ms`), timeoutMs)
    options.signal?.addEventListener('abort', abort, { once: true })
    if (options.signal?.aborted) abort()
    child.stdout.on('data', (data: Buffer) => { stdout = appendBounded(stdout, data, maxOutputBytes) })
    child.stderr.on('data', (data: Buffer) => { stderr = appendBounded(stderr, data, maxOutputBytes) })
    child.on('error', error => { terminationReason = String(error); finish(-1) })
    child.on('close', code => finish(code ?? -1))
  })
}

export function runDvbfixer(
  command: string,
  inputFile: string,
  outputFile: string,
  extraArgs: string[],
): Promise<{ code: number; stdout: string; stderr: string }> {
  return runDvbfixerArgs(command, [inputFile, '-o', outputFile, ...extraArgs])
}
