import { spawn } from 'node:child_process'

export function runDvbfixerArgs(
  command: string,
  args: string[],
  cwd = process.cwd(),
): Promise<{ code: number; stdout: string; stderr: string }> {
  const cmd = process.env.DVBFIXER_CMD || 'dvbfixer'
  const parts = cmd.split(/\s+/).filter(Boolean)
  const exe = parts[0]
  const commandArgs = [...parts.slice(1), command, ...args]
  return new Promise((resolve) => {
    const child = spawn(exe, commandArgs, { cwd })
    let stdout = ''
    let stderr = ''
    child.stdout.on('data', (data) => { stdout += data.toString() })
    child.stderr.on('data', (data) => { stderr += data.toString() })
    child.on('error', (error) => resolve({ code: -1, stdout, stderr: `${stderr}\n${String(error)}` }))
    child.on('close', (code) => resolve({ code: code ?? -1, stdout, stderr }))
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
