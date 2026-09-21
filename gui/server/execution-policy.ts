import { COMMANDS } from './dvbfixer-spec'

export type ApiExecutionPolicy = 'synchronous-transform' | 'managed-workflow'

/**
 * Public execution policy for every CLI operation.
 *
 * Naming is the only bounded synchronous V1 transform. All other commands use
 * persisted managed jobs, including quick report/rewriting commands, so the
 * HTTP server has one cancellation, provenance, and artifact-publication path.
 */
export const API_EXECUTION_POLICY: Readonly<Record<string, ApiExecutionPolicy>> = Object.freeze(
  Object.fromEntries(COMMANDS.map(command => [
    command.name,
    command.name === 'atom-names' ? 'synchronous-transform' : 'managed-workflow',
  ] satisfies [string, ApiExecutionPolicy])) as Record<string, ApiExecutionPolicy>,
)

export function executionPolicy(command: string): ApiExecutionPolicy | null {
  return API_EXECUTION_POLICY[command] || null
}
