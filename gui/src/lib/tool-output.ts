import type { ManagedJobRecord } from './managed-jobs'

export interface ToolReportPreview {
  file: string
  text: string
  truncated: boolean
  error?: string
  /** A synchronous response, not a downloadable workspace report file. */
  inline?: boolean
}

const TEXT_REPORT = /\.(txt|json|md|log|csv|tsv|fasta|fa|faa|pir|aln|dat|top|itp|mdp|ya?ml|toml|xml)$/i
const DIAGNOSTIC_REPORT = /(?:^|[/. _-])diagnose\.(json|txt|md)$/i

export function diagnosticLine(raw: string): { text: string; severity?: 'error' | 'warning' } {
  // eslint-disable-next-line no-control-regex -- Render terminal text without ANSI codes.
  const text = raw.replace(/\u001b\[[0-9;]*m/g, '')
  const normalized = text.trim().replace(/^!!! (?:ERROR|WARNING) !!!\s*/, '').replace(/^\[(?:ERROR|WARNING)\]\s*/, '')
  if (/^ERROR\b/.test(normalized) && !/^ERROR:\s*0\b/.test(normalized)) return { text, severity: 'error' }
  // Historical runs already declare scientific rejection, but labelled it WARNING.
  if (/^(?:WARNING:\s*)?Diffusion candidate (?:REJECTED:|was published but failed independent DVBFixer validation)/.test(normalized)) return { text, severity: 'error' }
  if (/^WARNING\b/.test(normalized)) return { text, severity: 'warning' }
  return { text }
}

/** Presentation only: surface explicit Python errors, never infer geometry. */
export function outputErrors(...streams: string[]): string[] {
  const errors = new Set<string>()
  for (const stream of streams) {
    for (const raw of stream.split('\n')) {
      const line = diagnosticLine(raw)
      if (line.severity !== 'error') continue
      errors.add(line.text.trim().replace(/^!!! ERROR !!!\s*/, '').replace(/^\[ERROR\]\s*/, '').slice(0, 4096))
      if (errors.size >= 100) return [...errors]
    }
  }
  return [...errors]
}

/** Present Python's explicit notices; do not calculate or reinterpret geometry. */
export function anchorRefinementNotices(...streams: string[]): string[] {
  const notices = new Set<string>()
  for (const stream of streams) {
    for (const raw of stream.split('\n')) {
      // eslint-disable-next-line no-control-regex -- ANSI color sequences intentionally begin with ESC.
      const line = raw.replace(/\u001b\[[0-9;]*m/g, '')
      const start = line.search(/Observed anchor .+ coordinates modified:|Observed-anchor B weighting uses uniform fallback:/)
      if (start >= 0) notices.add(line.slice(start, start + 4096).trim())
      if (notices.size >= 1000) return [...notices]
    }
  }
  return [...notices]
}

export function toolReportFiles(job: ManagedJobRecord): string[] {
  const primary = job.outputFile && TEXT_REPORT.test(job.outputFile) ? [job.outputFile] : []
  const diagnostics = (job.provenance?.outputs || []).map(output => output.file)
    .filter(file => DIAGNOSTIC_REPORT.test(file))
  return [...new Set([...primary, ...diagnostics])].slice(0, 3)
}

/** Read only a small preview; never download a large topology/report just to summarize it. */
async function readPreview(response: Response, maxBytes = 24 * 1024): Promise<{ text: string; truncated: boolean }> {
  if (!response.body) return { text: '', truncated: false }
  const reader = response.body.getReader()
  const decoder = new TextDecoder()
  let text = ''
  let remaining = maxBytes
  let truncated = false
  try {
    while (true) {
      const { value, done } = await reader.read()
      if (done) break
      const kept = value.subarray(0, remaining)
      text += decoder.decode(kept, { stream: true })
      remaining -= kept.length
      if (value.length > kept.length || !remaining) {
        truncated = true
        await reader.cancel()
        break
      }
    }
    if (!truncated) text += decoder.decode()
  } finally { reader.releaseLock() }
  // JSON formatting is presentation only; do not reinterpret diagnostic findings.
  if (!truncated && (text.trimStart().startsWith('{') || text.trimStart().startsWith('['))) {
    try { text = JSON.stringify(JSON.parse(text), null, 2) } catch { /* retain original text */ }
  }
  const lines = text.split('\n')
  return { text: lines.slice(0, 45).join('\n'), truncated: truncated || lines.length > 45 }
}

export async function loadToolReports(
  job: ManagedJobRecord,
  fetchFile: (file: string) => Promise<Response>,
): Promise<ToolReportPreview[]> {
  return Promise.all(toolReportFiles(job).map(async file => {
    try {
      const response = await fetchFile(file)
      if (!response.ok) throw new Error(`HTTP ${response.status}`)
      return { file, ...await readPreview(response) }
    } catch (error) {
      if (error instanceof Error && error.name === 'AbortError') throw error
      return { file, text: '', truncated: false, error: error instanceof Error ? error.message : String(error) }
    }
  }))
}
