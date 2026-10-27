/** Transport only: selection policy and PDB materialization belong to Python. */
export type AltlocPolicy = 'manual' | 'occupancy' | 'first' | 'prefer'
export interface AltlocRow {
  key: string
  model: number
  chain: string
  residueNumber: string
  insertionCode: string
  occurrence: number
  states: Array<{ label: string; residueName: string; occupancy: number | null; atomNames: string[]; serials: number[]; atomCount: number }>
  selected: string | null
  method: string
  sharedSerials: number[]
  warnings: string[]
}
export interface AltlocReport {
  schemaVersion: 1
  sourceSha256: string
  normalizedSha256: string
  policy: AltlocPolicy
  prefer: string | null
  fallback: 'error' | 'occupancy' | 'first'
  overrides: Record<string, string>
  residues: AltlocRow[]
  warnings: string[]
  unresolved: string[]
  inputAtoms: number
  removedAtoms: number | null
  outputSha256: string | null
  outputWritten: boolean
  normalizedPdb?: string
  toolVersion: string
  chainMap: Record<string, string>
}
export interface AltlocRequest {
  inputArtifactId: string
  action: 'analyze' | 'apply'
  policy: AltlocPolicy
  prefer?: string
  fallback: 'error' | 'occupancy' | 'first'
  plan?: { schemaVersion: 1; sourceSha256: string; selections: Record<string, string> }
}
export interface AltlocResponse {
  report: AltlocReport
  reportFile?: string
  outputArtifact: { id: string; file: string; name: string; kind: 'structure' } | null
}
