import { create } from 'zustand'
import type { AltlocPolicy, AltlocReport, AltlocRequest } from '../lib/altloc-types'
import { chooseAltloc } from '../lib/altloc-selection'

export interface AltlocSession {
  input: string
  inputMode: 'auto' | 'manual'
  policy: AltlocPolicy
  prefer: string
  fallback: AltlocRequest['fallback']
  report: AltlocReport | null
  choices: Record<string, string>
  activeKey: string
  focusLabel: string
  mode: 'all' | 'selected' | 'compare'
  viewerRequest: number
  busy: boolean
  error: string
  result: string
  reportFile: string
  generation: number
}
export const EMPTY_ALTLOC_SESSION: AltlocSession = {
  input: '', inputMode: 'auto', policy: 'occupancy', prefer: 'A', fallback: 'error', report: null, choices: {},
  activeKey: '', focusLabel: '', mode: 'compare', viewerRequest: 0,
  busy: false, error: '', result: '', reportFile: '', generation: 0,
}
interface AltlocStore {
  sessions: Record<string, AltlocSession>
  patch: (workspaceId: string, patch: Partial<AltlocSession>) => void
  configure: (workspaceId: string, patch: Partial<Pick<AltlocSession, 'input' | 'inputMode' | 'policy' | 'prefer' | 'fallback'>>) => void
  choose: (workspaceId: string, key: string, label: string) => void
  focus: (workspaceId: string, key: string, label?: string) => void
}
export const useAltlocStore = create<AltlocStore>((set, get) => ({
  sessions: {},
  patch: (id, patch) => set(state => ({ sessions: { ...state.sessions, [id]: { ...(state.sessions[id] || EMPTY_ALTLOC_SESSION), ...patch } } })),
  configure: (id, patch) => {
    const current = get().sessions[id] || EMPTY_ALTLOC_SESSION
    if (current.busy) return
    get().patch(id, { ...patch, report: null, choices: {}, activeKey: '', focusLabel: '',
      error: '', result: '', reportFile: '', viewerRequest: 0, generation: current.generation + 1 })
  },
  choose: (id, key, label) => {
    const session = get().sessions[id]
    if (!session?.report || session.busy) return
    get().patch(id, { choices: chooseAltloc(session.report.residues, session.choices, key, label),
      viewerRequest: session.viewerRequest + 1 })
  },
  focus: (id, key, label) => {
    const session = get().sessions[id]
    const row = session?.report?.residues.find(row => row.key === key)
    if (!session || !row) return
    get().patch(id, { activeKey: key, focusLabel: label || (session.activeKey === key && session.focusLabel) || row.states[0]?.label || '', viewerRequest: session.viewerRequest + 1 })
  },
}))
