import { useEffect, useRef } from 'react'
import { MolScriptBuilder as MS } from 'molstar/lib/mol-script/language/builder'
import { Structure } from 'molstar/lib/mol-model/structure'
import { setSubtreeVisibility } from 'molstar/lib/mol-plugin/behavior/static/state'
import type { PluginUIContext } from 'molstar/lib/mol-plugin-ui/context'
import { useStructureStore } from '../stores/structureStore'
import { useWorkspaceStore } from '../stores/workspaceStore'
import { useAltlocStore } from '../stores/altlocStore'
import { altlocSerialExpression, altlocStateOpacity, displayedAltloc } from '../lib/altloc-preview'
import { structureMetaFromArtifact } from '../lib/workspace-metadata'

interface Overlay {
  plugin: PluginUIContext
  refs: string[]
  originals: Array<{ ref: string; hidden: boolean }>
}

async function clearOverlay(overlay: Overlay | null): Promise<void> {
  if (!overlay) return
  const { plugin } = overlay
  const update = plugin.build()
  for (const ref of overlay.refs) if (plugin.state.data.cells.has(ref)) update.delete(ref)
  await update.commit()
  for (const original of overlay.originals) if (plugin.state.data.cells.has(original.ref)) {
    setSubtreeVisibility(plugin.state.data, original.ref, original.hidden)
  }
}

/** Owns temporary components in the ordinary A viewer; never creates another viewer. */
export function useAltlocViewer() {
  const plugin = useStructureStore(state => state.plugin)
  const workspaceId = useWorkspaceStore(state => state.active?.id || '')
  const report = useAltlocStore(state => state.sessions[workspaceId]?.report)
  const request = useAltlocStore(state => state.sessions[workspaceId]?.viewerRequest || 0)
  const queue = useRef<Promise<unknown>>(Promise.resolve())
  const generation = useRef(0)
  const overlay = useRef<Overlay | null>(null)
  const loaded = useRef<{ plugin: PluginUIContext; digest: string; model: number; ref: string } | null>(null)
  const baseline = useRef({ workspaceId, request: 0 })
  useEffect(() => {
    const counter = generation
    const token = ++counter.current
    if (baseline.current.workspaceId !== workspaceId) baseline.current = { workspaceId, request }
    const explicitRequest = request > baseline.current.request
    queue.current = queue.current.catch(() => {}).then(async () => {
      await clearOverlay(overlay.current)
      overlay.current = null
      if (!plugin || !report?.normalizedPdb || !explicitRequest || token !== generation.current) return
      const session = useAltlocStore.getState().sessions[workspaceId]
      const workspace = useWorkspaceStore.getState().active
      if (!session || session.viewerRequest !== request || workspace?.id !== workspaceId) return
      const row = report.residues.find(item => item.key === session.activeKey)
      const artifact = workspace.artifacts.find(item => item.id === session.input)
      if (!row || !artifact) return
      let ref = loaded.current?.ref
      if (loaded.current?.plugin !== plugin || loaded.current.digest !== report.normalizedSha256 ||
          loaded.current.model !== row.model || !ref || !plugin.state.data.cells.has(ref) ||
          useStructureStore.getState().fileName !== artifact.file) {
        await plugin.clear()
        if (token !== generation.current || useWorkspaceStore.getState().active?.id !== workspaceId) return
        const data = await plugin.builders.data.rawData({ data: report.normalizedPdb, label: artifact.name })
        const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb')
        const model = await plugin.builders.structure.createModel(trajectory, { modelIndex: row.model - 1 })
        const structure = await plugin.builders.structure.createStructure(model)
        await plugin.builders.structure.representation.applyPreset(structure, 'auto')
        ref = structure.ref
        loaded.current = { plugin, digest: report.normalizedSha256, model: row.model, ref }
        if (token !== generation.current || useWorkspaceStore.getState().active?.id !== workspaceId) return
        useStructureStore.getState().setFileName(artifact.file)
        useStructureStore.getState().setMeta(structureMetaFromArtifact(artifact))
        useWorkspaceStore.getState().update({ primaryFile: artifact.file })
      }
      const state = plugin.managers.structure.hierarchy.current.structures.find(item => item.cell.transform.ref === ref)
      if (!state || token !== generation.current) return
      const owner: Overlay = { plugin, refs: [], originals: state.components.map(item => ({ ref: item.cell.transform.ref, hidden: !!item.cell.state.isHidden })) }
      overlay.current = owner
      const focused = altlocSerialExpression([...row.sharedSerials, ...row.states.flatMap(item => item.serials)])
      const make = async (expression: ReturnType<typeof altlocSerialExpression>, key: string) => {
        const component = await plugin.builders.structure.tryCreateComponentFromExpression(ref!, expression, `altloc-${key}`)
        if (component) owner.refs.push(component.ref)
        return component
      }
      const context = await make(MS.struct.modifier.exceptBy({ 0: MS.struct.generator.all(), by: focused }), 'context')
      if (context) await plugin.builders.structure.representation.addRepresentation(context, { type: 'cartoon', color: 'chain-id', typeParams: { alpha: .6 } })
      const surroundings = await make(MS.struct.modifier.exceptBy({
        0: MS.struct.modifier.includeSurroundings({ 0: focused, radius: 5, 'as-whole-residues': true }), by: focused,
      }), 'surroundings')
      if (surroundings) await plugin.builders.structure.representation.addRepresentation(surroundings, { type: 'ball-and-stick', color: 'element-symbol', typeParams: { alpha: .2 } })
      const shared = await make(altlocSerialExpression(row.sharedSerials), 'shared')
      if (shared) await plugin.builders.structure.representation.addRepresentation(shared, { type: 'ball-and-stick', color: 'element-symbol' })
      const colors = [0x36a2eb, 0xff9f40, 0x9966ff, 0x4bc0c0, 0xff6384]
      const tapped = displayedAltloc(row, session.focusLabel)
      let focus = shared
      for (const [index, alternative] of row.states.entries()) {
        if (session.mode === 'selected' && alternative.label !== tapped) continue
        const component = await make(altlocSerialExpression(alternative.serials), alternative.label)
        if (!component) continue
        await plugin.builders.structure.representation.addRepresentation(component, { type: 'ball-and-stick', color: 'uniform',
          colorParams: { value: colors[index % colors.length] },
          typeParams: { alpha: altlocStateOpacity(session.mode, alternative.label, tapped) } })
        if (alternative.label === tapped || !focus) focus = component
      }
      if (token !== generation.current) return
      for (const original of owner.originals) setSubtreeVisibility(plugin.state.data, original.ref, true)
      if (focus?.obj) plugin.managers.camera.focusLoci(Structure.toStructureElementLoci(focus.obj.data), { minRadius: 8, extraRadius: 5 })
    }).catch(async cause => {
      await clearOverlay(overlay.current).catch(() => {})
      overlay.current = null
      if (token === generation.current) useAltlocStore.getState().patch(workspaceId, {
        error: `3D altloc display failed: ${cause instanceof Error ? cause.message : 'unknown error'}. Choices are unchanged.`,
      })
    })
    return () => { counter.current++ }
  }, [plugin, workspaceId, report, request])
  useEffect(() => () => { void queue.current.finally(() => clearOverlay(overlay.current)).catch(() => {}) }, [])
}
