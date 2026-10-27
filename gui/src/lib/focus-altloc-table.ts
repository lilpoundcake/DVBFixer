import { Actions, DockLocation, Model, TabNode, TabSetNode } from 'flexlayout-react'
import { PINNED_ALTLOC_TAB } from '../app-panels'

/** Put the editor beside the ordinary viewer, not inside the launcher tabset. */
export function focusAltlocTable(model: Model): void {
  let table: string | undefined
  let destination: string | undefined
  let first: string | undefined
  let clash: TabNode | undefined
  let tableNode: TabNode | undefined
  model.visitNodes(node => {
    if (node instanceof TabSetNode) first ??= node.getId()
    if (!(node instanceof TabNode)) return
    if (node.getComponent() === 'altlocs') { table ??= node.getId(); tableNode ??= node }
    if (node.getComponent() === 'clashes') clash ??= node
    if (node.getComponent() === 'elements' || node.getComponent() === 'interactions') destination = node.getParent()?.getId()
  })
  if (clash?.getParent()) {
    const parent = clash.getParent()!
    destination = parent.getId()
    const siblings = parent.getChildren()
    if (table && tableNode && siblings[siblings.indexOf(clash) + 1] !== tableNode) {
      // FlexLayout adjusts the insertion index itself when moving within a tabset.
      model.doAction(Actions.moveNode(table, destination, DockLocation.CENTER, siblings.indexOf(clash) + 1))
    }
    if (!table) model.doAction(Actions.addNode(PINNED_ALTLOC_TAB, destination, DockLocation.CENTER, siblings.indexOf(clash) + 1, true))
  } else if (!table && (destination || first)) {
    model.doAction(Actions.addNode(PINNED_ALTLOC_TAB, destination || first!, DockLocation.CENTER, -1, true))
  }
  if (table) {
    model.doAction(Actions.updateNodeAttributes(table, { enableDrag: false, enableClose: false }))
    model.doAction(Actions.selectTab(table))
  }
}

export function focusAltlocViewer(model: Model): void {
  let viewer: string | undefined
  let destination: string | undefined
  model.visitNodes(node => {
    if (node instanceof TabNode && node.getComponent() === 'viewer') viewer ??= node.getId()
    if (node instanceof TabNode && node.getComponent() === 'dvbfixer') destination = node.getParent()?.getId()
  })
  if (viewer) model.doAction(Actions.selectTab(viewer))
  else if (destination) model.doAction(Actions.addNode(
    { type: 'tab', name: '3D Structure', component: 'viewer' }, destination, DockLocation.CENTER, -1, true,
  ))
}
