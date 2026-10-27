import { Model, TabNode, TabSetNode } from 'flexlayout-react'
import { expect, it } from 'vitest'
import { focusAltlocTable, focusAltlocViewer } from './focus-altloc-table'

it('places one table in the data tabset and reveals the ordinary viewer separately', () => {
  const model = Model.fromJson({ global: {}, borders: [], layout: { type: 'row', children: [
    { type: 'tabset', id: 'main', children: [{ type: 'tab', component: 'viewer', id: 'viewer-a' }, { type: 'tab', component: 'dvbfixer' }] },
    { type: 'tabset', id: 'data', children: [{ type: 'tab', component: 'elements' }, { type: 'tab', component: 'clashes', id: 'clashes' }] },
  ] } })
  focusAltlocTable(model); focusAltlocTable(model); focusAltlocViewer(model)
  const tables: TabNode[] = []
  model.visitNodes(node => { if (node instanceof TabNode && node.getComponent() === 'altlocs') tables.push(node) })
  expect(tables).toHaveLength(1)
  expect(tables[0].getParent()?.getId()).toBe('data')
  const siblings = tables[0].getParent()!.getChildren()
  expect(siblings[siblings.indexOf(tables[0]) - 1].getId()).toBe('clashes')
  expect(tables[0].isEnableDrag()).toBe(false)
  expect(tables[0].isEnableClose()).toBe(false)
  expect((model.getNodeById('main') as TabSetNode).getSelectedNode()?.getId()).toBe('viewer-a')
  expect((model.getNodeById('data') as TabSetNode).getSelectedNode()).toBe(tables[0])
})

it('moves a pre-existing editor beside Clashes instead of creating another tab', () => {
  const model = Model.fromJson({ global: {}, borders: [], layout: { type: 'row', children: [
    { type: 'tabset', id: 'main', children: [{ type: 'tab', component: 'altlocs', id: 'existing' }, { type: 'tab', component: 'dvbfixer' }] },
    { type: 'tabset', id: 'data', children: [{ type: 'tab', component: 'clashes', id: 'clashes' }, { type: 'tab', component: 'elements' }] },
  ] } })
  focusAltlocTable(model)
  const siblings = model.getNodeById('data')!.getChildren()
  expect(siblings.map(node => node.getId()).slice(0, 2)).toEqual(['clashes', 'existing'])
})

it.each([true, false])('reorders the existing tab in the same tabset (before Clashes=%s)', before => {
  const table = { type: 'tab' as const, component: 'altlocs', id: 'existing' }
  const clash = { type: 'tab' as const, component: 'clashes', id: 'clashes' }
  const elements = { type: 'tab' as const, component: 'elements', id: 'elements' }
  const model = Model.fromJson({ global: {}, borders: [], layout: { type: 'row', children: [
    { type: 'tabset', id: 'data', children: before ? [table, clash, elements] : [clash, elements, table] },
  ] } })
  focusAltlocTable(model)
  focusAltlocTable(model)
  expect(model.getNodeById('data')!.getChildren().map(node => node.getId())).toEqual(['clashes', 'existing', 'elements'])
})
