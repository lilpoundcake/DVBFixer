import { renderToStaticMarkup } from 'react-dom/server'
import { expect, it } from 'vitest'
import { ToolGroupMenu } from './ToolGroupMenu'

it('visibly identifies the group control as a dropdown', () => {
  const html = renderToStaticMarkup(<ToolGroupMenu commands={[]} category="Utilities" onSelect={() => {}} />)
  expect(html).toContain('Workflow group: Utilities')
  expect(html).toContain('ArrowDropDownIcon')
  expect(html).toContain('aria-haspopup="menu"')
  expect(html).toContain('aria-expanded="false"')
})
