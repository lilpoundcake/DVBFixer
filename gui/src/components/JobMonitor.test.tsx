import { renderToStaticMarkup } from 'react-dom/server'
import { expect, it } from 'vitest'
import { JobMonitor } from './JobMonitor'

it('is an independently dockable monitor and never promises fabricated progress or scientific acceptance', () => {
  const html = renderToStaticMarkup(<JobMonitor />)
  expect(html).toContain('Job Monitor')
  expect(html).toContain('Progress is indeterminate')
  expect(html).toContain('Completion is not scientific acceptance')
  expect(html).toContain('Loading tasks')
  expect(html).toContain('<header')
  expect(html).toMatch(/<section[^>]*aria-label="Tasks"/)
  expect(html.replace(/<style[^>]*>[\s\S]*?<\/style>/g, '')).not.toContain('100%')
})
