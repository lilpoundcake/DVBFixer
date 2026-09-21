import { renderToStaticMarkup } from 'react-dom/server'
import { describe, expect, it } from 'vitest'
import { AuthGate } from './AuthGate'

describe('AuthGate', () => {
  it('does not render protected children before authentication completes', () => {
    const markup = renderToStaticMarkup(
      <AuthGate showLogout={false}>
        <div>protected workspace</div>
      </AuthGate>,
    )

    expect(markup).not.toContain('protected workspace')
    expect(markup).toContain('DVBfixer workspace')
  })
})
