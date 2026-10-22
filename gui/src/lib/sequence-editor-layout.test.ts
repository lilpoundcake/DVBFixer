import { describe, expect, it } from 'vitest'
import { syncSequenceEditorLayout } from './sequence-editor-layout'

describe('FASTA editing viewport', () => {
  it('uses the native client width rather than the box width including a scrollbar', () => {
    const editor = { clientWidth: 284, clientHeight: 51, scrollTop: 84, scrollLeft: 0 }
    const mirror = { style: {} as CSSStyleDeclaration, scrollTop: 0, scrollLeft: 0 }
    syncSequenceEditorLayout(editor, mirror)
    expect(mirror.style.width).toBe('284px')
    expect(mirror.style.height).toBe('51px')
    expect(mirror.scrollTop).toBe(84)
  })

  it('resets the mirror after deleting a prefix clamps scrolling and removes the scrollbar', () => {
    const editor = { clientWidth: 284, clientHeight: 51, scrollTop: 168, scrollLeft: 12 }
    const mirror = { style: {} as CSSStyleDeclaration, scrollTop: 0, scrollLeft: 0 }
    syncSequenceEditorLayout(editor, mirror)
    Object.assign(editor, { clientWidth: 300, scrollTop: 0, scrollLeft: 0 })
    syncSequenceEditorLayout(editor, mirror)
    expect(mirror.style.width).toBe('300px')
    expect(mirror.scrollTop).toBe(0)
    expect(mirror.scrollLeft).toBe(0)
  })

  it('tracks wrapping changes when the dock is resized', () => {
    const editor = { clientWidth: 300, clientHeight: 51, scrollTop: 0, scrollLeft: 0 }
    const mirror = { style: {} as CSSStyleDeclaration, scrollTop: 0, scrollLeft: 0 }
    syncSequenceEditorLayout(editor, mirror)
    Object.assign(editor, { clientWidth: 180, clientHeight: 68, scrollTop: 34 })
    syncSequenceEditorLayout(editor, mirror)
    expect(mirror.style.width).toBe('180px')
    expect(mirror.style.height).toBe('68px')
    expect(mirror.scrollTop).toBe(34)
  })
})
