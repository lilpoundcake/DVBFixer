/** Match the mirror to the native editing viewport, excluding its scrollbars. */
export function syncSequenceEditorLayout(
  editor: Pick<HTMLTextAreaElement, 'clientWidth' | 'clientHeight' | 'scrollTop' | 'scrollLeft'>,
  mirror: Pick<HTMLDivElement, 'style' | 'scrollTop' | 'scrollLeft'>,
): void {
  mirror.style.width = `${editor.clientWidth}px`
  mirror.style.height = `${editor.clientHeight}px`
  mirror.scrollTop = editor.scrollTop
  mirror.scrollLeft = editor.scrollLeft
}
