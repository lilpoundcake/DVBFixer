export function openedAltlocArtifact(
  artifacts: readonly { id: string; file: string; kind: string; hidden?: boolean }[],
  openedFile: string | null, current: string, inputMode: 'auto' | 'manual',
): string | null {
  if (inputMode === 'manual' && artifacts.some(item => item.id === current && !item.hidden)) return null
  const opened = artifacts.find(item => item.file === openedFile && item.kind === 'structure' && !item.hidden && /\.(pdb|ent|cif|mmcif)$/i.test(item.file))
  return opened && opened.id !== current ? opened.id : null
}
