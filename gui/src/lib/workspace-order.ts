export interface OrderableWorkspaceArtifact {
  id: string
}

/**
 * Move selected visible artifacts by one row without disturbing artifacts
 * excluded by the current view (for example hidden or non-PDB files).
 */
export function reorderVisibleArtifacts<T extends OrderableWorkspaceArtifact>(
  artifacts: readonly T[],
  visibleIds: readonly string[],
  selectedIds: ReadonlySet<string>,
  direction: -1 | 1,
): T[] {
  const visibleSet = new Set(visibleIds)
  const visibleIndexes: number[] = []
  const visible: T[] = []

  artifacts.forEach((artifact, index) => {
    if (!visibleSet.has(artifact.id)) return
    visibleIndexes.push(index)
    visible.push(artifact)
  })

  if (direction < 0) {
    for (let index = 1; index < visible.length; index++) {
      if (selectedIds.has(visible[index].id) && !selectedIds.has(visible[index - 1].id)) {
        ;[visible[index - 1], visible[index]] = [visible[index], visible[index - 1]]
      }
    }
  } else {
    for (let index = visible.length - 2; index >= 0; index--) {
      if (selectedIds.has(visible[index].id) && !selectedIds.has(visible[index + 1].id)) {
        ;[visible[index], visible[index + 1]] = [visible[index + 1], visible[index]]
      }
    }
  }

  const reordered = [...artifacts]
  visibleIndexes.forEach((artifactIndex, visibleIndex) => {
    reordered[artifactIndex] = visible[visibleIndex]
  })
  return reordered
}
