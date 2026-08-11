import type { StructureMeta } from '../stores/structureStore'
import type { WorkspaceArtifact } from '../stores/workspaceStore'

/** Build a complete metadata snapshot so loading a file cannot leak fields from the prior structure. */
export function structureMetaFromArtifact(artifact: WorkspaceArtifact): StructureMeta {
  return {
    name: artifact.name,
    organism: artifact.organism || '',
    method: artifact.method || '',
    resolution: artifact.resolution === undefined || artifact.resolution === null ? '' : String(artifact.resolution),
    description: artifact.description || '',
    iggSubtype: artifact.iggSubtype || '',
    allotype: artifact.allotype || '',
    equivalentChains: artifact.equivalentChains ?? undefined,
  }
}

export function workspaceArtifactNameMap(artifacts: readonly WorkspaceArtifact[]): Map<string, string> {
  return new Map(artifacts.map(artifact => [artifact.file, artifact.name]))
}
