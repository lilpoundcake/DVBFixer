import { Type, type Static } from '@sinclair/typebox'
import { Value } from '@sinclair/typebox/value'

const fields = {
  version: Type.Literal(1),
  stage: Type.Union((['preparation', 'runner', 'sampling', 'refinement', 'validation', 'publication'] as const).map(stage => Type.Literal(stage))),
  seed: Type.Optional(Type.Integer({ minimum: 0, maximum: Number.MAX_SAFE_INTEGER })),
}
const count = () => Type.Integer({ minimum: 0, maximum: Number.MAX_SAFE_INTEGER })
export const JobProgressSchema = Type.Union([
  Type.Object(fields, { additionalProperties: false }),
  Type.Object({ ...fields, completed: count(), total: Type.Integer({ minimum: 1, maximum: Number.MAX_SAFE_INTEGER }) }, { additionalProperties: false }),
])
export type JobProgress = Static<typeof JobProgressSchema>

/** Bounded NDJSON decoder on a dedicated pipe, never free-text log matching. */
export function createJobProgressDecoder(onProgress: (event: JobProgress) => void): (chunk: Buffer) => void {
  let pending = ''
  let discarding = false
  let delivered = 0
  return chunk => {
    if (delivered >= 4096) return
    // ASCII-only protocol. A malformed/oversized frame is ignored through newline.
    for (const byte of chunk) {
      if (byte === 10) {
        if (!discarding && pending) {
          try {
            const event: unknown = JSON.parse(pending)
            if (Value.Check(JobProgressSchema, event) &&
              (!('completed' in event) || event.completed <= event.total) && delivered < 4096) {
              delivered += 1
              onProgress(event)
            }
          } catch { /* Invalid telemetry cannot fail scientific execution. */ }
        }
        pending = ''
        discarding = false
      } else if (!discarding) {
        if (byte > 127 || pending.length >= 1024) { pending = ''; discarding = true }
        else pending += String.fromCharCode(byte)
      }
    }
  }
}
