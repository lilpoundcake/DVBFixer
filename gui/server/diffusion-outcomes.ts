import { Type, type Static } from '@sinclair/typebox'

export const DiffusionOutcomesSchema = Type.Object({
  report: Type.String(),
  outcomes: Type.Array(Type.Object({
    seed: Type.Integer({ minimum: 0, maximum: Number.MAX_SAFE_INTEGER }),
    candidateId: Type.Union([Type.String(), Type.Null()]),
    status: Type.Union([Type.Literal('accepted'), Type.Literal('rejected'), Type.Literal('not_returned')]),
    findings: Type.Array(Type.String()),
    hardGateFailures: Type.Array(Type.String()),
    omittedFindings: Type.Integer({ minimum: 0 }),
  }, { additionalProperties: false }), { maxItems: 1000 }),
}, { additionalProperties: false })

export type DiffusionOutcomes = Static<typeof DiffusionOutcomesSchema>

/** Make retained links relative to the workspace after native staging is copied. */
export function relocateDiffusionOutcomes(outcomes: DiffusionOutcomes, relativeReport: string): DiffusionOutcomes {
  return { ...outcomes, report: relativeReport }
}
