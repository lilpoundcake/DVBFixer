import { Type } from '@sinclair/typebox'

const text = () => Type.String({ maxLength: 2 * 1024 * 1024 })
const identity = () => Type.String({ minLength: 1, maxLength: 256 })
const span = Type.Object({ start: Type.Integer({ minimum: 0 }), end: Type.Integer({ minimum: 1 }) }, { additionalProperties: false })

/** Transport bounds only. Python still owns alignment/mask/scientific admission. */
export const HomologyJobRequestSchema = Type.Object({
  project: Type.Object({
    version: Type.Literal(1),
    id: Type.String({ minLength: 1, maxLength: 256, pattern: '^[a-zA-Z0-9_-]+$' }),
    name: Type.String({ maxLength: 256 }), targetFasta: text(),
    templates: Type.Array(Type.Object({
      id: identity(), file: Type.String({ minLength: 1, maxLength: 4096 }), chain: Type.String({ maxLength: 128 }), targetChain: identity(),
      label: Type.Optional(text()), fittedFile: Type.Optional(text()),
    }, { additionalProperties: false }), { maxItems: 1000 }),
    engine: Type.Union([Type.Literal('mafft'), Type.Literal('muscle'), Type.Literal('clustalo')]),
    alignmentGroups: Type.Array(Type.Object({
      chainId: identity(),
      rows: Type.Array(Type.Object({ id: identity(), kind: Type.Union([Type.Literal('target'), Type.Literal('template')]),
        sequence: text(), templateId: Type.Optional(identity()),
      }, { additionalProperties: false }), { maxItems: 1000 }),
      masks: Type.Record(Type.String(), Type.Array(span, { maxItems: 100_000 })),
      maskModes: Type.Optional(Type.Record(Type.String(), Type.Union([Type.Literal('all'), Type.Literal('ranges'), Type.Literal('none')]))),
    }, { additionalProperties: false }), { maxItems: 1000 }),
    structuralAlignment: Type.Optional(Type.Union([text(), Type.Record(Type.String(), text())])),
    modelOptions: Type.Record(Type.String(), Type.Unknown()),
    createdAt: identity(), updatedAt: identity(),
  }, { additionalProperties: false }),
}, { additionalProperties: false, $id: 'HomologyJobRequest' })
