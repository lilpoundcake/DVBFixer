import { Type, type Static } from '@sinclair/typebox'

const InputValueSchema = Type.Union([
  Type.String({ minLength: 1 }),
  Type.Array(Type.String({ minLength: 1 }), { minItems: 1 }),
])

export const ManagedJobRequestSchema = Type.Object({
  workspaceId: Type.Optional(Type.String({ minLength: 1 })),
  command: Type.String({ minLength: 1 }),
  inputFile: Type.Optional(Type.String({ minLength: 1 })),
  inputs: Type.Optional(Type.Record(Type.String(), InputValueSchema)),
  values: Type.Optional(Type.Record(Type.String(), Type.Unknown())),
  fastaContent: Type.Optional(Type.String({ maxLength: 10 * 1024 * 1024 })),
}, { additionalProperties: false, $id: 'ManagedJobRequest' })

export const ManagedJobRecordSchema = Type.Object({
  version: Type.Literal(1),
  id: Type.String({ format: 'uuid' }),
  workspaceId: Type.String(),
  status: Type.Union([
    Type.Literal('queued'), Type.Literal('running'), Type.Literal('succeeded'),
    Type.Literal('failed'), Type.Literal('cancelled'),
  ]),
  createdAt: Type.String(),
  startedAt: Type.Union([Type.String(), Type.Null()]),
  finishedAt: Type.Union([Type.String(), Type.Null()]),
  command: Type.String(),
  args: Type.Array(Type.String()),
  exitCode: Type.Union([Type.Integer(), Type.Null()]),
  stdoutLog: Type.String(),
  stderrLog: Type.String(),
  outputDir: Type.String(),
  outputFile: Type.Union([Type.String(), Type.Null()]),
  requestedByPrincipalId: Type.String(),
  error: Type.Optional(Type.String()),
  provenance: Type.Optional(Type.Object({
    apiVersion: Type.Literal(1),
    command: Type.String(),
    serviceVersion: Type.String(),
    dvbfixerVersion: Type.String(),
    request: Type.Object({
      inputs: Type.Record(Type.String(), InputValueSchema),
      values: Type.Record(Type.String(), Type.Unknown()),
      fastaSha256: Type.Union([Type.String({ pattern: '^[a-f0-9]{64}$' }), Type.Null()]),
    }, { additionalProperties: false }),
    resolvedInputs: Type.Array(Type.Object({
      artifactId: Type.Union([Type.String(), Type.Null()]),
      file: Type.String(), sha256: Type.String({ pattern: '^[a-f0-9]{64}$' }),
    }, { additionalProperties: false })),
    outputs: Type.Array(Type.Object({
      artifactId: Type.String(), file: Type.String(),
      sha256: Type.String({ pattern: '^[a-f0-9]{64}$' }),
    }, { additionalProperties: false })),
  }, { additionalProperties: false })),
}, { additionalProperties: false, $id: 'ManagedJobRecord' })

export type ManagedJobRequestPayload = Static<typeof ManagedJobRequestSchema>
