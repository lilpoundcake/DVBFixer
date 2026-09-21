import type { DvbfixerProcessAdmissionSnapshot } from './dvbfixer-runner'
import type { HttpObservation, HttpObservationStart } from './request-observability'

export interface MetricsConfig {
  enabled: boolean
}

const BUCKETS = [0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1, 2.5, 5, 10, 30, 60, 300, 1800]

export function parseMetricsConfig(environment: NodeJS.ProcessEnv = process.env): MetricsConfig {
  const value = environment.DVBFIXER_METRICS ?? 'off'
  if (value !== 'off' && value !== 'prometheus') throw new Error('DVBFIXER_METRICS must be off or prometheus')
  return { enabled: value === 'prometheus' }
}

function labels(values: Record<string, string>): string {
  return `{${Object.entries(values).sort(([a], [b]) => a.localeCompare(b))
    .map(([key, value]) => `${key}="${value.replace(/\\/g, '\\\\').replace(/"/g, '\\"').replace(/\n/g, '\\n')}"`)
    .join(',')}}`
}

export class ApiMetrics {
  private readonly requests = new Map<string, { labels: Record<string, string>; value: number }>()
  private readonly durations = new Map<string, { route: string; count: number; sum: number; buckets: number[] }>()
  private readonly active = new Map<string, { labels: Record<string, string>; value: number }>()

  start(observation: HttpObservationStart): void {
    if (observation.route === '/api/metrics') return
    const key = `${observation.method}\0${observation.route}`
    const entry = this.active.get(key) || { labels: { method: observation.method, route: observation.route }, value: 0 }
    entry.value += 1
    this.active.set(key, entry)
  }

  finish(observation: HttpObservation): void {
    if (observation.route === '/api/metrics') return
    const activeKey = `${observation.method}\0${observation.route}`
    const active = this.active.get(activeKey)
    if (active) active.value = Math.max(0, active.value - 1)
    const requestLabels = {
      method: observation.method, outcome: observation.outcome,
      route: observation.route, status: observation.status === null ? 'none' : String(observation.status),
    }
    const requestKey = Object.values(requestLabels).join('\0')
    const request = this.requests.get(requestKey) || { labels: requestLabels, value: 0 }
    request.value += 1
    this.requests.set(requestKey, request)
    const duration = this.durations.get(observation.route) || {
      route: observation.route, count: 0, sum: 0, buckets: BUCKETS.map(() => 0),
    }
    duration.count += 1
    duration.sum += observation.durationSeconds
    BUCKETS.forEach((bucket, index) => {
      if (observation.durationSeconds <= bucket) duration.buckets[index] += 1
    })
    this.durations.set(observation.route, duration)
  }

  render(admission: DvbfixerProcessAdmissionSnapshot): string {
    const lines = [
      '# HELP dvbfixer_http_requests_total Completed API requests.',
      '# TYPE dvbfixer_http_requests_total counter',
    ]
    for (const entry of [...this.requests.values()].sort((a, b) => labels(a.labels).localeCompare(labels(b.labels)))) {
      lines.push(`dvbfixer_http_requests_total${labels(entry.labels)} ${entry.value}`)
    }
    lines.push('# HELP dvbfixer_http_request_duration_seconds API request duration.', '# TYPE dvbfixer_http_request_duration_seconds histogram')
    for (const entry of [...this.durations.values()].sort((a, b) => a.route.localeCompare(b.route))) {
      BUCKETS.forEach((bucket, index) => lines.push(
        `dvbfixer_http_request_duration_seconds_bucket${labels({ route: entry.route, le: String(bucket) })} ${entry.buckets[index]}`,
      ))
      lines.push(`dvbfixer_http_request_duration_seconds_bucket${labels({ route: entry.route, le: '+Inf' })} ${entry.count}`)
      lines.push(`dvbfixer_http_request_duration_seconds_sum${labels({ route: entry.route })} ${entry.sum}`)
      lines.push(`dvbfixer_http_request_duration_seconds_count${labels({ route: entry.route })} ${entry.count}`)
    }
    lines.push('# HELP dvbfixer_http_active_requests In-flight API requests.', '# TYPE dvbfixer_http_active_requests gauge')
    for (const entry of [...this.active.values()].sort((a, b) => labels(a.labels).localeCompare(labels(b.labels)))) {
      lines.push(`dvbfixer_http_active_requests${labels(entry.labels)} ${entry.value}`)
    }
    lines.push(
      '# HELP dvbfixer_child_processes_active Active DVBFixer child-process permits.',
      '# TYPE dvbfixer_child_processes_active gauge',
      `dvbfixer_child_processes_active ${admission.active}`,
      '# HELP dvbfixer_child_processes_queued Queued DVBFixer child-process requests.',
      '# TYPE dvbfixer_child_processes_queued gauge',
      `dvbfixer_child_processes_queued ${admission.queued}`,
      '# HELP dvbfixer_child_process_limit Configured child-process concurrency limit.',
      '# TYPE dvbfixer_child_process_limit gauge',
      `dvbfixer_child_process_limit ${admission.limit}`,
      '# HELP dvbfixer_child_process_admission_accepting Whether child-process admission is open.',
      '# TYPE dvbfixer_child_process_admission_accepting gauge',
      `dvbfixer_child_process_admission_accepting ${admission.accepting ? 1 : 0}`,
    )
    return `${lines.join('\n')}\n`
  }
}
