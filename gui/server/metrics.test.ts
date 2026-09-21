import { describe, expect, it } from 'vitest'
import { ApiMetrics, parseMetricsConfig } from './metrics'

describe('API metrics', () => {
  it('validates opt-in configuration', () => {
    expect(parseMetricsConfig({})).toEqual({ enabled: false })
    expect(parseMetricsConfig({ DVBFIXER_METRICS: 'prometheus' })).toEqual({ enabled: true })
    expect(() => parseMetricsConfig({ DVBFIXER_METRICS: 'json' })).toThrow(/off or prometheus/)
  })

  it('renders bounded request, duration, active, and process-admission metrics', () => {
    const metrics = new ApiMetrics()
    metrics.start({ method: 'GET', route: '/api/health' })
    metrics.finish({
      method: 'GET', route: '/api/health', status: 200, durationSeconds: 0.02, outcome: 'completed',
    })
    const output = metrics.render({ active: 1, queued: 2, limit: 3, queueLimit: 4, accepting: true })
    expect(output).toContain('dvbfixer_http_requests_total{method="GET",outcome="completed",route="/api/health",status="200"} 1')
    expect(output).toContain('dvbfixer_http_request_duration_seconds_bucket{le="0.025",route="/api/health"} 1')
    expect(output).toContain('dvbfixer_http_active_requests{method="GET",route="/api/health"} 0')
    expect(output).toContain('dvbfixer_child_processes_active 1')
    expect(output).toContain('dvbfixer_child_processes_queued 2')
    expect(output).toContain('dvbfixer_child_process_queue_limit 4')
    expect(output.endsWith('\n')).toBe(true)
  })

  it('excludes metrics scrapes from HTTP accounting', () => {
    const metrics = new ApiMetrics()
    metrics.start({ method: 'GET', route: '/api/metrics' })
    metrics.finish({
      method: 'GET', route: '/api/metrics', status: 200, durationSeconds: 1, outcome: 'completed',
    })
    expect(metrics.render({ active: 0, queued: 0, limit: 1, queueLimit: 16, accepting: true }))
      .not.toContain('route="/api/metrics"')
  })
})
