/** Read-only Chrome acceptance against the actual running scientific cohort.
 * Never submits, cancels, retries or changes a scientific job.
 */
import assert from 'node:assert/strict'
import { spawn } from 'node:child_process'
import fs from 'node:fs'
import path from 'node:path'

const root = path.resolve(process.argv[2] || '../tmp/real-complex-cc94-2026-11-06')
const service = JSON.parse(fs.readFileSync(path.join(root, 'managed-execution/service.json'), 'utf8'))
assert.equal(new URL(service.origin).hostname, '127.0.0.1')
const output = fs.mkdtempSync(path.join(root, 'chrome-live-'))
const chrome = process.env.CHROME_EXECUTABLE || '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome'
const child = spawn(chrome, ['--headless=new', '--no-first-run', '--no-default-browser-check', '--disable-gpu',
  '--remote-debugging-port=0', `--user-data-dir=${path.join(output, 'profile')}`, 'about:blank'],
{ stdio: ['ignore', 'pipe', 'pipe'] })
let stderr = '', stdout = '', socket
child.stderr.on('data', chunk => { stderr += String(chunk) })
child.stdout.on('data', chunk => { stdout += String(chunk) })
const delay = ms => new Promise(resolve => setTimeout(resolve, ms))
async function until(check, label, seconds = 90) {
  const deadline = Date.now() + seconds * 1000
  while (Date.now() < deadline) {
    const result = await check()
    if (result) return result
    await delay(250)
  }
  throw Error(label)
}
try {
  const port = await until(() => stderr.match(/DevTools listening on ws:\/\/[^:]+:(\d+)/)?.[1], 'Chrome startup failed', 20)
  const target = await (await fetch(`http://127.0.0.1:${port}/json/new?${encodeURIComponent(service.origin)}`,
    { method: 'PUT' })).json()
  socket = new WebSocket(target.webSocketDebuggerUrl)
  await new Promise((resolve, reject) => {
    socket.addEventListener('open', resolve, { once: true })
    socket.addEventListener('error', reject, { once: true })
  })
  let sequence = 0
  const pending = new Map()
  socket.addEventListener('message', event => {
    const message = JSON.parse(String(event.data)), request = pending.get(message.id)
    if (!request) return
    pending.delete(message.id)
    clearTimeout(request.timer)
    if (message.error) request.reject(Error(JSON.stringify(message.error)))
    else request.resolve(message.result)
  })
  const call = (method, params = {}) => new Promise((resolve, reject) => {
    const id = ++sequence
    const timer = setTimeout(() => { pending.delete(id); reject(Error(`CDP timeout: ${method}`)) }, 10000)
    pending.set(id, { resolve, reject, timer })
    socket.send(JSON.stringify({ id, method, params }))
  })
  const evaluate = async expression => {
    const result = await call('Runtime.evaluate', { expression, returnByValue: true, awaitPromise: true })
    if (result.exceptionDetails) throw Error(JSON.stringify(result.exceptionDetails))
    return result.result.value
  }
  const openMonitor = () => until(() => evaluate(`(()=>{const tab=[...document.querySelectorAll('.flexlayout__tab_button')]
    .find(el=>el.textContent.trim()==='Job Monitor');if(!tab)return false;tab.click();return true})()`), 'Job Monitor unavailable')
  await call('Emulation.setDeviceMetricsOverride', { width: 1440, height: 1000, deviceScaleFactor: 1, mobile: false })
  await openMonitor()
  const active = await until(async () => {
    const response = await fetch(`${service.origin}/api/v1/jobs`)
    const page = await response.json()
    return page.activeJobs.find(job => job.progress?.seed !== undefined && ['sampling', 'refinement'].includes(job.progress.stage))
  }, 'No measured sampling/refinement stage available from an actual checkpoint-backed job')
  await until(() => evaluate(`document.body.innerText.includes(${JSON.stringify(active.id)}) &&
    document.body.innerText.includes('Seed ${active.progress.seed}')`), 'Live job/seed missing from monitor')
  fs.writeFileSync(path.join(output, 'active.png'), Buffer.from((await call('Page.captureScreenshot', { format: 'png' })).data, 'base64'))
  await call('Page.reload')
  await openMonitor()
  await until(() => evaluate(`document.body.innerText.includes(${JSON.stringify(active.id)})`), 'Real scientific job missing after reload')
  const response = await fetch(`${service.origin}/api/v1/workspaces/${active.workspaceId}/jobs/${active.id}`)
  const restored = await response.json()
  assert.ok(['running', 'succeeded', 'failed'].includes(restored.status))
  assert.notEqual(restored.status, 'cancelled')
  fs.writeFileSync(path.join(output, 'reloaded.png'), Buffer.from((await call('Page.captureScreenshot', { format: 'png' })).data, 'base64'))
  fs.writeFileSync(path.join(output, 'summary.json'), JSON.stringify({ passed: true,
    scope: 'real cc94 checkpoint-backed scientific job; read-only browser monitoring and reload, not final geometry acceptance',
    checks: ['actual measured stage and seed visible', 'job restored after browser reload without cancellation'],
    active, restored }, null, 2))
  console.log(JSON.stringify({ passed: true, output, job: active.id, progress: active.progress }, null, 2))
} catch (error) {
  fs.writeFileSync(path.join(output, 'summary.json'), JSON.stringify({ passed: false, error: String(error) }, null, 2))
  throw error
} finally {
  socket?.close()
  const exited = new Promise(resolve => child.once('exit', resolve))
  child.kill('SIGTERM')
  let timer
  await Promise.race([exited, new Promise(resolve => {
    timer = setTimeout(() => { child.kill('SIGKILL'); resolve() }, 5000)
  })]).finally(() => clearTimeout(timer))
  fs.writeFileSync(path.join(output, 'chrome.stderr.log'), stderr)
  fs.writeFileSync(path.join(output, 'chrome.stdout.log'), stdout)
}
