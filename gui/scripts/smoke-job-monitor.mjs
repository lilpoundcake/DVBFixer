/** Optional native Chrome smoke. Uses synthetic CLI work, never model weights. */
import assert from 'node:assert/strict'
import { spawn } from 'node:child_process'
import fs from 'node:fs'
import path from 'node:path'

const gui = path.resolve(import.meta.dirname, '..')
const outputParent = path.resolve(gui, '../tmp')
fs.mkdirSync(outputParent, { recursive: true })
const output = fs.mkdtempSync(path.join(outputParent, 'job-monitor-browser-'))
const chrome = process.env.CHROME_EXECUTABLE || '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome'
assert.ok(fs.existsSync(chrome), 'Set CHROME_EXECUTABLE to a native Chrome executable')
const fixture = path.join(output, 'synthetic-cli.cjs')
fs.writeFileSync(fixture, `
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto');
const argv=process.argv.slice(2),command=argv[0];
if(command==='--version'){console.log('dvbfixer 0.0.0 (synthetic browser smoke)');process.exit(0)}
const emit=(stage)=>fs.writeSync(Number(process.env.DVBFIXER_JOB_PROGRESS_FD),JSON.stringify({version:1,stage,seed:19})+'\\n');
emit('sampling');console.log('Synthetic work; no checkpoint sampling.');
if(command==='doctor'){
  emit('refinement');setInterval(()=>{if(fs.existsSync('finish'))process.exit(0)},20);
}else if(command==='homology'){
  const plan=JSON.parse(fs.readFileSync(argv[argv.indexOf('--template-plan')+1],'utf8'));
  if(!plan.templates[0].path.startsWith(process.cwd()+path.sep+'inputs'+path.sep))throw Error('Homology templates were not frozen');
  const bytes=fs.readFileSync(plan.templates[0].path),out=argv[argv.indexOf('-o')+1];
  emit('refinement');setInterval(()=>{if(fs.existsSync('finish')){
    fs.writeFileSync(out+'.pdb',bytes);fs.writeFileSync(out+'.dat','{}');process.exit(0);
  }},20);
}else if(command==='model'){
  const seeds=argv.flatMap((arg,i)=>arg==='--diffusion-seed'?[Number(argv[i+1])]:[]);
  if(JSON.stringify(seeds)!=='[7,19,101]')throw Error('incorrect repeated seed transport');
  const out=argv[argv.indexOf('-o')+1];fs.mkdirSync(out,{recursive:true});
  function write(name,data){const file=path.join(out,name);fs.mkdirSync(path.dirname(file),{recursive:true});
    fs.writeFileSync(file,data);return {path:name,sha256:crypto.createHash('sha256').update(data).digest('hex')};}
  const refs={};for(const label of ['pdb','dat','provenance','sampler_trace'])refs[label]=write('accepted.'+label,
    label==='provenance'?JSON.stringify({request_sha256:'a'.repeat(64)}):label);
  const report={version:1,request_sha256:'a'.repeat(64),requested_seeds:seeds,outcomes:[
    {seed:7,candidate_id:'accepted',status:'accepted',validation:{passed:true,hard_gate_failures:[],warnings:[]},coordinate:null},
    {seed:19,candidate_id:'rejected-19',status:'rejected',validation:{passed:false,hard_gate_failures:['synthetic-angle'],
      warnings:['Final validation: d/ASP28A: synthetic angle = 79.1 degrees; allowed 95-125 degrees.']},coordinate:write('rejected/seed-19.pdb','SYNTHETIC REJECTED COORDINATES')},
    {seed:101,candidate_id:null,status:'not_returned',validation:null,coordinate:null}]};
  emit('validation');const reportRef=write('seed-outcomes.json',JSON.stringify(report));emit('publication');
  write('bundle.json',JSON.stringify({schema_version:2,status:'success',seed_outcomes:reportRef,
    candidates:[{candidate_id:'accepted',seed:7,validation_passed:true,...refs}]}));
  console.log('Per-seed validation report: '+path.join(out,'seed-outcomes.json'));
}else throw Error('unsupported synthetic command: '+command);
`)

const delay = ms => new Promise(resolve => setTimeout(resolve, ms))
async function until(check, message, timeoutMs = 20_000) {
  const deadline = Date.now() + timeoutMs
  while (Date.now() < deadline) {
    const result = await check()
    if (result) return result
    await delay(100)
  }
  throw new Error(message)
}

function launch(executable, args, env, label) {
  const child = spawn(executable, args, { cwd: gui, env, stdio: ['ignore', 'pipe', 'pipe'] })
  const streams = { stdout: '', stderr: '' }
  for (const stream of ['stdout', 'stderr']) child[stream].on('data', chunk => { streams[stream] += String(chunk) })
  const exited = new Promise(resolve => { child.once('exit', resolve); child.once('error', error => resolve(String(error))) })
  return { child, streams, async stop() {
    if (child.exitCode === null && child.signalCode === null) {
      child.kill('SIGTERM')
      let timer
      await Promise.race([exited, new Promise(resolve => {
        timer = setTimeout(() => { child.kill('SIGKILL'); resolve() }, 5000)
      })]).finally(() => clearTimeout(timer))
    }
    for (const stream of ['stdout', 'stderr']) fs.writeFileSync(path.join(output, `${label}.${stream}.log`), streams[stream])
  } }
}

function connect(url) {
  const socket = new WebSocket(url)
  let sequence = 0
  const pending = new Map()
  socket.addEventListener('message', event => {
    const message = JSON.parse(String(event.data))
    const request = pending.get(message.id)
    if (!request) return
    pending.delete(message.id); clearTimeout(request.timer)
    if (message.error) request.reject(new Error(JSON.stringify(message.error)))
    else request.resolve(message.result)
  })
  return { socket, ready: new Promise((resolve, reject) => {
    socket.addEventListener('open', resolve, { once: true })
    socket.addEventListener('error', reject, { once: true })
  }), call(method, params = {}) {
    return new Promise((resolve, reject) => {
      const id = ++sequence
      const timer = setTimeout(() => { pending.delete(id); reject(new Error(`Chrome timed out: ${method}`)) }, 10_000)
      pending.set(id, { resolve, reject, timer }); socket.send(JSON.stringify({ id, method, params }))
    })
  } }
}

const dataRoot = path.join(output, 'data')
const host = launch(process.execPath, ['dist-server/server.js'], { ...process.env,
  DVBFIXER_HOST: '127.0.0.1', DVBFIXER_PORT: '0', DVBFIXER_GUI_DATA_DIR: dataRoot,
  DVBFIXER_EXECUTABLE: process.execPath, DVBFIXER_ARGS: JSON.stringify([fixture]),
  DVBFIXER_MAX_CONCURRENT_PROCESSES: '2',
}, 'host')
let browser, cdp
const checks = []
try {
  const port = await until(() => host.streams.stdout.match(/listening on http:\/\/[^:]+:(\d+)/)?.[1], 'host did not start')
  const origin = `http://127.0.0.1:${port}`
  const api = async (route, method = 'GET', body) => {
    const response = await fetch(origin + route, { method, headers: { 'Content-Type': 'application/json' },
      ...(body === undefined ? {} : { body: JSON.stringify(body) }) })
    const data = await response.json()
    assert.ok(response.ok, `${method} ${route}: ${JSON.stringify(data)}`)
    return data
  }
  const firstWorkspace = (await api('/api/workspaces'))[0]
  const secondWorkspace = await api('/api/workspaces', 'POST', { name: 'Synthetic second workspace' })
  const jobs = []
  for (const workspace of [firstWorkspace, secondWorkspace]) {
    jobs.push(await api(`/api/v1/workspaces/${workspace.id}/jobs`, 'POST', { command: 'doctor' }))
  }
  browser = launch(chrome, ['--headless=new', '--no-first-run', '--no-default-browser-check',
    '--remote-debugging-port=0', `--user-data-dir=${path.join(output, 'chrome-profile')}`, 'about:blank'], process.env, 'chrome')
  const debugPort = await until(() => browser.streams.stderr.match(/DevTools listening on ws:\/\/[^:]+:(\d+)/)?.[1], 'Chrome did not start')
  const target = await (await fetch(`http://127.0.0.1:${debugPort}/json/new?${encodeURIComponent(origin)}`, { method: 'PUT' })).json()
  cdp = connect(target.webSocketDebuggerUrl)
  await cdp.ready
  const evaluate = async expression => {
    const result = await cdp.call('Runtime.evaluate', { expression, returnByValue: true, awaitPromise: true })
    if (result.exceptionDetails) throw new Error(JSON.stringify(result.exceptionDetails))
    return result.result.value
  }
  const text = () => evaluate('document.body.innerText')
  const openMonitor = () => until(() => evaluate(`(()=>{const tab=[...document.querySelectorAll('.flexlayout__tab_button')]
    .find(el=>el.textContent.trim()==='Job Monitor');if(!tab)return false;tab.click();return true})()`), 'Job Monitor tab unavailable')
  await openMonitor()
  await until(async () => {
    const body = await text()
    return jobs.every(job => body.includes(job.id)) && body.includes('Refinement · Seed 19')
  }, 'simultaneous jobs or authoritative stage missing from monitor')
  checks.push('two real child processes visible across workspaces with authoritative stages')
  await evaluate(`(()=>{const button=[...document.querySelectorAll('button')].find(el=>el.textContent==='Stop'&&
    el.parentElement.parentElement.textContent.includes(${JSON.stringify(jobs[0].id)}));if(!button)throw Error('Stop not found');button.click()})()`)
  await until(async () => (await api(`/api/v1/workspaces/${jobs[0].workspaceId}/jobs/${jobs[0].id}`)).status === 'cancelled', 'first job did not stop')
  assert.equal((await api(`/api/v1/workspaces/${jobs[1].workspaceId}/jobs/${jobs[1].id}`)).status, 'running')
  checks.push('browser Stop cancels one child without stopping the other')
  fs.writeFileSync(path.join(dataRoot, 'projects', jobs[1].workspaceId, jobs[1].outputDir, 'finish'), 'done')
  await until(async () => (await api(`/api/v1/workspaces/${jobs[1].workspaceId}/jobs/${jobs[1].id}`)).status === 'succeeded', 'second job did not complete')
  fs.writeFileSync(path.join(dataRoot, 'projects', firstWorkspace.id, 'input.pdb'), 'END\n')
  const mixed = await api(`/api/v1/workspaces/${firstWorkspace.id}/jobs`, 'POST', {
    command: 'model', inputFile: 'input.pdb', values: { '--backend': 'diffusion', '--diffusion-seed': '7, 19 101' },
  })
  const finished = await until(async () => {
    const job = await api(`/api/v1/workspaces/${mixed.workspaceId}/jobs/${mixed.id}`)
    if (job.status === 'failed') throw new Error(job.error)
    return job.status === 'succeeded' ? job : null
  }, 'mixed synthetic outcome did not finish')
  assert.deepEqual(finished.diffusion.outcomes.map(row => row.status), ['accepted', 'rejected', 'not_returned'])
  assert.ok(finished.outputFile.endsWith('/accepted.pdb'))
  const manifest = await api(`/api/workspaces/${mixed.workspaceId}`)
  const rejected = manifest.artifacts.find(artifact => artifact.file.endsWith('/rejected/seed-19.pdb'))
  assert.equal(rejected.kind, 'artifact')
  assert.ok(rejected.name.startsWith('[rejected]'))
  const checkOutcomes = async () => {
    const body = await text()
    return body.includes('Seed 7: Accepted') && body.includes('Seed 19: Rejected') &&
      body.includes('Seed 101: Not returned for final validation') && body.includes('d/ASP28A')
  }
  await until(checkOutcomes, 'mixed seed outcomes not displayed')
  checks.push('all synthetic seed outcomes displayed without promoting rejected structures')
  await cdp.call('Emulation.setDeviceMetricsOverride', { width: 900, height: 900, deviceScaleFactor: 1, mobile: false })
  await until(checkOutcomes, 'outcomes unavailable after viewport resize')
  await cdp.call('Page.reload')
  await openMonitor()
  await until(checkOutcomes, 'persisted outcomes unavailable after page reload')
  checks.push('monitor restores persisted outcomes after viewport resize and page reload')
  await evaluate(`(()=>{const button=document.querySelector('button[aria-label="Open per-seed validation report"]');
    if(!button)throw Error('report link unavailable');button.click()})()`)
  await until(async () => (await text()).includes('"requested_seeds"'), 'retained report link did not open report text')
  checks.push('retained report opens in Text Files after native staging relocation')
  const screenshot = await cdp.call('Page.captureScreenshot', { format: 'png' })
  fs.writeFileSync(path.join(output, 'monitor-report.png'), Buffer.from(screenshot.data, 'base64'))

  // Real GUI submission/restoration, still with the explicitly synthetic CLI.
  fs.writeFileSync(path.join(dataRoot, 'projects', secondWorkspace.id, 'template.pdb'), 'END\n')
  const createdProject = await api('/api/homology/projects', 'POST', { workspaceId: secondWorkspace.id })
  const project = { ...createdProject, name: 'Persistent Homology smoke', targetFasta: '>H\nMA\n',
    templates: [{ id: 'template-1', file: 'template.pdb', chain: 'H', targetChain: 'H' }],
    alignmentGroups: [{ chainId: 'H', rows: [{ id: 'H', kind: 'target', sequence: 'MA' }], masks: {} }],
    modelOptions: { '--num-models': 1, '--md-level': 'none' },
  }
  await api(`/api/homology/projects/${project.id}`, 'PUT', { ...project, workspaceId: secondWorkspace.id })
  const openHomology = async () => {
    await until(() => evaluate(`(()=>{const workspace=[...document.querySelectorAll('[role="button"].MuiListItemButton-root')]
      .find(el=>el.textContent.includes('Synthetic second workspace'));if(!workspace)return false;workspace.click();return true})()`), 'workspace selector unavailable')
    await until(() => evaluate(`(()=>{const tab=[...document.querySelectorAll('.flexlayout__tab_button')]
      .find(el=>el.textContent.trim()==='Homology');if(!tab)return false;tab.click();return true})()`), 'Homology tab unavailable')
    await until(async () => (await text()).includes('4. Model'), 'Homology project did not load')
  }
  const submitHomology = async previousId => {
    await evaluate(`(()=>{const tab=[...document.querySelectorAll('button[role="tab"]')].find(el=>el.textContent==='4. Model');tab.click()})()`)
    await until(() => evaluate(`(()=>{const button=[...document.querySelectorAll('button')].find(el=>el.textContent==='Run Modeller');
      if(!button||button.disabled)return false;button.click();return true})()`), 'Run Modeller unavailable')
    return until(async () => {
      const all = await api(`/api/v1/workspaces/${secondWorkspace.id}/jobs`)
      return all.find(job => job.command === 'homology' && job.workflow?.projectId === project.id && job.id !== previousId && job.status === 'running')
    }, 'GUI did not create a persistent Homology job')
  }
  await openHomology()
  const homology = await submitHomology(null)
  await cdp.call('Page.navigate', { url: 'about:blank' })
  assert.equal((await api(`/api/v1/workspaces/${secondWorkspace.id}/jobs/${homology.id}`)).status, 'running')
  await cdp.call('Page.navigate', { url: origin })
  await openHomology()
  await until(async () => (await text()).includes('Continues after closing or reloading this tab.'), 'Homology did not restore its active job')
  const activeScreenshot = await cdp.call('Page.captureScreenshot', { format: 'png' })
  fs.writeFileSync(path.join(output, 'homology-restored.png'), Buffer.from(activeScreenshot.data, 'base64'))
  checks.push('GUI-submitted Homology survives navigation away and restores the active project job')
  await evaluate(`(()=>{const button=[...document.querySelectorAll('button')].find(el=>el.textContent==='Stop Homology');
    if(!button||button.disabled)throw Error('Stop Homology unavailable');button.click()})()`)
  await until(async () => (await api(`/api/v1/workspaces/${secondWorkspace.id}/jobs/${homology.id}`)).status === 'cancelled', 'Homology Stop did not cancel the child')
  await until(async () => (await text()).includes('Homology cancelled; diagnostics retained.'), 'Homology cancellation not confirmed in the panel')
  checks.push('restored Homology can be stopped explicitly without relying on browser disconnect')
  const completedHomology = await submitHomology(homology.id)
  await cdp.call('Page.reload')
  await openHomology()
  await until(async () => (await text()).includes('Continues after closing or reloading this tab.'), 'second Homology did not restore after reload')
  fs.writeFileSync(path.join(dataRoot, 'projects', secondWorkspace.id, completedHomology.outputDir, 'finish'), 'done')
  await until(async () => (await api(`/api/v1/workspaces/${secondWorkspace.id}/jobs/${completedHomology.id}`)).status === 'succeeded', 'Homology completion failed')
  await until(async () => (await text()).includes('Model complete:'), 'Homology panel did not show completed output')
  checks.push('Homology reload/reconnect completes with frozen-template output and retained diagnostics')

  // Exercise the approved group rename with historical persisted panel state.
  await cdp.call('Page.navigate', { url: 'about:blank' })
  const savedWorkspace = await api(`/api/workspaces/${secondWorkspace.id}`)
  const oldPanel = savedWorkspace.toolState.dvbfixer || {}
  await api(`/api/workspaces/${secondWorkspace.id}`, 'PATCH', { revision: savedWorkspace.revision,
    toolState: { ...savedWorkspace.toolState, dvbfixer: { ...oldPanel, category: 'Glycoprotein preparation',
      values: { ...(oldPanel.values || {}), convert: { '--navigation-smoke-marker': 'keep' } },
    } },
  })
  await cdp.call('Page.navigate', { url: origin })
  await until(() => evaluate(`(()=>{const tab=[...document.querySelectorAll('.flexlayout__tab_button')]
    .find(el=>el.textContent.trim()==='DVBFixer');if(!tab)return false;tab.click();return true})()`), 'DVBFixer tab unavailable')
  await until(() => evaluate(`Boolean(document.querySelector('button[aria-label="Workflow group: Glycans"]'))`), 'saved glycan group was not migrated')
  const key = async (name, code) => {
    await cdp.call('Input.dispatchKeyEvent', { type: 'keyDown', key: name, code: name, windowsVirtualKeyCode: code,
      ...(name === 'Enter' ? { text: '\r' } : {}) })
    await cdp.call('Input.dispatchKeyEvent', { type: 'keyUp', key: name, code: name, windowsVirtualKeyCode: code })
  }
  const activeTool = name => evaluate(`[...document.querySelectorAll('button[role="tab"][aria-selected="true"]')]
    .some(el=>el.textContent.trim()===${JSON.stringify(name)})`)
  await evaluate(`document.querySelector('button[aria-label="Workflow group: Glycans"]').focus()`)
  await key('Enter', 13)
  await until(() => evaluate(`document.activeElement?.textContent.trim()==='Glycans'`), 'group menu did not focus the selected group')
  const groupLabels = await evaluate(`[...document.querySelectorAll('#workflow-groups [role="menuitem"]')].map(el=>el.textContent.trim())`)
  assert.ok(!groupLabels.includes('Analysis') && !groupLabels.includes('Glycoprotein preparation'))
  await key('ArrowRight', 39)
  await until(() => evaluate(`Boolean(document.querySelector('[aria-label="Glycans tools"]'))`), 'keyboard submenu unavailable')
  assert.deepEqual(await evaluate(`[...document.querySelectorAll('[aria-label="Glycans tools"] [role="menuitem"]')].map(el=>el.textContent.trim())`), ['Convert', 'Transplant', 'Cluster'])
  await key('ArrowDown', 40)
  await key('Enter', 13)
  await until(() => activeTool('Transplant'), 'keyboard did not select Transplant')
  checks.push('historical Glycoprotein preparation state restores to Glycans with keyboard cascading-menu selection')
  const center = async (selector, text) => until(() => evaluate(`(()=>{const element=[...document.querySelectorAll(${JSON.stringify(selector)})]
    .find(el=>${text === undefined ? 'true' : `el.textContent.trim()===${JSON.stringify(text)}`});if(!element)return null;
    const rect=element.getBoundingClientRect();if(!rect.width||!rect.height)return null;return {x:rect.x+rect.width/2,y:rect.y+rect.height/2}})()`), 'menu target unavailable')
  const click = async point => {
    await cdp.call('Input.dispatchMouseEvent', { type: 'mousePressed', button: 'left', clickCount: 1, ...point })
    await cdp.call('Input.dispatchMouseEvent', { type: 'mouseReleased', button: 'left', clickCount: 1, ...point })
  }
  await click(await center('button[aria-label="Workflow group: Glycans"]'))
  await cdp.call('Input.dispatchMouseEvent', { type: 'mouseMoved', ...await center('#workflow-groups [role="menuitem"]', 'Glycans') })
  const cluster = await center('[aria-label="Glycans tools"] [role="menuitem"]', 'Cluster')
  await cdp.call('Input.dispatchMouseEvent', { type: 'mouseMoved', ...cluster })
  await click(cluster)
  await until(() => activeTool('Cluster'), 'pointer transition did not reach Cluster')
  checks.push('pointer hover and transition into the submenu opens Cluster directly')
  await cdp.call('Emulation.setDeviceMetricsOverride', { width: 560, height: 900, deviceScaleFactor: 1, mobile: true })
  await cdp.call('Emulation.setTouchEmulationEnabled', { enabled: true })
  const tap = async point => {
    await cdp.call('Input.dispatchTouchEvent', { type: 'touchStart', touchPoints: [{ ...point, id: 1 }] })
    await cdp.call('Input.dispatchTouchEvent', { type: 'touchEnd', touchPoints: [] })
  }
  await tap(await center('button[aria-label="Workflow group: Glycans"]'))
  await tap(await center('#workflow-groups [role="menuitem"]', 'Glycans'))
  await until(() => evaluate(`(()=>{const menu=document.querySelector('[aria-label="Glycans tools"]');if(!menu)return false;
    const r=menu.getBoundingClientRect();return r.left>=0&&r.top>=0&&r.right<=innerWidth&&r.bottom<=innerHeight})()`), 'narrow touch submenu exceeds the viewport')
  await tap(await center('[aria-label="Glycans tools"] [role="menuitem"]', 'Convert'))
  await until(() => activeTool('Convert'), 'touch selection did not open Convert')
  const panelAfter = await until(async () => {
    const state = (await api(`/api/workspaces/${secondWorkspace.id}`)).toolState.dvbfixer
    return state.category === 'Glycans' ? state : null
  }, 'renamed category was not persisted')
  assert.equal(panelAfter.values.convert['--navigation-smoke-marker'], 'keep')
  checks.push('touch selection works in a narrow viewport and preserves per-tool saved values')
  const navigationScreenshot = await cdp.call('Page.captureScreenshot', { format: 'png' })
  fs.writeFileSync(path.join(output, 'glycans-navigation.png'), Buffer.from(navigationScreenshot.data, 'base64'))
  fs.writeFileSync(path.join(output, 'summary.json'), JSON.stringify({ passed: true, scope: 'real Chrome and child-process transport with synthetic CLI outputs; no scientific sampling', checks }, null, 2))
  console.log(JSON.stringify({ passed: true, output, checks }, null, 2))
} catch (error) {
  fs.writeFileSync(path.join(output, 'summary.json'), JSON.stringify({ passed: false, checks, error: String(error) }, null, 2))
  console.error(`Browser smoke failed; retained evidence: ${output}`)
  throw error
} finally {
  cdp?.socket.close()
  await browser?.stop()
  await host.stop()
}
