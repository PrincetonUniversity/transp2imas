'use strict';
/* The Curate view: choose a CDF and the evidence, policy and strictness, preview the exact command, run the read-only
   screening, follow its log and open the report. The form is built from curate.py's own option schema. Curation is a
   screening, never a scientific certification, and an exit status that found nothing is not an acceptance. */
(() => {
  const { h, badge, dl, table, card, api, post, state, button, field, fileName, views } = window.GUI;

  const form = { cdf: '', options: { trdat_mode: 'auto', trdat: '', evidence_kind: 'none', evidence: '', policy: '', fail_on: 'fail', require_complete: false, validate_only: false } };
  const ui = { schema: null, cdfs: [], files: { evidence: [], pilot: [], policy: [] }, plan: null, error: null, job: null, jobs: [], offset: 0, timer: null, log: '', notice: '' };

  function request() {
    const options = {};
    for (const [key, value] of Object.entries(form.options)) {
      if (value === '' || value === null) continue;
      if (key === 'trdat' && form.options.trdat_mode !== 'path') continue;
      if (key === 'evidence' && form.options.evidence_kind === 'none') continue;
      options[key] = value;
    }
    return { cdf: form.cdf || undefined, options };
  }
  let planTimer = null;
  const schedulePlan = () => { clearTimeout(planTimer); planTimer = setTimeout(refreshPlan, 350); };
  async function refreshPlan() {
    if (!form.cdf && !form.options.validate_only) { ui.plan = null; ui.error = { message: 'Choose a CDF to screen, or tick "Validate the files only".', errors: [] }; drawPlan(); return; }
    try { ui.plan = (await post('curate/plan', request())).plan; ui.error = null; }
    catch (error) { ui.plan = null; ui.error = { message: error.message, errors: (error.body && error.body.errors) || [] }; }
    applyFieldErrors();
    drawPlan();
  }
  const errorFor = (name) => ui.error && ui.error.errors.find((e) => e.field === name);
  function applyFieldErrors() {
    if (!ui.formcol) return;
    for (const wrapper of ui.formcol.querySelectorAll('.fieldwrap[data-name]')) {
      const err = errorFor(wrapper.dataset.name);
      wrapper.classList.toggle('invalid', !!err);
      let box = wrapper.querySelector(':scope > .ferror');
      if (err) { if (!box) { box = h('div', { class: 'ferror' }); wrapper.append(box); } box.textContent = err.message; }
      else if (box) box.remove();
    }
  }
  const help = (dest) => { const o = ui.schema && ui.schema.options.find((x) => x.dest === dest); return o ? o.help : ''; };
  function wrap(name, label, control, helpText) {
    const err = errorFor(name);
    return h('div', { class: `fieldwrap${err ? ' invalid' : ''}`, 'data-name': name }, field(label, control, helpText), err ? h('div', { class: 'ferror', text: err.message }) : null);
  }
  function select(name, label, choices, get, set, helpText) {
    const el = h('select', { 'aria-label': label }, choices.map(([value, text]) => h('option', { value, text, selected: String(get()) === String(value) })));
    el.addEventListener('change', () => { set(el.value); schedulePlan(); });
    return wrap(name, label, el, helpText);
  }
  function check(label, get, set, helpText) {
    const el = h('input', { type: 'checkbox', checked: !!get(), 'aria-label': label });
    el.addEventListener('change', () => { set(el.checked); schedulePlan(); });
    return h('div', { class: 'fieldwrap' }, h('label', { class: 'check' }, el, h('span', { text: label })), helpText ? h('span', { class: 'fhelp', text: helpText }) : null);
  }
  /* A picker over the files found below the input roots, plus a box for any other path below them. */
  function filePicker(name, label, files, get, set, helpText) {
    const known = files.some((f) => f.path === get());
    const picker = h('select', { 'aria-label': label }, h('option', { value: '', text: files.length ? 'choose a file…' : 'no matching file found below the input roots' }),
      files.map((f) => h('option', { value: f.path, selected: f.path === get(), text: `${f.name} · ${f.path.split('/').slice(-3, -1).join('/')}` })));
    picker.addEventListener('change', () => { set(picker.value); schedulePlan(); });
    const manual = h('input', { type: 'text', value: known ? '' : get(), placeholder: 'or a path below the input roots', 'aria-label': `${label} path` });
    manual.addEventListener('input', () => { set(manual.value.trim()); schedulePlan(); });
    return h('div', {}, wrap(name, label, picker, helpText), h('div', { class: 'fieldwrap' }, manual));
  }

  function sourceSection() {
    const known = ui.cdfs.some((c) => c.path === form.cdf), only = form.options.validate_only;
    const picker = h('select', { 'aria-label': 'CDF file', disabled: false }, h('option', { value: '', text: ui.cdfs.length ? 'choose a CDF…' : 'no <runid>.CDF files under the input roots' }),
      ui.cdfs.map((c) => h('option', { value: c.path, selected: c.path === form.cdf, text: `${fileName(c.path)} · ${(c.size / 1e6).toFixed(1)} MB · TR.DAT ${c.trdat ? 'found' : 'not found'} · ${c.path.split('/').slice(-3, -1).join('/')}` })));
    picker.addEventListener('change', () => { form.cdf = picker.value; schedulePlan(); });
    const manual = h('input', { type: 'text', value: known ? '' : form.cdf, placeholder: 'or a path below the input roots', 'aria-label': 'CDF path' });
    manual.addEventListener('input', () => { form.cdf = manual.value.trim(); schedulePlan(); });
    return h('section', { class: 'fsection' }, h('h3', { text: '1. What to screen' }),
      wrap('cdf', 'TRANSP output (CDF)', picker, only ? 'Optional when only validating the files.' : 'The CDF is only read, never modified.'), h('div', { class: 'fieldwrap' }, manual),
      check('Validate the files only (no CDF screening)', () => form.options.validate_only, (v) => { form.options.validate_only = v; }, help('validate_only')));
  }
  function inputsSection() {
    const o = form.options, kind = o.evidence_kind;
    const trdatChoices = [['auto', 'the TR.DAT beside the CDF, if there is one'], ['none', 'none (the timeline is reported as unknown)'], ['path', 'a file I choose']];
    return h('section', { class: 'fsection' }, h('h3', { text: '2. Evidence and policy' }),
      select('trdat_mode', 'TR.DAT', trdatChoices, () => o.trdat_mode, (v) => { o.trdat_mode = v; drawForm(); }, help('trdat')),
      o.trdat_mode === 'path' ? filePicker('trdat', 'TR.DAT file', [], () => o.trdat, (v) => { o.trdat = v; }) : null,
      select('evidence_kind', 'Evidence', [['none', 'none'], ['evidence', 'an evidence sidecar'], ['pilot', 'pilot facts']], () => o.evidence_kind, (v) => { o.evidence_kind = v; o.evidence = ''; drawForm(); },
        'Without either, rules that need facts you did not supply are reported as not_assessed, never as passed. ' + help('evidence')),
      kind !== 'none' ? filePicker('evidence', kind === 'pilot' ? 'Pilot facts file' : 'Evidence sidecar', ui.files[kind] || [], () => o.evidence, (v) => { o.evidence = v; }, 'Recognised by the schema_version inside the file.') : null,
      filePicker('policy', 'Policy', ui.files.policy || [], () => o.policy, (v) => { o.policy = v; }, help('policy')));
  }
  function strictnessSection() {
    const o = form.options;
    return h('section', { class: 'fsection' }, h('h3', { text: '3. How strict' }),
      select('fail_on', 'Exit status is "findings" at', [['fail', 'fail or error (default)'], ['error', 'error only'], ['warning', 'warning or worse']], () => o.fail_on, (v) => { o.fail_on = v; }, help('fail_on')),
      check('Also stop short when a rule did not run', () => o.require_complete, (v) => { o.require_complete = v; }, help('require_complete')),
      h('p', { class: 'fhelp', text: 'Plugins (--plugin) are not offered here: a plugin is a Python module that would run inside the screening. Use curate.py on the command line for that.' }));
  }

  /* ---- plan, run and results ---- */
  const KINDS = { succeeded: 'pass', findings: 'warning', failed: 'fail', cancelled: 'warning', lost: 'warning', running: 'not_assessed' };
  const WORDS = { succeeded: 'completed: nothing at the threshold', findings: 'completed: findings at the threshold' };
  function jobBadge(job) {
    return h('span', {}, badge('status', KINDS[job.state] || 'not_assessed'), ` ${WORDS[job.state] || job.state}${job.returncode !== null && job.returncode !== undefined ? ` (exit ${job.returncode})` : ''}`);
  }
  function planCard() {
    if (ui.schema && ui.schema.available === false) return card('Plan', h('p', { class: 'notice', text: ui.schema.unavailable_reason }));
    if (ui.error && !ui.plan) {
      const detail = ui.error.errors.length ? h('ul', {}, ui.error.errors.map((e) => h('li', { text: `${e.field}: ${e.message}` }))) : h('p', { text: ui.error.message });
      return card('Plan', h('p', { class: 'notice', text: 'This request cannot run yet.' }), detail);
    }
    const plan = ui.plan;
    if (!plan) return card('Plan', h('p', { class: 'empty', text: 'Choose a CDF to see the command.' }));
    const running = ui.job && ui.job.state === 'running' && ui.job.kind === 'curate';
    const run = h('button', { id: 'run-button', type: 'button', class: 'primary', text: running ? 'Running…' : plan.validate_only ? 'Validate files' : 'Run curation', disabled: !!running });
    run.addEventListener('click', start);
    return card('Plan', h('p', { class: 'muted', text: 'The exact command that will run (an argument list, no shell):' }),
      h('pre', { class: 'command', text: plan.display }), h('div', { class: 'toolbar' }, button('Copy command', () => navigator.clipboard && navigator.clipboard.writeText(plan.display))),
      dl([['CDF', plan.cdf], ['TR.DAT', plan.trdat || 'none'], ['evidence', plan.evidence ? `${plan.evidence.kind}: ${plan.evidence.path}` : 'none'], ['policy', plan.policy || 'default (every rule)'],
        ['exit status is "findings" at', plan.fail_on], ['report', plan.report ? h('code', { text: plan.report }) : 'none (validation only)']]),
      plan.warnings.length ? h('ul', { class: 'warnings' }, plan.warnings.map((w) => h('li', { text: w }))) : null,
      h('div', { class: 'toolbar' }, run, running ? button('Cancel', cancel) : null));
  }
  async function start() {
    try { const result = await post('curate/start', request()); ui.job = result.job; ui.offset = 0; ui.log = ''; ui.notice = ''; pollJob(); drawPlan(); }
    catch (error) { ui.notice = error.message; drawPlan(); }
  }
  async function cancel() { try { await post(`jobs/${ui.job.id}/cancel`, {}); } catch (error) { ui.notice = error.message; drawPlan(); } }
  async function pollJob() {
    clearTimeout(ui.timer);
    if (!ui.job) return;
    try {
      const chunk = await api(`jobs/${ui.job.id}/log`, { offset: ui.offset });
      ui.offset = chunk.next; ui.log += chunk.text;
      ui.job = (await api(`jobs/${ui.job.id}`)).job;
      updateLog();
      if (ui.job.state === 'running') ui.timer = setTimeout(pollJob, 800);
      else { ui.jobs = (await api('jobs')).jobs; drawPlan(); }
    } catch (error) { ui.notice = error.message; drawPlan(); }
  }
  function updateLog() {
    const logEl = document.getElementById('job-log'), status = document.getElementById('job-state');
    if (logEl) { logEl.textContent = ui.log || ''; logEl.scrollTop = logEl.scrollHeight; }
    if (status && ui.job) status.replaceChildren(jobBadge(ui.job));
  }
  function outcomeBlock(job) {
    const o = job.outcome;
    if (!o) return null;
    return h('div', { id: 'job-outcome' }, h('p', { class: 'muted', text: 'Findings by status (a screening record, not a certification):' }),
      h('div', { class: 'toolbar' }, Object.entries(o.counts).map(([key, n]) => h('span', {}, badge('status', key), ` ${n}  `))),
      o.scientifically_certified ? null : h('p', { class: 'fhelp', text: 'Scientific certification: NOT performed.' }));
  }
  function jobCard() {
    const job = ui.job;
    if (!job) return card('Run', h('p', { class: 'empty', text: 'No curation has been started in this session. Past jobs are listed below.' }));
    const actions = [];
    if (job.state !== 'running' && job.report_written) actions.push(button('Open curation report', () => window.GUI.openReport(job.report), { id: 'open-report' }));
    return card('Run', h('div', { class: 'toolbar' }, h('span', { id: 'job-state' }, jobBadge(job)), h('code', { id: 'job-id', text: job.id })),
      outcomeBlock(job), h('pre', { id: 'job-log', class: 'joblog', text: ui.log || '' }), h('div', { class: 'toolbar' }, actions),
      job.note ? h('p', { class: 'notice', text: job.note }) : null);
  }
  async function showJob(summary) { clearTimeout(ui.timer); ui.job = (await api(`jobs/${summary.id}`)).job; ui.offset = 0; ui.log = ''; pollJob(); drawPlan(); }
  function historyCard() {
    const rows = ui.jobs.filter((j) => j.kind === 'curate');
    return card('Recent curation jobs', table([{ label: 'Job', render: (j) => h('code', { text: j.id }) }, { label: 'State', render: (j) => jobBadge(j) },
      { label: 'CDF', render: (j) => fileName((j.plan && j.plan.cdf) || '') || '(validation only)' },
      { label: 'Findings', render: (j) => (j.outcome ? Object.entries(j.outcome.counts).filter(([, n]) => n).map(([k, n]) => `${k} ${n}`).join(', ') || 'none' : '') },
      { label: 'Started', render: (j) => new Date(j.started * 1000).toLocaleString() }, { label: '', render: (j) => button('Show', () => showJob(j)) }], rows, 'No curation jobs yet.'));
  }

  function keepFocus(container, rebuild) {
    const active = document.activeElement, label = active && container.contains(active) && active.getAttribute && active.getAttribute('aria-label');
    const start = active && active.selectionStart, end = active && active.selectionEnd;
    rebuild();
    if (!label) return;
    const again = container.querySelector(`[aria-label="${label}"]`);
    if (again && again.focus) { again.focus(); try { if (start !== null && start !== undefined) again.setSelectionRange(start, end); } catch (error) { /* no caret API */ } }
  }
  function drawForm() {
    if (!ui.formcol || state.view !== 'curate') return;
    keepFocus(ui.formcol, () => ui.formcol.replaceChildren(sourceSection(), inputsSection(), strictnessSection()));
    applyFieldErrors();
  }
  function drawPlan() {
    if (!ui.plancol || state.view !== 'curate') return;
    keepFocus(ui.plancol, () => ui.plancol.replaceChildren(...[ui.notice ? h('div', { class: 'notice', text: ui.notice }) : null, planCard(), jobCard(), historyCard()].filter(Boolean)));
    updateLog();
  }

  async function mount(main, argument) {
    ui.root = main; ui.formcol = h('div', { class: 'formcol' }); ui.plancol = h('div', { class: 'plancol' });
    main.replaceChildren(h('p', { class: 'empty', text: 'Loading the curation options…' }));
    try {
      const [schema, inputs, jobs] = await Promise.all([api('curate/schema'), api('curate/inputs'), api('jobs')]);
      ui.schema = schema; ui.cdfs = inputs.cdfs; ui.files = { evidence: [], pilot: [], policy: [], ...inputs.files }; ui.jobs = jobs.jobs;
      if (!form.cdf && ui.cdfs.length === 1) form.cdf = ui.cdfs[0].path;
      if (argument && argument.cdf) form.cdf = argument.cdf;
      if (state.view !== 'curate' || ui.root !== main) return;      // the user switched views while this was loading
      main.replaceChildren(h('div', { class: 'convert' }, ui.formcol, ui.plancol));
      drawForm();
      await refreshPlan();
      if (ui.job && ui.job.state === 'running') pollJob();
    } catch (error) { main.replaceChildren(h('p', { class: 'notice', text: `Cannot load the curation interface: ${error.message}` })); }
  }
  views.curate = { mount };
})();
