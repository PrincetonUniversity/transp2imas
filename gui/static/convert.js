'use strict';
/* The Convert view: choose what to convert and how, preview the exact command and the target entry, run it, follow the
   log, and open the report or the converted entry. The form is built from the converter's own option schema. */
(() => {
  const { h, text, badge, pre, dl, table, card, api, post, state, button, field, fileName, views } = window.GUI;

  const CORE = ['core_profiles', 'equilibrium', 'summary', 'core_sources', 'core_transport', 'radiation', 'workflow'];
  const HARDWARE = ['wall', 'tf', 'pf_active', 'pf_passive', 'ec_launchers', 'ic_antennas', 'lh_antennas', 'nbi'];
  const form = { cdf: '', options: { backend: 'HDF5', dd_version: '4.1.1', new_db: true, ids: null, input_format: 'netcdf', pln_time_selection: 'all', ip_direction: 'auto', bt_direction: 'auto', moments: 'no',
    ec_nml_launch: 'no', nzones: 0, verbose: 0, log: '', run: '', occur: 0, machine: '', pf_passive_geometry: '', pf_active_geometry: '',
    lenient_mapping: false, debug: false, debug_max: false }, database: { mode: 'scratch', path: '' }, confirm: '' };
  const ui = { schema: null, inputs: [], plan: null, error: null, job: null, jobs: [], poll: null, offset: 0, timer: null, root: null };

  function clean(value) { return value === '' || value === null || value === undefined ? undefined : value; }
  function request() {
    const options = {};
    for (const [key, value] of Object.entries(form.options)) {
      if (['run', 'occur', 'nzones', 'verbose'].includes(key)) {
        if (value === '' || value === null) continue;
        const n = Number(value);
        if (!Number.isNaN(n)) options[key] = n;
        continue;
      }
      if (value === '') continue;
      options[key] = value;
    }
    if (options.run === 0) delete options.run;       // 0 means "derive the run number from the file name"
    return { cdf: form.cdf, options, database: { mode: form.database.mode, path: clean(form.database.path) }, confirm_overwrite: clean(form.confirm) };
  }

  let planTimer = null;
  function schedulePlan() { clearTimeout(planTimer); planTimer = setTimeout(refreshPlan, 350); }
  async function refreshPlan() {
    if (!form.cdf) { ui.plan = null; ui.error = { message: 'Choose a TRANSP input file.', errors: [] }; draw(); return; }
    try { ui.plan = (await post('convert/plan', request())).plan; ui.error = null; }
    catch (error) { ui.plan = null; ui.error = { message: error.message, errors: (error.body && error.body.errors) || [] }; }
    applyFieldErrors();
    drawPlan();
  }
  /* Mark invalid fields in place: the form is never rebuilt while someone is typing in it. */
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

  /* ---- form controls (each updates the form and schedules a new plan) ---- */
  const errorFor = (name) => (ui.error && ui.error.errors.find((e) => e.field === name || e.field === `options.${name}`));
  function help(dest) { const o = ui.schema && ui.schema.options.find((x) => x.dest === dest); return o ? o.help : ''; }
  function wrap(name, label, control, helpText) {
    const err = errorFor(name);
    return h('div', { class: `fieldwrap${err ? ' invalid' : ''}`, 'data-name': name }, field(label, control, helpText), err ? h('div', { class: 'ferror', text: err.message }) : null);
  }
  function select(name, label, choices, get, set, helpText) {
    const el = h('select', { 'aria-label': label }, choices.map((c) => h('option', { value: c, text: c === '' ? '(default)' : c, selected: String(get()) === String(c) })));
    el.addEventListener('change', () => { set(el.value); schedulePlan(); });
    return wrap(name, label, el, helpText);
  }
  function input(name, label, get, set, type, helpText, placeholder) {
    const el = h('input', { type: type || 'text', value: get() === undefined || get() === null ? '' : get(), placeholder, 'aria-label': label, ...(type === 'number' ? { min: 0 } : {}) });
    el.addEventListener('input', () => { set(type === 'number' ? el.value : el.value); schedulePlan(); });
    return wrap(name, label, el, helpText);
  }
  function check(name, label, get, set, helpText) {
    const el = h('input', { type: 'checkbox', checked: !!get(), 'aria-label': label });
    el.addEventListener('change', () => { set(el.checked); schedulePlan(); });
    return h('div', { class: 'fieldwrap' }, h('label', { class: 'check' }, el, h('span', { text: label })), helpText ? h('span', { class: 'fhelp', text: helpText }) : null);
  }

  function inputSection() {
    const cdfs = ui.inputs, known = cdfs.some((c) => c.path === form.cdf);
    const picker = h('select', { 'aria-label': 'TRANSP input file' }, h('option', { value: '', text: cdfs.length ? 'choose a TRANSP file…' : 'no CDF or XF.PLN files under the input roots' }),
      cdfs.map((c) => h('option', { value: c.path, selected: c.path === form.cdf, text: `${fileName(c.path)} · ${c.format === 'pln' ? 'native PLN' : 'NetCDF'} · ${(c.size / 1e6).toFixed(1)} MB · ${c.path.split('/').slice(-3, -1).join('/')}` })));
    picker.addEventListener('change', () => {
      form.cdf = picker.value;
      const selected = cdfs.find((c) => c.path === form.cdf);
      if (selected) form.options.input_format = selected.format;
      schedulePlan(); drawForm();
    });
    const manual = h('input', { type: 'text', value: known ? '' : form.cdf, placeholder: 'or a path below the input roots', 'aria-label': 'TRANSP input path' });
    manual.addEventListener('input', () => {
      form.cdf = manual.value.trim();
      if (form.cdf.toUpperCase().endsWith('XF.PLN')) form.options.input_format = 'pln';
      else if (form.cdf.toUpperCase().endsWith('.CDF')) form.options.input_format = 'netcdf';
      schedulePlan(); drawForm();
    });
    const fmt = select('input_format', 'Input format', ['netcdf', 'pln'], () => form.options.input_format, (v) => { form.options.input_format = v; }, 'NetCDF CDF remains the default. Native PLN uses XF.PLN and YF.PLN; TF.PLN supplies labels, or a matching CDF can supply the schema.');
    const plnTime = select('pln_time_selection', 'Native PLN times', ['all', 'last'], () => form.options.pln_time_selection, (v) => { form.options.pln_time_selection = v; }, 'Use all samples for the history or only the final sample for a time slice.');
    return h('section', { class: 'fsection' }, h('h3', { text: '1. What to convert' }), fmt,
      form.options.input_format === 'pln' ? plnTime : null,
      wrap('cdf', 'TRANSP output file', picker, 'CDF entries use NetCDF. For native input select the run XF.PLN; the matching YF.PLN and either TF.PLN or CDF schema must be beside it.'),
      h('div', { class: 'fieldwrap' }, manual));
  }

  function targetSection() {
    const o = form.options, dd = ui.schema.dd_versions, custom = !dd.includes(o.dd_version);
    const ddSelect = h('select', { 'aria-label': 'DD version' }, [...dd, 'custom…'].map((v) => h('option', { value: v, text: v, selected: v === (custom ? 'custom…' : o.dd_version) })));
    const ddText = h('input', { type: 'text', value: custom ? o.dd_version : '', placeholder: 'e.g. 4.0.0', hidden: !custom, 'aria-label': 'custom DD version' });
    ddSelect.addEventListener('change', () => { if (ddSelect.value === 'custom…') { ddText.hidden = false; o.dd_version = ddText.value || ''; } else { ddText.hidden = true; o.dd_version = ddSelect.value; } schedulePlan(); });
    ddText.addEventListener('input', () => { o.dd_version = ddText.value.trim(); schedulePlan(); });
    return h('section', { class: 'fsection' }, h('h3', { text: '2. How: target entry' }),
      select('backend', 'IMAS backend', ui.schema.backends, () => o.backend, (v) => { o.backend = v; }, 'HDF5 entries can be plotted here; MDSplus entries cannot (and their existence cannot be checked).'),
      wrap('dd_version', 'Data Dictionary version', h('div', { class: 'inline' }, ddSelect, ddText), help('dd_version')),
      input('machine', 'Machine / database name', () => o.machine, (v) => { o.machine = v; }, 'text', 'CDF default: the device named in the CDF. Native PLN needs this value.', 'machine name'),
      input('run', 'Run number', () => o.run, (v) => { o.run = v; }, 'number', help('run'), 'derived from the file name'),
      input('occur', 'Occurrence', () => o.occur, (v) => { o.occur = v; }, 'number', help('occur')),
      check('new_db', 'Create a new entry', () => o.new_db, (v) => { o.new_db = v; }, 'Replaces an existing entry of the same machine, DD major, shot and run; the next section shows whether it exists.'));
  }

  function databaseSection() {
    const config = state.config, mode = form.database.mode;
    const radio = (value, label, note, disabled) => {
      const el = h('input', { type: 'radio', name: 'dbmode', value, checked: mode === value, disabled: !!disabled, 'aria-label': label });
      el.addEventListener('change', () => { form.database.mode = value; form.confirm = ''; schedulePlan(); drawForm(); });
      return h('label', { class: `radio${disabled ? ' off' : ''}` }, el, h('span', {}, h('strong', { text: label }), h('span', { class: 'fhelp', text: note })));
    };
    return h('section', { class: 'fsection' }, h('h3', { text: '3. Where to write' }),
      radio('scratch', 'Scratch database (recommended)', config.scratch_database, false),
      radio('custom', 'A directory I choose', 'must be inside an output root: ' + config.output_roots.join(', '), false),
      mode === 'custom' ? input('database.path', 'Database directory', () => form.database.path, (v) => { form.database.path = v; }, 'text', '', config.output_roots[0]) : null,
      radio('default', 'My default IMAS database', config.allow_default_database ? config.default_database : 'not enabled on this server (--allow-default-database)', !config.allow_default_database));
  }

  function idsSection() {
    const all = ui.schema.ids, o = form.options, chosen = new Set(o.ids || all);
    const apply = (list) => { o.ids = list.length === all.length ? null : [...list]; schedulePlan(); drawForm(); };
    const boxes = all.map((name) => {
      const el = h('input', { type: 'checkbox', checked: chosen.has(name), 'aria-label': name });
      el.addEventListener('change', () => { const next = new Set(o.ids || all); el.checked ? next.add(name) : next.delete(name); apply(all.filter((n) => next.has(n))); });
      return h('label', { class: 'idsbox' }, el, ` ${name}`);
    });
    return h('section', { class: 'fsection' }, h('h3', { text: '4. Which IDSs' }),
      h('div', { class: 'toolbar' }, button('All', () => apply(all)), button('None', () => apply([])), button('Core plasma', () => apply(all.filter((n) => CORE.includes(n)))), button('Hardware', () => apply(all.filter((n) => HARDWARE.includes(n))))),
      h('div', { class: `fieldwrap${errorFor('ids') ? ' invalid' : ''}`, 'data-name': 'ids' }, h('div', { class: 'idsgrid' }, boxes), errorFor('ids') ? h('div', { class: 'ferror', text: errorFor('ids').message }) : null),
      h('p', { class: 'fhelp', text: `${chosen.size} of ${all.length} selected${o.ids ? ' — a partial entry' : ''}. Writers run in a fixed dependency-safe order.` }));
  }

  function mappingSection() {
    const o = form.options, sch = ui.schema.options;
    const choices = (dest) => { const x = sch.find((s) => s.dest === dest); return (x && x.choices) || []; };
    return h('section', { class: 'fsection' }, h('h3', { text: '5. Mapping and physics options' }),
      check('lenient_mapping', 'Skip variables that cannot be mapped', () => o.lenient_mapping, (v) => { o.lenient_mapping = v; },
        'A variable that exists but cannot be mapped (bad shape, unreadable) is recorded with its reason in the conversion report instead of stopping the conversion. Default: stop at the first one. Converter flag: --lenient-mapping.'),
      select('ip_direction', 'Plasma current direction', choices('ip_direction'), () => o.ip_direction, (v) => { o.ip_direction = v; }, help('ip_direction')),
      select('bt_direction', 'Toroidal field direction', choices('bt_direction'), () => o.bt_direction, (v) => { o.bt_direction = v; }, help('bt_direction')),
      input('pf_passive_geometry', 'Passive structure geometry (JSON)', () => o.pf_passive_geometry, (v) => { o.pf_passive_geometry = v; }, 'text', help('pf_passive_geometry'), 'optional path'),
      input('pf_active_geometry', 'Active coil geometry (JSON)', () => o.pf_active_geometry, (v) => { o.pf_active_geometry = v; }, 'text', help('pf_active_geometry'), 'optional path'),
      select('moments', 'Equilibrium moments', choices('moments'), () => o.moments, (v) => { o.moments = v; }, help('moments')),
      input('nzones', 'Poloidal points per surface', () => o.nzones, (v) => { o.nzones = v; }, 'number', help('nzones')),
      select('ec_nml_launch', 'EC launch from namelist', choices('ec_nml_launch'), () => o.ec_nml_launch, (v) => { o.ec_nml_launch = v; }, help('ec_nml_launch')));
  }

  function diagnosticsSection() {
    const o = form.options, sch = ui.schema.options, levels = ['', ...((sch.find((s) => s.dest === 'log') || {}).choices || [])];
    return h('section', { class: 'fsection' }, h('h3', { text: '6. Diagnostics' }),
      select('log', 'Log level', levels, () => o.log, (v) => { o.log = v; }, help('log')),
      input('verbose', 'Verbosity (0–3)', () => o.verbose, (v) => { o.verbose = v; }, 'number', help('verbose')),
      check('debug', 'Extra debug prints', () => o.debug, (v) => { o.debug = v; }),
      check('debug_max', 'All debug output (including plots)', () => o.debug_max, (v) => { o.debug_max = v; }));
  }

  /* ---- plan, run and results ---- */
  function planCard() {
    if (ui.error && !ui.plan) {
      const general = ui.error.errors.length ? h('ul', {}, ui.error.errors.map((e) => h('li', { text: `${e.field}: ${e.message}` }))) : h('p', { text: ui.error.message });
      return card('Plan', h('p', { class: 'notice', text: 'This request cannot run yet.' }), general);
    }
    const plan = ui.plan;
    if (!plan) return card('Plan', h('p', { class: 'empty', text: 'Choose a CDF to see the command and the target entry.' }));
    const t = plan.target, config = state.config;
    const phraseNeeded = plan.needs_confirmation, confirmed = !phraseNeeded || form.confirm === plan.confirm_phrase;
    const running = ui.job && ui.job.state === 'running';
    const setConfirm = (value, input) => { form.confirm = value; if (input) input.value = value; const run = document.getElementById('run-button'); if (run) run.disabled = !(value === plan.confirm_phrase && config.run_enabled && !running); };
    const confirmInput = h('input', { type: 'text', value: form.confirm, placeholder: plan.confirm_phrase, 'aria-label': 'confirm entry path' });
    confirmInput.addEventListener('input', () => setConfirm(confirmInput.value));
    const alternatives = [];
    if (t.exists && t.next_free_run) alternatives.push(button(`Use run ${t.next_free_run} instead (nothing exists there)`, () => { form.options.run = t.next_free_run; form.confirm = ''; drawForm(); refreshPlan(); }));
    const confirmBox = phraseNeeded ? h('div', { class: 'confirm' },
      h('p', { class: 'notice', text: t.exists ? `This entry already exists. Running will ${form.options.new_db ? 'DELETE its contents and write only the selected IDSs' : 'write the selected IDSs into it, replacing those IDSs'}. To go ahead, confirm below; to keep it, choose another run number.` : 'Whether this target exists cannot be checked for this backend, so it is treated as possibly existing.' }),
      h('div', { class: 'toolbar' }, ...alternatives),
      h('label', { class: 'field' }, h('span', { class: 'flabel', text: 'To overwrite: enter the entry path (or use the button) and press Run' }), confirmInput),
      h('div', { class: 'toolbar' }, button('Fill in the path', () => setConfirm(plan.confirm_phrase, confirmInput), { 'aria-label': 'fill in the entry path' }))) : null;
    const run = h('button', { id: 'run-button', type: 'button', class: 'primary', text: running ? 'Running…' : 'Run conversion', disabled: !(config.run_enabled && confirmed && !running) });
    run.addEventListener('click', start);
    const reasons = [];
    if (!config.run_requested) reasons.push('Running is switched off on this server, so the Run button is disabled. Restart it with --enable-run.');
    if (!config.imas_python) reasons.push(`Conversion needs an interpreter with IMAS-Python${config.imas_note ? ` (${config.imas_note})` : ''}. Restart with --imas-python /path/to/python, or the virtual-environment directory that has IMAS.`);
    const why = reasons.length ? h('div', { class: 'notice' }, reasons.map((r) => h('p', { text: r }))) : null;
    return card('Plan', h('p', { class: 'muted', text: 'The exact command that will run (an argument list, no shell):' }),
      h('pre', { class: 'command', text: plan.display }), h('div', { class: 'toolbar' }, button('Copy command', () => navigator.clipboard && navigator.clipboard.writeText(plan.display))),
      dl([['machine', t.machine], ['shot', t.shot], ['run', t.run], ['backend', t.backend], ['DD version', t.dd_version], ['occurrence', t.occurrence], ['database', `${t.database_mode}: ${t.database_root}`],
        ['entry', h('code', { text: t.primary })], ['entry exists', t.exists_unknown ? 'unknown (MDSplus)' : t.exists ? 'yes' : 'no'], ['IDSs', `${plan.selected_ids.length} selected`]]),
      plan.warnings.length ? h('ul', { class: 'warnings' }, plan.warnings.map((w) => h('li', { text: w }))) : null, confirmBox,
      h('div', { class: 'toolbar' }, run, ui.job && ui.job.state === 'running' ? button('Cancel', cancel) : null), why);
  }

  async function start() {
    try {
      const result = await post('convert/start', request());
      ui.job = result.job; ui.offset = 0; ui.log = ''; form.confirm = ''; pollJob(); draw();
    } catch (error) {
      ui.error = { message: error.message, errors: (error.body && error.body.errors) || [] };
      if (error.status === 409 && error.body && error.body.confirm_phrase) ui.error = null;
      ui.notice = error.message; draw();
    }
  }
  async function cancel() { try { await post(`jobs/${ui.job.id}/cancel`, {}); } catch (error) { ui.notice = error.message; draw(); } }

  async function pollJob() {
    clearTimeout(ui.timer);
    if (!ui.job) return;
    try {
      const chunk = await api(`jobs/${ui.job.id}/log`, { offset: ui.offset });
      ui.offset = chunk.next; ui.log = (ui.log || '') + chunk.text;
      ui.job = (await api(`jobs/${ui.job.id}`)).job;
      updateLog();
      if (ui.job.state === 'running') ui.timer = setTimeout(pollJob, 1000);
      else { ui.jobs = (await api('jobs')).jobs; draw(); }
    } catch (error) { ui.notice = error.message; draw(); }
  }
  function updateLog() {
    const logEl = document.getElementById('job-log'), status = document.getElementById('job-state');
    if (logEl) { logEl.textContent = ui.log || ''; logEl.scrollTop = logEl.scrollHeight; }
    if (status && ui.job) status.replaceChildren(jobBadge(ui.job));
  }
  function jobBadge(job) {
    const kind = { succeeded: 'pass', failed: 'fail', cancelled: 'warning', lost: 'warning', running: 'not_assessed' }[job.state] || 'not_assessed';
    return h('span', {}, badge('status', kind), ` ${job.state}${job.returncode !== null && job.returncode !== undefined ? ` (exit ${job.returncode})` : ''}`);
  }

  function jobCard() {
    const job = ui.job;
    if (!job) return card('Run', h('p', { class: 'empty', text: 'No conversion has been started in this session. Past jobs are listed below.' }));
    const entry = job.plan && job.plan.target && job.plan.target.primary;
    const actions = [];
    if (job.state !== 'running' && job.report_written) actions.push(button('Open conversion report', () => window.GUI.openReport(job.report)));
    if (job.state === 'succeeded' && entry && job.plan.target.backend.startsWith('HDF')) actions.push(button('Explore converted entry', () => window.GUI.showView('explore', { entry })));
    return card('Run', h('div', { class: 'toolbar' }, h('span', { id: 'job-state' }, jobBadge(job)), h('code', { id: 'job-id', text: job.id })),
      h('pre', { id: 'job-log', class: 'joblog', text: ui.log || '' }), h('div', { class: 'toolbar' }, actions),
      job.note ? h('p', { class: 'notice', text: job.note }) : null);
  }

  async function showJob(summary) {
    clearTimeout(ui.timer);
    ui.job = (await api(`jobs/${summary.id}`)).job; ui.offset = 0; ui.log = '';
    pollJob(); draw();
  }
  function historyCard() {
    return card('Recent jobs', table([{ label: 'Job', render: (j) => h('code', { text: j.id }) }, { label: 'State', render: (j) => jobBadge(j) },
      { label: 'Input', render: (j) => fileName((j.plan && j.plan.cdf) || '') }, { label: 'Entry', render: (j) => h('code', { text: (j.plan && j.plan.target && j.plan.target.primary) || '' }) },
      { label: 'Started', render: (j) => new Date(j.started * 1000).toLocaleString() }, { label: '', render: (j) => button('Show', () => showJob(j)) }], ui.jobs, 'No jobs yet.'));
  }

  /* The form column is rebuilt only for structural changes (database mode, IDS choice). Everything else (plan, job, history,
     notices) lives in the plan column, which is rebuilt on its own so that typing in the form is never disturbed. */
  function keepFocus(container, rebuild) {
    const active = document.activeElement, label = active && container.contains(active) && active.getAttribute && active.getAttribute('aria-label');
    const start = active && active.selectionStart, end = active && active.selectionEnd;
    rebuild();
    if (!label) return;
    const again = container.querySelector(`[aria-label="${label}"]`);
    if (again && again.focus) { again.focus(); try { if (start !== null && start !== undefined) again.setSelectionRange(start, end); } catch (error) { /* number inputs have no caret API */ } }
  }
  function drawForm() {
    if (!ui.formcol || state.view !== 'convert') return;
    keepFocus(ui.formcol, () => ui.formcol.replaceChildren(inputSection(), targetSection(), databaseSection(), idsSection(), mappingSection(), diagnosticsSection()));
  }
  function drawPlan() {
    if (!ui.plancol || state.view !== 'convert') return;
    keepFocus(ui.plancol, () => ui.plancol.replaceChildren(...[ui.notice ? h('div', { class: 'notice', text: ui.notice }) : null, planCard(), jobCard(), historyCard()].filter(Boolean)));
    updateLog();
  }
  const draw = drawPlan;

  async function mount(main, argument) {
    ui.root = main;
    ui.formcol = h('div', { class: 'formcol' });
    ui.plancol = h('div', { class: 'plancol' });
    main.replaceChildren(h('p', { class: 'empty', text: 'Loading the converter options…' }));
    try {
      [ui.schema, ui.inputs, ui.jobs] = await Promise.all([api('convert/schema'), api('inputs').then((r) => r.cdfs), api('jobs').then((r) => r.jobs)]);
      if (!form.cdf && ui.inputs.length === 1) {
        form.cdf = ui.inputs[0].path;
        form.options.input_format = ui.inputs[0].format || 'netcdf';
      }
      if (argument && argument.cdf) {
        form.cdf = argument.cdf;
        if (form.cdf.toUpperCase().endsWith('XF.PLN')) form.options.input_format = 'pln';
      }
      if (state.view !== 'convert' || ui.root !== main) return;      // the user switched views while this was loading
      main.replaceChildren(h('div', { class: 'convert' }, ui.formcol, ui.plancol));
      drawForm();
      await refreshPlan();
      if (ui.job && ui.job.state === 'running') pollJob();
    } catch (error) { main.replaceChildren(h('p', { class: 'notice', text: `Cannot load the conversion interface: ${error.message}` })); }
  }
  views.convert = { mount };
})();
