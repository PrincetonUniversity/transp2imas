'use strict';
/* The Explore view: plot converted IMAS data. An entry is read by a worker in the IMAS environment (read-only); the page
   shows curated recipes and any numeric node, scrubs time locally for profiles, and can overlay times, compare two
   entries, switch to a log axis and export PNG, SVG or CSV. */
(() => {
  const { h, text, fmt, badge, dl, table, card, api, state, button, fileName, download, views } = window.GUI;

  const ui = { df: null, root: null, entries: [], entry: '', describe: null, nodes: null, plot: null, compare: '', busy: '', error: '', chart: null, timer: null };
  const query = (params) => Object.fromEntries(Object.entries(params).filter(([, v]) => v !== undefined && v !== null && v !== ''));

  async function inspect(path) {
    clearInterval(ui.timer); ui.timer = null;
    ui.entry = path; ui.describe = null; ui.nodes = null; ui.plot = null; ui.df = null; ui.error = ''; ui.busy = 'Reading the entry…'; draw();
    try {
      ui.describe = await api('imas/describe', { entry: path });
      if (ui.describe.ids.some((i) => i.ids === 'distributions' && !i.error)) await loadDf();
    } catch (error) { ui.error = error.message; }
    ui.busy = ''; draw();
  }

  async function loadDf() {
    try {
      const info = await api('imas/df', { entry: ui.entry, what: 'info' });
      ui.df = { info, dist: 0, zone: '', e_min: '', e_max: '', p_min: '', p_max: '' };
    } catch (error) { ui.df = { error: error.message }; }
  }
  const signature = (d) => JSON.stringify((d && d.ids || []).map((i) => [i.ids, i.occurrence, i.bytes]));
  /* Coming back to an entry that may have been overwritten or extended meanwhile: re-read it, and drop the old plot and node
     list only if what is on disk changed. */
  async function refreshEntry() {
    ui.refreshing = true;
    draw();
    try {
      const fresh = await api('imas/describe', { entry: ui.entry });
      if (signature(fresh) !== signature(ui.describe)) { ui.plot = null; ui.nodes = null; ui.df = null; clearInterval(ui.timer); if (fresh.ids.some((i) => i.ids === 'distributions' && !i.error)) await loadDf(); }
      ui.describe = fresh; ui.error = '';
    } catch (error) { ui.describe = null; ui.plot = null; ui.error = error.message; }
    ui.refreshing = false;
    draw();
  }

  async function fetchPlot(entry, spec, timeIndex) {
    if (spec.df) return api('imas/df', query({ entry, ...spec.df, time_index: timeIndex || 0 }));
    const params = { entry, ...spec.query, time_index: timeIndex || 0 };
    if (spec.kind === 'profile') params.frames = '1';
    return api('imas/data', query(params));
  }

  async function plot(spec) {
    clearInterval(ui.timer);
    ui.busy = `Reading ${spec.title}…`; ui.error = ''; draw();
    try {
      const payload = await fetchPlot(ui.entry, spec, 0);
      ui.plot = { spec, payload, frame: 0, time: 0, overlay: false, log: false, compare: null, compareName: '', playing: false, fill: true, hide: new Set(), zone: spec.df && spec.df.zone !== '' && spec.df.zone !== undefined ? Number(spec.df.zone) : undefined, contours: spec.query && spec.query.recipe === 'eq_psi_map' ? 10 : 0 };
      if (ui.compare) await loadCompare();
    } catch (error) { ui.error = error.message; ui.plot = null; }
    ui.busy = ''; draw();
  }
  async function loadCompare() {
    const p = ui.plot;
    if (!p || !ui.compare || !['profile', 'timetrace'].includes(p.payload.kind) || p.spec.df) { if (p) p.compare = null; return; }
    try { p.compare = await fetchPlot(ui.compare, p.spec, 0); p.compareName = fileName(ui.compare.replace(/\/$/, '')); }
    catch (error) { p.compare = null; ui.error = `Comparison entry: ${error.message}`; }
  }

  /* ---- chart building from a payload ---- */
  function nearestFrame(frames, time) {
    let best = 0, gap = Infinity;
    frames.forEach((f, i) => { const d = Math.abs((f.time ?? 0) - time); if (d < gap) { gap = d; best = i; } });
    return best;
  }
  function seriesFor(p) {
    const payload = p.payload, out = [];
    if (payload.kind === 'timetrace') {
      payload.series.forEach((sr) => out.push({ name: sr.name, x: sr.x, y: sr.y }));
      if (p.compare) p.compare.series.forEach((sr) => out.push({ name: `${p.compareName}: ${sr.name}`, x: sr.x, y: sr.y, dash: true }));
    } else if (payload.kind === 'profile') {
      const frames = payload.frames || [], chosen = p.overlay && frames.length ? sampleFrames(frames) : [p.frame];
      // overlaid times: one colour per quantity, fading in with time (the legend lists each quantity once, at the latest time)
      chosen.forEach((index, order) => (frames[index] ? frames[index].series : payload.series).forEach((sr, k) => out.push(chosen.length > 1
        ? { name: `${sr.name} (latest t = ${fmt(frames[index].time, 4)} s)`, x: sr.x, y: sr.y, colour: k, opacity: 0.25 + 0.75 * (order / Math.max(1, chosen.length - 1)), legend: order === chosen.length - 1 }
        : { name: sr.name, x: sr.x, y: sr.y })));
      if (p.compare && p.compare.frames) {
        const base = frames[p.frame] ? frames[p.frame].time : 0, index = nearestFrame(p.compare.frames, base);
        p.compare.frames[index].series.forEach((sr) => out.push({ name: `${p.compareName}: ${sr.name} @ ${fmt(p.compare.frames[index].time, 4)} s`, x: sr.x, y: sr.y, dash: true }));
      }
    }
    return out;
  }
  function sampleFrames(frames) { const n = Math.min(8, frames.length); return [...new Set(Array.from({ length: n }, (_, i) => Math.round((i * (frames.length - 1)) / Math.max(1, n - 1))))]; }

  function plotPanel() {
    const p = ui.plot;
    if (!p) return card('Plot', h('p', { class: 'empty', text: ui.busy || 'Pick a recipe or select nodes to draw.' }), ui.error ? h('p', { class: 'notice', text: ui.error }) : null);
    const payload = p.payload, isMap = payload.kind === 'map', isShape = payload.kind === 'rz' || payload.kind === 'polys';
    const title = `${payload.title} — ${payload.ids}`;
    let chart;
    if (isShape) {
      chart = window.GUI.charts.shapeMap({ title, onZone: (zone) => { if (ui.df) { ui.df.zone = String(zone); } p.zone = zone; chart.set(payload, zone); draw(); } });
      chart.set(payload, p.zone);
    } else if (isMap) {
      chart = window.GUI.charts.heatmap({ title }); chart.set(payload);
      chart.configure({ fill: p.fill, contours: p.contours, hide: [...(p.hide || [])].map((n) => [n, true]) });
    } else {
      const first = payload.series[0] || {};
      chart = window.GUI.charts.lineChart({ title, xLabel: payload.x_label, xUnits: payload.x_units, yLabel: payload.series.length === 1 ? payload.series[0].name : 'value', yUnits: [...new Set(payload.series.map((x) => x.units).filter(Boolean))].join(', ') || first.units });
      chart.set(seriesFor(p)); chart.setLog(p.log);
    }
    ui.chart = chart;
    const controls = [];
    const time = payload.time;
    if (payload.kind === 'profile' && payload.frames && payload.frames.length) {
      const frames = payload.frames, slider = h('input', { type: 'range', min: 0, max: frames.length - 1, value: p.frame, 'aria-label': 'time' });
      const label = h('span', { class: 'timelabel', text: `t = ${fmt(frames[p.frame].time, 5)} s (slice ${frames[p.frame].index}${payload.time.stride > 1 || frames.length < payload.time.count ? ` of ${payload.time.count}` : ''})` });
      const move = (value) => { p.frame = Math.max(0, Math.min(frames.length - 1, value)); slider.value = p.frame; label.textContent = `t = ${fmt(frames[p.frame].time, 5)} s (slice ${frames[p.frame].index} of ${payload.time.count})`; chart.set(seriesFor(p)); };
      slider.addEventListener('input', () => move(Number(slider.value)));
      const play = button(p.playing ? 'Pause' : 'Play', () => {
        p.playing = !p.playing; play.textContent = p.playing ? 'Pause' : 'Play'; clearInterval(ui.timer);
        if (p.playing) ui.timer = setInterval(() => move(p.frame + 1 >= frames.length ? 0 : p.frame + 1), 120);
      });
      const overlay = h('input', { type: 'checkbox', checked: p.overlay, 'aria-label': 'overlay several times' });
      overlay.addEventListener('change', () => { p.overlay = overlay.checked; chart.set(seriesFor(p)); });
      controls.push(h('div', { class: 'toolbar' }, button('◀', () => move(p.frame - 1), { 'aria-label': 'previous time' }), play, button('▶', () => move(p.frame + 1), { 'aria-label': 'next time' }), slider, label,
        h('label', { class: 'check' }, overlay, h('span', { text: 'overlay 8 times' }))));
    }
    if ((isMap || isShape || (payload.kind === 'profile' && !(payload.frames && payload.frames.length))) && time && time.count > 1) {
      const slider = h('input', { type: 'range', min: 0, max: time.count - 1, value: p.time, 'aria-label': 'time' });
      const label = h('span', { class: 'timelabel', text: `t = ${fmt(time.value, 5)} s (slice ${time.index} of ${time.count})` });
      const loadSlice = async (index) => {
        if (p.loadingTimeSlice) return;
        const target = Math.max(0, Math.min(time.count - 1, index));
        slider.value = target;
        p.loadingTimeSlice = true;
        try {
          p.payload = await fetchPlot(ui.entry, p.spec, target);
          p.time = target;
          ui.error = '';
        } catch (error) {
          ui.error = error.message;
          p.playing = false;
          clearInterval(ui.timer);
          ui.timer = null;
        } finally {
          p.loadingTimeSlice = false;
          draw();
        }
      };
      slider.addEventListener('input', () => { label.textContent = `slice ${slider.value} of ${time.count} (release to load)`; });
      slider.addEventListener('change', () => loadSlice(Number(slider.value)));
      const row = [button('◀', () => loadSlice(p.time - 1), { 'aria-label': 'previous time' })];
      if (isMap) {
        const play = button(p.playing ? 'Pause' : 'Play', () => {
          p.playing = !p.playing;
          play.textContent = p.playing ? 'Pause' : 'Play';
          clearInterval(ui.timer);
          ui.timer = null;
          if (p.playing) ui.timer = setInterval(() => {
            if (!p.loadingTimeSlice) loadSlice(p.time + 1 >= time.count ? 0 : p.time + 1);
          }, 350);
        });
        row.push(play);
      }
      row.push(button('▶', () => loadSlice(p.time + 1), { 'aria-label': 'next time' }), slider, label);
      controls.push(h('div', { class: 'toolbar' }, ...row));
    }
    if (isMap) {
      const fill = h('input', { type: 'checkbox', checked: p.fill, 'aria-label': 'colour the map' });
      fill.addEventListener('change', () => { p.fill = fill.checked; chart.configure({ fill: p.fill }); });
      const count = h('input', { type: 'number', min: 0, max: 60, step: 1, value: p.contours, class: 'narrow', 'aria-label': 'number of contour lines' });
      count.addEventListener('input', () => { p.contours = Math.max(0, Math.min(60, Math.round(Number(count.value) || 0))); chart.configure({ contours: p.contours }); });
      const row = [h('label', { class: 'check' }, fill, h('span', { text: 'colour map' })),
        h('label', { class: 'field inline' }, h('span', { class: 'flabel', text: payload.flux ? 'Flux surfaces' : 'Contour lines' }), count)];
      [...new Set((payload.overlays || []).map((o) => o.name.split(' ')[0]))].forEach((name) => {
        const box = h('input', { type: 'checkbox', checked: !p.hide.has(name), 'aria-label': `show ${name}` });
        box.addEventListener('change', () => { box.checked ? p.hide.delete(name) : p.hide.add(name); chart.configure({ hide: [[name, !box.checked]] }); });
        row.push(h('label', { class: 'check' }, box, h('span', { text: name })));
      });
      controls.push(h('div', { class: 'toolbar' }, ...row));
    }
    if (!isMap && !isShape) {
      const log = h('input', { type: 'checkbox', checked: p.log, 'aria-label': 'logarithmic y axis' });
      log.addEventListener('change', () => { p.log = log.checked; chart.setLog(p.log); });
      controls.push(h('div', { class: 'toolbar' }, h('label', { class: 'check' }, log, h('span', { text: 'log y axis' }))));
    }
    const noCompare = isMap || isShape || !!p.spec.df;
    const compareSelect = h('select', { 'aria-label': 'compare with entry' }, h('option', { value: '', text: noCompare ? 'comparison is available for profiles and time traces' : 'none' }),
      ...ui.entries.filter((e) => e.path !== ui.entry && !noCompare).map((e) => h('option', { value: e.path, selected: e.path === ui.compare, text: entryLabel(e) })));
    compareSelect.disabled = noCompare;
    compareSelect.addEventListener('change', async () => { ui.compare = compareSelect.value; ui.busy = 'Reading the comparison entry…'; draw(); await loadCompare(); ui.busy = ''; draw(); });
    controls.push(h('div', { class: 'toolbar' }, h('label', { class: 'field inline' }, h('span', { class: 'flabel', text: 'Compare with' }), compareSelect)));
    const stem = `${payload.ids}-${(payload.title || 'plot').replace(/[^A-Za-z0-9]+/g, '_')}`;
    controls.push(h('div', { class: 'toolbar' }, h('span', { class: 'muted', text: 'Save:' }), button('PNG', async () => download(`${stem}.png`, 'image/png', await chart.exportPng())),
      ...(isMap || isShape ? [] : [button('SVG', () => download(`${stem}.svg`, 'image/svg+xml', chart.exportSvg()))]), button('CSV', () => download(`${stem}.csv`, 'text/csv', chart.csv()))));
    const notes = [];
    if (payload.missing && payload.missing.length) notes.push(h('p', { class: 'notice', text: `Not present in this entry (not drawn): ${[...new Set(payload.missing)].join(', ')}` }));
    if (payload.series && payload.series.some((x) => x.stride > 1)) notes.push(h('p', { class: 'muted', text: 'Long series are thinned for drawing (every ' + Math.max(...payload.series.map((x) => x.stride || 1)) + 'th sample); the data on disk is unchanged.' }));
    if (p.compare && payload.time && p.compare.time) notes.push(h('p', { class: 'muted', text: `Time bases differ: this entry spans ${fmt(payload.time.values[0], 4)} to ${fmt(payload.time.values[payload.time.values.length - 1], 4)} s, the comparison entry ${fmt(p.compare.time.values[0], 4)} to ${fmt(p.compare.time.values[p.compare.time.values.length - 1], 4)} s.` }));
    if (payload.kind === 'map') notes.push(h('p', { class: 'muted', text: 'dim1 and dim2 follow the grid type stored with the entry (for a rectangular grid they are R and Z in metres).' }));
    if (payload.kind === 'map' && payload.stride && (payload.stride[0] > 1 || payload.stride[1] > 1)) notes.push(h('p', { class: 'muted', text: `The map is thinned by ${payload.stride.join(' x ')} for drawing.` }));
    return card('Plot', ui.busy ? h('p', { class: 'muted', text: ui.busy }) : null, ui.error ? h('p', { class: 'notice', text: ui.error }) : null, chart.el, ...controls, ...notes);
  }

  /* ---- entry, recipes and node browser ---- */
  const entryLabel = (e) => `${e.machine || '?'} ${e.shot || ''}/${e.run || ''} (DD ${e.dd_major || '?'}) · ${e.path.split('/').slice(-5, -4)[0] || ''} · ${(e.bytes / 1e6).toFixed(1)} MB`;

  function entryCard() {
    const select = h('select', { 'aria-label': 'entry' }, h('option', { value: '', text: ui.entries.length ? 'choose a converted entry…' : 'no HDF5 entries found' }),
      ui.entries.map((e) => h('option', { value: e.path, selected: e.path === ui.entry, text: entryLabel(e) })));
    select.addEventListener('change', () => { if (select.value) inspect(select.value); });
    const manual = h('input', { type: 'text', placeholder: 'or the path of an entry directory (…/<machine>/<DD>/<shot>/<run>)', value: ui.entries.some((e) => e.path === ui.entry) ? '' : ui.entry, 'aria-label': 'entry path' });
    manual.addEventListener('keydown', (ev) => { if (ev.key === 'Enter' && manual.value.trim()) inspect(manual.value.trim()); });
    return card('Entry', h('div', { class: 'toolbar' }, select, button('Rescan', refreshEntries)), h('div', { class: 'toolbar' }, manual, button('Open', () => manual.value.trim() && inspect(manual.value.trim()))),
      state.config.read_database_roots.length ? h('p', { class: 'fhelp', text: 'Readable databases: ' + state.config.read_database_roots.join(', ') }) : h('p', { class: 'fhelp', text: 'No database directory exists yet; convert something first.' }));
  }

  function contentsCard() {
    const d = ui.describe;
    if (!d) return ui.error ? card('Entry contents', h('p', { class: 'notice', text: ui.error })) : card('Entry contents', h('p', { class: 'empty', text: ui.busy || 'Choose an entry.' }));
    const byIds = {};
    d.recipes.forEach((r) => (byIds[r.ids] = byIds[r.ids] || []).push(r));
    return card(`Entry contents — DD ${d.dd_version}`, ui.refreshing ? h('p', { class: 'muted', id: 'refreshing', text: 'Checking the entry for changes…' }) : null, ...d.ids.map((r) => h('div', { class: 'idsblock' },
      h('div', { class: 'idshead' }, h('strong', { text: r.ids + (r.occurrence ? ` [${r.occurrence}]` : '') }),
        h('span', { class: 'muted', text: ` ${(r.bytes / 1e6).toFixed(2)} MB` + (r.n_times ? ` · ${r.n_times} slices · ${fmt(r.t_first, 5)} to ${fmt(r.t_last, 5)} s` : ' · no time axis') }),
        button('Browse nodes', () => loadNodes(r.ids, r.occurrence), { class: 'small' })),
      r.error ? h('p', { class: 'notice', text: r.error }) : null,
      ...recipeGroups(byIds[r.ids] || [], r),
      r.ids === 'distributions' ? dfPanel() : null)));
  }

  function recipeGroups(list, r) {
    const groups = new Map();
    list.forEach((rc) => (groups.get(rc.group || 'Other') || groups.set(rc.group || 'Other', []).get(rc.group || 'Other')).push(rc));
    return [...groups].map(([name, items]) => h('div', { class: 'recipegroup' }, groups.size > 1 ? h('span', { class: 'grouplabel', text: name }) : null,
      h('div', { class: 'recipes' }, items.map((rc) => button(rc.title, () => plot({ title: rc.title, kind: rc.kind, query: { recipe: rc.id, occurrence: r.occurrence } }))))));
  }

  /* The fast-ion distribution function f(zone, energy, pitch): its own panel because it is a 3-D array on a general grid. */
  function dfPanel() {
    const df = ui.df;
    if (!df) return null;
    if (df.error) return h('p', { class: 'notice', text: `Distribution function: ${df.error}` });
    const info = df.info;
    const number = (key, label, help) => {
      const input = h('input', { type: 'number', class: 'narrow', value: df[key], placeholder: help, 'aria-label': label, step: 'any' });
      input.addEventListener('input', () => { df[key] = input.value; });
      return h('label', { class: 'field inline' }, h('span', { class: 'flabel', text: label }), input);
    };
    const which = h('select', { 'aria-label': 'distribution' }, info.distributions.map((d) => h('option', { value: d.index, selected: d.index === df.dist, text: d.label })));
    which.addEventListener('change', () => { df.dist = Number(which.value); });
    const request = (what) => ({ what, dist: df.dist, zone: df.zone, e_min: df.e_min, e_max: df.e_max, p_min: df.p_min, p_max: df.p_max, occurrence: 0 });
    const go = (title, what) => button(title, () => plot({ title, kind: what === 'zones' ? 'polys' : what === 'ep' ? 'map' : 'profile', df: request(what) }), { id: `df-${what}` });
    return h('div', { class: 'dfpanel' }, h('div', { class: 'grouplabel', text: 'Fast-ion distribution function' }),
      h('p', { class: 'fhelp', text: `${info.zones} zones × ${info.energy_bins} energy × ${info.pitch_bins} pitch bins, energy ${fmt(info.energy_range[0], 4)}–${fmt(info.energy_range[1], 5)} eV, pitch ${info.pitch_range[0]} to ${info.pitch_range[1]}. `
        + 'With no zone chosen the plots sum over the plasma volume (particle numbers); with a zone they show the stored value per m³ in that zone. Click a zone on the R,Z map to choose it.' }),
      h('div', { class: 'toolbar' }, which, number('zone', 'Zone', 'all'), number('e_min', 'Energy from [eV]', 'none'), number('e_max', 'to', 'none'), number('p_min', 'Pitch from', 'none'), number('p_max', 'to', 'none')),
      h('div', { class: 'recipes' }, go('Energy–pitch map', 'ep'), go('Energy spectrum', 'energy'), go('Pitch distribution', 'pitch'), go('Density by zone in R,Z', 'zones')));
  }

  async function loadNodes(ids, occurrence) {
    ui.busy = `Listing the nodes of ${ids}…`; ui.error = ''; draw();
    try { ui.nodes = { ...(await api('imas/nodes', { entry: ui.entry, ids, occurrence })), selected: new Set(), filter: '', kind: 'all' }; }
    catch (error) { ui.error = error.message; }
    ui.busy = ''; draw();
  }
  function nodesCard() {
    const n = ui.nodes;
    if (!n) return null;
    const filter = h('input', { type: 'text', placeholder: 'filter by path', value: n.filter, 'aria-label': 'filter nodes' });
    const kind = h('select', { 'aria-label': 'node kind' }, ['all', 'profile', 'timetrace', 'map'].map((k) => h('option', { value: k, selected: k === n.kind, text: k })));
    const list = h('div', { class: 'nodelist' });
    const redraw = () => {
      const rows = n.nodes.filter((x) => (n.kind === 'all' || x.kind === n.kind) && x.path.toLowerCase().includes(n.filter.toLowerCase()));
      list.replaceChildren(h('p', { class: 'muted', text: `${rows.length} of ${n.nodes.length} nodes` }), ...rows.map((x) => {
        const box = h('input', { type: 'checkbox', checked: n.selected.has(x.path), 'aria-label': x.path });
        box.addEventListener('change', () => { box.checked ? n.selected.add(x.path) : n.selected.delete(x.path); });
        return h('label', { class: 'node' }, box, h('code', { text: x.path }), h('span', { class: 'muted', text: ` ${x.kind} ${x.shape.length ? '[' + x.shape.join('×') + ']' : ''} ${x.units || ''}` }));
      }));
    };
    filter.addEventListener('input', () => { n.filter = filter.value; redraw(); });
    kind.addEventListener('change', () => { n.kind = kind.value; redraw(); });
    redraw();
    const draw1 = button('Plot selected', () => {
      const chosen = n.nodes.filter((x) => n.selected.has(x.path));
      if (!chosen.length) { ui.error = 'Select at least one node.'; draw(); return; }
      const kinds = new Set(chosen.map((x) => x.kind));
      if (kinds.size > 1) { ui.error = 'Select nodes of one kind (profile, time trace or map) to draw together.'; draw(); return; }
      plot({ title: chosen.length === 1 ? chosen[0].label : `${chosen.length} nodes of ${n.ids}`, kind: [...kinds][0], query: { ids: n.ids, occurrence: n.occurrence, paths: JSON.stringify(chosen.map((x) => x.path)) } });
    }, { class: 'primary' });
    return card(`Nodes of ${n.ids}`, h('p', { class: 'fhelp', text: 'Every numeric node found in the first time slice. [:] is the time index; draw nodes of one kind together.' }), h('div', { class: 'toolbar' }, filter, kind, draw1), list);
  }

  function draw() {
    if (!ui.root || state.view !== 'explore') return;
    ui.root.replaceChildren(h('div', { class: 'explore' }, h('div', { class: 'explorecol' }, entryCard(), contentsCard(), nodesCard()), h('div', { class: 'plotcol' }, plotPanel())));
  }
  async function refreshEntries() {
    try { ui.entries = (await api('entries')).entries; } catch (error) { ui.error = error.message; }
    draw();
  }
  async function mount(main, argument) {
    ui.root = main;
    main.replaceChildren(h('p', { class: 'empty', text: 'Looking for converted entries…' }));
    if (!state.config.imas_python) { main.replaceChildren(h('p', { class: 'notice', text: 'Plotting needs an interpreter with IMAS-Python. Start the viewer with --imas-python PATH (or run it in an IMAS environment).' })); return; }
    await refreshEntries();
    if (state.view !== 'explore' || ui.root !== main) return;        // the user switched views while this was loading
    const wanted = argument && argument.entry;
    if (wanted) inspect(wanted);
    else if (ui.entry) refreshEntry();
    else if (ui.entries.length === 1) inspect(ui.entries[0].path);
    else draw();
  }
  views.explore = { mount };
})();
