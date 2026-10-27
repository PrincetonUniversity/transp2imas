'use strict';
/* Report tabs and the view switcher. Shared helpers come from core.js. */
const { h, s, text, fmt, badge, pre, dl, table, card, api, post, showError, state, views } = window.GUI;

/* ---------- curation tabs ---------- */
function overview(v) {
  const counts = v.summary.counts, total = Object.values(counts).reduce((a, b) => a + b, 0) || 1;
  const src = v.header.source || {};
  const bar = h('div', { class: 'bar', role: 'img', 'aria-label': 'finding status counts' });
  for (const status of v.status_order) {
    if (!counts[status]) continue;
    const seg = h('span', { class: `seg-${status}`, title: `${status}: ${counts[status]}` });
    seg.style.width = `${(100 * counts[status]) / total}%`;
    bar.append(seg);
  }
  const cov = v.summary.rule_coverage;
  return h('div', { class: 'grid' },
    card('Findings', bar, dl(v.status_order.filter((k) => counts[k] !== undefined).map((k) => [k, counts[k]]))),
    card('Applications (preliminary)', table([{ label: 'Application', key: 'name' }, { label: 'Status', render: (r) => h('span', { text: r.status }) }],
      Object.entries(v.summary.applications).map(([name, status]) => ({ name, status })), 'No applications reported.'),
      v.summary.excused_findings && Object.values(v.summary.excused_findings).some((l) => l.length)
        ? h('p', { class: 'muted', text: 'Some findings are excused for an application by a declared event; they stay flagged in Findings.' }) : null),
    card('Source', dl([['path', src.path || 'in-memory dataset'], ['backend', src.backend], ['sha256', h('code', { text: src.sha256 || 'not available' })],
      ['evidence file', src.evidence_file && src.evidence_file.path], ['evidence kind', src.evidence_file && src.evidence_file.kind]])),
    card('Run', dl([['curation', v.header.curation_version], ['report schema', v.header.schema_version], ['generated', v.header.generated_at_utc],
      ['policy', v.header.policy_name], ['device', v.run_identity && v.run_identity.effective && v.run_identity.effective.device],
      ['shot', v.run_identity && v.run_identity.effective && v.run_identity.effective.shot],
      ['TRANSP version', v.run_identity && v.run_identity.transp_version]])),
    card('Fingerprints', v.header.fingerprints ? dl([['evidence', h('code', { text: v.header.fingerprints.evidence_sha256 })],
      ['policy', h('code', { text: v.header.fingerprints.policy_sha256 })], ['assessment', h('code', { text: v.header.fingerprints.assessment_sha256 })]])
      : h('p', { class: 'empty', text: 'Not in this report version.' })),
    card('Rule coverage', cov ? dl([['required rules executed', `${cov.required.executed} of ${cov.required.total}`],
      ['optional rules executed', `${cov.optional.executed} of ${cov.optional.total}`],
      ['not executed (unassessed)', Object.keys({ ...cov.required.not_executed, ...cov.optional.not_executed }).join(', ') || 'none']])
      : h('p', { class: 'empty', text: 'Not in this report version.' })),
    v.migration ? card('Evidence migration', dl([['declared', v.migration.from], ['raw sha256', h('code', { text: v.migration.raw_sha256 })]]),
      h('ul', {}, (v.migration.steps || []).map((step) => h('li', { text: step })))) : null);
}

function findings(v) {
  const statuses = v.status_order.filter((st) => v.findings.some((f) => f.status === st));
  const active = new Set(statuses.filter((st) => st !== 'pass'));
  const rules = [...new Set(v.findings.map((f) => f.rule).filter(Boolean))].sort();
  const search = h('input', { type: 'text', placeholder: 'search id or message', 'aria-label': 'search findings' });
  const rule = h('select', { 'aria-label': 'rule' }, h('option', { value: '', text: 'all rules' }), rules.map((r) => h('option', { value: r, text: r })));
  const list = h('div', {});
  const draw = () => {
    const q = search.value.toLowerCase();
    const shown = v.findings.filter((f) => active.has(f.status) && (!rule.value || f.rule === rule.value)
      && (!q || f.id.toLowerCase().includes(q) || (f.message || '').toLowerCase().includes(q)));
    shown.sort((a, b) => v.status_order.indexOf(a.status) - v.status_order.indexOf(b.status) || a.id.localeCompare(b.id));
    list.replaceChildren(h('p', { class: 'muted', text: `${shown.length} of ${v.findings.length} findings` }), ...shown.map(findingRow));
  };
  const toggles = statuses.map((st) => {
    const box = h('input', { type: 'checkbox', checked: active.has(st), 'aria-label': st, on: { change: () => { box.checked ? active.add(st) : active.delete(st); draw(); } } });
    return h('label', {}, box, badge('status', st), ` ${v.findings.filter((f) => f.status === st).length}`);
  });
  search.addEventListener('input', draw); rule.addEventListener('change', draw);
  draw();
  return h('div', {}, h('div', { class: 'toolbar' }, toggles, search, rule), list);
}
function findingRow(f) {
  const counts = f.counts || {};
  return h('details', { class: 'finding' },
    h('summary', {}, badge('status', f.status), h('code', { text: f.id }), h('span', { text: f.message })),
    h('div', { class: 'body' },
      dl([['rule', f.rule ? `${f.rule} v${f.rule_version}` : ''], ['variables', (f.variables || []).join(', ')], ['applications', (f.applications || []).join(', ')],
        ['samples', counts.evaluated_count !== undefined ? `${counts.evaluated_count} evaluated, ${counts.flagged_count} flagged` : ''],
        ['recommendation', f.recommendation]]),
      f.policy_exclusion ? card('Policy exclusion (raw finding unchanged)', pre(f.policy_exclusion)) : null,
      f.metrics && Object.keys(f.metrics).length ? card('Metrics', pre(f.metrics)) : null,
      f.evidence && Object.keys(f.evidence).length ? card('Evidence', pre(f.evidence)) : null,
      f.scope && Object.keys(f.scope).length ? card('Scope', pre(f.scope)) : null));
}

function rules(v) {
  const ledger = v.ledger || [];
  return h('div', {},
    table([{ label: 'Rule', render: (r) => h('code', { text: r.rule }) }, { label: 'Version', key: 'version' }, { label: 'Required', render: (r) => (r.required ? 'yes' : 'no') },
      { label: 'State', render: (r) => badge('state', r.state) },
      { label: 'Findings', render: (r) => Object.entries(r.findings || {}).map(([k, n]) => `${k} ${n}`).join(', ') },
      { label: 'Note', render: (r) => (r.error ? `${r.error.type}: ${r.error.message}` : r.missing ? `missing ${r.missing.join(', ')}` : r.note || '') }],
    ledger, 'This report has no execution ledger (older report format).'),
    h('p', { class: 'muted', text: 'A rule that did not execute leaves its checks unassessed: it is not a pass.' }),
    v.registry && v.registry.plugins && v.registry.plugins.length ? card('Plugins (trusted code)', table(
      [{ label: 'Module', key: 'module' }, { label: 'Version', key: 'version' }, { label: 'Rules', render: (r) => r.rules.join(', ') }, { label: 'sha256', render: (r) => h('code', { text: r.sha256 }) }], v.registry.plugins)) : null);
}

function evidence(v) {
  const g = v.evidence.graph;
  const detail = h('div', { class: 'card' }, h('p', { class: 'muted', text: 'Select a node.' }));
  let graphEl = h('p', { class: 'empty', text: 'No evidence graph in this report.' });
  if (g && g.nodes.length) {
    const pad = 16, w = 180, hh = 52;
    const svg = s('svg', { width: g.width + w + 2 * pad, height: g.height + hh + 2 * pad, role: 'img', 'aria-label': 'evidence lineage graph' },
      s('defs', {}, s('marker', { id: 'arrow', viewBox: '0 0 8 8', refX: 7, refY: 4, markerWidth: 7, markerHeight: 7, orient: 'auto' }, s('path', { d: 'M0,0 L8,4 L0,8 z', fill: 'var(--muted)' }))));
    const at = Object.fromEntries(g.nodes.map((n) => [n.id, n]));
    for (const e of g.edges) {
      const a = at[e.from], b = at[e.to];
      svg.append(s('path', { class: 'edge', 'marker-end': 'url(#arrow)', d: `M${a.x + pad + w},${a.y + pad + hh / 2} C${a.x + pad + w + 40},${a.y + pad + hh / 2} ${b.x + pad - 40},${b.y + pad + hh / 2} ${b.x + pad},${b.y + pad + hh / 2}` }));
    }
    for (const n of g.nodes) {
      const cls = ['node', n.usage.includes('independent_validation') ? 'validation' : n.usage.length ? 'constraint' : '', n.scope === false ? 'inactive' : ''].join(' ');
      const group = s('g', { class: cls, transform: `translate(${n.x + pad},${n.y + pad})`, tabindex: 0, role: 'button' },
        s('rect', { width: w, height: hh, rx: 6 }), s('text', { x: 8, y: 20 }, n.id.length > 24 ? `${n.id.slice(0, 23)}…` : n.id),
        s('text', { x: 8, y: 38, class: 'tick' }, n.usage.join(', ') || n.kind));
      const choose = () => detail.replaceChildren(h('h3', { text: n.id }), dl([['kind', n.kind], ['role', n.role], ['usage', n.usage.join(', ')], ['source', n.source], ['sha256', h('code', { text: n.sha256 || '' })],
        ['selector', n.selector && JSON.stringify(n.selector)], ['time scope', n.scope === false ? 'acts nowhere in the assessed interval' : n.scope ? `${n.scope[0]} to ${n.scope[1]}` : 'unscoped'],
        ['parents', n.parents.join(', ')], ['description', n.description]]));
      group.addEventListener('click', choose);
      group.addEventListener('keydown', (ev) => { if (ev.key === 'Enter') choose(); });
      svg.append(group);
    }
    graphEl = h('div', { style: false }, svg);
    graphEl.style.overflow = 'auto';
  }
  const unc = Object.entries(v.uncertainty || {});
  return h('div', {}, h('h3', { text: 'Lineage' }), h('p', { class: 'muted', text: 'Edges run from a source to what was derived from it. Shared lineage warns; no shared lineage is not proof of independence.' }), graphEl, detail,
    h('h3', { text: 'Time coverage of evidence' }),
    table([{ label: 'Finding', render: (r) => h('code', { text: r.id }) }, { label: 'Status', render: (r) => badge('status', r.status) }, { label: 'Coverage', key: 'state' }, { label: 'Clock', key: 'clock' },
      { label: 'Support', render: (r) => (r.support || []).map((x) => fmt(x)).join(' to ') }, { label: 'Record interval', render: (r) => (r.record_interval || []).map((x) => (x === null ? 'open' : fmt(x))).join(' to ') },
      { label: 'Unassessed', render: (r) => (r.unassessed_intervals || []).map((i) => i.map((x) => fmt(x)).join('–')).join(', ') }], v.evidence.coverage, 'No evidence record declared a time applicability.'),
    unc.length ? h('div', {}, h('h3', { text: 'Propagated uncertainty (supplied covariance)' }), unc.map(([id, u]) => card(id, dl([['scope', u.scope], ['cross-time', u.cross_time], ['quantity', u.quantity],
      ['slices', `${u.slices_evaluated} of ${u.slices_total}`], ['sigma range', u.std_gradient_range && u.std_gradient_range.map((x) => fmt(x)).join(' to ')]])))) : null);
}

function events(v) {
  const ev = v.events, declared = ev.declared, detected = ev.detected;
  const all = [...declared.map((e) => ({ ...e, kind: 'declared' })), ...detected.map((e) => ({ ...e, id: `${e.detector}`, kind: 'candidate' }))];
  let timeline = h('p', { class: 'empty', text: 'No events were declared and no detector produced candidates. That is not evidence of an event-free plasma.' });
  if (all.length) {
    const first = Math.min(...all.map((e) => e.start)), last = Math.max(...all.map((e) => e.end)), margin = (last - first || 1) * 0.15;
    const lo = first - margin, hi = last + margin, span = hi - lo, W = 760, rowH = 22;
    const x = (t) => 10 + ((t - lo) / span) * (W - 20);
    const svg = s('svg', { width: W, height: all.length * rowH + 34, role: 'img', 'aria-label': 'event timeline' });
    all.forEach((e, i) => {
      const cls = e.kind === 'candidate' ? 'ev-candidate' : e.acts_as_exclusion === false ? 'ev-inactive' : 'ev-declared';
      svg.append(s('rect', { class: cls, x: x(e.start), y: i * rowH + 4, width: Math.max(3, x(e.end) - x(e.start)), height: rowH - 8, rx: 3 }),
        s('text', { x: x(e.end) + 6 < W - 190 ? x(e.end) + 6 : x(e.start) + 6, y: i * rowH + 16 }, `${e.id} (${e.kind})`));
    });
    svg.append(s('line', { class: 'axis', x1: 10, x2: W - 10, y1: all.length * rowH + 8, y2: all.length * rowH + 8 }),
      s('text', { class: 'tick', x: 10, y: all.length * rowH + 24 }, `${fmt(lo)} s`), s('text', { class: 'tick', x: W - 70, y: all.length * rowH + 24 }, `${fmt(hi)} s`));
    timeline = h('div', {}, svg, h('div', { class: 'legend' }, h('span', {}, h('i', { class: 'ev-declared' }), 'declared (may exclude)'), h('span', {}, h('i', { class: 'ev-inactive' }), 'declared (cannot exclude)'), h('span', {}, h('i', { class: 'ev-candidate' }), 'automatic candidate')));
  }
  return h('div', {}, h('h3', { text: 'Timeline' }), timeline,
    h('h3', { text: 'Declared events' }),
    table([{ label: 'Event', render: (e) => h('code', { text: e.id }) }, { label: 'Type', key: 'type' }, { label: 'Origin', key: 'origin' },
      { label: 'Interval (s)', render: (e) => `${fmt(e.start)} to ${fmt(e.end)}` }, { label: 'Padding', key: 'padding' }, { label: 'Review', key: 'review_state' },
      { label: 'Acts as exclusion', render: (e) => (e.acts_as_exclusion ? 'yes' : 'no') }, { label: 'Experimental evidence', render: (e) => (e.independent_experimental_evidence ? 'yes' : 'no') },
      { label: 'Matched findings', render: (e) => (e.matched_findings || []).length }, { label: 'Not waivable', render: (e) => (e.mandatory_rules_not_waivable || []).join(', ') },
      { label: 'Notes', render: (e) => [e.outside_source_time ? 'outside source time' : '', e.unused ? 'matched nothing it can act on' : ''].filter(Boolean).join('; ') }],
      (v.findings.find((f) => f.id === 'EVENTS.INVENTORY') || { evidence: { events: [] } }).evidence.events || [], 'No events declared.'),
    h('h3', { text: 'Detected candidates' }),
    table([{ label: 'Detector', key: 'detector' }, { label: 'Diagnostic', key: 'diagnostic' }, { label: 'Interval (s)', render: (e) => `${fmt(e.start)} to ${fmt(e.end)}` },
      { label: 'Score (uncalibrated)', render: (e) => fmt(e.score) }, { label: 'Relative drop', render: (e) => fmt(e.relative_drop) }, { label: 'Review', key: 'review_state' }], detected, 'No candidates. This is not evidence of an event-free plasma.'),
    h('h3', { text: 'Policy exclusions (raw findings are unchanged)' }),
    table([{ label: 'Rule', render: (r) => h('code', { text: r.rule }) }, { label: 'Events', render: (r) => r.events.join(', ') }, { label: 'Origins', render: (r) => Object.values(r.origins || {}).join(', ') },
      { label: 'Excluded', key: 'excluded_count' }, { label: 'Eligible', key: 'eligible_count' }, { label: 'Unassessed', key: 'unassessed_count' }, { label: 'Flagged excluded', key: 'flagged_excluded_count' },
      { label: 'Covers all flagged', render: (r) => (r.covers_all_flagged ? 'yes' : 'no') }, { label: 'Policy status', key: 'policy_status' }, { label: 'Support', render: (r) => Object.values(r.support_basis || {}).join(', ') }],
      Object.entries(ev.policy_exclusions || {}).map(([rule, r]) => ({ rule, ...r })), 'No exclusions applied.'));
}

const COLORS = ['--c0', '--c1', '--c2', '--c3', '--c4'];
function masks(v) {
  const withMasks = v.findings.filter((f) => f.has_masks);
  if (!withMasks.length) return h('p', { class: 'empty', text: 'No mask file is stored beside this report.' });
  const select = h('select', { 'aria-label': 'finding' }, withMasks.map((f) => h('option', { value: f.mask_prefix, text: `${f.id} (${f.status})` })));
  const out = h('div', {});
  const load = async () => {
    try {
      const m = await api('masks', { path: state.current.path, prefix: select.value });
      out.replaceChildren(drawMask(m));
    } catch (error) { out.replaceChildren(h('p', { class: 'notice', text: error.message })); }
  };
  select.addEventListener('change', load);
  load();
  return h('div', {}, h('div', { class: 'toolbar' }, h('label', {}, 'Finding ', select)), out);
}
function drawMask(m) {
  const rows = m.grid.length, cols = m.grid[0].length, cell = Math.max(3, Math.min(60, Math.floor(640 / cols)));
  const canvas = h('canvas', { width: cols * cell, height: rows * cell, role: 'img', 'aria-label': `mask of ${m.prefix}` });
  const style = getComputedStyle(document.documentElement), ctx = canvas.getContext('2d');
  m.grid.forEach((row, i) => row.forEach((c, j) => { ctx.fillStyle = style.getPropertyValue(COLORS[c]).trim(); ctx.fillRect(j * cell, i * cell, cell, cell); }));
  const tip = h('div', { class: 'muted', text: 'Hover a cell.' });
  canvas.addEventListener('mousemove', (ev) => {
    const box = canvas.getBoundingClientRect(), j = Math.floor(((ev.clientX - box.left) / box.width) * cols), i = Math.floor(((ev.clientY - box.top) / box.height) * rows);
    if (i < 0 || j < 0 || i >= rows || j >= cols) return;
    const [a, b] = m.axes, first = (i * m.block[0]), second = (j * m.block[1]);
    const val = (axis, idx) => (axis && axis.values && axis.values[idx] !== undefined && axis.values[idx] !== null ? fmt(axis.values[idx]) : `index ${idx}`);
    tip.textContent = m.grid.length === 1 || m.shape.length === 1
      ? `${m.dimensions[0] || 'axis'} ${val(a, second)}: ${m.categories[m.grid[i][j]]}`
      : `${m.dimensions[0]} ${val(a, first)}, ${m.dimensions[1]} ${val(b, second)}: ${m.categories[m.grid[i][j]]}`;
  });
  return h('div', {}, canvas, tip, h('div', { class: 'legend' }, m.categories.map((name, k) => { const sw = h('i', {}); sw.style.background = `var(${COLORS[k]})`; return h('span', {}, sw, `${name} (${m.counts[name]})`); })),
    h('p', { class: 'muted', text: `Shape ${m.shape.join(' x ')}; each drawn cell covers ${m.block.join(' x ')} samples and shows the most severe category. Rows: ${m.dimensions[0]}${m.dimensions[1] ? `, columns: ${m.dimensions[1]}` : ''}.` }));
}

function inputs(v) {
  const t = v.timeline;
  if (!t.available) return h('p', { class: 'empty', text: 'No TR.DAT was supplied; controls and defaults remain unknown.' });
  const files = (t.files && t.files.entries) || [];
  return h('div', { class: 'grid' },
    card('Parser', dl([['kind', t.parser && t.parser.kind], ['version', t.parser && t.parser.version], ['states', t.states], ['update times valid', t.temporal_valid], ['assignments recorded', `${t.assignments_total}${t.assignments_truncated ? ' (truncated)' : ''}`], ['source sha256', h('code', { text: (t.source && t.source.sha256) || '' })]]),
      h('p', { class: 'muted', text: 'This is the declared state in file order, not the resolved TRANSP configuration.' })),
    card('Files followed', table([{ label: 'Role', key: 'role' }, { label: 'Path', render: (r) => h('code', { text: r.path || r.resolved_path || '' }) }, { label: 'sha256', render: (r) => h('code', { text: (r.sha256 || '').slice(0, 16) }) }], files, 'No referenced files.')),
    card('Unresolved constructs', table([{ label: 'Kind', key: 'kind' }, { label: 'Line', key: 'line' }, { label: 'Reason', key: 'reason' }], t.unresolved, 'None reported.')),
    card('Parser gaps', table([{ label: 'Kind', key: 'kind' }, { label: 'Line', key: 'line' }, { label: 'Name', key: 'name' }, { label: 'Reason', key: 'reason' }], t.gaps, 'None.')),
    v.code_identity ? card('Code identity', dl([['modules hashed', Object.keys(v.code_identity.modules || {}).length], ['dependencies', Object.entries(v.code_identity.dependencies || {}).map(([k, x]) => `${k} ${x || 'absent'}`).join(', ')]])) : null);
}

/* ---------- conversion tabs ---------- */
function conversionSummary(v) {
  const hd = v.header, p = v.publication;
  return h('div', { class: 'grid' },
    card('Conversion', dl([['status', badge('status', hd.status === 'PASS' ? 'pass' : 'fail')], ['machine', hd.machine], ['shot', hd.shot], ['run', hd.run], ['occurrence', hd.occurrence], ['backend', hd.backend], ['DD version', hd.dd_version], ['mapping policy', hd.mapping_policy], ['source CDF', hd.source_cdf], ['output', hd.output_directory], ['error', hd.error]])),
    card('Publication to the backend', dl([['status', p.status], ['accepted writes', p.accepted.length], ['failed', p.failed && JSON.stringify(p.failed)], ['unattempted', p.unattempted.length], ['failure kind', p.kind]]),
      h('p', { class: 'muted', text: 'Writes are sequential, not transactional: accepted writes stay accepted after a failure.' })),
    card('IDSs', dl([['selected', v.selected_ids.join(', ')], ['written', v.written_ids.map((w) => `${w.ids}[${w.occurrence}]`).join(', ')]])));
}
function writers(v) {
  return table([{ label: 'IDS', key: 'ids' }, { label: 'Occurrences', render: (r) => (r.occurrences || []).join(', ') }, { label: 'Status', key: 'status' }, { label: 'Seconds', render: (r) => fmt(r.seconds, 3) },
    { label: 'Writes', key: 'write_count' }, { label: 'Variables read', key: 'variables_read' }, { label: 'Consumed', key: 'consumed' }, { label: 'Mapping entries', key: 'mapping_entries' }, { label: 'Error', key: 'error' }], v.writers);
}
function journal(v) {
  return table([{ label: 'IDS', key: 'ids' }, { label: 'Kind', key: 'kind' }, { label: 'Target', key: 'target' }, { label: 'Count', key: 'count' }, { label: 'Reason', key: 'reason' }], v.mapping_journal,
    'No skipped, missing-optional or adjusted fields were recorded.');
}

const CURATION_TABS = [['Overview', overview], ['Findings', findings], ['Rules', rules], ['Evidence', evidence], ['Events', events], ['Masks', masks], ['Inputs', inputs]];
const CONVERSION_TABS = [['Summary', conversionSummary], ['Writers', writers], ['Mapping journal', journal]];

function render() {
  const cur = state.current, view = document.getElementById('view'), tabs = document.getElementById('tabs');
  if (state.view !== 'reports') return;
  if (!cur) { tabs.replaceChildren(); view.replaceChildren(h('p', { class: 'empty', text: 'Choose a report above.' })); return; }
  const defs = cur.view.kind === 'curation' ? CURATION_TABS : CONVERSION_TABS;
  if (!defs.some(([name]) => name === state.tab)) state.tab = defs[0][0];
  tabs.replaceChildren(...defs.map(([name]) => h('button', { class: 'tab', role: 'tab', type: 'button', 'aria-selected': name === state.tab, text: name,
    on: { click: () => { state.tab = name; render(); } } })));
  document.getElementById('notices').replaceChildren(...cur.view.notices.map((n) => h('div', { class: 'notice', text: n })));
  const [, draw] = defs.find(([name]) => name === state.tab);
  view.replaceChildren(draw(cur.view));
}

async function open(report) {
  try {
    const result = await api(report.kind, { path: report.path });
    state.current = { path: report.path, view: result.view };
    render();
  } catch (error) { showError(error); }
}
async function openPath(path) {
  for (const kind of ['curation', 'conversion']) {
    try { const result = await api(kind, { path }); state.current = { path, view: result.view }; render(); return; } catch (error) { if (kind === 'conversion') showError(error); }
  }
}
async function scan() {
  const select = document.getElementById('report-select');
  try {
    state.reports = (await api('scan')).reports;
    select.replaceChildren(h('option', { value: '', text: state.reports.length ? 'choose a report…' : 'no reports under the roots' }),
      ...state.reports.map((r, i) => h('option', { value: i, text: `${r.kind}: ${r.path}` })));
  } catch (error) { showError(error); }
}

const VIEW_NAMES = [['convert', 'Convert'], ['explore', 'Explore'], ['reports', 'Reports']];
function showView(name, argument) {
  state.view = name;
  for (const button of document.querySelectorAll('#views .view-button')) button.setAttribute('aria-selected', button.dataset.view === name);
  const reports = name === 'reports';
  for (const id of ['reportbar', 'tabs', 'notices']) document.getElementById(id).hidden = !reports;
  const main = document.getElementById('view');
  main.replaceChildren();
  if (reports) { render(); if (!state.current && state.reports.length === 1) open(state.reports[0]); return; }
  document.getElementById('notices').replaceChildren();
  views[name].mount(main, argument);
}
window.GUI.showView = showView;
window.GUI.openReport = async (path) => { showView('reports'); await openPath(path); };

async function start() {
  state.config = await api('config');
  const config = state.config;
  document.getElementById('roots').textContent = `Roots: ${config.roots.join(', ')} · viewer ${config.version}`;
  const badge = document.getElementById('mode-badge');
  badge.textContent = config.run_enabled ? 'conversion enabled' : 'read-only';
  badge.setAttribute('class', `pill ${config.run_enabled ? 'enabled' : 'readonly'}`);
  badge.setAttribute('title', config.run_enabled ? 'This server can start conversions into the scratch IMAS database' : 'No conversion can be started from this server');
  document.getElementById('views').replaceChildren(...VIEW_NAMES.map(([name, label]) => {
    const button = h('button', { class: 'view-button', type: 'button', 'data-view': name, text: label, 'aria-selected': false });
    button.addEventListener('click', () => showView(name));
    return button;
  }));
  document.getElementById('report-select').addEventListener('change', (ev) => { if (ev.target.value !== '') open(state.reports[Number(ev.target.value)]); });
  document.getElementById('rescan').addEventListener('click', scan);
  document.getElementById('open-form').addEventListener('submit', (ev) => { ev.preventDefault(); openPath(document.getElementById('open-path').value.trim()); });
  await scan();
  // Land on Convert (the first tab) unless the URL asks for a view or a report.
  const params = new URLSearchParams(location.search), wanted = params.get('open');
  const requested = params.get('view'), initial = VIEW_NAMES.some(([name]) => name === requested) ? requested : wanted ? 'reports' : 'convert';
  showView(initial);
  if (initial === 'reports' && wanted) openPath(wanted);
}
window.addEventListener('load', () => start().catch(showError));
