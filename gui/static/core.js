'use strict';
/* Shared helpers for every view. Every value from a report, an entry or a form is placed with textContent, a text
   node, or setAttribute on a fixed attribute name: such text never reaches an HTML parser, there are no inline handlers
   and no remote assets. */
(() => {

const SVG_NS = 'http://www.w3.org/2000/svg';
const state = { reports: [], current: null, tab: null, view: 'reports', config: null };

function h(tag, props, ...children) {
  const node = document.createElement(tag);
  fill(node, props, children);
  return node;
}
function s(tag, props, ...children) {
  const node = document.createElementNS(SVG_NS, tag);
  fill(node, props, children);
  return node;
}
function fill(node, props, children) {
  for (const [key, value] of Object.entries(props || {})) {
    if (value === null || value === undefined || value === false) continue;
    if (key === 'class') node.setAttribute('class', value);
    else if (key === 'text') node.textContent = value;
    else if (key === 'on') for (const [event, fn] of Object.entries(value)) node.addEventListener(event, fn);
    else node.setAttribute(key, value === true ? '' : String(value));
  }
  for (const child of children.flat()) {
    if (child === null || child === undefined || child === false) continue;
    node.append(child instanceof Node ? child : document.createTextNode(String(child)));
  }
}
const text = (value) => (value === null || value === undefined ? '' : typeof value === 'object' ? JSON.stringify(value) : String(value));
const fmt = (value, digits = 4) => (typeof value === 'number' ? Number(value.toPrecision(digits)).toString() : text(value));

function badge(kind, value) { return h('span', { class: `badge ${kind}-${value}`, text: value }); }
function pre(value) { return h('pre', { text: JSON.stringify(value, null, 2) }); }
function dl(pairs) {
  return h('dl', { class: 'kv' }, pairs.filter(([, v]) => v !== undefined && v !== null && v !== '')
    .flatMap(([k, v]) => [h('dt', { text: k }), h('dd', {}, v instanceof Node ? v : h('span', { text: text(v) }))]));
}
function table(columns, rows, empty) {
  if (!rows.length) return h('p', { class: 'empty', text: empty || 'Nothing to show.' });
  return h('table', {}, h('thead', {}, h('tr', {}, columns.map((c) => h('th', { text: c.label })))),
    h('tbody', {}, rows.map((row) => h('tr', {}, columns.map((c) => h('td', {}, c.render ? c.render(row) : text(row[c.key])))))));
}
function card(title, ...body) { return h('div', { class: 'card' }, h('h3', { text: title }), ...body); }

async function api(path, params) {
  const url = `/api/${path}?${new URLSearchParams(params || {})}`;
  const response = await fetch(url, { credentials: 'same-origin' });
  const body = await response.json();
  if (!response.ok) { const error = new Error(body.error || response.statusText); error.body = body; error.status = response.status; throw error; }
  return body;
}
/* A POST carries the token in a header (a cookie alone is refused), as JSON; the page learns the token from /api/config. */
async function post(path, payload) {
  const response = await fetch(`/api/${path}`, { method: 'POST', credentials: 'same-origin',
    headers: { 'Content-Type': 'application/json', 'X-GUI-Token': state.config.token }, body: JSON.stringify(payload || {}) });
  const body = await response.json();
  if (!response.ok) { const error = new Error(body.error || response.statusText); error.body = body; error.status = response.status; throw error; }
  return body;
}
function showError(error) {
  document.getElementById('view').replaceChildren(h('p', { class: 'notice', text: `Cannot show this: ${error.message}` }));
}
const fileName = (path) => String(path).split('/').pop();
function button(label, onClick, props) { return h('button', { type: 'button', ...(props || {}), on: { click: onClick }, text: label }); }
function field(label, control, help) {
  return h('label', { class: 'field' }, h('span', { class: 'flabel', text: label }), control, help ? h('span', { class: 'fhelp', text: help }) : null);
}
function download(name, mime, data) {
  const blob = data instanceof Blob ? data : new Blob([data], { type: mime });
  const url = URL.createObjectURL(blob);
  const link = h('a', { href: url, download: name });
  document.body.append(link); link.click(); link.remove();
  setTimeout(() => URL.revokeObjectURL(url), 2000);
}
const views = {};   // name -> { mount(container) }, registered by each view script

window.GUI = { h, s, text, fmt, badge, pre, dl, table, card, api, post, showError, state, fileName, button, field, download, views, fill };
})();
