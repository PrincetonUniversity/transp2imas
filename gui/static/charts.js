'use strict';
/* Plain charts, no library: an SVG line chart (hover readout, legend toggles, log axis, gaps for missing values) and a
   canvas heat map with overlays. A nonfinite sample arrives as null and is drawn as a gap, never as zero. */
(() => {
  const { h, s, fmt, download } = window.GUI;
  const PALETTE = 8;

  function niceStep(span, target) {
    const raw = span / target, magnitude = 10 ** Math.floor(Math.log10(raw)), unit = raw / magnitude;
    return (unit < 1.5 ? 1 : unit < 3 ? 2 : unit < 7 ? 5 : 10) * magnitude;
  }
  function ticks(lo, hi, target = 6) {
    if (!(hi > lo)) return [lo];
    const step = niceStep(hi - lo, target), start = Math.ceil(lo / step) * step, out = [];
    for (let v = start; v <= hi + step * 1e-9; v += step) out.push(Math.abs(v) < step * 1e-9 ? 0 : v);
    return out;
  }
  function logTicks(lo, hi) {
    const out = [];
    for (let e = Math.floor(Math.log10(lo)); e <= Math.ceil(Math.log10(hi)); e += 1) out.push(10 ** e);
    return out;
  }
  const tickLabel = (v) => (v === 0 ? '0' : Math.abs(v) >= 1e5 || Math.abs(v) < 1e-3 ? v.toExponential(1) : Number(v.toPrecision(4)).toString());

  function lineChart(options) {
    const W = 780, H = 400, m = { l: 74, r: 18, t: 30, b: 50 }, uid = `c${Math.random().toString(36).slice(2, 8)}`;
    const model = { series: [], hidden: new Set(), log: false, extra: null };
    const svg = s('svg', { viewBox: `0 0 ${W} ${H}`, class: 'plot', role: 'img', 'aria-label': options.title || 'line chart' });
    const tip = h('div', { class: 'tip', hidden: true });
    const legend = h('div', { class: 'legend chart-legend' });
    const root = h('div', { class: 'chart' }, svg, tip, legend);

    function visible() { return model.series.filter((_, i) => !model.hidden.has(i)); }
    function extent() {
      let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
      for (const sr of visible()) {
        sr.x.forEach((x, i) => { const y = sr.y[i]; if (x === null || y === null || (model.log && y <= 0)) return; if (x < x0) x0 = x; if (x > x1) x1 = x; if (y < y0) y0 = y; if (y > y1) y1 = y; });
      }
      if (!isFinite(x0)) return null;
      if (x0 === x1) { x0 -= 0.5; x1 += 0.5; }
      if (y0 === y1) { const pad = Math.abs(y0) * 0.1 || 1; y0 -= pad; y1 += pad; }
      else if (!model.log) { const pad = (y1 - y0) * 0.05; y0 -= pad; y1 += pad; }
      return { x0, x1, y0, y1 };
    }
    function draw() {
      svg.replaceChildren();
      const box = extent(), pw = W - m.l - m.r, ph = H - m.t - m.b;
      svg.append(s('defs', {}, s('clipPath', { id: uid }, s('rect', { x: m.l, y: m.t, width: pw, height: ph }))));
      svg.append(s('text', { x: m.l, y: 18, class: 'title' }, options.title || ''));
      if (!box) { svg.append(s('text', { x: W / 2, y: H / 2, 'text-anchor': 'middle', class: 'tick' }, 'No finite data to draw')); legendDraw(); return; }
      const ly = (v) => Math.log10(v), yv = (v) => (model.log ? ly(v) : v);
      const yLo = yv(box.y0 > 0 || !model.log ? box.y0 : 1e-12), yHi = yv(box.y1);
      const px = (x) => m.l + ((x - box.x0) / (box.x1 - box.x0)) * pw, py = (y) => m.t + ph - ((yv(y) - yLo) / (yHi - yLo)) * ph;
      for (const t of ticks(box.x0, box.x1)) {
        svg.append(s('line', { class: 'grid', x1: px(t), x2: px(t), y1: m.t, y2: m.t + ph }), s('text', { class: 'tick', x: px(t), y: m.t + ph + 16, 'text-anchor': 'middle' }, tickLabel(t)));
      }
      for (const t of model.log ? logTicks(box.y0, box.y1) : ticks(box.y0, box.y1)) {
        if (t < box.y0 * (model.log ? 0.999 : 1) - 1e-300 || t > box.y1 * (model.log ? 1.001 : 1)) continue;
        svg.append(s('line', { class: 'grid', x1: m.l, x2: m.l + pw, y1: py(t), y2: py(t) }), s('text', { class: 'tick', x: m.l - 6, y: py(t) + 4, 'text-anchor': 'end' }, tickLabel(t)));
      }
      svg.append(s('rect', { class: 'frame', x: m.l, y: m.t, width: pw, height: ph }));
      svg.append(s('text', { class: 'label', x: m.l + pw / 2, y: H - 8, 'text-anchor': 'middle' }, `${options.xLabel || ''}${options.xUnits ? ` [${options.xUnits}]` : ''}`));
      svg.append(s('text', { class: 'label', transform: `translate(14,${m.t + ph / 2}) rotate(-90)`, 'text-anchor': 'middle' }, `${options.yLabel || ''}${options.yUnits ? ` [${options.yUnits}]` : ''}${model.log ? ' (log)' : ''}`));
      const group = s('g', { 'clip-path': `url(#${uid})` });
      model.series.forEach((sr, i) => {
        if (model.hidden.has(i)) return;
        let d = '', pen = false;
        sr.x.forEach((x, k) => {
          const y = sr.y[k];
          if (x === null || y === null || (model.log && y <= 0)) { pen = false; return; }
          d += `${pen ? 'L' : 'M'}${px(x).toFixed(1)},${py(y).toFixed(1)}`; pen = true;
        });
        const path = s('path', { class: `line s${(sr.colour ?? i) % PALETTE}${sr.dash ? ' dash' : ''}`, d });
        if (sr.opacity !== undefined) path.setAttribute('opacity', sr.opacity);
        group.append(path);
      });
      svg.append(group);
      const cross = s('line', { class: 'cross', y1: m.t, y2: m.t + ph, visibility: 'hidden' });
      const dots = s('g', {});
      const hit = s('rect', { class: 'hit', x: m.l, y: m.t, width: pw, height: ph });
      hit.addEventListener('mousemove', (ev) => {
        const bounds = svg.getBoundingClientRect(), x = box.x0 + (((ev.clientX - bounds.left) / bounds.width) * W - m.l) / pw * (box.x1 - box.x0);
        const rows = [], found = [];
        model.series.forEach((sr, i) => {
          if (model.hidden.has(i)) return;
          let best = -1, gap = Infinity;
          sr.x.forEach((xv, k) => { if (xv === null || sr.y[k] === null) return; const dx = Math.abs(xv - x); if (dx < gap) { gap = dx; best = k; } });
          if (best >= 0) { rows.push([sr.name, sr.x[best], sr.y[best], i]); found.push([px(sr.x[best]), py(sr.y[best]), i]); }
        });
        if (!rows.length) return;
        cross.setAttribute('visibility', 'visible'); cross.setAttribute('x1', found[0][0]); cross.setAttribute('x2', found[0][0]);
        dots.replaceChildren(...found.filter((f) => isFinite(f[1])).map((f) => s('circle', { class: `dot s${f[2] % PALETTE}`, cx: f[0], cy: f[1], r: 3.5 })));
        tip.hidden = false;
        tip.replaceChildren(h('div', { class: 'muted', text: `${options.xLabel || 'x'} = ${fmt(rows[0][1], 5)}` }), ...rows.map(([name, , y, i]) => h('div', {}, h('i', { class: `swatch s${i % PALETTE}` }), ` ${name}: ${fmt(y, 5)}${options.yUnits ? ' ' + options.yUnits : ''}`)));
        tip.style.left = `${Math.min(ev.clientX - root.getBoundingClientRect().left + 14, root.clientWidth - 220)}px`;
        tip.style.top = `${ev.clientY - root.getBoundingClientRect().top + 12}px`;
      });
      hit.addEventListener('mouseleave', () => { tip.hidden = true; cross.setAttribute('visibility', 'hidden'); dots.replaceChildren(); });
      svg.append(cross, dots, hit);
      legendDraw();
    }
    function legendDraw() {
      legend.replaceChildren(...model.series.map((sr, i) => {
        if (sr.legend === false) return null;
        const box = h('input', { type: 'checkbox', checked: !model.hidden.has(i), 'aria-label': `show ${sr.name}` });
        box.addEventListener('change', () => { box.checked ? model.hidden.delete(i) : model.hidden.add(i); draw(); });
        return h('label', {}, box, h('i', { class: `swatch s${(sr.colour ?? i) % PALETTE}${sr.dash ? ' dashed' : ''}` }), ` ${sr.name}`);
      }).filter(Boolean));
    }
    const api = {
      el: root, svg,
      set(series) { model.series = series; model.hidden = new Set([...model.hidden].filter((i) => i < series.length)); draw(); },
      setLog(on) { model.log = !!on; draw(); },
      csv() {
        const rows = ['series,x,y'];
        for (const sr of model.series) sr.x.forEach((x, i) => rows.push(`"${sr.name.replaceAll('"', '""')}",${x === null ? '' : x},${sr.y[i] === null ? '' : sr.y[i]}`));
        return rows.join('\n') + '\n';
      },
      exportSvg() { return serialize(svg); },
      async exportPng() { return rasterize(svg, W, H); },
    };
    draw();
    return api;
  }

  /* Inline the resolved colours and fonts so a saved image does not depend on this page's stylesheet. */
  function serialize(svg) {
    const clone = svg.cloneNode(true), sources = [svg, ...svg.querySelectorAll('*')], targets = [clone, ...clone.querySelectorAll('*')];
    sources.forEach((node, i) => {
      const style = getComputedStyle(node), target = targets[i];
      for (const prop of ['stroke', 'fill', 'stroke-width', 'stroke-dasharray', 'font-size', 'font-family', 'font-weight', 'opacity']) target.setAttribute(prop, style.getPropertyValue(prop));
      if (target.getAttribute('class') && /\bhit\b/.test(target.getAttribute('class'))) target.setAttribute('fill', 'none');
    });
    clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
    const background = s('rect', { x: 0, y: 0, width: '100%', height: '100%', fill: getComputedStyle(document.body).backgroundColor });
    clone.insertBefore(background, clone.firstChild);
    return new XMLSerializer().serializeToString(clone);
  }
  function rasterize(svg, width, height) {
    return new Promise((resolve, reject) => {
      const image = new Image();
      image.onload = () => {
        const canvas = h('canvas', { width: width * 2, height: height * 2 }), ctx = canvas.getContext('2d');
        ctx.scale(2, 2); ctx.drawImage(image, 0, 0, width, height); canvas.toBlob((blob) => resolve(blob), 'image/png');
      };
      image.onerror = () => reject(new Error('could not render the image'));
      image.src = `data:image/svg+xml;base64,${btoa(unescape(encodeURIComponent(serialize(svg))))}`;
    });
  }

  const VIRIDIS = ['#440154', '#482878', '#3e4989', '#31688e', '#26828e', '#1f9e89', '#35b779', '#6ece58', '#b5de2b', '#fde725'];
  function colour(t) {
    const x = Math.min(0.999999, Math.max(0, t)) * (VIRIDIS.length - 1), i = Math.floor(x), f = x - i;
    const a = VIRIDIS[i], b = VIRIDIS[i + 1], rgb = (c) => [1, 3, 5].map((k) => parseInt(c.slice(k, k + 2), 16));
    return rgb(a).map((v, k) => Math.round(v + (rgb(b)[k] - v) * f));
  }

  function heatmap(options) {
    const canvas = h('canvas', { class: 'heat', role: 'img', 'aria-label': options.title || 'map' });
    const bar = h('canvas', { class: 'colorbar', width: 22, height: 240, role: 'img', 'aria-label': 'colour scale' });
    const info = h('div', { class: 'muted', text: 'Hover the map for values.' });
    const root = h('div', { class: 'chart' }, h('div', { class: 'heatrow' }, canvas, bar), info);
    let current = null;
    const view = { fill: true, contours: 0, hidden: new Set() };
    const margin = { l: 78, r: 8, t: 8, b: 40 };

    function draw(payload) {
      current = payload;
      const xs = payload.x.filter((v) => v !== null), ys = payload.y.filter((v) => v !== null);
      const x0 = Math.min(...xs), x1 = Math.max(...xs), y0 = Math.min(...ys), y1 = Math.max(...ys);
      const zs = payload.z.flat().filter((v) => v !== null), z0 = Math.min(...zs), z1 = Math.max(...zs);
      const pw = 520, ph = Math.max(220, Math.min(620, Math.round(pw * ((y1 - y0) / (x1 - x0 || 1)) * (payload.x_units === payload.y_units ? 1 : 0.8) || 360)));
      canvas.width = pw + margin.l + margin.r; canvas.height = ph + margin.t + margin.b;
      const ctx = canvas.getContext('2d');
      ctx.clearRect(0, 0, canvas.width, canvas.height);
      const rows = payload.z.length, cols = payload.z[0].length, image = new ImageData(cols, rows);
      payload.z.forEach((row, r) => row.forEach((v, c) => {
        const at = ((rows - 1 - r) * cols + c) * 4, rgb = v === null ? [128, 128, 128] : colour((v - z0) / (z1 - z0 || 1));
        image.data[at] = rgb[0]; image.data[at + 1] = rgb[1]; image.data[at + 2] = rgb[2]; image.data[at + 3] = v === null ? 40 : 255;
      }));
      const tmp = h('canvas', { width: cols, height: rows });
      tmp.getContext('2d').putImageData(image, 0, 0);
      ctx.imageSmoothingEnabled = false;
      if (view.fill) ctx.drawImage(tmp, margin.l, margin.t, pw, ph);
      const gx = (x) => margin.l + ((x - x0) / (x1 - x0 || 1)) * pw, gy = (y) => margin.t + ph - ((y - y0) / (y1 - y0 || 1)) * ph;
      const ink = getComputedStyle(document.body).color;
      ctx.strokeStyle = ink; ctx.fillStyle = ink; ctx.font = '11px system-ui, sans-serif'; ctx.lineWidth = 1;
      ctx.strokeRect(margin.l, margin.t, pw, ph);
      for (const t of ticks(x0, x1, 6)) { ctx.fillText(tickLabel(t), gx(t) - 12, margin.t + ph + 14); ctx.beginPath(); ctx.moveTo(gx(t), margin.t + ph); ctx.lineTo(gx(t), margin.t + ph + 4); ctx.stroke(); }
      ctx.textAlign = 'right';
      for (const t of ticks(y0, y1, 6)) { ctx.fillText(tickLabel(t), margin.l - 8, gy(t) + 4); ctx.beginPath(); ctx.moveTo(margin.l - 4, gy(t)); ctx.lineTo(margin.l, gy(t)); ctx.stroke(); }
      ctx.textAlign = 'left';
      ctx.fillText(`${payload.x_label}${payload.x_units ? ` [${payload.x_units}]` : ''}`, margin.l + pw / 2 - 30, canvas.height - 6);
      ctx.save(); ctx.translate(12, margin.t + ph / 2 + 30); ctx.rotate(-Math.PI / 2); ctx.fillText(`${payload.y_label}${payload.y_units ? ` [${payload.y_units}]` : ''}`, 0, 0); ctx.restore();
      ctx.save(); ctx.beginPath(); ctx.rect(margin.l, margin.t, pw, ph); ctx.clip();
      if (view.contours > 0) {
        current.levels = contourLevels(payload, view.contours, z0, z1);
        ctx.lineJoin = 'round';
        current.levels.forEach((level, k) => {
          ctx.beginPath();
          for (const [a, b] of contourSegments(payload, level)) { ctx.moveTo(gx(a[0]), gy(a[1])); ctx.lineTo(gx(b[0]), gy(b[1])); }
          if (view.fill) { ctx.strokeStyle = 'rgba(255,255,255,0.9)'; ctx.lineWidth = 1.1; } else { ctx.strokeStyle = ink; ctx.lineWidth = 1.1; }
          ctx.stroke();
        });
      } else current.levels = [];
      for (const line of (payload.overlays || []).filter((o) => !view.hidden.has(o.name.split(' ')[0]))) {
        ctx.beginPath(); let pen = false;
        line.x.forEach((x, i) => { const y = line.y[i]; if (x === null || y === null) { pen = false; return; } pen ? ctx.lineTo(gx(x), gy(y)) : ctx.moveTo(gx(x), gy(y)); pen = true; });
        ctx.strokeStyle = '#ffffff'; ctx.lineWidth = 3; ctx.stroke(); ctx.strokeStyle = '#000000'; ctx.lineWidth = 1.2; ctx.stroke();
      }
      for (const point of payload.points || []) {
        ctx.strokeStyle = '#ff2a2a'; ctx.lineWidth = 2; const px = gx(point.x), py = gy(point.y);
        ctx.beginPath(); ctx.moveTo(px - 6, py); ctx.lineTo(px + 6, py); ctx.moveTo(px, py - 6); ctx.lineTo(px, py + 6); ctx.stroke();
      }
      ctx.restore();
      bar.width = 78; bar.height = ph + margin.t + margin.b;
      const bctx = bar.getContext('2d');
      for (let y = 0; y < ph; y += 1) { const rgb = colour(1 - y / ph); bctx.fillStyle = `rgb(${rgb[0]},${rgb[1]},${rgb[2]})`; bctx.fillRect(0, margin.t + y, 22, 1); }
      bctx.fillStyle = ink; bctx.font = '11px system-ui, sans-serif';
      for (const t of ticks(z0, z1, 5)) bctx.fillText(tickLabel(t), 26, margin.t + ph - ((t - z0) / (z1 - z0 || 1)) * ph + 4);
      current.geometry = { x0, x1, y0, y1, pw, ph, z0, z1 };
    }
    canvas.addEventListener('mousemove', (ev) => {
      if (!current) return;
      const g = current.geometry, box = canvas.getBoundingClientRect();
      const px = ((ev.clientX - box.left) / box.width) * canvas.width - margin.l, py = ((ev.clientY - box.top) / box.height) * canvas.height - margin.t;
      if (px < 0 || py < 0 || px > g.pw || py > g.ph) return;
      const x = g.x0 + (px / g.pw) * (g.x1 - g.x0), y = g.y1 - (py / g.ph) * (g.y1 - g.y0);
      const c = Math.min(current.x.length - 1, Math.round((px / g.pw) * (current.x.length - 1))), r = Math.min(current.y.length - 1, Math.round((1 - py / g.ph) * (current.y.length - 1)));
      const v = current.z[r] ? current.z[r][c] : null;
      info.textContent = `${current.x_label} ${fmt(x, 4)}, ${current.y_label} ${fmt(y, 4)}: ${current.z_label} = ${v === null || v === undefined ? 'no data' : fmt(v, 5)} ${current.z_units || ''}`;
    });
    return {
      el: root, canvas,
      set(payload) { draw(payload); },
      /* fill: colour the map; contours: number of contour lines (for a flux map: surfaces evenly spaced in flux from the axis to the boundary);
         hide(name, on): hide the overlays whose name starts with ``name`` (wall, boundary). */
      configure(options) { if (options.fill !== undefined) view.fill = !!options.fill; if (options.contours !== undefined) view.contours = Math.max(0, Math.min(60, Math.round(options.contours))); if (options.hide) options.hide.forEach(([name, on]) => (on ? view.hidden.add(name) : view.hidden.delete(name))); if (current) draw(current); },
      levels() { return current ? current.levels || [] : []; },
      csv() { const rows = ['x,y,z']; current.z.forEach((row, r) => row.forEach((v, c) => rows.push(`${current.x[c]},${current.y[r]},${v === null ? '' : v}`))); return rows.join('\n') + '\n'; },
      exportPng() { return new Promise((resolve) => canvas.toBlob((blob) => resolve(blob), 'image/png')); },
    };
  }

  /* Contour levels: for a flux map with known axis and boundary values, n surfaces evenly spaced from the axis to the boundary
     (the last is the boundary itself); otherwise n levels strictly inside the data range. */
  function contourLevels(payload, n, z0, z1) {
    if (payload.flux && Number.isFinite(payload.flux.axis) && Number.isFinite(payload.flux.boundary)) {
      return Array.from({ length: n }, (_, i) => payload.flux.axis + ((payload.flux.boundary - payload.flux.axis) * (i + 1)) / n);
    }
    return Array.from({ length: n }, (_, i) => z0 + ((z1 - z0) * (i + 1)) / (n + 1));
  }
  /* Marching squares on the payload grid (z[row][col], x by column, y by row): line segments of the level, linearly interpolated.
     A cell with a missing corner is skipped; an ambiguous saddle cell pairs its crossings in edge order. */
  function contourSegments(payload, level) {
    const { x, y, z } = payload, segments = [];
    for (let r = 0; r + 1 < z.length; r += 1) {
      for (let c = 0; c + 1 < x.length; c += 1) {
        const v = [z[r][c], z[r][c + 1], z[r + 1][c + 1], z[r + 1][c]];
        if (v.some((q) => q === null || q === undefined)) continue;
        const p = [[x[c], y[r]], [x[c + 1], y[r]], [x[c + 1], y[r + 1]], [x[c], y[r + 1]]], hits = [];
        for (let k = 0; k < 4; k += 1) {
          const a = v[k], b = v[(k + 1) % 4];
          if ((a < level) !== (b < level)) { const t = (level - a) / (b - a), u = p[k], w = p[(k + 1) % 4]; hits.push([u[0] + (w[0] - u[0]) * t, u[1] + (w[1] - u[1]) * t]); }
        }
        for (let k = 0; k + 1 < hits.length; k += 2) segments.push([hits[k], hits[k + 1]]);
      }
    }
    return segments;
  }
  function pointInPolygon(px, py, xy) {
    let inside = false;
    for (let i = 0, j = xy.length - 1; i < xy.length; j = i, i += 1) {
      const [xi, yi] = xy[i], [xj, yj] = xy[j];
      if ((yi > py) !== (yj > py) && px < ((xj - xi) * (py - yi)) / (yj - yi) + xi) inside = !inside;
    }
    return inside;
  }

  /* An equal-aspect picture in the R,Z plane: filled polygons coloured by value (the fast-ion zones), and outlines (the wall).
     Hovering a polygon names its zone and value; clicking calls options.onZone(zone). */
  function shapeMap(options) {
    const canvas = h('canvas', { class: 'heat shapes', role: 'img', 'aria-label': options.title || 'R,Z picture' });
    const bar = h('canvas', { class: 'colorbar', width: 78, height: 240, role: 'img', 'aria-label': 'colour scale' });
    const info = h('div', { class: 'muted', text: 'Hover for values.' });
    const root = h('div', { class: 'chart' }, h('div', { class: 'heatrow' }, canvas, bar), info);
    const margin = { l: 62, r: 8, t: 8, b: 40 };
    let current = null, selected = null;
    function draw(payload) {
      current = payload;
      const lines = [...(payload.overlays || []), ...(payload.curves || [])], polys = payload.polys || [];
      const xs = [], ys = [];
      [...lines.map((l) => [l.x, l.y]), ...polys.map((q) => [q.xy.map((p) => p[0]), q.xy.map((p) => p[1])])].forEach(([a, b]) => { a.forEach((v) => v !== null && xs.push(v)); b.forEach((v) => v !== null && ys.push(v)); });
      const ctx = canvas.getContext('2d');
      if (!xs.length) { canvas.width = 400; canvas.height = 80; ctx.fillText('Nothing to draw', 10, 40); return; }
      let x0 = Math.min(...xs), x1 = Math.max(...xs), y0 = Math.min(...ys), y1 = Math.max(...ys);
      const padx = (x1 - x0) * 0.03 || 0.1, pady = (y1 - y0) * 0.03 || 0.1;
      x0 -= padx; x1 += padx; y0 -= pady; y1 += pady;
      const scale = Math.min(520 / (x1 - x0), 640 / (y1 - y0)), pw = Math.round((x1 - x0) * scale), ph = Math.round((y1 - y0) * scale);
      canvas.width = pw + margin.l + margin.r; canvas.height = ph + margin.t + margin.b;
      const gx = (x) => margin.l + (x - x0) * scale, gy = (y) => margin.t + ph - (y - y0) * scale;
      const values = polys.map((q) => q.value), v0 = values.length ? Math.min(...values) : 0, v1 = values.length ? Math.max(...values) : 1;
      const ink = getComputedStyle(document.body).color;
      ctx.clearRect(0, 0, canvas.width, canvas.height);
      ctx.lineJoin = 'round';
      polys.forEach((q) => {
        ctx.beginPath(); q.xy.forEach((p, i) => (i ? ctx.lineTo(gx(p[0]), gy(p[1])) : ctx.moveTo(gx(p[0]), gy(p[1])))); ctx.closePath();
        const rgb = colour((q.value - v0) / (v1 - v0 || 1)); ctx.fillStyle = `rgb(${rgb[0]},${rgb[1]},${rgb[2]})`; ctx.fill();
        // Match the plotting script's edgecolors="face": a white outline here
        // produces bright seams between otherwise adjacent MC-zone polygons.
        ctx.strokeStyle = ctx.fillStyle; ctx.lineWidth = 0.8; ctx.stroke();
      });
      for (const line of lines) {
        ctx.beginPath(); let pen = false;
        line.x.forEach((x, i) => { const y = line.y[i]; if (x === null || y === null) { pen = false; return; } pen ? ctx.lineTo(gx(x), gy(y)) : ctx.moveTo(gx(x), gy(y)); pen = true; });
        if (polys.length) { ctx.strokeStyle = '#ffffff'; ctx.lineWidth = 3; ctx.stroke(); ctx.strokeStyle = '#000000'; ctx.lineWidth = 1.2; } else { ctx.strokeStyle = ink; ctx.lineWidth = 1.8; }
        ctx.stroke();
      }
      const chosen = polys.find((q) => q.zone === selected);
      if (chosen) { ctx.beginPath(); chosen.xy.forEach((p, i) => (i ? ctx.lineTo(gx(p[0]), gy(p[1])) : ctx.moveTo(gx(p[0]), gy(p[1])))); ctx.closePath(); ctx.strokeStyle = '#ff2a2a'; ctx.lineWidth = 2.5; ctx.stroke(); }
      ctx.strokeStyle = ink; ctx.fillStyle = ink; ctx.font = '11px system-ui, sans-serif'; ctx.lineWidth = 1;
      ctx.strokeRect(margin.l, margin.t, pw, ph);
      for (const t of ticks(x0, x1, 6)) { ctx.fillText(tickLabel(t), gx(t) - 12, margin.t + ph + 14); ctx.beginPath(); ctx.moveTo(gx(t), margin.t + ph); ctx.lineTo(gx(t), margin.t + ph + 4); ctx.stroke(); }
      ctx.textAlign = 'right';
      for (const t of ticks(y0, y1, 6)) { ctx.fillText(tickLabel(t), margin.l - 8, gy(t) + 4); ctx.beginPath(); ctx.moveTo(margin.l - 4, gy(t)); ctx.lineTo(margin.l, gy(t)); ctx.stroke(); }
      ctx.textAlign = 'left';
      ctx.fillText(`${payload.x_label} [${payload.x_units}]`, margin.l + pw / 2 - 20, canvas.height - 6);
      ctx.save(); ctx.translate(12, margin.t + ph / 2 + 20); ctx.rotate(-Math.PI / 2); ctx.fillText(`${payload.y_label} [${payload.y_units}]`, 0, 0); ctx.restore();
      bar.hidden = !polys.length;
      if (polys.length) {
        bar.width = 90; bar.height = ph + margin.t + margin.b;
        const bctx = bar.getContext('2d');
        for (let y = 0; y < ph; y += 1) { const rgb = colour(1 - y / ph); bctx.fillStyle = `rgb(${rgb[0]},${rgb[1]},${rgb[2]})`; bctx.fillRect(0, margin.t + y, 22, 1); }
        bctx.fillStyle = ink; bctx.font = '11px system-ui, sans-serif';
        for (const t of ticks(v0, v1, 5)) bctx.fillText(tickLabel(t), 26, margin.t + ph - ((t - v0) / (v1 - v0 || 1)) * ph + 4);
        bctx.fillText(payload.z_units || '', 26, margin.t - 0);
      }
      current.geometry = { x0, y1, scale, pw, ph };
    }
    function locate(ev) {
      const g = current.geometry, box = canvas.getBoundingClientRect();
      const px = ((ev.clientX - box.left) / box.width) * canvas.width - margin.l, py = ((ev.clientY - box.top) / box.height) * canvas.height - margin.t;
      if (px < 0 || py < 0 || px > g.pw || py > g.ph) return null;
      return { x: g.x0 + px / g.scale, y: g.y1 - py / g.scale };
    }
    canvas.addEventListener('mousemove', (ev) => {
      if (!current) return;
      const at = locate(ev);
      if (!at) return;
      const hit = (current.polys || []).find((q) => pointInPolygon(at.x, at.y, q.xy));
      info.textContent = `R ${fmt(at.x, 4)}, Z ${fmt(at.y, 4)}` + (hit ? `: zone ${hit.zone}, ${current.z_label} = ${fmt(hit.value, 5)} ${current.z_units || ''} (click to plot this zone)` : '');
    });
    canvas.addEventListener('click', (ev) => {
      if (!current || !options.onZone) return;
      const at = locate(ev), hit = at && (current.polys || []).find((q) => pointInPolygon(at.x, at.y, q.xy));
      if (hit) options.onZone(hit.zone);
    });
    return {
      el: root, canvas,
      set(payload, zone) { selected = zone === undefined ? null : zone; draw(payload); },
      csv() {
        const rows = ['kind,name_or_zone,x_or_value,y'];
        (current.curves || []).forEach((c) => c.x.forEach((x, i) => rows.push(`outline,"${c.name}",${x === null ? '' : x},${c.y[i] === null ? '' : c.y[i]}`)));
        (current.polys || []).forEach((q) => rows.push(`zone,${q.zone},${q.value},`));
        return rows.join('\n') + '\n';
      },
      exportPng() { return new Promise((resolve) => canvas.toBlob((blob) => resolve(blob), 'image/png')); },
    };
  }
  window.GUI.charts = { lineChart, heatmap, shapeMap, ticks, contourSegments, contourLevels, pointInPolygon };
})();
