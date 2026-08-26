/*
 * latticePanel.js — the trading lattice, docked in screen space.
 *
 * WHY THIS IS NOT A GLOBE LAYER
 * The geolocator already decided this question. `_endpoint()` in api.py tags
 * each band endpoint with a coordinate frame, and the `fuel` frame ships
 * geometry:null on purpose, because "a correlation between two price series has
 * no location" and anchoring it anywhere "invents geography the statistic does
 * not have". A trading symbol has no coordinates either. So it lives beside the
 * globe, not on it, and the join between them is the measured band.
 *
 * FIDELITY TO THE BOARD ON :4199
 * hash32 / ruleBits / caStep / buildAutomaton below are the board's own
 * functions. Its seed string is twelve fields, the last five of which come from
 * an open position -- and when a cell has no position those five evaluate to "".
 * This panel always joins five trailing empty fields, so every UNHELD cell seeds
 * bit-identically to the board.
 *
 * Held symbols are the deliberate divergence: position side/vol/leverage/entry/
 * liq are redacted at the boundary, so those rows seed differently and rule 30
 * (the board's open-short rule) never appears here. This page cannot even tell
 * you which rows those are, because has_position is redacted too. The panel
 * says so rather than implying a parity it does not have.
 */

import { T, FEED } from '../tokens.js';
import { fmtAge } from '../layerManager.js';

const COLS = 18;
const ROWS = 6;
const MAX_ROWS = 60;          // the panel is a readout, not the full 1,102-row atlas
const FUNCTOR = ['sense', 'context', 'infer', 'gate', 'size', 'commit'];

/* ── the board's automaton, verbatim in behaviour ─────────────────────────── */
function hash32(str) {
  str = String(str == null ? '' : str);
  let h = 2166136261 >>> 0;
  for (let i = 0; i < str.length; i++) {
    h ^= str.charCodeAt(i);
    h = Math.imul(h, 16777619);
  }
  return h >>> 0;
}

function ruleBits(rule) {
  const bits = [];
  for (let i = 0; i < 8; i++) bits.push((rule >> i) & 1);
  return bits;
}

function caStep(prev, rule) {
  const bits = ruleBits(rule);
  const next = [];
  for (let i = 0; i < prev.length; i++) {
    const l = prev[(i - 1 + prev.length) % prev.length];
    const c = prev[i];
    const r = prev[(i + 1) % prev.length];
    next.push(bits[(l << 2) | (c << 1) | r] ? 1 : 0);
  }
  return next;
}

/** Rule 30 is unreachable here: it is the board's open-SHORT rule, and side is
 *  redacted. stock 110 / etf 184 / crypto 90 are reproduced exactly. */
function ruleForCell(c) {
  if (c.class === 'stock') return 110;
  if (c.class === 'etf') return 184;
  return 90;
}

function buildAutomaton(c, cols = COLS, rows = ROWS) {
  const seedText = [
    c.symbol, c.quote, c.class, c.rank,
    c.signal ? c.signal.direction : '',
    c.signal ? c.signal.score : '',
    c.signal ? c.signal.z_4h : '',
    '', '', '', '', '',            // the five redacted position fields
  ].join('|');
  const seed = hash32(seedText);
  const rule = ruleForCell(c);
  const row = new Array(cols);
  for (let i = 0; i < cols; i++) {
    row[i] = ((seed >> (i % 24)) ^ (seed >>> ((i + 7) % 24)) ^ (i * 2654435761 >>> 0)) & 1;
  }
  if (!row.some((v) => v)) row[seed % cols] = 1;
  const grid = [row];
  for (let r = 1; r < rows; r++) grid.push(caStep(grid[r - 1], rule));
  return { rule, seed, grid };
}

/* ── panel ────────────────────────────────────────────────────────────────── */
export function createLatticePanel(root) {
  root.innerHTML = `
    <div class="lat-head">
      <span class="lat-eyebrow">MEXC · STREAM</span>
      <span class="lat-chip" id="lat-chip">—</span>
    </div>
    <div class="lat-mission" id="lat-mission"></div>
    <canvas id="lat-canvas"></canvas>
    <div class="lat-functor" id="lat-functor"></div>
    <div class="lat-foot" id="lat-foot"></div>`;

  const canvas = root.querySelector('#lat-canvas');
  const ctx = canvas.getContext('2d');
  let doc = null;
  let hoverIndex = -1;
  const listeners = [];

  function rowsFor(d) {
    return (d.cells || [])
      .slice()
      .sort((a, b) => {
        const sa = a.signal ? a.signal.score : -Infinity;
        const sb = b.signal ? b.signal.score : -Infinity;
        if (sa !== sb) return sb - sa;
        return String(a.symbol).localeCompare(String(b.symbol));
      })
      .slice(0, MAX_ROWS);
  }

  function draw() {
    if (!doc) return;
    const rows = rowsFor(doc);
    const dpr = Math.min(devicePixelRatio || 1, 2);
    const cell = 3;
    const rowH = ROWS * cell + 7;
    const w = root.clientWidth - 24;
    const h = rows.length * rowH;
    canvas.width = w * dpr;
    canvas.height = h * dpr;
    canvas.style.width = w + 'px';
    canvas.style.height = h + 'px';
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, h);
    ctx.font = `10px ${T.mono}`;
    ctx.textBaseline = 'top';

    rows.forEach((c, i) => {
      const y = i * rowH;
      const auto = buildAutomaton(c);
      const long = c.signal && c.signal.direction === 'LONG';
      const on = c.signal ? (long ? '#6be8ff' : '#e07b39') : 'rgba(147,161,173,0.45)';
      if (i === hoverIndex) {
        ctx.fillStyle = 'rgba(107,232,255,0.08)';
        ctx.fillRect(0, y - 1, w, rowH);
      }
      for (let r = 0; r < ROWS; r++) {
        for (let x = 0; x < COLS; x++) {
          if (!auto.grid[r][x]) continue;
          ctx.fillStyle = on;
          ctx.globalAlpha = 0.35 + 0.65 * (r / ROWS);
          ctx.fillRect(x * cell, y + r * cell, cell - 0.5, cell - 0.5);
        }
      }
      ctx.globalAlpha = 1;
      const tx = COLS * cell + 8;
      ctx.fillStyle = c.signal ? T.ink : T.inkSoft;
      ctx.fillText(String(c.symbol).replace(/_USDT$|_USDC$/, ''), tx, y + 2);
      if (c.signal) {
        ctx.fillStyle = on;
        ctx.textAlign = 'right';
        ctx.fillText(`${long ? 'L' : 'S'} ${c.signal.score.toFixed(2)}`, w, y + 2);
        ctx.textAlign = 'left';
      }
      ctx.fillStyle = 'rgba(147,161,173,0.35)';
      ctx.fillText(`r${auto.rule}`, tx, y + 13);
    });
  }

  canvas.addEventListener('mousemove', (e) => {
    const rect = canvas.getBoundingClientRect();
    const rowH = ROWS * 3 + 7;
    const i = Math.floor((e.clientY - rect.top) / rowH);
    if (i !== hoverIndex) { hoverIndex = i; draw(); }
    const rows = rowsFor(doc || { cells: [] });
    listeners.forEach((fn) => fn(rows[i] || null));
  });
  canvas.addEventListener('mouseleave', () => {
    hoverIndex = -1; draw(); listeners.forEach((fn) => fn(null));
  });

  function update(d) {
    doc = d;
    const feed = FEED[d.status] || FEED.loading;
    const chip = root.querySelector('#lat-chip');
    chip.textContent = feed.label;
    chip.style.color = feed.color;
    chip.style.borderColor = feed.color;

    const m = d.mission || {};
    root.querySelector('#lat-mission').innerHTML = d.status === 'unavailable'
      ? `<span class="k">board</span><span class="v" style="color:${FEED.unavailable.color}">unreachable</span>`
      : [
        ['phase', m.phase ?? '—'],
        ['entries', m.entries_allowed == null ? '—' : (m.entries_allowed ? 'allowed' : 'held')],
        ['age', fmtAge(d.age_s)],
      ].map(([k, v]) => `<span class="k">${k}</span><span class="v">${v}</span>`).join('');

    const gated = m.entries_allowed === false;
    root.querySelector('#lat-functor').innerHTML = FUNCTOR
      .map((s, i) => {
        const dim = gated && i >= 3;
        return `<span class="fn${dim ? ' dim' : ''}">${s}</span>`;
      })
      .join('<span class="fn-arr">→</span>');

    const withSignal = (d.cells || []).filter((c) => c.signal).length;
    root.querySelector('#lat-foot').innerHTML =
      `<div>${(d.cells || []).length.toLocaleString()} cells · ${withSignal} signals · top ${MAX_ROWS} shown</div>` +
      `<div class="caveat">Position fields are redacted at the boundary. Rows for held ` +
      `symbols seed differently from the board and rule 30 never appears here — and ` +
      `this page cannot tell you which rows those are.</div>`;
    draw();
  }

  return { update, redraw: draw, onHover: (fn) => listeners.push(fn) };
}

/* ── the data layer that feeds it ─────────────────────────────────────────── */
export function createLatticeLayer({ panel, manager }) {
  let stats = { status: 'loading' };

  async function update() {
    try {
      const r = await fetch('/api/lattice', { cache: 'no-store' });
      const d = await r.json();
      stats = { status: d.status, count: (d.cells || []).length,
                age_s: d.age_s, authority: d.authority,
                source_timestamp: d.generated, error: d.error };
      panel.update(d);
    } catch (err) {
      stats = { status: 'unavailable', error: String(err) };
      panel.update({ status: 'unavailable', cells: [], mission: {} });
    }
    manager.refresh();
  }

  return {
    id: 'lattice',
    name: 'Trading Lattice',
    icon: '▤',
    layerClass: 'live',
    source: 'lattice board :4199 (redacted)',
    enabledByDefault: true,
    refreshInterval: 5000,      // the board itself is rewritten roughly every 60s
    enable: () => { document.getElementById('lattice').classList.add('on'); return update(); },
    disable: () => { document.getElementById('lattice').classList.remove('on'); },
    update,
    getStats: () => stats,
  };
}
