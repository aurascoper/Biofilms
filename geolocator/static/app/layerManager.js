/*
 * layerManager.js — a layer registry with a freshness model.
 *
 * The contract is reimplemented from God's Eye View's DataLayerManager rather
 * than copied: upstream's manager.js is 90 KB of lifecycle, QA seams and share-
 * link serialisation, and this page needs the shape, not the machinery.
 *
 *   { id, name, icon, layerClass, init, enable, disable, update, getStats }
 *
 * The part worth having is the six-state feed model, because the pre-fusion page
 * had NO freshness signal of any kind -- no timestamps, no last-updated, no
 * error states -- and its one data-integrity indicator was a truncation chip.
 *
 * Freshness is never computed here. Each layer reports what the SERVER resolved,
 * and the server resolves it from a timestamp inside the data (see
 * geolocator/freshness.py). A layer that has no such timestamp reports
 * `fallback`, and the chip says so rather than dressing an mtime up as one.
 */

import { FEED } from './tokens.js';

export function createLayerManager() {
  const layers = new Map();
  let sealed = false;
  let panel = null;

  function register(mod) {
    if (sealed) throw new Error('registrations are finalized');
    if (!mod?.id) throw new Error('layer needs a stable id');
    if (layers.has(mod.id)) throw new Error(`duplicate layer id: ${mod.id}`);
    layers.set(mod.id, { mod, enabled: !!mod.enabledByDefault, timer: null });
  }

  function finalize() { sealed = true; }

  async function setEnabled(id, on) {
    const rec = layers.get(id);
    if (!rec || rec.enabled === on) return;
    rec.enabled = on;
    try {
      await (on ? rec.mod.enable?.() : rec.mod.disable?.());
    } catch (err) {
      console.warn(`[layer] ${id} ${on ? 'enable' : 'disable'} failed:`, err);
    }
    if (on && rec.mod.refreshInterval) {
      rec.timer = setInterval(() => rec.mod.update?.(), rec.mod.refreshInterval);
    } else if (rec.timer) {
      clearInterval(rec.timer);
      rec.timer = null;
    }
    renderPanel();
  }

  const isEnabled = (id) => !!layers.get(id)?.enabled;
  const enabledIds = () => [...layers].filter(([, r]) => r.enabled).map(([id]) => id);

  /** What the chip shows. `off` beats every other state: a layer you turned off
   *  is not stale, it is off, and colouring it amber would be a false alarm. */
  function feedState(id) {
    const rec = layers.get(id);
    if (!rec) return 'unavailable';
    if (!rec.enabled) return 'off';
    const s = rec.mod.getStats?.() || {};
    return FEED[s.status] ? s.status : 'loading';
  }

  function buildTogglePanel(container) {
    panel = container;
    renderPanel();
  }

  function renderPanel() {
    if (!panel) return;
    panel.innerHTML = '';
    for (const [id, rec] of layers) {
      if (rec.mod.showInTogglePanel === false) continue;
      const s = rec.mod.getStats?.() || {};
      const state = feedState(id);
      const feed = FEED[state];

      const row = document.createElement('button');
      row.type = 'button';
      row.className = 'layer-row' + (rec.enabled ? ' on' : '');
      row.title = describe(rec.mod, s);
      row.innerHTML =
        `<span class="l-icon">${rec.mod.icon || '•'}</span>` +
        `<span class="l-name">${rec.mod.name}</span>` +
        `<span class="l-count">${s.count == null ? '' : Number(s.count).toLocaleString()}</span>` +
        `<span class="l-chip" style="color:${feed.color};border-color:${feed.color}">${feed.label}</span>`;
      row.addEventListener('click', () => setEnabled(id, !rec.enabled));
      panel.appendChild(row);
    }
  }

  /** The tooltip is where the honesty lives: which timestamp was believed, and
   *  when there wasn't one, that there wasn't one. */
  function describe(mod, s) {
    const out = [`${mod.name} — ${mod.layerClass || 'layer'}`];
    if (mod.source) out.push(`source: ${mod.source}`);
    if (s.authority === 'filesystem-mtime') {
      out.push('no source timestamp exists for this dataset');
      out.push(`filesystem mtime: ${s.source_timestamp || '—'} (not a source timestamp)`);
    } else if (s.authority === 'dataset-vintage') {
      out.push(`vintage: ${s.vintage || s.source_timestamp || '—'}`);
      out.push('a versioned dataset has a vintage, not an age');
    } else if (s.source_timestamp) {
      out.push(`generated: ${s.source_timestamp}`);
      if (s.age_s != null) out.push(`age: ${fmtAge(s.age_s)}`);
    }
    if (s.error) out.push(`error: ${s.error}`);
    return out.join('\n');
  }

  return {
    register, finalize, setEnabled, isEnabled, enabledIds, feedState,
    buildTogglePanel, refresh: renderPanel,
    get: (id) => layers.get(id)?.mod,
    ids: () => [...layers.keys()],
  };
}

export function fmtAge(s) {
  if (s == null) return '—';
  if (s < 90) return `${Math.round(s)}s`;
  if (s < 5400) return `${Math.round(s / 60)}m`;
  if (s < 172800) return `${(s / 3600).toFixed(1)}h`;
  return `${(s / 86400).toFixed(1)}d`;
}
