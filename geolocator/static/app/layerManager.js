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

  const rowNodes = new Map();

  /* Rows are built ONCE and then updated in place.
   *
   * The first version rebuilt panel.innerHTML on every refresh, and refresh is
   * called by the lattice layer every 5 seconds -- so eleven buttons were being
   * destroyed and recreated forever, throwing away hover and focus state and
   * dropping any click that landed mid-rebuild. Nothing about a status chip
   * changing requires new nodes. */
  function renderPanel() {
    if (!panel) return;
    for (const [id, rec] of layers) {
      if (rec.mod.showInTogglePanel === false) continue;
      const s = rec.mod.getStats?.() || {};
      const feed = FEED[feedState(id)];

      let node = rowNodes.get(id);
      if (!node) {
        const row = document.createElement('button');
        row.type = 'button';
        row.className = 'layer-row';
        row.innerHTML =
          `<span class="l-icon">${rec.mod.icon || '•'}</span>` +
          `<span class="l-name">${rec.mod.name}</span>` +
          `<span class="l-count"></span><span class="l-chip"></span>`;
        row.addEventListener('click', () => setEnabled(id, !layers.get(id).enabled));
        panel.appendChild(row);
        node = {
          row,
          count: row.querySelector('.l-count'),
          chip: row.querySelector('.l-chip'),
        };
        rowNodes.set(id, node);
      }

      const count = s.count == null ? '' : Number(s.count).toLocaleString();
      if (node.count.textContent !== count) node.count.textContent = count;
      if (node.chip.textContent !== feed.label) node.chip.textContent = feed.label;
      if (node.chip.style.color !== feed.color) {
        node.chip.style.color = feed.color;
        node.chip.style.borderColor = feed.color;
      }
      node.row.classList.toggle('on', rec.enabled);
      const tip = describe(rec.mod, s);
      if (node.row.title !== tip) node.row.title = tip;
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
