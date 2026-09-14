/*
 * layers.js — the six site layers and the three epistemic band layers.
 *
 * Rendering behaviour is carried over from the pre-fusion page unchanged, and
 * the comments that record earlier fixes are carried with it. The six site
 * layers deliberately share one renderer rather than owning their own meshes,
 * because size normalisation is cross-layer: a 4-degree grid cell and a single
 * mine are not comparable, so maxCap is computed per layer across the set that
 * is actually on screen.
 *
 * What is new is that each layer now reports a status, and the status comes
 * from the server's per-layer freshness resolution rather than being invented
 * here.
 */

import { BAND_STATUS } from './tokens.js';
import { latLonToVec3, raDecToVec3, R_EARTH, R_CELESTIAL } from './globe.js';

export const PLANT_LIMIT = 5000;

const SITE_LAYERS = [
  { id: 'power',     name: 'Power Plants',   icon: '⚡' },
  { id: 'worldgrid', name: 'World Grid',     icon: '▦' },
  { id: 'nuclear',   name: 'Nuclear Cycle',  icon: '☢' },
  { id: 'battery',   name: 'Battery Cycle',  icon: '▮' },
  { id: 'gridcoin',  name: 'Gridcoin/BOINC', icon: '◈' },
  { id: 'stars',     name: 'Star Systems',   icon: '✦' },
  // mw: false -- this layer carries no MW figure (capacity_mw is null for every county),
  // so the Min MW filter does not apply to it and the HUD must not report it in MW.
  { id: 'agri_overlay', name: 'Agri Overlay', icon: '🌾', mw: false },
];

/** Whether a layer's items carry a capacity in MW. Non-MW layers are neither filtered
 *  nor summed by the MW controls: an MW predicate over a null is a statement about a
 *  quantity the layer does not have. */
export const hasCapacity = (id) => SITE_LAYERS.find((l) => l.id === id)?.mw !== false;

// Marker size is normalised per layer, and the magnitude is layer-specific: the agri
// overlay sizes from the NDVI anomaly instead of being drawn at the minimum radius.
function magnitude(id, p) {
  if (id === 'agri_overlay') return Math.abs(p.extra?.ndvi?.z ?? 0);
  return p.capacity_mw || 0;
}

export function createSiteSystem({ markerRoot, manager, onRender }) {
  const groups = {};
  const features = {};
  const health = { freshness: {} };
  let markers = [];

  SITE_LAYERS.forEach(({ id }) => {
    groups[id] = new THREE.Group();
    markerRoot.add(groups[id]);
  });

  /** Drop a layer's cached features when the poll says its source changed, so the next
   *  render refetches. `ensureFetched` caches on first fetch; without this a layer
   *  fetched while its source was missing stayed empty, and a routine export update
   *  stayed invisible, until a page reload. Only ids the previous poll had already
   *  reported are compared, so the first poll does not refetch what was just loaded. */
  function invalidateChanged(before, after) {
    return Object.keys(features).filter((id) => {
      const a = before[id]; const b = after[id] || {};
      return a && (a.fingerprint !== b.fingerprint || a.status !== b.status);
    }).map((id) => { delete features[id]; return id; });
  }

  async function pollHealth() {
    const before = { ...(health.freshness || {}) };
    try {
      const r = await fetch('/api/health', { cache: 'no-store' });
      Object.assign(health, await r.json());
    } catch { /* the panel will show what it last knew */ }
    const after = health.freshness || {};
    const changed = invalidateChanged(before, after);
    manager.refresh();
    if (changed.some((id) => manager.isEnabled(id))) await refreshAndRender();
    // The HUD reads freshness only during a render. A poll that changed freshness but
    // invalidated nothing (the first poll, whose ids the previous poll had not reported)
    // must still re-render, or an enabled layer that is unavailable from page load reads
    // 0 MW until the next poll or a user action.
    else if (JSON.stringify(before) !== JSON.stringify(after)) render();
  }

  async function ensureFetched(id) {
    if (id in features) return;
    const r = await fetch(`/api/plants?layer=${id}&limit=${PLANT_LIMIT}`, { cache: 'no-store' });
    if (!r.ok) { features[id] = []; return; }
    features[id] = (await r.json()).features;
  }

  const enabled = () => SITE_LAYERS.map((l) => l.id).filter((id) => manager.isEnabled(id));

  function activeFeatures() {
    return enabled().flatMap((id) => (features[id] || []).map((f) => [id, f]));
  }

  /** Tears down and rebuilds every mesh. Geometries and materials are disposed
   *  on the way out -- the pre-fusion version dropped them on the floor, which
   *  leaked GPU memory every time the 34,936-plant layer was toggled. */
  function clearGroup(g) {
    while (g.children.length) {
      const m = g.children[0];
      g.remove(m);
      m.geometry?.dispose();
      m.material?.dispose();
    }
  }

  function render() {
    SITE_LAYERS.forEach(({ id }) => clearGroup(groups[id]));
    markers = [];

    const fuel = document.getElementById('fuel-filter')?.value || '';
    const minCap = parseFloat(document.getElementById('min-cap')?.value) || 0;
    const shown = activeFeatures().filter(([id, f]) =>
      (!fuel || f.properties.color_key === fuel) &&
      (!hasCapacity(id) || (f.properties.capacity_mw || 0) >= minCap));

    // Normalise per layer: a 4-degree grid cell and a single mine are not comparable.
    const maxCap = {};
    shown.forEach(([id, f]) => {
      maxCap[id] = Math.max(maxCap[id] || 1, magnitude(id, f.properties));
    });

    shown.forEach(([id, f]) => {
      const p = f.properties;
      const isStar = id === 'stars';
      const size = isStar ? 0.05 : (0.004 + 0.03 * Math.sqrt(magnitude(id, p) / maxCap[id]));
      const mesh = new THREE.Mesh(
        new THREE.SphereGeometry(size, 8, 8),
        new THREE.MeshBasicMaterial({ color: new THREE.Color(p.color || '#BDC3C7') }));
      if (isStar) {
        // stars ship coordinates [0,0]; their real position is a STRING in extra.
        const e = p.extra || {};
        mesh.position.copy(raDecToVec3(parseFloat(e.ra) || 0, parseFloat(e.dec) || 0, R_CELESTIAL));
      } else {
        mesh.position.copy(
          latLonToVec3(f.geometry.coordinates[1], f.geometry.coordinates[0], R_EARTH));
      }
      mesh.userData = { ...p, _layer: id };
      groups[id].add(mesh);
      markers.push(mesh);
    });

    onRender?.(shown, features, enabled(), health.freshness || {});
  }

  async function refreshAndRender() {
    const wanted = enabled().filter((id) => !(id in features));
    const loading = document.getElementById('loading');
    if (wanted.length && loading) loading.style.display = 'flex';
    await Promise.all(wanted.map(ensureFetched));
    if (loading) loading.style.display = 'none';
    populateFuelFilter();
    render();
  }

  function populateFuelFilter() {
    const sel = document.getElementById('fuel-filter');
    if (!sel) return;
    const keep = sel.value;
    const keys = [...new Set(activeFeatures().map(([, f]) => f.properties.color_key))].sort();
    sel.innerHTML = '<option value="">All</option>' + keys.map((k) => `<option>${k}</option>`).join('');
    if (keys.includes(keep)) sel.value = keep;
  }

  const modules = SITE_LAYERS.map((l) => ({
    ...l,
    layerClass: 'dataset',
    enabledByDefault: l.id === 'power',
    enable: refreshAndRender,
    disable: () => { clearGroup(groups[l.id]); render(); },
    getStats: () => health.freshness?.[l.id] || { status: 'loading' },
  }));

  return {
    modules,
    groups,
    markers: () => markers,
    render,
    refreshAndRender,
    pollHealth,
    pickTargets: () => enabled().flatMap((id) => groups[id].children),
  };
}

/* ── bands ──────────────────────────────────────────────────────────────── */

export function createBandSystem({ bandRoot, manager, onMeasured }) {
  const groups = {};
  let bands = [];
  let measuredFreshness = { status: 'loading' };

  Object.keys(BAND_STATUS).forEach((st) => {
    groups[st] = new THREE.Group();
    bandRoot.add(groups[st]);
  });

  function resolveEnd(e) {
    if (e.frame === 'geo') return latLonToVec3(e.lat, e.lon, R_EARTH);
    if (e.frame === 'celestial') return raDecToVec3(e.ra, e.dec, R_CELESTIAL);
    return null;   // fuel frame has no position, and must not be given one
  }

  async function load() {
    const r = await fetch('/api/links', { cache: 'no-store' });
    const doc = await r.json();
    bands = doc.features || [];
    measuredFreshness = doc.measured_freshness || { status: 'loading' };
    draw();
    onMeasured?.(measured(), measuredFreshness);
    manager.refresh();
  }

  const measured = () => bands.filter((b) => b.properties.status === 'measured');

  function draw() {
    Object.values(groups).forEach((g) => {
      while (g.children.length) {
        const c = g.children[0];
        g.remove(c);
        c.geometry?.dispose();
        c.material?.dispose();
      }
    });
    bands.forEach((b) => {
      const p = b.properties;
      const a = resolveEnd(p.a);
      const c = resolveEnd(p.b);
      if (!a || !c) return;                     // measured bands render in the HUD
      const spec = p.status === 'speculative';
      // Bulge past the far endpoint so an Earth-to-sky band reads as an escape
      // rather than a chord through the planet.
      const bulge = Math.max(a.length(), c.length()) * (spec ? 1.12 : 1.28);
      const mid = a.clone().add(c).multiplyScalar(0.5).normalize().multiplyScalar(bulge);
      const geo = new THREE.BufferGeometry()
        .setFromPoints(new THREE.QuadraticBezierCurve3(a, mid, c).getPoints(48));
      let line;
      if (spec) {
        // Dashes must be scaled to arc length; these arcs are ~0.3-2.5 world units,
        // so a dash size of 1 would render solid.
        line = new THREE.Line(geo, new THREE.LineDashedMaterial({
          color: p.color, transparent: true, opacity: 0.5, dashSize: 0.035, gapSize: 0.025 }));
        line.computeLineDistances();
      } else {
        line = new THREE.Line(geo, new THREE.LineBasicMaterial({
          color: p.color, transparent: true, opacity: 0.32 }));
      }
      line.userData = p;
      groups[p.status].add(line);
    });
    applyVisibility();
  }

  function applyVisibility() {
    Object.entries(groups).forEach(([st, g]) => { g.visible = manager.isEnabled(`band:${st}`); });
  }

  const modules = Object.entries(BAND_STATUS).map(([st, meta]) => ({
    id: `band:${st}`,
    name: meta.label,
    icon: st === 'measured' ? '∿' : (st === 'speculative' ? '┄' : '─'),
    layerClass: st === 'measured' ? 'live' : 'derived',
    source: meta.hint,
    enabledByDefault: true,
    enable: () => { applyVisibility(); onMeasured?.(measured(), measuredFreshness); },
    disable: () => { applyVisibility(); onMeasured?.(measured(), measuredFreshness); },
    // Only the measured band has a freshness authority of its own: it is only as
    // current as the newest bar behind its statistic.
    // The measured band's freshness comes from the bar cache, but its COUNT is
    // the number of bands -- not the 560,758 rows behind them.
    getStats: () => (st === 'measured'
      ? { ...measuredFreshness, count: countOf(st) }
      : { status: 'nominal', count: countOf(st) }),
  }));

  const countOf = (st) => bands.filter((b) => b.properties.status === st).length;

  return {
    modules, groups, load, applyVisibility,
    bands: () => bands,
    measured,
    pickTargets: () => Object.entries(groups)
      .filter(([st]) => manager.isEnabled(`band:${st}`))
      .flatMap(([, g]) => g.children),
  };
}
