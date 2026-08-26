/*
 * main.js — wiring only. Every behaviour lives in the module that owns it.
 */

import { installTokens, BAND_STATUS, FEED } from './tokens.js';
import { createGlobe } from './globe.js';
import { createImageryLayer, imageryCredit } from './imagery.js';
import { createLayerManager, fmtAge } from './layerManager.js';
import { createSiteSystem, createBandSystem, PLANT_LIMIT } from './layers.js';
import { createLatticePanel, createLatticeLayer } from './hud/latticePanel.js';
import { createDetectionOverlay } from './hud/detectionOverlay.js';
import { createStyleChain } from './postprocess.js';

installTokens();

const $ = (id) => document.getElementById(id);
const globe = createGlobe($('globe'));
const manager = createLayerManager();

/* ── base imagery ─────────────────────────────────────────────────────────── */
/* Registered first so it heads the panel as the base layer, and so a provider
 * outage is visible before anything drawn on top of it. */
let imageryMeta = null;
const imageryLayer = createImageryLayer({
  globe, manager,
  onState: (meta) => { imageryMeta = meta; renderProvenance(); },
});
manager.register(imageryLayer);

/* ── site layers ──────────────────────────────────────────────────────────── */
const sites = createSiteSystem({
  markerRoot: globe.markerRoot,
  manager,
  onRender: (shown, features, enabled) => { updateStats(shown, features, enabled); renderLegend(shown); },
});
sites.modules.forEach((m) => manager.register(m));

/* ── epistemic bands ──────────────────────────────────────────────────────── */
const bandsys = createBandSystem({
  bandRoot: globe.bandRoot,
  manager,
  onMeasured: renderMeasured,
});
bandsys.modules.forEach((m) => manager.register(m));

/* ── the trading lattice ──────────────────────────────────────────────────── */
const latticePanel = createLatticePanel($('lattice-body'));
const latticeLayer = createLatticeLayer({ panel: latticePanel, manager });
manager.register(latticeLayer);
manager.finalize();
manager.buildTogglePanel($('layer-panel'));

/* ── HUD readouts ─────────────────────────────────────────────────────────── */
function updateStats(shown, features, enabled) {
  const cap = shown.reduce((a, [, f]) => a + (f.properties.capacity_mw || 0), 0);
  $('stat-count').textContent = shown.length.toLocaleString();
  $('stat-cap').textContent = Math.round(cap).toLocaleString();
  // The fetch limit used to truncate power from 34,936 to 5,000 with no notice.
  const trunc = enabled.filter((id) => (features[id] || []).length >= PLANT_LIMIT);
  const el = $('stat-trunc');
  el.style.display = trunc.length ? 'block' : 'none';
  if (trunc.length) {
    el.textContent = `Showing first ${PLANT_LIMIT.toLocaleString()} of each: ${trunc.join(', ')}`;
  }
}

function renderLegend(shown) {
  const by = {};
  shown.forEach(([, f]) => {
    const k = f.properties.color_key;
    (by[k] = by[k] || { n: 0, color: f.properties.color }).n++;
  });
  $('legend').innerHTML = Object.entries(by)
    .sort((a, b) => b[1].n - a[1].n)
    .map(([k, v]) =>
      `<div class="row"><span class="swatch" style="background:${v.color}"></span>${k} <b>${v.n.toLocaleString()}</b></div>`)
    .join('');
}

function renderMeasured(ms, freshness) {
  const el = $('measured');
  if (!ms.length || !manager.isEnabled('band:measured')) { el.style.display = 'none'; return; }
  el.style.display = 'block';
  const feed = FEED[freshness?.status] || FEED.loading;
  el.innerHTML =
    `<div class="stat" style="margin-bottom:4px"><b>Measured relationships</b>` +
    `<span class="lat-chip" style="color:${feed.color};border-color:${feed.color};margin-left:6px">${feed.label}</span></div>` +
    ms.map((b) => {
      const p = b.properties, s = p.stats;
      return `<div class="row" style="display:block;margin-bottom:6px">
        <span class="swatch" style="background:${p.color};display:inline-block"></span>
        <b>${p.a.name} and ${p.b.name}</b><br>
        r = ${s.r_partial} (raw ${s.r_raw}) · n = ${s.n} · p = ${s.p}<br>
        ${s.basis}, control ${s.control}</div>`;
    }).join('') +
    `<div class="row" style="display:block">No position on the globe: this is a relationship ` +
    `between price series, not places.</div>` +
    (freshness?.source_timestamp
      ? `<div class="row" style="display:block">Newest bar ${freshness.source_timestamp.slice(0, 10)} · ${fmtAge(freshness.age_s)} old</div>`
      : '');
}

/* ── hover ────────────────────────────────────────────────────────────────── */
const mouse = new THREE.Vector2();
globe.renderer.domElement.addEventListener('mousemove', (e) => {
  const rect = globe.renderer.domElement.getBoundingClientRect();
  mouse.x = ((e.clientX - rect.left) / rect.width) * 2 - 1;
  mouse.y = -((e.clientY - rect.top) / rect.height) * 2 + 1;
  globe.raycaster.setFromCamera(mouse, globe.camera);
  const hits = globe.raycaster.intersectObjects([...sites.pickTargets(), ...bandsys.pickTargets()]);
  const tip = $('tooltip');
  if (!hits.length) { tip.style.display = 'none'; return; }
  const d = hits[0].object.userData;
  tip.style.display = 'block';
  tip.style.left = (e.clientX + 14) + 'px';
  tip.style.top = (e.clientY + 14) + 'px';
  tip.innerHTML = d.band ? bandTip(d) : siteTip(d);
});

function bandTip(p) {
  const meta = BAND_STATUS[p.status];
  let html = `<div class="t-name">${p.a.name} to ${p.b.name}</div>
    <div class="t-row">${p.label}</div>
    <div class="t-row"><b>${meta.label}</b> · ${meta.hint}</div>`;
  if (p.status === 'structural')
    html += '<div class="t-row">Shows which stages connect, not that these two sites trade.</div>';
  if (p.status === 'speculative')
    html += '<div class="t-row">Generated for this layer. Not observed.</div>';
  return html;
}

function siteTip(p) {
  let html = `<div class="t-name">${p.name}</div>
    <div class="t-row">${p.color_key}${p.country ? ' · ' + p.country : ''}</div>`;
  if (p.capacity_mw) html += `<div class="t-row">${p.capacity_mw.toLocaleString()} MW</div>`;
  if (p.extra && p.extra.distance_ly)
    html += `<div class="t-row">${p.extra.distance_ly} ly · ${p.extra.star_type || ''}</div>`;
  if (p.note) html += `<div class="t-row">${p.note}</div>`;
  return html;
}

/* ── sensor styles + detection ────────────────────────────────────────────── */
const detection = createDetectionOverlay($('detection-root'), {
  camera: globe.camera,
  targets: () => sites.pickTargets(),
});

const styles = createStyleChain({
  renderer: globe.renderer,
  scene: globe.scene,
  camera: globe.camera,
  onChange: (style, broken) => {
    $('style-name').textContent = style.label;
    $('style-keys').innerHTML = styles.styles
      .map((s) => `<span class="sk${s === style ? ' on' : ''}${broken.includes(s.id) ? ' bad' : ''}">${s.key}</span>`)
      .join('');
  },
});
styles.select(0);
styles.resize(innerWidth, innerHeight);

window.addEventListener('keydown', (e) => {
  if (e.target.matches?.('input, select, textarea')) return;
  if (e.key === 'h' || e.key === 'H') $('hud').classList.toggle('hidden');
  if (e.key === 'd' || e.key === 'D') $('det-state').textContent = detection.toggle() ? 'ON' : 'OFF';
});

window.addEventListener('resize', () => {
  globe.camera.aspect = innerWidth / innerHeight;
  globe.camera.updateProjectionMatrix();
  globe.renderer.setSize(innerWidth, innerHeight);
  styles.resize(innerWidth, innerHeight);
  latticePanel.redraw();
});

/* ── loop ─────────────────────────────────────────────────────────────────── */
let last = performance.now();
(function animate(now) {
  requestAnimationFrame(animate);
  const dt = ((now || performance.now()) - last) / 1000;
  last = now || last;
  globe.controls.update();
  if (!styles.render(dt)) globe.renderer.render(globe.scene, globe.camera);
})();

/* ── boot ─────────────────────────────────────────────────────────────────── */
$('fuel-filter').addEventListener('change', sites.render);
$('min-cap').addEventListener('change', sites.render);

sites.refreshAndRender();
sites.pollHealth();
setInterval(sites.pollHealth, 15000);
bandsys.load();
latticeLayer.update();
setInterval(latticeLayer.update, latticeLayer.refreshInterval);

let powerProvenance = '';

/* Attribution is a licence condition, not decoration: it is rendered from the
 * metadata the server actually served, and it disappears when imagery does. */
function renderProvenance() {
  $('provenance').textContent = [imageryCredit(imageryMeta), powerProvenance]
    .filter(Boolean).join('   ·   ');
}

fetch('/api/layers', { cache: 'no-store' })
  .then((r) => r.json())
  .then((d) => {
    const power = d.layers.find((l) => l.id === 'power') || {};
    powerProvenance =
      `WRI Global Power Plant Database v1.3.0 · CC BY 4.0 · retrieved ${(power.vintage || '').slice(0, 10) || '—'} · ` +
      `${(power.count || 0).toLocaleString()} plants`;
    renderProvenance();
  })
  .catch(() => {});

imageryLayer.update();
