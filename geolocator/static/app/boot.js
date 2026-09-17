/*
 * boot.js — the boot order, as a function so `node --test` can execute it.
 *
 * Everything that does not depend on site freshness starts at once. The site fetch runs
 * INSIDE the first health poll, after that poll has stored the freshness baseline, so
 * the ordering the per-layer generation guard relies on is kept; only the poll interval
 * is registered behind the first response. A top-level `await` on the health request
 * used to hold bands, lattice, attribution and imagery behind the imagery provider's
 * 20-second timeout, or indefinitely when the request hung.
 */
export function boot({ sites, bandsys, latticeLayer, imageryLayer, loadProvenance,
                       setInterval = globalThis.setInterval }) {
  bandsys.load();
  latticeLayer.update();
  setInterval(latticeLayer.update, latticeLayer.refreshInterval);
  loadProvenance();
  imageryLayer.update();
  return sites.pollHealth().then(() => setInterval(sites.pollHealth, 15000));
}
