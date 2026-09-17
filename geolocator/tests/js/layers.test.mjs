// The site cache across health polls, executed against a scripted fetch. THREE and the
// document are read by createSiteSystem only at call time, so stubs suffice; the real
// layer manager is DOM-bound, so a fake stands in for it.
import test from 'node:test';
import assert from 'node:assert/strict';

class Vector3 {
  constructor(x = 0, y = 0, z = 0) { this.x = x; this.y = y; this.z = z; }
  copy(v) { this.x = v.x; this.y = v.y; this.z = v.z; return this; }
  set(x, y, z) { this.x = x; this.y = y; this.z = z; return this; }
}
globalThis.THREE = {
  Group: class { constructor() { this.children = []; } add(c) { this.children.push(c); }
                 remove(c) { const i = this.children.indexOf(c); if (i >= 0) this.children.splice(i, 1); } },
  Mesh: class { constructor(g, m) { this.geometry = g; this.material = m; this.position = new Vector3(); } },
  SphereGeometry: class { dispose() {} },
  MeshBasicMaterial: class { dispose() {} },
  Color: class {},
  Vector3,
};
globalThis.document = { getElementById: () => null };
const { createSiteSystem } = await import('../../static/app/layers.js');

const deferred = () => { let resolve, reject; const p = new Promise((a, b) => { resolve = a; reject = b; }); return { p, resolve, reject }; };
const tick = () => new Promise((r) => setImmediate(r));
const feat = (name) => ({ type: 'Feature', properties: { name, color_key: 'x', capacity_mw: 1 }, geometry: { coordinates: [0, 0] } });
const healthOk = (freshness) => () => Promise.resolve({ ok: true, json: async () => ({ status: 'ok', freshness }) });
const healthDown = () => Promise.reject(new Error('health unreachable'));
const plantsOk = (names) => ({ ok: true, json: async () => ({ features: names.map(feat) }) });

function harness() {
  const calls = [];
  const health = [];    // responders for /api/health, consumed in order; none => never answers
  const plants = [];    // one deferred per /api/plants request, in order
  globalThis.fetch = (url) => {
    calls.push(url);
    if (url.startsWith('/api/health')) { const next = health.shift(); return next ? next() : new Promise(() => {}); }
    if (url.startsWith('/api/plants')) { const d = deferred(); plants.push(d); return d.p; }
    throw new Error(`unexpected fetch ${url}`);
  };
  const sites = createSiteSystem({
    markerRoot: new THREE.Group(),
    manager: { isEnabled: (id) => id === 'power', refresh() {} },
    onRender: () => {},
  });
  return { calls, health, plants, sites, plantsCalls: () => calls.filter((u) => u.startsWith('/api/plants')).length };
}

test('the site fetch stays behind the first health response', async () => {
  const h = harness();                       // no health responder: the request never answers
  const pending = h.sites.pollHealth();
  await tick(); await tick();
  assert.equal(h.plantsCalls(), 0);
  const state = await Promise.race([pending.then(() => 'settled'), tick().then(() => 'pending')]);
  assert.equal(state, 'pending');
});

test('repeated failed polls preserve the loaded layer: one plants fetch, no marker rebuild', async () => {
  const h = harness();
  h.health.push(healthDown);
  const first = h.sites.pollHealth();        // health fails; power is enabled and absent, so it is fetched
  await tick(); await tick();
  assert.equal(h.plantsCalls(), 1);
  h.plants[0].resolve(plantsOk(['A']));
  await first;
  const built = h.sites.markers();
  assert.equal(built.length, 1);
  for (let i = 0; i < 3; i++) { h.health.push(healthDown); await h.sites.pollHealth(); }
  // absent -> absent on every side of every failed poll: nothing observed to have changed.
  assert.equal(h.plantsCalls(), 1);
  assert.equal(h.sites.markers(), built);   // the same array: render() was not called
});

test('a real fingerprint or status change still invalidates and refetches', async () => {
  const h = harness();
  h.health.push(healthOk({ power: { status: 'ok', fingerprint: '1' } }));   // absent -> present: the initial change
  const p1 = h.sites.pollHealth(); await tick(); await tick();
  assert.equal(h.plantsCalls(), 1);
  h.plants[0].resolve(plantsOk(['A'])); await p1;
  const m1 = h.sites.markers();
  h.health.push(healthOk({ power: { status: 'ok', fingerprint: '2' } }));   // fingerprint moved
  const p2 = h.sites.pollHealth(); await tick(); await tick();
  assert.equal(h.plantsCalls(), 2);
  h.plants[1].resolve(plantsOk(['B'])); await p2;
  assert.notEqual(h.sites.markers(), m1);
  assert.equal(h.sites.markers()[0].userData.name, 'B');
  h.health.push(healthOk({ power: { status: 'unavailable', fingerprint: '2' } }));   // status moved
  const p3 = h.sites.pollHealth(); await tick(); await tick();
  assert.equal(h.plantsCalls(), 3);
  h.plants[2].resolve(plantsOk(['C'])); await p3;
  h.health.push(healthOk({ power: { status: 'unavailable', fingerprint: '2' } }));   // present -> present, unchanged
  await h.sites.pollHealth();
  assert.equal(h.plantsCalls(), 3);
});

test('a plants response that started before an invalidation is discarded; the retry lands', async () => {
  const h = harness();
  h.health.push(healthOk({ power: { status: 'ok', fingerprint: '1' } }));
  const p1 = h.sites.pollHealth(); await tick(); await tick();
  assert.equal(h.plantsCalls(), 1);                                   // in flight, unanswered
  h.health.push(healthOk({ power: { status: 'ok', fingerprint: '2' } }));
  await h.sites.pollHealth();                                          // invalidates while the fetch is in flight
  assert.equal(h.plantsCalls(), 1);                                   // no second request beside one in flight
  h.plants[0].resolve(plantsOk(['OLD'])); await p1;
  assert.equal(h.sites.markers().length, 0);                          // the pre-invalidation response was rejected
  h.health.push(healthOk({ power: { status: 'ok', fingerprint: '2' } }));
  const p3 = h.sites.pollHealth(); await tick(); await tick();
  assert.equal(h.plantsCalls(), 2);
  h.plants[1].resolve(plantsOk(['NEW'])); await p3;
  assert.equal(h.sites.markers()[0].userData.name, 'NEW');
});
