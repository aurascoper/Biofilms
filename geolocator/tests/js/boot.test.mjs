// The boot order, executed. A health request that never answers must not hold the
// panels that do not depend on it; only the poll interval waits on the first response.
import test from 'node:test';
import assert from 'node:assert/strict';
import { boot } from '../../static/app/boot.js';

test('independent panels start before the first health poll resolves; only its interval waits', async () => {
  const started = [];
  const intervals = [];
  let resolveHealth;
  const health = new Promise((resolve) => { resolveHealth = resolve; });
  const done = boot({
    sites: { pollHealth: () => { started.push('health'); return health; } },
    bandsys: { load: () => started.push('bands') },
    latticeLayer: { update: () => started.push('lattice'), refreshInterval: 5000 },
    imageryLayer: { update: () => started.push('imagery') },
    loadProvenance: () => started.push('provenance'),
    setInterval: (fn, ms) => intervals.push(ms),
  });
  // Synchronously, before anything could have answered: every independent start is
  // recorded, the health poll has been issued, and only the lattice interval exists.
  assert.deepEqual(started, ['bands', 'lattice', 'provenance', 'imagery', 'health']);
  assert.deepEqual(intervals, [5000]);
  resolveHealth();
  await done;
  assert.deepEqual(intervals, [5000, 15000]);
});
