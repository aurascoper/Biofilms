/*
 * detectionOverlay.js — screen-space detection boxes, bounded by construction.
 *
 * The scene is mesh-per-feature with no instancing and no LOD, and the power
 * layer alone is 34,936 features (5,000 after the fetch limit). An overlay that
 * projected every marker every frame would be the thing that finally killed the
 * render loop, so the bounds are structural rather than an optimisation to add
 * later:
 *
 *   - a fixed POOL of DOM nodes, reused; never rebuilt
 *   - at most MAX_BOXES targets, chosen by camera distance
 *   - frustum-culled before sorting, so off-screen markers cost nothing
 *   - driven by a 10 Hz timer, not by requestAnimationFrame
 *
 * Plate model from God's Eye View: a DARK plate at ~0.5-0.62 alpha painted with
 * normal blending is self-adapting. Over the dark globe it nearly vanishes; over
 * a bright cluster of markers it does the full darkening job. That is why the
 * boxes do not need to know what is behind them.
 */

import { T } from '../tokens.js';

const MAX_BOXES = 100;
const HZ = 10;

export function createDetectionOverlay(root, { camera, targets }) {
  const pool = [];
  let on = false;
  let timer = null;

  function node(i) {
    if (!pool[i]) {
      const el = document.createElement('div');
      el.className = 'det-box';
      el.innerHTML = '<span class="det-id"></span>';
      root.appendChild(el);
      pool[i] = el;
    }
    return pool[i];
  }

  const v = new THREE.Vector3();
  const frustum = new THREE.Frustum();
  const mat = new THREE.Matrix4();

  function tick() {
    if (!on) return;
    camera.updateMatrixWorld();
    mat.multiplyMatrices(camera.projectionMatrix, camera.matrixWorldInverse);
    frustum.setFromProjectionMatrix(mat);

    const w = innerWidth, h = innerHeight;
    const cand = [];
    for (const m of targets()) {
      m.getWorldPosition(v);
      if (!frustum.containsPoint(v)) continue;         // cull before sorting
      cand.push([v.distanceTo(camera.position), m]);
      if (cand.length > MAX_BOXES * 8) break;          // hard ceiling on the scan
    }
    cand.sort((a, b) => a[0] - b[0]);
    const shown = cand.slice(0, MAX_BOXES);

    shown.forEach(([dist, m], i) => {
      m.getWorldPosition(v).project(camera);
      const x = (v.x * 0.5 + 0.5) * w;
      const y = (-v.y * 0.5 + 0.5) * h;
      const size = Math.max(10, Math.min(46, 900 / Math.max(dist, 0.3)));
      const el = node(i);
      const d = m.userData || {};
      el.style.cssText =
        `display:block;left:${(x - size / 2).toFixed(1)}px;top:${(y - size / 2).toFixed(1)}px;` +
        `width:${size.toFixed(0)}px;height:${size.toFixed(0)}px;border-color:${d.color || T.accent}`;
      const id = el.firstChild;
      const label = `${(d.color_key || d._layer || '').toString().slice(0, 14)}`;
      if (id.textContent !== label) id.textContent = label;
    });
    for (let i = shown.length; i < pool.length; i++) pool[i].style.display = 'none';
  }

  return {
    enabled: () => on,
    toggle(next) {
      on = next == null ? !on : next;
      if (on && !timer) timer = setInterval(tick, 1000 / HZ);
      if (!on) {
        clearInterval(timer); timer = null;
        pool.forEach((el) => { el.style.display = 'none'; });
      }
      return on;
    },
  };
}
