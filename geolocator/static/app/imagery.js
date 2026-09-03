/*
 * imagery.js — NASA Blue Marble as a registered reference layer.
 *
 * Imagery is not a globe material here. It arrives with a provider, a dataset
 * vintage, a licence and a required credit, exactly like every other source, so
 * it registers like every other source and reports the same shape. Bolting it
 * onto the sphere would put the one dataset with a legal attribution obligation
 * outside the model that tracks obligations.
 *
 * The layer is OPTIONAL by contract. If /api/imagery/blue-marble 503s -- NASA
 * unreachable and no cache -- the globe keeps its graticule and the chip reads
 * UNAVAILABLE. It never reads STALE: a 2004 composite is not out of date, and
 * not having the bytes is a different failure from having old ones.
 */

import { FEED } from './tokens.js';

export function createImageryLayer({ globe, manager, onState }) {
  let texture = null;
  let stats = { status: 'loading' };
  let meta = null;

  const mat = () => globe.earthMesh.material;

  /** Disposal is the whole reason this is centralised: a replaced GPU texture
   *  that nobody disposes is a leak with no visible symptom. */
  function detach() {
    if (texture) { mat().map = null; texture.dispose(); texture = null; }
    mat().color.setHex(globe.BARE_EARTH);
    globe.graticuleMaterial.opacity = 0.5;
    mat().needsUpdate = true;
    onState?.(null);
  }

  function attach(tex) {
    const previous = texture;
    texture = tex;
    // A texture on a near-black Phong material renders as near-black, so the
    // base colour goes white while imagery is mounted and is restored on detach.
    mat().map = tex;
    mat().color.setHex(0xffffff);
    mat().needsUpdate = true;
    // The graticule earns its contrast over a dark sphere; over true-colour
    // imagery at full strength it is just a cage.
    globe.graticuleMaterial.opacity = 0.16;
    if (previous) previous.dispose();
    onState?.(meta);
  }

  async function update() {
    try {
      const r = await fetch('/api/imagery/blue-marble/meta', { cache: 'no-store' });
      meta = await r.json();
      stats = meta;
    } catch (err) {
      stats = { status: 'unavailable', error: String(err) };
      meta = null;
    }
    manager.refresh();
    if (!manager.isEnabled('blue_marble')) return;
    if (stats.status !== 'nominal') { detach(); return; }

    await new Promise((resolve) => {
      new THREE.TextureLoader().load(
        '/api/imagery/blue-marble',
        (tex) => {
          // Equirectangular 2:1 maps straight onto SphereGeometry's default UVs.
          tex.colorSpace = THREE.SRGBColorSpace ?? undefined;
          if (THREE.sRGBEncoding !== undefined && tex.encoding !== undefined) {
            tex.encoding = THREE.sRGBEncoding;   // r147 path
          }
          tex.minFilter = THREE.LinearMipmapLinearFilter;
          tex.anisotropy = Math.min(8, globe.renderer.capabilities.getMaxAnisotropy());
          attach(tex);
          resolve();
        },
        undefined,
        () => {
          stats = { ...stats, status: 'unavailable', error: 'texture decode failed' };
          detach();
          manager.refresh();
          resolve();
        },
      );
    });
  }

  return {
    id: 'blue_marble',
    name: 'Blue Marble',
    icon: '◉',
    layerClass: 'reference',
    kind: 'imagery',
    source: 'NASA GIBS · BlueMarble_NextGeneration',
    enabledByDefault: true,
    enable: update,
    disable: detach,
    update,
    getStats: () => stats,
    currentMeta: () => meta,
  };
}

/** The credit is a licence condition, so it is rendered from the served
 *  metadata rather than hardcoded, and it disappears with the imagery. */
export function imageryCredit(meta) {
  if (!meta || meta.status !== 'nominal') return '';
  const bits = [meta.attribution, `${meta.dataset} (${meta.vintage})`];
  if (meta.availability === 'cached') bits.push('cached');
  return bits.filter(Boolean).join(' · ');
}
