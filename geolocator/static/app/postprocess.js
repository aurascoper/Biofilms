/*
 * postprocess.js — the sensor-style chain.
 *
 * Keys 1-7 cycle NORMAL plus the six vendored looks. NORMAL is not a style: it
 * is the composer switched off entirely, which is what makes it the safe floor.
 *
 * Every style is compiled behind a guard. A shader that fails to compile marks
 * itself broken, logs once, and drops the view back to NORMAL -- one bad look
 * must never take the globe with it.
 */

import { adapt } from './styles/adapter.js';
import { thermalShader } from './styles/thermal.js';
import { nightVisionShader } from './styles/surveillance.js';
import { retroShader } from './styles/retro.js';
import { noirShader } from './styles/noir.js';
import { snowShader } from './styles/snow.js';
import { animeShader } from './styles/anime.js';

const STYLES = [
  { key: '1', id: 'normal', label: 'NORMAL', shader: null },
  { key: '2', id: 'surveillance', label: 'NVG', shader: nightVisionShader },
  { key: '3', id: 'thermal', label: 'FLIR', shader: thermalShader },
  { key: '4', id: 'retro', label: 'CRT', shader: retroShader },
  { key: '5', id: 'noir', label: 'NOIR', shader: noirShader },
  { key: '6', id: 'snow', label: 'SNOW', shader: snowShader },
  { key: '7', id: 'anime', label: 'ANIME', shader: animeShader },
];

export function createStyleChain({ renderer, scene, camera, onChange }) {
  const available = typeof THREE.EffectComposer === 'function';
  let composer = null;
  let pass = null;
  let active = 0;
  const broken = new Set();

  if (available) {
    composer = new THREE.EffectComposer(renderer);
    composer.addPass(new THREE.RenderPass(scene, camera));
  }

  function build(style) {
    const spec = adapt(style.shader);
    const p = new THREE.ShaderPass({
      uniforms: spec.uniforms,
      vertexShader: spec.vertexShader,
      fragmentShader: spec.fragmentShader,
    });
    // r147's ShaderPass builds the ShaderMaterial itself and does not forward
    // glslVersion, so it has to be set after construction.
    p.material.glslVersion = THREE.GLSL3;
    p.renderToScreen = true;
    return p;
  }

  function select(i) {
    const style = STYLES[i];
    if (!style || broken.has(style.id)) return;
    if (pass) { composer.removePass(pass); pass.dispose?.(); pass = null; }
    active = i;
    if (style.shader && composer) {
      try {
        pass = build(style);
        composer.addPass(pass);
      } catch (err) {
        broken.add(style.id);
        console.warn(`[style] ${style.id} disabled:`, err);
        return select(0);
      }
    }
    onChange?.(style, [...broken]);
  }

  function cycle(dir) {
    const usable = STYLES.filter((s) => !broken.has(s.id));
    const here = usable.findIndex((s) => s === STYLES[active]);
    const next = usable[(here + dir + usable.length) % usable.length];
    select(STYLES.indexOf(next));
  }

  window.addEventListener('keydown', (e) => {
    if (e.target.matches?.('input, select, textarea')) return;
    const i = STYLES.findIndex((s) => s.key === e.key);
    if (i >= 0) select(i);
    else if (e.key === '[') cycle(-1);
    else if (e.key === ']') cycle(1);
  });

  return {
    styles: STYLES,
    activeStyle: () => STYLES[active],
    select,
    resize(w, h) {
      composer?.setSize(w, h);
      if (pass) pass.uniforms.tDiffuseSize.value = [w, h];
    },
    /** @returns {boolean} true if the composer drew; false means render normally. */
    render(dtSeconds) {
      if (!pass || !composer) return false;
      pass.uniforms.time.value += dtSeconds;
      try {
        composer.render();
        return true;
      } catch (err) {
        broken.add(STYLES[active].id);
        console.warn(`[style] ${STYLES[active].id} failed at draw:`, err);
        select(0);
        return false;
      }
    },
  };
}
