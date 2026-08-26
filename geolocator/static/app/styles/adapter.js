/*
 * adapter.js — run God's Eye View's Cesium post-process shaders under three.js.
 *
 * The vendored shader bodies in this directory are byte-verbatim from upstream so
 * they stay diffable when upstream fixes one. Everything Cesium-specific about
 * them is NAMING, not API, so the whole port is a handful of substitutions done
 * here at runtime rather than edits smeared across six files:
 *
 *   colorTextureDimensions  ->  tDiffuseSize   (must go FIRST: `colorTexture`
 *                                               is a prefix of it)
 *   colorTexture            ->  tDiffuse       (what three.js ShaderPass binds)
 *   v_textureCoordinates    ->  vUv            (what three.js ShaderPass provides)
 *
 * The shaders are GLSL 300 es (`in`, `texture()`, `out_FragColor`). three.js can
 * run that directly via THREE.GLSL3 — but under GLSL3 three.js declares no output
 * variable at all (it only injects `pc_fragColor` for GLSL1), and Cesium used to
 * supply `out_FragColor` implicitly. So we declare it.
 */

const VERTEX = /* glsl */ `
varying vec2 vUv;
void main() {
  vUv = uv;
  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
}
`;

function translate(src) {
  let out = src
    .replace(/\bcolorTextureDimensions\b/g, 'tDiffuseSize')
    .replace(/\bcolorTexture\b/g, 'tDiffuse')
    .replace(/\bv_textureCoordinates\b/g, 'vUv');
  if (/\bout_FragColor\b/.test(out) && !/out\s+vec4\s+out_FragColor/.test(out)) {
    out = 'out vec4 out_FragColor;\n' + out;
  }
  return out;
}

/**
 * @param {{name:string, uniforms:Object, fragmentShader:string}} shader upstream module
 * @returns {{name, uniforms, vertexShader, fragmentShader, controls}} three.js-ready
 */
export function adapt(shader) {
  const uniforms = {
    tDiffuse: { value: null },
    tDiffuseSize: { value: [1, 1] },
    intensity: { value: 1.0 },
    time: { value: 0 },
  };
  const controls = [];
  for (const [key, meta] of Object.entries(shader.uniforms || {})) {
    uniforms[key] = { value: meta.default };
    controls.push({ key, ...meta });
  }
  return {
    name: shader.name,
    uniforms,
    controls,
    vertexShader: VERTEX,
    fragmentShader: translate(shader.fragmentShader),
  };
}
