/*
 * tokens.js — one source of truth for colour and type.
 *
 * Values adapted from God's Eye View, src/overlays/worldOverlayTokens.js
 * (https://github.com/bilawalsidhu/gods-eye-view, MIT, (c) 2026 Bilawal Sidhu).
 * See THIRD_PARTY_NOTICES.md.
 *
 * Two deliberate departures from upstream:
 *
 *  - Type is the SYSTEM mono stack, not JetBrains Mono. Upstream vendors its
 *    fonts; this page has no build step and already dies without two CDN
 *    scripts, and a webfont would be a third. The system stack also matches the
 *    lattice board on :4199, which is the other half of this fused surface.
 *
 *  - Tokens are exported as JS AND written to :root as custom properties from
 *    the same object. The canvas draws from the object, the panels style from
 *    the properties, and they cannot drift apart.
 */

export const T = {
  bg: '#040c10',
  panel: 'rgba(4, 12, 16, 0.82)',
  panelSelected: 'rgba(5, 18, 24, 0.94)',
  line: 'rgba(190, 232, 242, 0.18)',
  lineSelected: 'rgba(107, 232, 255, 0.72)',
  ink: 'rgba(232, 240, 244, 0.96)',
  inkSoft: 'rgba(147, 161, 173, 0.92)',
  leader: 'rgba(147, 213, 228, 0.58)',
  accent: '#6be8ff',
  radius: '4px',
  mono: 'ui-monospace, "SF Mono", SFMono-Regular, Menlo, Consolas, monospace',
};

/* The six-state feed model, ported from upstream's layerFeedState. FALLBACK is
 * amber rather than red on purpose: "no source timestamp exists" is a fact about
 * the dataset, not a fault to alarm about. */
export const FEED = {
  nominal:     { label: 'ON',          color: '#6be8ff' },
  loading:     { label: 'LOADING',     color: '#9aa7c0' },
  degraded:    { label: 'DEGRADED',    color: '#f4d03f' },
  stale:       { label: 'STALE',       color: '#e07b39' },
  fallback:    { label: 'FALLBACK',    color: '#c9a227' },
  unavailable: { label: 'UNAVAILABLE', color: '#e05a50' },
  off:         { label: 'OFF',         color: 'rgba(147,161,173,0.55)' },
};

export const BAND_STATUS = {
  structural:  { label: 'Supply chain', hint: 'stage adjacency in the data' },
  measured:    { label: 'Measured',     hint: 'carries a statistic' },
  speculative: { label: 'Speculative',  hint: 'generated, not observed' },
};

export function installTokens() {
  const r = document.documentElement.style;
  r.setProperty('--bg', T.bg);
  r.setProperty('--panel', T.panel);
  r.setProperty('--panel-sel', T.panelSelected);
  r.setProperty('--line', T.line);
  r.setProperty('--line-sel', T.lineSelected);
  r.setProperty('--ink', T.ink);
  r.setProperty('--ink-soft', T.inkSoft);
  r.setProperty('--accent', T.accent);
  r.setProperty('--radius', T.radius);
  r.setProperty('--mono', T.mono);
  for (const [k, v] of Object.entries(FEED)) r.setProperty(`--feed-${k}`, v.color);
}
