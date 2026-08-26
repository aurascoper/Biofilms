# Third-party notices

## God's Eye View

<https://github.com/bilawalsidhu/gods-eye-view> — MIT License, Copyright (c) 2026 Bilawal Sidhu.

The geolocator vendors and adapts source from this project. MIT permits the reuse and requires the
copyright and permission notice be retained, which is what this file does.

**Vendored verbatim** (bodies unmodified so they stay diffable against upstream; Cesium-flavoured
GLSL names are rewritten for three.js at runtime by `geolocator/static/app/styles/adapter.js`):

| File here | Upstream |
|---|---|
| `geolocator/static/app/styles/thermal.js` | `src/styles/thermal.js` |
| `geolocator/static/app/styles/surveillance.js` | `src/styles/surveillance.js` |
| `geolocator/static/app/styles/retro.js` | `src/styles/retro.js` |
| `geolocator/static/app/styles/noir.js` | `src/styles/noir.js` |
| `geolocator/static/app/styles/snow.js` | `src/styles/snow.js` |
| `geolocator/static/app/styles/anime.js` | `src/styles/anime.js` |

**Adapted, not copied** — reimplemented against upstream's design, with its structure credited in
each file's header:

- `geolocator/static/app/tokens.js` — colour and type values from `src/overlays/worldOverlayTokens.js`.
  Type is the system mono stack rather than upstream's JetBrains Mono; this page has no build step.
- `geolocator/static/app/layerManager.js` — the layer contract and six-state feed model from
  `src/data/manager.js` (`layerFeedState`). Upstream's manager is 90 KB of lifecycle, QA seams and
  share-link serialisation; only the shape is reused.
- `geolocator/static/app/hud/detectionOverlay.js` — the dark self-adapting plate model described in
  upstream's `DETECTION_PLATE_BAND` / `CARD_PLATE_ALPHA` commentary.

**Deliberately NOT copied.** Upstream's MIT licence covers source only. Its bundled datasets under
`src/data/local_data/` carry their own terms — TeleGeography submarine cables are CC BY-NC-SA 3.0
(NonCommercial + ShareAlike), and the datacenter/dam extracts are ODbL 1.0. None of that data is
used here.

---

## three.js

<https://github.com/mrdoob/three.js> — MIT License, Copyright (c) 2010-2026 three.js authors.
Loaded from CDN at the pinned version 0.147.0 (`build/three.min.js`, `examples/js/controls/OrbitControls.js`,
`examples/js/postprocessing/{EffectComposer,RenderPass,ShaderPass}.js`, `examples/js/shaders/CopyShader.js`).
The pin is deliberate: 0.160 dropped the non-module `OrbitControls` path this page depends on.

---

## Data

| Dataset | Source | Licence |
|---|---|---|
| Global Power Plant Database v1.3.0 (34,936 plants) | World Resources Institute | CC BY 4.0 |
| Energy world grid (`energy_world.json`, 4° bins) | Derived from the above | CC BY 4.0 |
| Nuclear / battery fuel-cycle sites | Hand-authored from IAEA, WNA, DOE, USGS and company disclosures | see `source` column per row |
| BOINC project server locations | Public project pages | — |
| Star systems | Real exoplanet hosts; `energy_class` is fiction and is labelled SPECULATIVE in the UI | — |
| MEXC daily futures bars | MEXC public futures kline API | — |
| Trading lattice cells | Local lattice board on :4199, redacted at the boundary | not published |

If NASA Blue Marble Next Generation imagery is added later, it is public domain and NASA Earth
Observatory credit belongs in the provenance footer.
