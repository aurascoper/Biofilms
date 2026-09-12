# `biofilms-4d-viewer.html` — external dependency record
Recorded 2026-09-11 on macOS 26.6.2 / M4. Offline operation is a **separate check, not yet run**.
- File: `~/Downloads/biofilms-4d-viewer.html`
- Size: 974,077 bytes
- SHA-256: `df4dae79391e89ea57e2703cfe2408bd351dbce9a2fd394411b70f4d00a8464c`
## Structure
Outer document is a wrapper whose entire body is one `<iframe sandbox="allow-scripts" srcdoc="...">`. The real page is HTML-entity-escaped inside that attribute, so naive tag-counting on the outer file reports zero scripts. Decoded inner document is 952,471 chars with 7 `<script>` blocks and a single base64 payload of ~864 KB carrying `<u4`/`<f8`/`u1` typed arrays plus a `signal_offsets` ragged index (sparse per-frame voxel storage).
## Subresources actually requested (require network)
- `https://unpkg.com/@floating-ui/core@1.7.3/dist/floating-ui.core.umd.min.js`
- `https://unpkg.com/@floating-ui/dom@1.7.4/dist/floating-ui.dom.umd.min.js`
- `https://unpkg.com/lucide@1.17.0/dist/umd/lucide.js`

All three are UI chrome — positioning and icons — not the renderer, which is inline. Expected offline behaviour is therefore **degraded, not broken**: missing icons and tooltips, working voxel view. **This expectation is untested.**

## CSP allowlist (hosts permitted, not hosts used)
- `https://cdn.jsdelivr.net`
- `https://cdnjs.cloudflare.com`
- `https://esm.sh`
- `https://fonts.bunny.net`
- `https://fonts.googleapis.com`
- `https://fonts.gstatic.com`
- `https://unpkg.com`

Inner CSP sets `connect-src blob: data:` and `default-src 'none'`, so the page cannot make XHR/fetch/WebSocket calls to any origin. Combined with `sandbox="allow-scripts"` (no `allow-same-origin`) and `referrerpolicy="no-referrer"`, the data it displays cannot leave the page.

## Offline check, when run
Disable the network, open the file, and record: does the voxel trajectory render; do the MCS controls step; which icons/tooltips are missing; any console errors. Until then this file's offline behaviour is inferred from which libraries it loads, not observed.
