# Inert signal diagnostic: executed validation

The delivered derived run is `inert_signal_seed42_v2`, driven by the immutable
`manuscript42_seed42_8cbb8ee` trajectory. Sources use occupancy at MCS k over
[k,k+1); masked walls absorb and outer cube faces have no flux.

- All 101 Float64 HDF5 signal frames and their VTI counterparts: 1,414 checks passed, including exact deterministic solver replay, axes, masks and endpoint counts.
- Numerical contract and production-source mutation controls: 38 passed. Tests include all seven distinct source rates, time and slab convergence, a negative field when stability substeps are removed, and a failed slab bound when decay is removed.
- Production input/configuration controls: 20 passed on the parent run; unknown acknowledgements, missing/duplicate MCS, wrong parcel configuration, altered hash and existing output are refused.
- Native viewer HDF5 reader and CLI: 107 passed, including all 101 frames and mismatched-time/altered-field controls.
- Full Julia suite: passed, including the existing basis-gate census, accepted-copy controls and the new numerical tests.
- Calibration tier: 428 passed, 5 expected skips for unavailable Dryad ND2 inputs.
- Browser: all 101 reconstructed label hashes match, all 14,281 accepted events replay, controls work, and light/dark narrow layouts have no runtime errors or horizontal overflow.
- Original trajectory full postflight: passed. This command rewrites receipts; the authoritative archive receipts were restored and their hashes reverified. Run it on a separate extraction for future checks.
- Staged manuscript: 40 pages; no undefined references or overfull boxes. The new section and figure pages were rendered and inspected.

At MCS 100, A >= 5 at 3,384 / 4,888 occupied interior voxels (69.2%). The first
saved occupied threshold crossing is MCS 8. These are conditional numerical
outputs, not biological QS-response evidence. The earlier prototype has a
separate convention and is not used for the delivered figure.

The native GLMakie window has not been executed here. Its reader and transfer
map resolve, but OpenGL rendering still needs a target Mac/Linux check. The
ParaView files were read back numerically; no local ParaView application was
available for a GUI test. The browser raster view was executed and inspected.

## Pinned files

| Artifact | SHA-256 |
|---|---|
| Parent run manifest | `f00e470afd982620a15def3a2e0624eb50f1ad7c09b3ef063b219c185d91c7fe` |
| Delivered derived manifest | `9768e0a9f4f6387315ccd863a75d198216c01458f546fa5b11c6ece81bac90c3` |
| Original Figure 5 PDF | `54fb219c6545f4d1536bd078f368613f6c05df82c4d2f2e1731ce5473df33db0` |
| New signal figure PDF | `c83b37e2fb9e88dea5e04b3186eb730d8c3fde4df370e8a9a3f7ddcfe4c09773` |
| Staged manuscript PDF | `dea95758b2c583ab816b82fe2088a19f50864f81ae6999051772e6bf1cf52d2c` |

The archive copy of the original manuscript remains authoritative for the
original evidence. The staged PDF includes a new signal section and figure.
PR #25, Figure 1, and the original Figure 5 are unchanged.
