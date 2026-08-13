# Branch report: feat/openmc-dose-coupling

**Start:** `ec8929d` (clean master baseline) · **End:** `9e63eb7` + this report ·
**10 commits, 44 files changed, ~3850 insertions** · 2026-08-13

Unrelated visualization WIP preserved separately on `feat/visualize-3d`.

## What was built (per the reviewer-amended sequence)

1. `dc5dcb5` audit + contract/deterministic tests
2. `d53f236` presentation defects (R sign dead code, neutral Fig-1 labels, 50 Gy literal removed)
3. `69eea21` serial refactor: clock, dose contract, lifecycle archive, manual `divide_cell!`, windowed API
4. `b674a72` transport snapshot + restart checkpoint schemas, orientation probes, bit-exact resume
5. `e57c3c6` pinned OpenMC stack, config contract, Python package, CI
6. `23af70c` exact fingerprints; transport_result / dose_attribution split
7. `7ca2f37` integration benchmarks (transmission, energy balance, feasibility)
8. `e1f1ead` offline sensitivity driver with explicit feedback gate
9. `86218dd` decision docs (online design + gates, ADDER, JACC port)
10. `9e63eb7` benchmark fixes after first live transport runs

## Physical assumptions supplied by configuration

**None.** Every physical value (voxel pitch, origin, cylinder/membrane geometry,
boundary conditions, seconds/MCS, source rate, spectrum, spatial/angular
distributions, material densities/compositions, per-consumer dose transforms)
ships REQUIRED-but-unset in `config/coupling_template.toml` and fails loudly.
Test fixtures use clearly-labeled synthetic demo values only.

## Assumptions intentionally left unresolved (stop conditions)

- All physical config values above — user-supplied when known.
- Membrane m-vs-P_eff constitutive choice (audit §6): STOP CONDITION — no OpenMC
  dose enters the radiodialysis model until declared.
- Membrane dose statistic (mass- vs area-weighted): config key exists, unset.
- Photon spectrum/activity for the thorium discussion: not derivable; REQUIRED input.

## Tests

| Suite | Result |
|---|---|
| Julia (`tests/runtests.jl`) | **135 passed, 0 failed** (contract CSV byte-identical; ΔH signs; multi-seed drift; lifecycle; windowed ≡ legacy; snapshot probes; bit-exact restart resume) |
| Python bare venv (no OpenMC) | **38 passed, 3 SKIPPED** (openmc-gated; skipif — never deselected, never passed) |
| Python conda tier (OpenMC 0.15.3 + ENDF/B-VIII.0) | **44 passed, 0 failed, 0 skipped** — all gated tiers executed |
| JACC `--selftest` | PASS against the refactored monolith |

## OpenMC execution

**Yes — real transport ran.** Pinned conda-forge `openmc 0.15.3`
(dagmc_nompi_py313h6424856_102), Python 3.13.15, HDF5 1.14.6 (docs/openmc_stack.md).
Photon interaction data verified: official ENDF/B-VIII.0 HDF5 library, archive
SHA-256 `200cc6b9…da9b970`, `cross_sections.xml` SHA-256 prefix `ba7f6d5a371b5a8d`.

Benchmark results (synthetic demo config, 250k histories):

- **Transmission**: uncollided partial current through 1/2/3 cm water is
  exponential (ln-steps agree within 15% tolerance; energy filter isolates the
  unscattered 1 MeV line).
- **Energy balance**: closed reflective water box deposits 1 000 000.0 eV per
  1 MeV source photon — relative error 5.2×10⁻⁹.
- **Feasibility (transport uncertainty & sparsity)**: at cell-scale voxels
  (10 µm demo pitch), 75% of occupied voxels scored no heating and median
  rel_err ≈ 1 at these statistics. Lineage-resolved dose will require large
  history counts, spatial coarsening, or both — exactly the coarsening
  decision the review flagged; recorded here, not hidden. All dose values are
  voxel-averaged absorbed dose under OpenMC's charged-particle local-deposition
  approximation, not single-cell microdosimetry.

## Offline feedback gate

**NOT EVALUATED** — by design: `biofilm-dose-scan` refuses to evaluate while
the physical config is unset (exit 2 with the full missing-key list). The gate
machinery itself is tested end-to-end with an injected fake runner (PASSED /
FAILED / NOT EVALUATED paths all covered).

## ADDER

Not in the runtime path (docs/adder_applicability.md): stock ADDER interfaces
MCNP5/6 + ORIGEN/CRAM; the repo's field is gamma-only and OpenMC photon
transport has no photoneutron path, so no depletion mechanism exists here;
biological/polymer updates are not depletion. Revisit triggers documented.

## Preprint

`preprint/**` untouched. Regenerating figures after commit 2 (neutral zone
labels, removed 50 Gy annotation) will differ from the committed versions —
recorded in audit §10, deliberately not regenerated on this branch.
