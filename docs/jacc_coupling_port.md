# JACC (GPU) coupling port — design only

No JACC coupling code exists on this branch, deliberately: the serial
implementation is the validation reference, and the JACC port's own contract
is *statistical* (checkerboard) rather than pathwise equivalence, so it must
never be the first causal comparison vehicle. This document specifies the
port so it can be implemented after the serial coupling is validated.

## What stays device-resident

The port (biofilms_potts_jacc.jl) keeps `lat`, `vols`, `spec`, parameter
vectors, and the ping-pong melanin/nutrient fields on device. New arrays for
coupling follow the same rule:

- `dose_rate_mean` / `dose_rate_sd` (Float32, N³): device-resident — read by
  the (future) dose-transformed Hamiltonian signal and melanin drive.
- `accumulated_dose` (Float32, N³): device-resident; per-MCS accrual is one
  elementwise kernel.
- `lineage_id` / `generation` per CELL (device vectors sized like `vols`):
  labels only. The per-voxel label arrays needed by the transport snapshot
  are derived on host at export time from `lat` + host mirrors — no
  device-side voxel label arrays needed.

## Host transfer points

- **Snapshot export** (every transport interval): `JACC.to_host(lat)` plus
  the accumulated-dose array. This pairs naturally with the existing
  `mcs % 10 == 1` host sync for radial biomass — transport intervals should
  be multiples of it.
- **Dose import**: one host→device copy per transport solve (fields are
  regenerated, not incrementally updated).
- HDF5 writing is host-only; reuse `export_checkpoint.jl`'s writers against
  host mirrors (`spec_h` already exists; add `lineage_h`, `generation_h`).

## Why exact serial/JACC trajectories are not expected

The 8-color checkerboard executes N³/8 site updates per color pass with a
role-inverted copy attempt and a counter-based RNG keyed by (site, mcs,
color). The update *schedule* differs from serial random-sequential dynamics
even though the single-flip acceptance physics is identical (the `--selftest`
asserts ΔH agreement site-by-site). Consequently only distributional
observables must match.

## The port is now executed (2026-08-16)

`--selftest` existed from the start and **nothing called it**. Neither
`tests/runtests.jl` nor `validate_serial.jl` touched the port — both load the
serial monolith through the split-marker trick — so it could have drifted or
broken outright while every suite stayed green. A radiodialysis substep guard
was added to it before any tier existed, and nothing would have caught a mistake.

`tests/jacc_port_tests.jl` closes that, from `runtests.jl`. It asserts the
kernel-level criterion this document already specifies, plus that the substep
rule is textually identical in both ports and is a no-op at the geometry every
call site uses — so adding it cannot have moved a published trajectory.

**It runs everywhere, which is the point of a portability layer.** JACC picks
its backend from `LocalPreferences.toml`, which is untracked; with no preference
file it defaults to threads, so a runner with no GPU executes these kernels for
real instead of skipping them. Verified on both:

| Backend | Result | Wall time |
|---|---|---|
| `ROCBackend` (gfx1150, Radeon 890M) | PASS | 13.9 s |
| `ThreadsBackend` (no preference file, as CI sees it) | PASS | 1.8 s |

The GPU being ~8x SLOWER here is not a defect and is worth reading correctly: at
N = 40 the lattice is ~6.4e4 sites, and the run is dominated by HIP kernel
compilation and launch overhead rather than by arithmetic. That matches the
local JACC 1.3.1 reduction benchmarks, where reusing a reducer beats the default
by 11-21x at n = 1e3-1e5 and the arms converge to within 1.07-1.23x by n = 1e7.
**The near-term work on this port is kernel structure, not more FLOPs**, and a
larger device would not change that.

## Ensemble observables that must remain statistically consistent

Across matched seed ensembles (the existing CSV contract, extended):

- per-species total volume, cell count, survival flags;
- mean melanin per species; membrane integrity `m`;
- NEW: per-lineage mass-weighted accumulated dose; division/death event
  counts per species; generation distributions.

Validation protocol stays the serial-vs-JACC CSV comparison across seeds,
now including the lifecycle/dose columns.

## Serializing genealogy events safely under checkerboard parallelism

- **ID allocation**: division needs new cell IDs while `vols`/`spec` are
  fixed-size device arrays. Pre-allocate capacity (2× initial cells is ample
  given division is a manual primitive) and allocate IDs with a device-side
  atomic counter bump.
- **Determinism**: thread order within a color pass is nondeterministic, so
  an atomic-append event log has nondeterministic ORDER. Key every event by
  `(site linear index, mcs, color)` — all available in-kernel — and sort
  host-side before logging/export.
- **Division itself cannot be an in-kernel operation**: it relabels many
  sites of one cell across all 8 colors. Detect candidates in-kernel (log
  only), then perform relabeling in a serialized host-side pass between MCS
  — same place the serial `divide_cell!` runs.
- **The stale-`vols` shortcut becomes load-bearing**: the documented
  `ponytail:` stale read in `delta_H` is bounded-unbiased for a Metropolis
  energy but NOT acceptable as a division trigger (double-fire/miss). Any
  volume-threshold check must use the atomic value or a post-pass reduction.
  (Moot while division stays manual — recorded for whenever a scheduler is
  scientifically justified.)
