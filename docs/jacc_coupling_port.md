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

## The decomposition is now measured (2026-08-29)

Kernel agreement above is blind to the one defect the checkerboard can have of its
own. A parity-correlated bias in accepted moves passes `jacc_port_tests.jl`
whenever **both** kernels carry it, and `tests/fixtures/serial_seed42.csv` pins the
serial stream, which has no sublattices at all. `tests/jacc_parity_tests.jl` is that
measurement. `cpm_color!` gained two write-only per-site arrays — `st`
(0 never proposed / 1 evaluated-rejected / 2 evaluated-accepted) and `dh` (ΔH of
evaluated proposals) — reduced per sweep and **never drawn**: `d404438` refused a
spatial map for the serial decisive-label tally at 6157 of 64000 voxels touched,
median 3, and a parity bias is a global count comparison regardless.

Three choices in that tier are not obvious and each avoids a false failure on the
first run:

**The table is conditioned on opportunity.** The kernel's early returns — wall,
same-σ, medium-into-medium, out-of-bounds — are geometry-dependent, so a uniform
null over eight classes would report the shape of the domain as a decomposition
artifact. `st` carries the denominator and the 2×8 accepted/rejected table asks
about the acceptance *rate* per class.

**The thresholds are effect sizes, not χ².** n is 1.3e5–3.1e5 evaluated proposals
per run and the cells are not independent — an accepted move changes the lattice
for every later pass — so χ² over-disperses from autocorrelation alone and a fixed
critical value would test "is there any asymmetry at all". Cramér's V and the max
per-class rate deviation are asserted; χ² and n are reported beside them.

**`color_order` separates the decomposition from the `vols` staleness.** The colour
loop is sequential and `vols` accumulates across passes, so the first pass evaluates
against sweep-start volumes and the last against volumes moved by seven passes —
a deterministic, parity-correlated difference with nothing to do with the
checkerboard, and `c` indexes both spatial class and sequence position. Reversal
would separate only a monotonic position effect; random permutations drop that
assumption. The RNG step key stays `mcs*8 + c`, keyed to the colour rather than to
its position, so permuting changes the pass order and nothing else.

Measured at N=20, 50 MCS, threads backend, seeds 42/43/44 × three orderings:
V 0.0050–0.0111, max per-class rate deviation 0.017–0.049. **No decomposition
artifact detected, and the pattern tracks neither spatial class nor pass position** —
the third disposition, which the tier reports as unresolved rather than attributing.

### Reproducibility has a narrower scope than this document implied

`delta_H` reads `vols` while `cpm_color!` mutates it in the same pass. The
`ponytail:` comment calls that staleness bounded and unbiased, which it is, but
bounded is not deterministic: **the port's trajectory is not reproducible across
thread counts.** Measured on the threads backend, seed 42, N=40, 100 MCS — one
thread gives the same lattice every run, four threads gave three different lattices
in three runs. Any lattice-equality check against this port is meaningful only
single-threaded, and there is deliberately no lattice fixture in `tests/fixtures/`
for that reason: it would sit in a compare-never-regenerate directory while being
hardware- and thread-count-dependent.

## A Metal backend would need explicit Float32 conversion at the kernel boundary

Not exercised today, and recorded here only as a forward-looking note: `Project.toml` has no
`Metal.jl` dependency and no backend-selection code in this repository names a Metal target —
only `ROCBackend` and `ThreadsBackend` are. If a Metal backend is ever added, it would need
more than a `JACC.set_backend` call. `CPMParams`'s host-side fields (`T_cpm`, `β_ion`, `λ_V`,
`I0`, `κ`, `D_M`, `dt_field`, `D_C`, `C_wall`, `α_M_species`, …) are declared `Float64`, and
Metal does not execute `Float64` kernels — Apple GPUs support only 32-bit floating point.
The device-resident arrays this port already keeps on device (`lat`, `vols`, `spec`, the
ping-pong melanin/nutrient fields, and the coupling arrays listed above) are already `Float32`
by the rule stated at the top of this document, so the array side is fine; the parameter side
is not. Any Metal port would need an explicit `Float32` conversion at the point where these
host `CPMParams` values cross into a kernel launch, not an implicit one — the same "producer
declares, consumer must not assume" discipline `AGENTS.md` rule 4 states elsewhere, applied to
a numeric type instead of a semantic sentinel.

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
