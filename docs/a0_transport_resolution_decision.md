# A0 transport resolution/history decision

**Date:** 2026-08-14 · **Data:** `data/a0_transport_convergence.csv` (100 runs + 1 reference) ·
**Driver:** `coupling/scripts/a0_sweep.py` · **Metrics:** `coupling/biofilm_openmc/convergence.py`

## What this sweep decides, and what it does not

For the frozen A0 geometry and spectrum (`config/reference_a0_water_phantom.toml`: a 15 cm
radius × 30 cm water cylinder, idealized 661.655 keV point source, vacuum boundaries), it
answers **how many histories and what tally resolution give a stable per-source transport
estimate**, and it demonstrates the method the biofilm sweep will reuse.

It does **not** set the biofilm's production mesh. That depends on the calibrated domain
dimensions, biofilm density and water fraction, membrane geometry, source placement, metal
loading, and the actual spatial gradients — none of which exist yet. A0 also has no cells, so it
cannot test cell-, lineage- or generation-dose stability at all; those enter in the second pass.

Only rows whose `sensitivity_domain` is `transport_numerical` in
`data/parameter_provenance.csv` may be ranked from this sweep. It varied the tally mesh and the
history count and nothing else.

## Design

| Axis | Values |
|---|---|
| Coarsening factor | 1, 2, 4, 8, 16 (over a 48³ base, so meshes are nested) |
| Histories | 2.5×10⁵, 10⁶, 4×10⁶, 1.6×10⁷ |
| Independent seeds | 5 at every point |
| Reference | 48³ at 6.4×10⁷ histories, seed 9999 — 4× the largest swept point, and a seed no swept point uses, so nothing coincides with it and reports a spurious zero difference |

Geometry, materials and spectrum held exactly fixed. Fixed history counts with explicit
zero-score accounting rather than tally triggers: a trigger that ignores unscored bins can fire
early precisely when sparsity is the thing being measured.

The staged design (cheap meshes first, then escalate) turned out to be unnecessary. Measured
throughput was ~5.4×10⁵ histories/s, so the full 100-run matrix plus the reference took **6m36s
wall**. It is recorded here because the biofilm sweep, on a far finer mesh, will not be so cheap.

## Two accounting rules that change the answer

**Sparsity is measured inside the cylinder only.** The mesh spans the circumscribing cube, so
1 − π/4 ≈ 21.5% of bins are void by construction and *must* score zero. Counting them would add
a fixed offset to every mesh's sparsity number.

**Two different field metrics, because one of them cannot see discretization.** Since the meshes
are nested, summing the reference's deposited energy onto a coarse grid reproduces that grid's
*exact* true answer — so `field_diff_vs_reference` has zero discretization error by construction
and measures only Monte Carlo noise. To see what coarsening destroys, `resolution_loss` broadcasts
the coarse dose *rate* (intensive) back onto the reference grid and compares there. Both are
reported; conflating them would have made every mesh look convergent.

## Results

Mean over 5 seeds. `noise` = `field_diff_vs_reference`, `resloss` = `resolution_loss_vs_reference`.

| f | bins | bin cm³ | across cyl. | histories | median rel err | p90 | zero % | noise | resloss | s | MB |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 110 592 | 0.244 | 48 | 2.5×10⁵ | 0.613 | 1.000 | 3.07 | 0.348 | 0.348 | 0.3 | 1.8 |
| 1 | | | | 1.6×10⁷ | 0.095 | 0.140 | 0.00 | 0.049 | 0.049 | 11.2 | 1.8 |
| 2 | 13 824 | 1.953 | 24 | 2.5×10⁵ | 0.242 | 0.434 | 0.00 | 0.125 | 0.250 | 0.3 | 0.2 |
| 2 | | | | 1.6×10⁷ | 0.034 | 0.050 | 0.00 | 0.018 | **0.177** | 9.5 | 0.2 |
| 4 | 1 728 | 15.6 | 12 | 1.6×10⁷ | 0.012 | 0.018 | 0.00 | 0.006 | **0.306** | 9.1 | 0.1 |
| 8 | 216 | 125 | 6 | 1.6×10⁷ | 0.005 | 0.007 | 0.00 | 0.003 | **0.502** | 9.1 | 0.0 |
| 16 | 27 | 1000 | 3 | 1.6×10⁷ | 0.002 | 0.003 | 0.00 | 0.001 | **0.620** | 9.5 | 0.0 |

Three things stand out.

**Sparsity is not a constraint at dosimetry scale.** The worst case in the entire matrix is 3.07%
of in-cylinder bins unscored, at the finest mesh and lowest history count; everywhere else it is
0.00%. This is the direct contrast with the foundation's biofilm feasibility run, where 75% of
occupied bins were unscored at a 10 µm pitch. Sparsity is a *scale* problem, not a method problem,
and the A0 geometry does not have it.

**Statistics behave, and are cheap.** Median relative error falls as 1/√N to within a few percent
at every mesh. Histories needed for a median relative error ≤ 10%, extrapolated from the 1.6×10⁷
points: 1.4×10⁷ at f1, 1.8×10⁶ at f2, 2.4×10⁵ at f4.

**Resolution loss is an N-independent floor.** It does not improve with more histories — f2 sits
at 0.177 whether run with 4×10⁶ or 1.6×10⁷ histories, f4 at 0.306, f8 at 0.502. That is the
signature of discretization rather than noise. The A0 point source makes this severe: a 1/r²
singularity at the centre means the near-field gradient is not representable at any coarsening,
and the mass-weighted metric is dominated by exactly those high-dose bins.

## The hard ceiling: coarsening ≤ 4

The mass-weighted mean dose — the aggregate the sensitivity gate actually compares — is stable
across the first three meshes and then breaks:

| f | mass-weighted mean (Gy/source) | drift vs next finer | ROI bins | in-cylinder fraction |
|---|---|---|---|---|
| 1 | 2.348906×10⁻¹⁵ | — | 86 592 | 0.783 (π/4 = 0.785) |
| 2 | 2.357055×10⁻¹⁵ | 0.35% | 10 752 | 0.778 |
| 4 | 2.347228×10⁻¹⁵ | 0.42% | 1 344 | 0.778 |
| 8 | 2.072563×10⁻¹⁵ | **11.70%** | 192 | 0.889 |
| 16 | 1.844198×10⁻¹⁵ | **11.02%** | 27 | **1.000** |

The cause is diagnosable and is not physics. At f16 the mesh is 3×3×3 and *every* bin centre
falls inside the cylinder, so the region of interest degenerates to the whole cube and void
corners are counted as water; at f8 the in-cylinder fraction is 0.889 against a true 0.785. The
drift is a region-of-interest and full-bin-mass discretization artifact — the documented
approximation in `phantom_mass_kg` and `voxel_mass_kg`, whose upgrade path is OpenMC's ray-traced
`mesh.material_volumes()`. Beyond f4 the cylinder is simply no longer resolved (6 and 3 bins
across it), and any aggregate over it stops meaning what it says.

Two independent checks confirm the physics itself is resolution-invariant, as it must be:

- **Absorbed fraction is identical to six significant figures at every mesh**: 0.469710. The
  tally is conservative and coarsening moves no energy.
- **Seed spread of the aggregate is 0.012–0.014% at every mesh**, from 5 independent replicates —
  so the 11% drifts above are far outside the statistical noise and are structural.

## Decision

For **Reference A0**, declared as engineering targets rather than physical laws:

```
transport.mesh.base_dimension   = [48, 48, 48]
transport.mesh.coarsening_factor = 2      # 24^3, 1.95 cm bins, 24 across the cylinder
transport.batches                = 80
transport.particles              = 50000  # 4e6 histories
```

`coarsening_factor = 2` is the recommendation because it resolves the region of interest properly
(in-cylinder fraction 0.778 against a true 0.785), holds the aggregate within 0.35% of the finest
mesh, reaches a median relative error of 6.7% and a p90 of 10.1% at 4×10⁶ histories in under 3
seconds, and costs 0.2 MB per statepoint against 1.8 MB at f1.

**Hard ceiling: `coarsening_factor ≤ 4`.** Beyond it the in-cylinder region is unresolved and
aggregates drift ~11% per step for reasons that have nothing to do with transport.

**Choose f1 instead when the consumer is a per-bin field rather than an aggregate.** f2 already
discards 17.7% of the field structure by the resolution-loss measure, and f4 discards 30.6%. Which
of these two criteria binds depends on what reads the dose field — and A0, having no cells, cannot
decide that. It is the question the biofilm sweep exists to answer, with cell-, lineage- and
generation-dose rank stability as its metrics.

## Caveats

- `peak_child_rss_kb_cumulative` in the CSV is a running high-water mark across all child
  processes, so only the largest value carries information. Per-run peak memory was not measured.
- Statepoint sizes are for a single heating tally; adding scores or filters scales them.
- Runtimes are wall-clock on one machine with the pinned OpenMC 0.15.3 build and ENDF/B-VIII.0
  data, sharing the machine with nothing else. They indicate scaling, not a portable budget.
- The reference is itself a Monte Carlo estimate (6.4×10⁷ histories, median relative error ~0.05
  at f1). Differences below roughly 1% at f1 are within its own noise.
