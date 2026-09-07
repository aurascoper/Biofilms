# 3-D lattice viewer and .vti export: GLMakie voxels, replacing feat/visualize-3d

Evidence discipline for this note. Every API claim cites a documentation URL and the version
consulted or, where the claim was checked by running code in this repository, the resolved
package version and the file that ran. Every "not present" claim states the scope searched.
No physical unit is attached to a lattice quantity anywhere below.

## Summary

Two tools, added on `feat/lattice-viewer` (2026-09-07), both reading the transport snapshot
that `export_checkpoint.jl` already writes and neither re-running the simulation:

1. `export_vti.jl`, a pure, CI-covered exporter from a transport snapshot to VTK ImageData
   (`.vti`) for ParaView, and from a directory of snapshots to a `.pvd` series keyed by MCS.
   Tests in `tests/vti_export_tests.jl`, run by `tests/runtests.jl`.
2. `viewer/visualize_lattice.jl`, an interactive viewer on GLMakie's `voxels` recipe, in its own
   Julia project so the main environment pulls no OpenGL dependency. Not run in CI (see
   "What CI covers").

Both write lattice units only. Spacing is 1.0 per site and the file says `units = lattice` in
its field data. No metre spacing and no MCS-to-seconds mapping is written, because neither is
measured. `D-PITCH` and `D-TIMESERIES` are `supplied_by = measured`, `status =
awaiting_measurement` in `data/calibration/reference_d_requirements.csv`. A physical spacing
can be written only by passing a declared pitch explicitly, and the file then says
`declared:` with the value, unit and source.

## What `feat/visualize-3d` was (read locally, 2026-09-07)

The scope searched was `git log origin/master..origin/feat/visualize-3d`, `git diff --stat`, and the file
itself, in a checkout of this repository on 2026-09-07.

One commit, 50228b8 (2026-08-13), on base ec8929d, 246 files behind master, without `tests/`,
the claims ledger, or the current preprint (its preprint is the pre-revision title). It added
`visualize_3d.jl` (80 lines) and nothing else. That script used CairoMakie `meshscatter` of
unit cubes on an `Axis3`, re-ran the simulation in-process with the `include_string` trick
from `validate_serial.jl`, never read a checkpoint, recorded an orbiting GIF to the repository
root, and titled every frame with the first word of `RETRACTED_IN_FIGURES`
(`calibration/tests/test_claims_ledger.py`). Its Project.toml deps were AMDGPU, CairoMakie and
JACC; it had no Manifest and no test. Nothing on it is a base. Two things were salvaged, the
seven-entry palette and labels (the serial script's `FIG_COLORS` / `FIG_LABELS`) and the shape
of the frame loop. Retiring the remote branch is Hunter's call; this note records the
decision when it is made.

## The snapshot both tools read

`export_checkpoint.jl` on master a3025e0 (lines 29 to 40 and 148 to 166) and
`docs/exchange_schema.md`. The transport snapshot carries:

| dataset | dtype | meaning |
|---|---|---|
| `lattice/cell_id` | Int32 (N,N,N) | 0 background, -1 wall, otherwise the cell id |
| `lattice/species_id`, `lattice/lineage_id`, `lattice/generation` | Int32 (N,N,N) | derived per site by `_label_arrays` |
| `lattice/interior_mask` | UInt8 (N,N,N) | the cylinder |
| `fields/radiation_cpm`, `fields/melanin` | Float64 (N,N,N) | CPM fields |
| `dose/accumulated_Gy` | Float64 (N,N,N) | physical Gy; zeros until a dose was imported |
| `orientation_probes` | Int64 (k,4) | 0-based x, y, z, value rows |

and the attributes `schema_version`, `logical_axis_order = "xyz"`, `dataset_axis_order_h5py =
"zyx"`, `coordinate_index_base = 0`, `cell_id_background = 0`, `cell_id_wall = -1`, `git_sha`,
`mcs`, `physical_time_s`. The tools read those attributes and carry them; they restate none of
them.

Two things the snapshot does not carry, with the dataset list written at those lines as the scope. One is
the nutrient field, which is only in the restart checkpoint (`fields/nutrient`), the other is the OpenMC
dose field. The per-source-particle field never reaches a file at the transport stage
(`coupling/biofilm_openmc/drivers.py:185` and `:328`); what `transport_result_*.h5` holds is
`mesh/dose_rate_mean_Gy_s`, in Gy/s with the source activity applied
(`coupling/biofilm_openmc/results.py:51`). The exporter takes both as optional inputs and
labels the dose array by that dataset name.

## GLMakie `voxels`

Claims from the documentation, as consulted for the draft of this note (2026-09-07), then what
running the viewer here established.

| Item | Value | Source | Version |
|---|---|---|---|
| `voxels` introduced | Makie v0.21, PR #3527 | https://makie.org/website/blogposts/v0.21/ ; https://docs.makie.org/v0.21/changelog | v0.21 |
| Signature | `voxels(chunk::Array{<:Real,3})`, `voxels(x, y, z, chunk)`; only the extrema of x, y, z are used | https://docs.makie.org/dev/reference/plots/voxels | dev |
| Representation | `Array{UInt8,3}`; `0x00` is always invisible air; ids 1 to 255 visible | same | dev |
| `color` per id | `color = [c1, c2, ...]` indexed by voxel id, skipping `0x00`; `colorrange` ignored for UInt8 input | same | dev |
| `is_air` | a predicate selecting values drawn as air | same | dev |
| Status | "experimental and may still see breaking changes in patch releases" | same | dev |
| Backends | dedicated implementation in GLMakie and WGLMakie; CairoMakie projects 3-D flat with no z-clipping | https://makie.org/website/blogposts/v0.21/ ; https://docs.makie.org/v0.21/explanations/backends/cairomakie | v0.21 |
| GPU | OpenGL 3.3 or higher | https://docs.makie.org/stable/explanations/backends/glmakie.html | stable |
| Headless | xvfb software rendering is how GLMakie's own tests run | https://docs.makie.org/dev/explanations/headless | dev |

Established by running `viewer/visualize_lattice.jl` on 2026-09-07 with the versions pinned
in `viewer/Manifest.toml` (GLMakie 0.13.14, Makie 0.24.14, HDF5 0.17.3), on a Radeon 890M at
OpenGL 4.6 with a display. A UInt8 grid with `0x00` for medium and wall, `color` as the
seven palette entries and `is_air = ==(0x00)` renders every occupied site of a 20-site
snapshot at MCS 20, with the legend built from the species present. A still and a
12-frame orbit were written to the scratch directory and not committed. `attributes` is
exported by both HDF5 and Makie, so the viewer qualifies it as `HDF5.attributes`.

The viewer's grid is a property of the file's own sentinels. A site is air when its
`cell_id` equals the snapshot's `cell_id_background` or `cell_id_wall` attribute, and
species id otherwise. That is rule 4 of AGENTS.md applied to a plot.

## The .vti exporter as built

`export_vti.jl`, functions plus a CLI guard, the `export_checkpoint.jl` pattern.

```
julia --project=. export_vti.jl <snapshot.h5> <out_stem> [--restart r.h5] [--dose d.h5]
julia --project=. export_vti.jl <snapshot_dir> <out_stem> [--restart-dir d] [--dose-dir d]
```

The grid is `vtk_grid(stem, 0:N, 0:N, 0:N; compress = false)`, which WriteVTK writes as ImageData
(`.vti`) with N cells per axis and spacing 1.0. A CPM label is a property of a site, so every
array is cell data. VTK x is Julia's first index. Both are fastest-varying, so `A[i,j,k]` lands
on VTK (x,y,z) with no permutation. That is a claim, so it is tested, twice (below).

| cell array | dtype | source |
|---|---|---|
| `species` | UInt8 | 0 where cell_id is background or wall; the species id elsewhere |
| `cell_id` | Int32 | as stored, sentinels kept |
| `lineage_id`, `generation` | Int32 | as stored |
| `interior_mask` | UInt8 | as stored |
| `radiation_cpm`, `melanin` | Float64 | as stored |
| `accumulated_dose_Gy` | Float64 | as stored; field data `accumulated_dose_Gy_units` says "Gy, physical, schema dose/accumulated_Gy: zeros until a dose was imported" |
| `nutrient` | Float64 | only with `--restart`, from `fields/nutrient` |
| `dose_rate_mean_Gy_s` | Float64 | only with `--dose`, from `mesh/dose_rate_mean_Gy_s`; refused unless its mesh is the lattice (the exporter resamples nothing; `viewer_bundle.h5` is where that happens) |

The field data holds `units` ("lattice", or a "declared" line with value, unit and source), `species_zero`,
the dose unit strings, and the snapshot's own attributes `schema_version`,
`coordinate_index_base`, `cell_id_background`, `cell_id_wall`, `mcs`, `physical_time_s` as
numbers and `logical_axis_order`, `dataset_axis_order_h5py`, `git_sha` as strings.
`physical_time_s` is carried as the attribute it is; it is not the `.pvd` time key and not a
spacing. The `.pvd` key is the `mcs` attribute, a step count. Two snapshots carrying one
`mcs` are refused.

The spacing rule is that `spacing != 1.0` throws unless `declared_pitch = (value, unit, source)` is
passed. That is the no-fabricated-metres rule as code rather than prose.

WriteVTK 1.22.0 and ReadVTK 0.2.6 are ordinary `[deps]` in `Project.toml`. CI's Julia job runs
`Pkg.instantiate()` and then `tests/runtests.jl` directly, never `Pkg.test`, so an `[extras]`
entry would not be installed there (`.github/workflows/coupling-tests.yml`, julia-tests).
The file is written uncompressed (`compress = false`) so its string field data can be read
back in the tests. ReadVTK reads numeric arrays only, and its own documentation calls it
incomplete (https://juliavtk.github.io/ReadVTK.jl/stable/). WriteVTK sources:
https://juliavtk.github.io/WriteVTK.jl/stable/ and
https://juliavtk.github.io/WriteVTK.jl/stable/grids/datasets/ (ImageData from ranges, cell
versus point placement, `VTKFieldData`, `paraview_collection`), consulted 2026-09-07.

## Tests and their controls (`tests/vti_export_tests.jl`, 44 assertions)

The fixture is an N = 12 snapshot with two cells per species written by `export_transport_snapshot`
after two MCS (N = 8 cannot be initialised; N = 10 places five of seven species). The
fixture's `basis_gate_ack` is enumerated in the census in `tests/radiodialysis_basis_gate.jl`.

- Round trip. Every site array reads back equal, element type included; `species` is 0 exactly
  where `cell_id <= 0` and the species id elsewhere; all seven species present; spacing 1.0;
  `units = lattice`; the dose unit string and `logical_axis_order` read back.
- Sentinels. Background 0 and wall -1 survive in `cell_id`.
- Orientation, twice. A 3x4x5 array valued 100i+10j+k reads back index for index, and the
  exporter's own `cell_id` output satisfies every row of the snapshot's `orientation_probes`.
- Refusal. Spacing 0.012 throws and writes nothing; with a declared pitch it writes 0.012 and
  `units` begins "declared:"; a snapshot whose `logical_axis_order` is not "xyz" is refused.
- Series. One `.vti` per snapshot, `.pvd` timesteps 2.0 and 5.0, a duplicate `mcs` refused.

Controls, each planted on committed state d49aac0, run through a narrow driver that includes
the test file alone, restored from a scratch copy and confirmed byte-identical, tree clean:

| planted in `export_vti.jl` | red | green after restore |
|---|---|---|
| `species` written as `permutedims(species, (3, 2, 1))` | 42 pass, 2 fail (round trip) | 44 |
| the spacing refusal replaced by `if false` | 42 pass, 2 fail (refusal) | 44 |
| one label corrupted, `species[2, 3, 4]` flipped | 43 pass, 1 fail (round trip) | 44 |

## What CI covers, and what it does not

The exporter and its tests run in the julia-tests job. The viewer does not, because GLMakie needs an
OpenGL 3.3 context, and the hosted runner has none. That is uncovered surface, stated here,
not a skip. An xvfb workflow is possible and was not added; if one is added it is opt-in.

One guard was missing and is now present. Planting the retracted word in the viewer's title
on 2026-09-07 left every suite green. The figure-vocabulary guard reads
`preprint/figures/*.txt` sidecars (`figure_sidecars()` in `test_claims_ledger.py`), and the
ceiling-vocabulary walk scans three other terms. `tests/manuscript_claims_tests.jl` now scans
every `.jl` under `viewer/` and `export_vti.jl` for that word with a synthetic control; the
same plant on committed state a44a792 fails it at the planted line, and restoring the file
returns 117 of 117.

Any image from either tool intended for the manuscript is registered with the figure
staleness guard first (FIG-01 to FIG-07 in `data/claims_ledger.csv`). Until then, viewer and
exporter output is exploratory.

## ParaView workflow

ParaView is not installed on the machine this was built on (the scope searched was `which paraview
pvpython` and `/opt`, 2026-09-07), so the steps below are from the ParaView documentation and
were not exercised. Open the `.vti`, or the `.pvd` for a series. Threshold on `species` to
drop 0 (air). Colour by `species` with "Interpret Values As Categories" and one name per id.
Clip or slice through the cylinder axis for a radial section. A Calculator expression
`sqrt((coordsX - x0)^2 + (coordsY - y0)^2)` gives a radial coordinate in sites. Colour by
`dose_rate_mean_Gy_s` or `accumulated_dose_Gy` for the dose, keeping the unit strings from
the field data in any caption. The Python side can read the same file with `pyvista.read`
(https://docs.pyvista.org/api/utilities/_autosummary/pyvista.read.html); pyvista, vtk and h5py
are in the coupling venv here, none of them in CI's tier.

## Open decisions

- Retiring `feat/visualize-3d` on origin is Hunter's call, recorded here when made.
- Per-source-particle dose. It never reaches a file at the transport stage, so the exporter
  cannot offer it. If it should, the writer in `drivers.py` is where it would be added.
- Nutrient in the transport snapshot. The exporter takes it from a restart file; adding it to
  the snapshot is a schema change with its own hash implications (`label_state_hash` excludes
  fields, so it may be safe; not checked).
- WGLMakie for a browser viewer was not tried.
- ParaView version to pin. None is installed; the draft's 6.1.0 (release notes,
  https://www.kitware.com/paraview-6-1-0-release-notes/) is unconfirmed here.
