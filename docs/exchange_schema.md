# HDF5 exchange schemas (CPM ↔ OpenMC)

Two deliberately distinct files. A **transport snapshot** is portable input for the
Python/OpenMC side: geometry-defining labels and reference fields only. A **restart
checkpoint** is the complete Julia simulation state, sufficient to resume the exact
stochastic trajectory. Conflating them was rejected in review: a snapshot cannot resume
the serial process (RNG position, evolving fields, radiodialysis state), and a restart
file is neither portable nor needed by transport.

Writer/reader: `export_checkpoint.jl` (Julia, HDF5.jl). Python readers live in
`coupling/biofilm_openmc/` (commit 5).

## Axis and index conventions (both files)

| Attr | Value | Meaning |
|---|---|---|
| `schema_version` | `1` | bump on breaking change |
| `logical_axis_order` | `"xyz"` | index semantics of every 3-D dataset, in Julia's native (column-major) layout |
| `dataset_axis_order_h5py` | `"zyx"` | what a C-order reader (h5py/numpy) sees: dims reversed. `arr.transpose(2,1,0)` in Python recovers logical (x,y,z). This is a storage-layout fact, not a coordinate transform |
| `coordinate_index_base` | `0` | orientation-probe and any stored coordinates are 0-based |
| `cell_id_background` | `0` | lattice value for medium |
| `cell_id_wall` | `-1` | lattice value for *outside the biological domain* — *not* a membrane material (audit §6) |
| `git_sha` | string | repo revision of the writer |
| `julia_version` | string | writer's Julia version |

ID conventions are separate attributes from coordinate conventions on purpose (review
amendment 9): one `index_base` must never describe both.

Coordinate mapping to OpenMC (Python side, commit 5): after recovering logical (x,y,z),
`RectLattice.universes` ordering is (z, y, x) **with increasing y-index at decreasing
physical y**, i.e. `arr.transpose(2,1,0)[:, ::-1, :]`. That mapping lives in exactly one
Python function and is pinned by the orientation probes end to end (Julia write → h5py →
RectLattice lookup → dose write → Julia read).

## `transport_snapshot.h5`

| Path | Type (logical shape) | Notes |
|---|---|---|
| `/` attrs | table above + `label_state_hash`, `mcs`, `physical_time_s` (NaN if unconfigured), `material_class_source` | `label_state_hash` = SHA-256 over cell/species/lineage/generation arrays |
| `/config_toml` | string | raw coupling config embedded verbatim; `""` when exported unconfigured. Python re-validates and fails loudly |
| `/lattice/cell_id` | Int32 (N,N,N) | 0 medium, −1 outside domain, >0 cell |
| `/lattice/species_id` | Int32 (N,N,N) | 0 where no cell; 1..7 authoritative registry |
| `/lattice/lineage_id` | Int32 (N,N,N) | 0 where no cell — **label array, never a material identity** |
| `/lattice/generation` | Int32 (N,N,N) | 0 where no cell |
| `/lattice/interior_mask` | UInt8 (N,N,N) | 1 = biological domain |
| `/fields/radiation_cpm` | Float64 (N,N,N) | dimensionless legacy Hamiltonian signal, reference only |
| `/fields/melanin` | Float64 (N,N,N) | dimensionless; input to any future config-defined classification |
| `/dose/accumulated_Gy` | Float64 (N,N,N) | physical; zeros until dose was ever imported |
| `/cells/{id,species,volume,lineage_id,parent_id,generation,birth_mcs}` | Int64 (n) | live-cell registry columns |
| `/cells/lifetime_dose_Gy` | Float64 (n) | |
| `/events/division/{...}` | Int64/Float64 columns | all 11 `DivisionEvent` fields |
| `/events/death/{...}` | columns | all `DeathEvent` fields |
| `/archive/{...}` | columns; `fate` as strings | archival cell table incl. extinct cells |
| `/orientation_probes` | Int64 (P, 4) | rows `[x, y, z, expected_cell_id]`, 0-based, ≥3 asymmetric points |

**No material composition is asserted here.** Material classes are a *config-defined*
mapping applied on the Python side from these label/field arrays; when the physical
config is absent the model builder fails loudly (`material_class_source = "absent"`).
This is what enforces "genealogy is a label, never a material" and keeps the class table
out of hard-coded Julia.

## `restart_checkpoint.h5`

Everything needed for `restore_restart_checkpoint` to produce a `CoupledSimulation`
whose continuation is bit-identical to an unbroken run (tested):

| Path | Notes |
|---|---|
| `/` attrs | conventions table + `mcs`, `physical_time_s`, `next_id`, `current_mcs` |
| `/lattice/cell_id`, `/lattice/interior_mask` | as above |
| `/fields/{radiation, melanin, melanin_drive, nutrient}` | all evolving/static CPM fields |
| `/dose/{rate_mean_Gy_s, rate_sd_Gy_s, increment_Gy, accumulated_Gy}` + attrs `dose_active`, `membrane_dose_rate_Gy_s` | full dose contract state |
| `/contaminant` | the 3-D projected contaminant field |
| `/cells/*`, `/archive/*`, `/events/*` | as above (complete tables) |
| `/rd/{r_grid, c, s}` + attrs `m`, `t` | radiodialysis state |
| `/rd_params/<field>` | every `RadiolysisParams` field as its own dataset (reconstructed via keyword constructor — no binary type serialization for ephemeral-module types) |
| `/cpm_params/<field>` | every `CPMParams` field likewise |
| `/rng/serialized` | `Serialization`-serialized `MersenneTwister` bytes (stdlib type: stable path). **Pinned to `julia_version`**; restoring under a different Julia version is refused unless `--allow-version-mismatch` |

Design note: `CPMParams`/`RadiolysisParams` are field-enumerated rather than
binary-serialized because they are defined inside the dynamically created sandbox module
(`load_serial()`), whose type paths are not stable across processes. The RNG is a
`Random` stdlib type, so binary serialization is stable there.

## Result files (commit 6, for reference)

`transport_result.h5` (physical dose field + run metadata + `transport_state_hash`) and
`dose_attribution.h5` (join of a possibly-cached field with *current* labels) are
separate files with separate hashes — a cached transport field may be legitimately
reused when only labels changed, without claiming the full checkpoint matched.

`transport_state_hash` covers what the OpenMC **run** depends on and deliberately excludes
`source.photons_per_second`: heating is tallied per source particle, so the activity cannot
change the transport result and must not invalidate it. A Gy/s field is keyed on
`fingerprint.dose_state_hash` (transport hash + activity) instead, so a scenario differing
only in activity reuses the run and still gets its own dose. Never reuse a `DoseResult`
across differing dose hashes.
