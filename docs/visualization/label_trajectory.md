# Verified manuscript-configuration label trajectory

`lattice_evidence.jl` is the dedicated entrypoint for N=40, six parcels for
each of seven species, seed 42, MersenneTwister and the existing coupled runner.
It writes MCS 0 before advancing and one snapshot per MCS through 100. The
two-per-species demo defaults in `export_checkpoint.jl` are unchanged.

Run from the repository root with Julia 1.12.6 and the project instantiated:

```sh
julia --project=. lattice_evidence.jl run manuscript42_seed42_v1
julia --project=. lattice_evidence.jl verify out/lattice_evidence/manuscript42_seed42_v1
python tools/render_label_trajectory.py out/lattice_evidence/manuscript42_seed42_v1 --install
```

The first command performs both the exact production export and its postflight.
The second is a separate verification entrypoint. The renderer requires numpy,
matplotlib, h5py, latexmk and Poppler. It stages the manuscript within the run
directory; `--install` copies only the new Figure 5 and generated numerical
fragment into `preprint`. The manuscript PDF is an isolated build artifact.
The renderer refuses an existing manuscript staging directory. Use a new run
ID for a fresh execution; existing run directories are always refused by the
producer. No open `snap.h5`, `snaps/`, PVD or ParaView state is touched.

`--config FILE` is optional embedded provenance, not a CPM parameter override.
The full restart checkpoints at MCS 0 and 100 preserve the RNG and all state.
Transport snapshots do not replace complete restart checkpoints. The underlying
initial restart's basis may be standalone; the diagnostic's separate
`c_s_analysis_blocked` declaration remains true at both endpoints.

The exact basis-gate acknowledgement census includes only the canonical factory
as the new opening. It permits stepping to record labels, never a quantitative
claim from mobile/sorbed fields. Their arrays exist only where required for a
complete restart. No analyser or renderer reads them. The postflight compares
every one-MCS label-state hash after multiplying the uptake parameters by ten;
the unit control also passes perturbed output through the real snapshot writer
and reader. Existing physical calibration and institutional gates are unchanged.

`on_accepted(event)` is called after the acceptance decision and before copying.
Its immutable value record includes no state or RNG object. `proposal_index`
counts all attempts, including skipped attempts. Sites are 1-based Julia linear
indices with x fastest. The columnar `accepted_copies.h5` file records the four
terms, their sum and the exact existing uniform draw, with NaN for downhill
moves that never draw. Mutation replay consumes the records in (MCS, proposal)
order and requires each saved label-state hash. Lifecycle mutation would need
an additional event schema before enabling it; this run has no automatic split
scheduler.

The fixed interior includes empty label 0. Walls are excluded using the mask,
not by a time-varying species threshold. Species-indicator means are sampled
occupancy; numeric species codes are never averaged. Identity statistics are
computed separately for parcel IDs and species. A completed return is a sampled
departure from the initial label, one or more different labels, then arrival at
the initial label. The initial observation is not a return. End-of-window
departures without a return remain incomplete. Persistence requires agreement
with the reference at every saved observation, while endpoint overlap permits
intervening departures.

The analyser exports maps and initial-label/initial-species/all-interior
denominators with persistence curves for cadences 1, 2, 5 and 10. All cadences
sample one trajectory over the same 0-100 window. A hidden reversal is one
site--coarse-interval pair with equal coarse endpoint labels and a visible
one-MCS change inside that interval. Multiple excursions inside that coarse
interval count once in this metric. No snapshot cadence resolves events that
reverse within a single MCS; only the accepted-copy log supplies that ordering.

`run_manifest.json` identifies the run, source commit and current tracked/source
bytes, Julia executable and project/manifest hashes, parameters, seed/RNG,
gate declarations, snapshots and all resulting artifacts. It enumerates every
run file except itself and its own receipt: `run_manifest.sha256` hashes the
manifest, avoiding an impossible self-referential checksum. Regenerating a
derived artifact requires refreshing the manifest through the declared build
workflow; changing a hash does not replace running verification.

Validation covers the 42-parcel initial inventory, all 101 MCS and masks,
registry/lattice correspondence, restart equivalence, RNG continuation,
instrumentation inertness, one-versus-100-MCS windows, and the legacy coupled
runner at 0, 20, 40, 60, 80 and 100. The label tests include equal occupancy with
different history, a hidden reversal, same-species parcel handoffs, a label-code
permutation, empty space and excluded walls. Corrupt snapshots and an undeclared
gate-opening source are rejected through the production readers/walker.

Figure 5 uses CPU rasterized voxel faces, a common camera and categorical species
colours. The transition and return maps use the lower central z plane (Julia
index 20 on the 40-cube). Figure 1 is hashed before execution and checked again
after verification and figure installation. The new panel describes CPM parcel
labels; it is not a membrane, EPS or biological reconstitution measurement.
