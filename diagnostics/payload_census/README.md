# Payload census: an exporter audit, and two representation dispositions

## What this found

The question was whether a sparse voxel DAG is worth adopting. Measuring to answer it turned
up something more useful: **more than half of every frame is written for no reason.**

```
101 frame(s), 64,000 cells each, 47 bytes/site

array                      type     B/site  %frame  distinct  const  static  derivable from
radiation_cpm              Float64       8   17.0%       198            yes
melanin                    Float64       8   17.0%     50201
accumulated_dose_Gy        Float64       8   17.0%         1    yes     yes
signal                     Float64       8   17.0%     50201
cell_id                    Int32         4    8.5%        44
lineage_id                 Int32         4    8.5%        43                 cell_id
generation                 Int32         4    8.5%         1    yes     yes
species                    UInt8         1    2.1%         8                 cell_id, lineage_id
interior_mask              UInt8         1    2.1%         2            yes  cell_id, radiation_cpm
occupied_above_threshold   UInt8         1    2.1%         2                 signal (weak)

  on disk                            290.14 MiB
  gzip -9                             72.40 MiB  (4.01x)
  static across frames                44.7% of every frame, written 101x, needed once
  derivable from another array        10.6%
  neither                             44.7%  <- the real payload
  binary occupancy, 1 bit/site       0.266%  <- the DAG ceiling (0.77 MiB)
```

`radiation_cpm` is **bit-identical across all 101 frames**: 49.3 MiB (8 × 64,000 × 101 =
51.7 MB) of a field that never changes. `accumulated_dose_Gy` is identically zero, another
49.3 MiB. `generation` is
identically zero because `divide_cell!` has no trigger. `interior_mask` never moves.

`lineage_id` equals `cell_id` on every occupied site and zero elsewhere, verified frame by
frame, and cannot diverge because divergence happens only at division. `species` is a
42-entry lookup on `cell_id`, consistent across the whole trajectory.

This redundancy is **semantic, not statistical**, and dropping it is lossless by
construction, which no compressor can claim.

### One correction to the derivable list

`occupied_above_threshold` is not `signal > 5.0`. Measured, it is
`(signal > 5.0) & (cell_id > 0)`, matching all 101 frames; the plain threshold does not.
The name says as much. `>` and `>=` both match, because no site sits exactly at 5.0, so the
data cannot distinguish them.

The census marks this one **(weak)** rather than derivable, on purpose. `signal` has 50,201
distinct values over 64,000 sites, so a pairwise dependence test finds almost everything
"determined" by it. A joint dependence on two arrays is real and invisible to a pairwise
search, which is why the flag says to check directly rather than reporting a saving.

## Disposition 1: sparse voxel DAG, no

A DAG compresses **binary geometry**. Binary geometry is 0.266% of this payload. Deleting it
outright, which no DAG can beat, removes 0.77 MiB of 290 MiB. `gzip -9` removes 217.7 MiB
today in a format the container already supports.

Three things make this unlikely to change. The fraction is **scale-invariant**, since every
array grows as N³. The paper leaves **material attributes and dynamic updates out of scope**,
and those are precisely the 68% that dominates. And the compression already exists where size
actually binds: the HTML viewer delta-encodes label changes and gzips them, reaching roughly
100× on the label class in about two lines.

Revisit only if the payload becomes geometry-dominated, if an attribute-carrying variant is
on the table (HashDAG supports edits and attributes; its implementation requires CUDA, so it
is an architectural reference rather than an Apple Silicon dependency), or if navigation is
ever *measured* as the bottleneck. Nothing in this repository measures it that way today.

## Disposition 2: signed distance field, not as storage, yes as an observable

Not as a representation. This lattice is a **42-way exact partition**, every site belonging
to exactly one parcel plus wall and medium sentinels; a distance field is intrinsically
two-phase. One field per parcel at Float32 is 10.9 MB per frame against the 256 KB the Int32
label array costs, and multiphase level sets buy vacuum-and-overlap pathologies at triple
junctions in exchange.

The exactness argument is stronger than the size one. The volume constraint is
`λ_V(V − V_target)²` on **integer** site counts, connectivity is checked exactly including the
disconnected-daughter refusal, and a Metropolis move is a label copy between neighbouring
sites. Reconstruction from a distance field is approximate exactly at boundary sites, which
are the only sites a copy attempt ever touches, so the round trip would be lossy precisely
where the model does all its work.

As a **derived observable**, it is worth computing. Two open questions here are both
approximating depth into the aggregate with poor proxies: distance from a centroid, which is
a bad proxy in a ragged multi-parcel aggregate, and the binding benchmark's 6-connected
occupied-neighbour fraction, which is a one-voxel-radius probe taking six distinct values
over 3,038 sites with 930 and 895 piled at 5/6 and 6/6. A true signed distance to the
occupied-set boundary is what both are groping toward: smooth, far more dynamic range, and it
would make the matched-total-capacity control considerably more discriminating. An exact
Euclidean transform is linear per axis and 64,000 sites is nothing.

**Keep it out of rendering.** A smooth isosurface implies sub-lattice spatial resolution, and
the pitch is a declared refusal here, D-PITCH blocked and unmeasured. Blocky voxel rendering
is honest about a resolution the calibration does not have.

## Does dropping the static arrays make anything faster? No.

`bench_reduction.py` re-emits the trajectory twice through the **same writer**, differing
only in array inventory, because comparing a reduced form against the original tier would
confound dropping arrays with swapping writers. The reduced variant reconstructs the full
ten-array state and is timed with that reconstruction included, since a smaller answer is
not a faster one.

```
  full       289.84 MiB   101 files
  reduced    124.62 MiB   103 files   (57.0% smaller)

timings, 15 samples each, warm cache (the repeat-viewing case)
  full trajectory, 10 arrays         median   54.7 ms   IQR  53.8- 65.9   min  53.0
  reduced + reconstruction           median   53.9 ms   IQR  50.7- 58.3   min  49.1
  median ratio 1.02x, but 7 of 15 reduced samples exceed the full median
  sample ranges OVERLAP -- these samples do not separate the two

  reconstruction exact on all 101 frames x 10 arrays: YES

  raw read, uncompressed             median   50.4 ms   IQR  49.8- 60.8   min  48.6
  gzip decompress, in memory         median  351.8 ms   IQR 346.8-361.9   min 345.0
```

The block above is the committed run's output with two labels corrected after review. The
second-to-last line was printed as "distributions separate: NO -- the difference is not
distinguishable from noise"; overlapping ranges from 15 samples show only that the samples
overlap, and the script now says exactly that. The last line was printed as "gzip decompress
on read ... compression makes a READ 7.0x slower": at that commit the compressed samples
decompressed buffers already in memory while the raw samples opened files, so the 7.0x
compared CPU against I/O. The script now reads the compressed files from disk in both
paths. **That comparison has not been rerun**: the 101-frame tier is not on the machine
that made this correction, so no read-against-read ratio is claimed here.

**The reduction is a storage argument and not a performance one.** 57% fewer bytes, and the
load time does not separate from noise in 15 samples. Reading 290 MiB from warm cache is
50 ms, so halving it saves nothing anyone can perceive.

**Compression is a per-read cost.** It is quoted as a 75% saving; decompressing 290 MiB
takes 352 ms of CPU against a 50 ms warm-cache read, and that is the floor of what a
compressed read costs, before its own file I/O. Write it once at 10.8 s if disk is the
constraint; do not pay it per read.

That is the honest ranking at this scale. Dropping the four static arrays is worth doing
because the bytes are meaningless, not because anything is slow. Nothing here is slow.

## Where the action lives

The exporter half of this is filed as **issue #34**, not fixed here: the change touches
`export_vti.jl` on `feat/lattice-viewer` (PR #25, on hold), and this diagnostic is on a
branch off `master`. The issue ranks three fixes by cost. The free one, with no change to
what a consumer must do, is to stop writing `generation` and `accumulated_dose_Gy` while
they are identically zero: 25.5% of every frame carrying no information at all.

## The constraint on all of it

The existing 101-frame tier **cannot be rewritten**, in either direction. 209 artifact hashes
are pinned in `render_manifest.json` and those pins are what the receipt chain rests on. Every
finding here applies to new runs only.

## Running it

```sh
python3 diagnostics/payload_census/census.py <file.vti | directory> [--json out.json]
python3 diagnostics/payload_census/test_census.py    # 14 tests, no data
```

numpy only; VTK deliberately absent, since an audit of the renderer's output should not need
the renderer's toolchain. About 42 s over 101 frames, dominated by the pairwise dependence
search.

`test_census.py` builds its own `.vti` files, header and appended binary both, so it takes no
data path. What it pins, each found by mutating `census.py` and watching for a green suite:
static detection must compare frames rather than one frame to itself, and a single frame is
never static; the weak-source threshold must survive; `PointData` sits inside `<Piece>` and
must not be counted per cell; **independent arrays must not be reported as dependent**,
without which a dependence test that answers "yes" unconditionally passes everything else
in the file; a lookup that changes between frames is not a lookup; `Int64` values past 2^53
are compared exactly, not through `Float64`; a mutually derivable pair is counted as
removable once, not twice; and a frame whose extent, dtype or array length disagrees with
the file's own header or with frame 0 is refused rather than divided by the wrong count.

`census.json` was written before the `removable` flag existed and does not carry it. Its
fractions are unaffected: `species` and `lineage_id` derive from `cell_id`, which is kept,
so nothing in that run was double-counted.

## Files

| | |
|---|---|
| `census.py` | per-array bytes, distinct values, constant, static, discovered dependences |
| `test_census.py` | 14 data-free tests |
| `census.json` | receipt of the run above |
| `bench_reduction.py` | does the reduction make anything faster? (no) |
