#!/usr/bin/env python3
"""Does dropping the static arrays actually make anything faster?

    python3 diagnostics/payload_census/bench_reduction.py <source .vti dir> <scratch dir> [--reps 5]

The census found 44.7% of every frame is bit-identical across the trajectory and another
10.6% is exactly derivable. That is a storage argument. Whether it is also a SPEED argument
is a separate question, and this measures it rather than assuming it.

Method, and the trap it avoids. Comparing a reduced trajectory against the ORIGINAL tier
would confound two changes: dropping arrays, and swapping the writer. So both variants are
re-emitted here by the same writer from the same in-memory arrays, and only the array
inventory differs. The source tier is read and never written.

  full/     all ten arrays, every frame
  reduced/  the three time-varying arrays per frame, plus one static.vti written once

`reduced` is only a fair comparison if it can reconstruct the full state, so it does, and
the reconstruction is timed as part of its cost. Anything else would be measuring a smaller
answer, not a faster one.
"""
import gzip, os, re, statistics, struct, sys, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from census import read_vti, DT

VTK = {np.dtype(v).name: k for k, v in DT.items()}
STATIC = ("generation", "interior_mask", "radiation_cpm", "accumulated_dose_Gy")
DERIVED = ("lineage_id", "species", "occupied_above_threshold")
KEPT = ("cell_id", "melanin", "signal")


def write_vti(path, cells, arrays):
    blob, offs = bytearray(), {}
    for n, a in arrays.items():
        offs[n] = len(blob)
        b = np.ascontiguousarray(a).tobytes()
        blob += struct.pack("<Q", len(b)) + b
    ext = f"0 {cells[0]} 0 {cells[1]} 0 {cells[2]}"
    decls = "".join(f'<DataArray type="{VTK[np.asarray(a).dtype.name]}" Name="{n}" '
                    f'format="appended" offset="{offs[n]}"/>' for n, a in arrays.items())
    head = (f'<?xml version="1.0"?>\n<VTKFile type="ImageData" header_type="UInt64">\n'
            f'<ImageData WholeExtent="{ext}">\n<Piece Extent="{ext}">'
            f'<CellData>{decls}</CellData></Piece>\n</ImageData>\n'
            f'<AppendedData encoding="raw">_')
    with open(path, "wb") as f:
        f.write(head.encode() + bytes(blob) + b'\n</AppendedData>\n</VTKFile>\n')


def reconstruct(kept, static, threshold):
    """Rebuild the three derivable arrays. Timed as part of the reduced variant's cost."""
    cid = kept["cell_id"]
    out = dict(kept)
    out.update(static)
    out["lineage_id"] = np.where(cid > 0, cid, 0).astype(np.int32)
    # species is a lookup on cell_id; the table rides in the static sidecar in a real
    # implementation, and is rebuilt here from the same information.
    out["species"] = static["_species_lut"][np.clip(cid, 0, None)].astype(np.uint8)
    out["occupied_above_threshold"] = ((kept["signal"] > threshold) & (cid > 0)).astype(np.uint8)
    return out


def build(src_files, scratch):
    full_d, red_d = os.path.join(scratch, "full"), os.path.join(scratch, "reduced")
    for d in (full_d, red_d):
        if os.path.exists(d):
            raise SystemExit(f"destination already exists: {d}")
        os.makedirs(d)
    cells = None
    lut = None
    for i, p in enumerate(src_files):
        arrays, fields, n, _, _ = read_vti(p)
        a = {k: v[0] for k, v in arrays.items()}
        if cells is None:
            side = round(n ** (1 / 3))
            cells = (side, side, side)
            # one lookup table for the whole trajectory, from the first frame
            lut = np.zeros(int(a["cell_id"].max()) + 1, dtype=np.uint8)
            occ = a["cell_id"] > 0
            lut[a["cell_id"][occ]] = a["species"][occ]
            write_vti(os.path.join(red_d, "static.vti"), cells,
                      {k: a[k] for k in STATIC})
            np.save(os.path.join(red_d, "species_lut.npy"), lut)
        else:
            o = a["cell_id"] > 0
            lut[a["cell_id"][o]] = a["species"][o]
        write_vti(os.path.join(full_d, f"f{i:06d}.vti"), cells, a)
        write_vti(os.path.join(red_d, f"f{i:06d}.vti"), cells, {k: a[k] for k in KEPT})
    np.save(os.path.join(red_d, "species_lut.npy"), lut)
    return full_d, red_d, cells


def dirbytes(d):
    return sum(os.path.getsize(os.path.join(d, f)) for f in os.listdir(d))


def load_full(d):
    fs = sorted(f for f in os.listdir(d) if f.endswith(".vti"))
    return [{k: v[0] for k, v in read_vti(os.path.join(d, f))[0].items()} for f in fs]


def load_reduced(d, threshold=5.0):
    static = {k: v[0] for k, v in read_vti(os.path.join(d, "static.vti"))[0].items()}
    static["_species_lut"] = np.load(os.path.join(d, "species_lut.npy"))
    fs = sorted(f for f in os.listdir(d) if f.endswith(".vti") and f != "static.vti")
    out = []
    for f in fs:
        kept = {k: v[0] for k, v in read_vti(os.path.join(d, f))[0].items()}
        out.append(reconstruct(kept, static, threshold))
    return out


def samples(fn, reps):
    """Timing samples, warmed. Returns them sorted, so the caller can see the spread.

    A median alone hides whether two distributions actually separate. At this scale they
    mostly do not, and reporting one number would manufacture a result out of noise.
    """
    r = fn()
    out = []
    for _ in range(reps):
        t0 = time.perf_counter()
        fn()
        out.append(time.perf_counter() - t0)
    return sorted(out), r


def show(name, s):
    q1, q3 = s[len(s) // 4], s[3 * len(s) // 4]
    print(f"  {name:34s} median {statistics.median(s)*1000:6.1f} ms   "
          f"IQR {q1*1000:5.1f}-{q3*1000:5.1f}   min {s[0]*1000:5.1f}")


def main(argv):
    src, scratch = argv[0], argv[1]
    reps = int(argv[argv.index("--reps") + 1]) if "--reps" in argv else 5
    files = sorted(os.path.join(src, f) for f in os.listdir(src) if f.endswith(".vti"))
    print(f"source: {len(files)} frames, {dirbytes(src)/2**20:.1f} MiB (read only, never written)")
    full_d, red_d, cells = build(files, scratch)

    fb, rb = dirbytes(full_d), dirbytes(red_d)
    print(f"\n  full     {fb/2**20:8.2f} MiB   {len(os.listdir(full_d))} files")
    print(f"  reduced  {rb/2**20:8.2f} MiB   {len(os.listdir(red_d))} files   "
          f"({100*(1-rb/fb):.1f}% smaller)")

    print(f"\ntimings, {reps} samples each, warm cache (the repeat-viewing case)")
    sf, af = samples(lambda: load_full(full_d), reps)
    sr, ar = samples(lambda: load_reduced(red_d), reps)
    show("full trajectory, 10 arrays", sf)
    show("reduced + reconstruction", sr)
    mf, mr = statistics.median(sf), statistics.median(sr)
    slower = sum(1 for x in sr if x > mf)
    sep = sr[-1] < sf[0] or sf[-1] < sr[0]
    print(f"  median ratio {mf/mr:.2f}x, but {slower} of {len(sr)} reduced samples exceed the "
          f"full median")
    print(f"  distributions separate: {'yes' if sep else 'NO -- the difference is not '
          'distinguishable from noise at this scale'}")

    # Is the reconstruction faithful? A faster wrong answer is not an answer.
    bad = []
    for i, (x, y) in enumerate(zip(af, ar)):
        for k in x:
            if k.startswith("_"):
                continue
            if not np.array_equal(x[k], y[k]):
                bad.append((i, k))
    print(f"\n  reconstruction exact on all {len(af)} frames x {len(af[0])} arrays: "
          f"{'YES' if not bad else 'NO -> ' + str(bad[:5])}")

    # What a READER pays. Compression is quoted as a saving; on read it is a cost, and
    # this is the comparison that says which dominates.
    blobs = [open(os.path.join(full_d, f), "rb").read() for f in sorted(os.listdir(full_d))]
    gz = [gzip.compress(b, 9) for b in blobs]
    sraw, _ = samples(lambda: [open(os.path.join(full_d, f), "rb").read()
                               for f in sorted(os.listdir(full_d))], reps)
    sdec, _ = samples(lambda: [gzip.decompress(b) for b in gz], reps)
    print()
    show("raw read, uncompressed", sraw)
    show("gzip decompress on read", sdec)
    ratio = statistics.median(sdec) / statistics.median(sraw)
    saved = 100 * (1 - sum(len(b) for b in gz) / sum(len(b) for b in blobs))
    print(f"  compression makes a READ {ratio:.1f}x slower, to save {saved:.0f}% of disk")

    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
