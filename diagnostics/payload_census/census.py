#!/usr/bin/env python3
"""An exporter audit for a .vti trajectory: what is written, and what need not be.

    python3 diagnostics/payload_census/census.py <file.vti | directory> [--json out.json]

This started as a compression argument and became something more useful. Byte shares alone
say which array is big; they do not say which array carries information. Three further
columns do: how many distinct values an array holds, whether it is constant inside a frame,
and whether it is bit-identical across frames. An array that is static across a trajectory
is being written once per frame and read once per frame for no reason, and no compressor
can fix a file that should not have been written.

The report also looks for functional dependence: if every site sharing a value of Y shares
the same value of X, then X is a lookup on Y and is derivable rather than measured. That
search is generic, so it finds relations without being told them.

numpy only. VTK is deliberately absent: an audit of the renderer's output should not need
the renderer's toolchain. The array reader below is minimal by design, not a copy of a
fuller one, because this file has to read what it audits.
"""
import gzip, hashlib, json, os, re, sys
import numpy as np

DT = {"Int8": np.int8, "UInt8": np.uint8, "Int16": np.int16, "UInt16": np.uint16,
      "Int32": np.int32, "UInt32": np.uint32, "Int64": np.int64, "UInt64": np.uint64,
      "Float32": np.float32, "Float64": np.float64}


def read_vti(path):
    """Cell arrays and field data from an uncompressed appended-binary .vti."""
    with open(path, "rb") as f:
        raw = f.read()
    cut = raw.find(b"<AppendedData")
    if cut < 0:
        raise ValueError(f"{path}: no AppendedData; this reader handles the appended form only")
    head = raw[:cut].decode("utf-8", "replace")
    if "compressor=" in head:
        raise ValueError(f"{path}: compressed .vti; this reader handles the raw form only")
    # The block-length prefix is read as 8 bytes little-endian below; a file declaring
    # anything else would be decoded misaligned rather than refused.
    if 'header_type="UInt64"' not in head:
        raise ValueError(f"{path}: header_type is not UInt64; this reader handles no other")
    if 'byte_order="BigEndian"' in head:
        raise ValueError(f"{path}: BigEndian; this reader handles little-endian only")
    m = re.search(r'WholeExtent="([^"]+)"', head)
    if not m:
        raise ValueError(f"{path}: no WholeExtent")
    e = [int(v) for v in m.group(1).split()]
    nx, ny, nz = e[1] - e[0], e[3] - e[2], e[5] - e[4]
    cells = nx * ny * nz
    if cells <= 0:
        raise ValueError(f"{path}: WholeExtent {e} encloses no cells")
    marker = raw.index(b"_", cut) + 1

    def grab(dtype, off):
        w = np.dtype(dtype).itemsize
        n = int(np.frombuffer(raw, "<u8", count=1, offset=marker + off)[0])
        # Flooring n // w would accept a truncated block and drop its tail silently.
        if n % w or marker + off + 8 + n > len(raw):
            raise ValueError(f"{path}: block at offset {off} declares {n} bytes, "
                             f"not a whole number of {w}-byte items inside the file")
        return np.frombuffer(raw, dtype, count=n // w, offset=marker + off + 8)

    # Only <CellData> is per-site. FieldData is a sibling of <Piece> and PointData sits
    # inside it; counting either per cell inflates bytes/site by a whole array's width.
    cd = re.search(r"<CellData.*?</CellData>", head, re.S)
    if not cd:
        raise ValueError(f"{path}: no CellData")
    arrays = {}
    for a in re.finditer(r'<DataArray type="(\w+)" Name="(\w+)"[^>]*offset="(\d+)"', cd.group(0)):
        t, n, off = a.group(1), a.group(2), int(a.group(3))
        if t not in DT:
            raise ValueError(f"{path}: unknown DataArray type {t}")
        arr = grab(DT[t], off)
        if len(arr) != cells:
            raise ValueError(f"{path}: cell array {n} holds {len(arr)} values for {cells} cells")
        arrays[n] = (arr, t, np.dtype(DT[t]).itemsize)
    if not arrays:
        raise ValueError(f"{path}: no CellData arrays")

    fields = {}
    fd = re.search(r"<FieldData>(.*?)</FieldData>", head, re.S)
    if fd:
        for a in re.finditer(r'<DataArray type="(\w+)" Name="(\w+)"[^>]*offset="(\d+)"', fd.group(1)):
            if a.group(1) in DT:
                fields[a.group(2)] = float(grab(DT[a.group(1)], int(a.group(3)))[0])
    return arrays, fields, cells, len(raw), raw, (nx, ny, nz)


def audit(files):
    frames, sizes, gz = [], [], []
    for p in files:
        arrays, fields, cells, nbytes, raw, dims = read_vti(p)
        frames.append((arrays, fields, cells, dims))
        sizes.append(nbytes)
        gz.append(len(gzip.compress(raw, 9)))
    names = list(frames[0][0])
    cells, dims = frames[0][2], frames[0][3]
    # Names alone let a frame with another extent or another dtype through, and every
    # bytes/site figure below would then be divided by the first frame's cell count. The
    # three dimensions are compared, not their product: 2x4x8 and 4x4x4 both hold 64
    # cells and are not site-for-site aligned when flattened.
    schema = [(n, frames[0][0][n][1]) for n in names]
    for i, (arrays, _, _, d) in enumerate(frames[1:], 1):
        if [(n, arrays[n][1]) for n in arrays] != schema or d != dims:
            raise ValueError(f"frame {i} ({files[i]}): inventory, types or extent "
                             f"differ from frame 0")
    frames = [(a, f, c) for a, f, c, _ in frames]

    rows = []
    for n in names:
        v0, t, w = frames[0][0][n]
        hs = {hashlib.sha256(np.ascontiguousarray(a[n][0]).tobytes()).hexdigest() for a, _, _ in frames}
        distinct = max(len(np.unique(a[n][0])) for a, _, _ in frames)
        rows.append({"name": n, "type": t, "bytes_per_site": w, "distinct_max": int(distinct),
                     "constant_in_frame": bool(distinct == 1),
                     # One frame has nothing to be static against; the flag was true for
                     # every array of a single-file audit and the whole frame read as waste.
                     "static_across_frames": len(frames) > 1 and len(hs) == 1})

    # Functional dependence: X is derivable from Y when no two sites sharing a Y value
    # disagree about X, across the WHOLE trajectory. The pairs of every frame are pooled
    # before counting, so a lookup that changes between frames is not a lookup; and the
    # pair keeps both native dtypes, since a Float64 cast merges Int64 values past 2^53.
    def determines(other, name):
        y = np.concatenate([np.ascontiguousarray(a[other][0]).ravel() for a, _, _ in frames])
        x = np.concatenate([np.ascontiguousarray(a[name][0]).ravel() for a, _, _ in frames])
        pair = np.empty(len(y), dtype=[("y", y.dtype), ("x", x.dtype)])
        pair["y"], pair["x"] = y, x
        return len(np.unique(pair)) == len(np.unique(y))

    for r in rows:
        if r["constant_in_frame"]:
            r["derivable_from"] = None
            continue
        src = []
        for other in names:
            if other == r["name"]:
                continue
            if determines(other, r["name"]):
                # A near-unique source determines everything trivially: if almost every
                # site has its own Y value, no two sites can disagree about X. Record the
                # source cardinality so a reader can tell a real lookup from an artefact,
                # and keep only the informative ones in the byte total.
                card = max(len(np.unique(a[other][0])) for a, _, _ in frames)
                src.append({"array": other, "distinct": int(card),
                            "informative": bool(card <= cells // 8)})
        r["derivable_from"] = src or None

    # What can actually be dropped. An array is removable only while an informative
    # source of it is still kept, largest first; otherwise a mutual pair (the fixture's
    # label and label_squared) was counted twice although one of them must stay.
    kept = set(names)
    for r in sorted(rows, key=lambda r: (-r["bytes_per_site"], r["name"])):
        srcs = [d["array"] for d in (r["derivable_from"] or []) if d["informative"]]
        r["removable"] = not r["static_across_frames"] and any(s in kept for s in srcs)
        if r["removable"]:
            kept.discard(r["name"])
    return rows, cells, sum(sizes), sum(gz), len(files)


def report(files, out_json=None):
    rows, cells, tot, tot_gz, n = audit(files)
    per = sum(r["bytes_per_site"] for r in rows)
    print(f"{n} frame(s), {cells:,} cells each, {per} bytes/site\n")
    print(f"{'array':26s} {'type':8s} {'B/site':>6s} {'%frame':>7s} {'distinct':>9s} "
          f"{'const':>6s} {'static':>7s}  derivable from")
    for r in sorted(rows, key=lambda r: -r["bytes_per_site"]):
        print(f"{r['name']:26s} {r['type']:8s} {r['bytes_per_site']:6d} "
              f"{100*r['bytes_per_site']/per:6.1f}% {r['distinct_max']:9d} "
              f"{'yes' if r['constant_in_frame'] else '':>6s} "
              f"{'yes' if r['static_across_frames'] else '':>7s}  "
              f"{', '.join(d['array'] + ('' if d['informative'] else ' (weak)') for d in r['derivable_from']) if r['derivable_from'] else ''}")

    static = sum(r["bytes_per_site"] for r in rows if r["static_across_frames"])
    deriv = sum(r["bytes_per_site"] for r in rows if r["removable"])
    occ = cells * n / 8
    print(f"\n  on disk                         {tot/2**20:9.2f} MiB")
    print(f"  gzip -9                         {tot_gz/2**20:9.2f} MiB  ({tot/tot_gz:.2f}x)")
    print(f"  static across frames            {100*static/per:8.1f}% of every frame, "
          f"written {n}x, needed once" + ("" if n > 1 else "  (one frame: not assessable)"))
    print(f"  derivable from another array    {100*deriv/per:8.1f}%")
    print(f"  neither                         {100*(per-static-deriv)/per:8.1f}%  <- the real payload")
    print(f"  binary occupancy, 1 bit/site    {100*(1/8)/per:8.3f}%  <- the DAG ceiling "
          f"({occ/2**20:.2f} MiB)")
    weak = [r["name"] for r in rows if r["derivable_from"]
            and not any(d["informative"] for d in r["derivable_from"])]
    if weak:
        print(f"\n  marked (weak): {', '.join(weak)} -- the only source that determines these")
        print(f"  is near-unique, so the pairwise test is uninformative. A JOINT dependence on")
        print(f"  two arrays is real but invisible to a pairwise search; check it directly.")
    doc = {"frames": n, "cells": cells, "bytes_per_site": per, "on_disk_bytes": tot,
           "gzip_bytes": tot_gz, "gzip_ratio": tot / tot_gz,
           "static_fraction": static / per, "derivable_fraction": deriv / per,
           "irreducible_fraction": (per - static - deriv) / per,
           "occupancy_fraction_of_frame": (1 / 8) / per, "arrays": rows}
    if out_json:
        with open(out_json, "w") as f:
            json.dump(doc, f, indent=2)
            f.write("\n")
        print(f"\nwrote {out_json}")
    return doc


def main(argv):
    if not argv:
        raise SystemExit("usage: census.py <file.vti | directory> [--json out.json]")
    target = argv[0]
    out = argv[argv.index("--json") + 1] if "--json" in argv else None
    files = (sorted(os.path.join(target, f) for f in os.listdir(target) if f.endswith(".vti"))
             if os.path.isdir(target) else [target])
    if not files:
        raise SystemExit(f"no .vti files under {target}")
    report(files, out)


if __name__ == "__main__":
    main(sys.argv[1:])
