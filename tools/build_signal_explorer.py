#!/usr/bin/env python3
"""Package verified labels and accepted events without inventing intermediate frames.

Signal colour values are quantized only for display (0..10 / 255); endpoint
numbers are the unrounded Julia HDF5 results. The full field stays in HDF5/VTI.
"""
import argparse
import base64
import gzip
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def build(parent, derived, output, template):
    parent, derived, output = map(Path, (parent, derived, output))
    pm = json.loads((parent / "run_manifest.json").read_text())
    dm = json.loads((derived / "derived_manifest.json").read_text())
    assert digest(parent / "run_manifest.json") == dm["parent_manifest_sha256"]
    for root, manifest in ((parent, pm), (derived, dm)):
        for name, expected in manifest["artifacts"].items():
            assert digest(root / name) == expected, name
    assert [r["mcs"] for r in pm["snapshots"]] == list(range(101))
    chunks, specs, cursor = [], {}, 0

    def add(name, arr, dtype):
        nonlocal cursor
        a = np.asarray(arr, dtype=dtype).ravel()
        b = a.tobytes()
        specs[name] = dict(offset=cursor, length=len(a), dtype=dtype)
        chunks.append(b)
        cursor += len(b)

    frames, signal, reference_mask = [], [], None
    for t in range(101):
        with h5py.File(parent / f"snapshots/snap_mcs{t:06}.h5") as f:
            assert f.attrs["dataset_axis_order_h5py"] in ("zyx", b"zyx")
            ids = f["lattice/cell_id"][:].ravel()
            mask = f["lattice/interior_mask"][:].ravel().astype(bool)
            if reference_mask is not None:
                assert np.array_equal(reference_mask, mask)
            reference_mask = mask
            assert f.attrs["mcs"] == t
            assert set(ids[ids > 0]) == set(range(1, 43))
            if t == 0:
                mapping = np.zeros(43, dtype="u1")
                mapping[f["cells/id"][:]] = f["cells/species"][:]
            frames.append(np.where(mask, ids, 255).astype("u1"))
        with h5py.File(derived / f"fields/signal_mcs{t:06}.h5") as f:
            assert f.attrs["mcs"] == t
            a = f["fields/signal"][:].ravel()
            assert np.min(a) >= 0 and np.max(a) <= 10
            signal.append(np.rint(a[(ids > 0) & mask] * 255 / 10).astype("u1"))
    add("initial", frames[0], "u1")
    add("species", mapping, "u1")
    changes, offsets = [], [0]
    for t in range(1, 101):
        sites = np.flatnonzero(frames[t] != frames[t-1])
        changes.extend(((sites.astype("u4") << 6) | frames[t][sites]).tolist())
        offsets.append(len(changes))
    add("changes", changes, "<u4")
    add("offsets", offsets, "<u4")
    add("signal", np.concatenate(signal), "u1")
    add("signal_offsets", np.cumsum([0] + [len(a) for a in signal]), "<u4")
    # Verify every ordered event against both pre-copy labels and every endpoint.
    with h5py.File(parent / "accepted_copies.h5") as f:
        assert f.attrs["linear_index_base"] == 1
        records = {k: f["accepted"][k][:] for k in f["accepted"]}
        floats = {"adh", "vol", "rad", "mel", "delta_h", "draw"}
        for key, values in records.items():
            add("event_" + key, values, "<f8" if key in floats else "<u4")
        replay, index = frames[0].copy(), 0
        event_offsets = [0, 0]
        for t in range(1, 101):
            last_proposal = 0
            while index < len(records["mcs"]) and records["mcs"][index] == t:
                source = records["donor_site"][index] - 1
                target = records["recipient_site"][index] - 1
                assert replay[source] == records["donor_id"][index]
                assert replay[target] == records["recipient_id"][index]
                assert records["proposal_index"][index] > last_proposal
                last_proposal = records["proposal_index"][index]
                replay[target] = replay[source]
                index += 1
            assert np.array_equal(replay, frames[t]), t
            event_offsets.append(index)
        assert index == 14281
        add("event_offsets", event_offsets, "<u4")
    metadata = dict(specs=specs, metrics=json.loads((derived / "metrics.json").read_text()),
                    parent_manifest_sha256=dm["parent_manifest_sha256"],
                    derived_manifest_sha256=digest(derived / "derived_manifest.json"),
                    initial_parcels=42, seed=42, frames=101, accepted_copies=14281,
                    display_signal_error_bound=10/510,
                    frame_sha256=[hashlib.sha256(a.tobytes()).hexdigest() for a in frames])
    packed = base64.b64encode(gzip.compress(b"".join(chunks), compresslevel=9, mtime=0)).decode()
    html = Path(template).read_text().replace("__METADATA__", json.dumps(metadata, separators=(",", ":")))
    html = html.replace("__PAYLOAD__", packed)
    assert "__PAYLOAD__" not in html and "__METADATA__" not in html
    assert len(html.encode()) < 1_000_000, len(html.encode())
    output.write_text(html)
    print(json.dumps(dict(bytes=output.stat().st_size, exact_frames=101, replayed_copies=index,
                          output=str(output), sha256=digest(output))))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("parent"); p.add_argument("derived"); p.add_argument("output")
    p.add_argument("--template", default=str(Path(__file__).parents[1] / "viewer/signal_explorer.html"))
    args = p.parse_args()
    build(args.parent, args.derived, args.output, args.template)
