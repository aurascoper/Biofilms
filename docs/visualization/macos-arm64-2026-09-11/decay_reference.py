"""Decay reference: the smallest defensible isotope addition.

Species labels beside a HYPOTHETICAL normalized activity distribution fading by the
documented decay law, on a FROZEN consortium snapshot.

    a(x, t) = a0(x) * 2^(-t / t_half),   t in DAYS

**Why the geometry is frozen.** There is no established MCS->seconds conversion in this
repository; `seconds_per_mcs` ships NaN and the exporter refuses a physical pitch
(`D-PITCH` is awaiting_measurement). Advancing biological rearrangement alongside isotope
decay would silently introduce that mapping as an assumption. Freezing the snapshot keeps
the isotope clock independent, in days, and answerable.

**What this does.** It demonstrates the decay law on this geometry.
**What this does not do.** It computes NO binding and NO dose. There is no source term,
no dose point kernel, and no transport. Nothing here establishes anything biological.

**a0(x) is an INPUT, never a result.** The default is a declared uniform distribution over
the occupied interior voxels of the frozen frame, normalized to sum to 1. It is a
hypothetical placeholder chosen for being obviously arbitrary. No measurement supports it.

**On naming.** The file, the fields and the receipt say "decay reference". Naming the
isotope is legitimate here and only here, because the only thing used is a pinned, citable
physical constant -- the half-life, decaying to a stable daughter -- and nothing else. The
repository's source-term gate leaves isotope identity unestablished and the material path
is elemental, not isotopic; this file does not disturb that.

usage: decay_reference.py <inert_signal_dir> <out_stem> [--mcs N] [--days D] [--step S]
                          [--receipt render_manifest.json]
"""
import sys, os, json, math, hashlib
import numpy as np
from vti_read import read_vti

# Pinned constant. The plan carried two values -- 6.6443 d in its sourced-constants
# section and 6.647 d in the decay equation. They are not interchangeable, and the plan's
# own decay constant settles it: ln2 / (6.6443 d * 86400 s/d) = 1.207431e-6 1/s, which
# matches the pinned lambda to six figures, while 6.647 d gives 1.206941e-6 and does not.
# 6.647 d is a stale figure. CONFIRM against LNHB/CEA or NNDC/ENSDF before any use that
# depends on the last digits; this arithmetic is a consistency check, not a source.
T_HALF_DAYS = 6.6443
T_HALF_UNCERTAINTY_DAYS = 0.0009
LAMBDA_PER_S = math.log(2) / (T_HALF_DAYS * 86400.0)
DAUGHTER = "stable"

NOTE = ("would be Lu-177 if the evidence audit and the source term land; "
        "the material path is elemental, not isotopic, and no isotope identity "
        "is established by this file")


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _vti_bytes(dims, cell_arrays, field_num, field_str):
    nx, ny, nz = dims
    ext = "0 %d 0 %d 0 %d" % (nx, ny, nz)
    blocks, decls_cell, decls_field, off = [], [], [], 0

    def add(name, arr, kind):
        nonlocal off
        b = arr.tobytes(order="F") if arr.ndim == 3 else arr.tobytes()
        blocks.append(np.uint64(len(b)).tobytes() + b)
        vtk_t = {"uint8": "UInt8", "int32": "Int32", "float64": "Float64",
                 "float32": "Float32"}[arr.dtype.name]
        # NumberOfTuples is mandatory on FieldData and only there: CellData infers its
        # tuple count from the extent, FieldData has no extent to infer from. Omitting it
        # makes VTK read zero tuples and drop the array -- ParaView then reports zero
        # field-data entries with no error, which is provenance that silently is not there.
        ntup = '' if kind == "cell" else 'NumberOfTuples="1" '
        decl = ('<DataArray type="%s" Name="%s" NumberOfComponents="1" %s'
                'format="appended" offset="%d"/>' % (vtk_t, name, ntup, off))
        (decls_cell if kind == "cell" else decls_field).append(decl)
        off += len(b) + 8

    for n, a in cell_arrays.items():
        add(n, a, "cell")
    for n, v in field_num.items():
        add(n, np.array([v], dtype=np.float64), "field")
    for n, v in field_str.items():
        # VTK stores a String array's payload NUL-terminated, and the UInt64 byte count
        # INCLUDES the terminator. Omitting it does not merely lose the string: it
        # desynchronises the appended section and VTK then reads zero cells from the
        # whole file, silently. Verified against the exporter's own output, where
        # units = b"lattice\x00" with a declared length of 8.
        b = v.encode("utf-8") + b"\x00"
        blocks.append(np.uint64(len(b)).tobytes() + b)
        decls_field.append('<Array type="String" Name="%s" NumberOfComponents="1" '
                           'NumberOfTuples="1" format="appended" offset="%d"/>' % (n, off))
        off += len(b) + 8

    head = ('<?xml version="1.0" encoding="utf-8"?>\n'
            '<VTKFile type="ImageData" version="1.0" byte_order="LittleEndian" '
            'header_type="UInt64">\n'
            '  <ImageData WholeExtent="%s" Origin="0.0 0.0 0.0" Spacing="1.0 1.0 1.0">\n'
            '    <Piece Extent="%s">\n'
            '      <CellData>\n        %s\n      </CellData>\n'
            '    </Piece>\n'
            '    <FieldData>\n      %s\n    </FieldData>\n'
            '  </ImageData>\n'
            '  <AppendedData encoding="raw">\n_'
            % (ext, ext, "\n        ".join(decls_cell), "\n      ".join(decls_field)))
    return head.encode("utf-8") + b"".join(blocks) + b"\n  </AppendedData>\n</VTKFile>\n"


def main(argv):
    root = os.path.abspath(argv[1]); stem = os.path.abspath(argv[2])
    getopt = lambda f, d: argv[argv.index(f) + 1] if f in argv else d
    frozen_mcs = int(getopt("--mcs", "0"))
    days = float(getopt("--days", "20")); step = float(getopt("--step", "0.5"))
    receipt_path = getopt("--receipt", None)
    if not math.isfinite(step) or step <= 0:
        raise SystemExit("FATAL: --step must be finite and positive, got %r" % step)
    if not math.isfinite(days) or days < 0:
        raise SystemExit("FATAL: --days must be finite and non-negative, got %r" % days)

    src = os.path.join(root, "paraview", "signal_mcs%06d.vti" % frozen_mcs)
    if not os.path.isfile(src):
        raise SystemExit("FATAL: no source frame at %s" % src)
    # `fields` used to be discarded as `_`, so the frame's own mcs, git_sha and
    # parent_manifest_sha256 were thrown away and the receipt recorded the COMMAND LINE
    # rather than a reading of the data. A frame misnamed by the exporter, or --mcs aimed
    # at a relabelled copy, propagated into every output without a murmur.
    a, fields, dims = read_vti(src)

    # --- identity: EXACT integral equality, not rounded ---
    embedded = fields.get("mcs")
    if embedded is None:
        raise SystemExit("FATAL: %s declares no `mcs` field; refusing to assert a frozen "
                         "MCS the artifact does not state" % src)
    emb = float(embedded)
    # `int(round(emb)) != frozen_mcs` accepted a frame declaring mcs=0.25 as MCS 0: round()
    # masked the disagreement it was supposed to detect. A frozen geometry must be frozen at
    # an actual MCS, so a non-integral or non-finite declaration is itself a refusal.
    if not math.isfinite(emb):
        raise SystemExit("FATAL: %s declares a non-finite mcs=%s" % (src, embedded))
    if emb != int(emb):
        raise SystemExit("FATAL: %s declares a non-integral mcs=%s; a frozen snapshot must "
                         "sit at an integer MCS" % (src, embedded))
    if int(emb) != frozen_mcs:
        raise SystemExit("FATAL: --mcs %d but %s declares mcs=%s" % (frozen_mcs, src, embedded))

    # --- geometry: compare what the writer will emit against what the source declares ---
    # The writer hard-codes Origin 0,0,0 and Spacing 1,1,1 (lattice sites). Anything else
    # must be refused rather than silently relabelled. Checking coordinate_index_base alone
    # did not check geometry at all, and the reader did not return these attributes.
    EMIT_ORIGIN, EMIT_SPACING = (0.0, 0.0, 0.0), (1.0, 1.0, 1.0)
    for name, got, want in (("Origin", fields.get("_vti_origin"), EMIT_ORIGIN),
                            ("Spacing", fields.get("_vti_spacing"), EMIT_SPACING)):
        if got is None:
            raise SystemExit("FATAL: %s declares no %s; refusing to guess it" % (src, name))
        if any(abs(g - w) > 1e-12 for g, w in zip(got, want)):
            raise SystemExit("FATAL: %s declares %s=%s but this writer emits %s. It will not "
                             "silently relabel the geometry; propagate it or pick another "
                             "source." % (src, name, got, want))
    for key, want in (("coordinate_index_base", 0.0),):
        got = fields.get(key)
        if got is not None and abs(float(got) - want) > 1e-12:
            raise SystemExit("FATAL: %s declares %s=%s; this writer emits %s and will not "
                             "silently relabel it" % (src, key, got, want))

    # --- the source bytes must match the bundle's manifest, and optionally a frozen pin ---
    # Recording the source's own sha256 captures what was read; it does not validate it.
    man_path = os.path.join(root, "derived_manifest.json")
    rel = os.path.relpath(src, root)
    manifest_ok = None
    if os.path.isfile(man_path):
        man = json.load(open(man_path))
        want_hash = man.get("artifacts", {}).get(rel)
        if want_hash is None:
            raise SystemExit("FATAL: %s is not a registered artifact in %s; refusing to "
                             "freeze a geometry the bundle does not account for"
                             % (rel, man_path))
        got_hash = _sha256(src)
        if got_hash != want_hash:
            raise SystemExit("FATAL: %s does not match its manifest hash\n  manifest %s\n  "
                             "actual   %s" % (rel, want_hash, got_hash))
        manifest_ok = want_hash
        # A manifest beside its own data is not an authority -- a consistently tampered
        # clone regenerates it for free. --receipt pins it to a frozen expectation.
        if receipt_path:
            if not os.path.isfile(receipt_path):
                raise SystemExit("FATAL: --receipt given but not readable: %s" % receipt_path)
            pin = json.load(open(receipt_path)).get("pinned", {})
            want_man = pin.get("derived_manifest_sha256")
            if want_man is None:
                raise SystemExit("FATAL: %s carries no pinned.derived_manifest_sha256"
                                 % receipt_path)
            got_man = _sha256(man_path)
            if got_man != want_man:
                raise SystemExit("FATAL: derived_manifest.json does not match the frozen "
                                 "receipt\n  pinned %s\n  actual %s" % (want_man, got_man))
    elif receipt_path:
        raise SystemExit("FATAL: --receipt supplied but no derived_manifest.json at %s"
                         % man_path)
    species, mask = a["species"], a["interior_mask"]
    occupied = (species > 0) & (mask == 1)
    n_occ = int(occupied.sum())
    if n_occ == 0:
        raise SystemExit("frozen frame MCS %d has no occupied interior voxels; "
                         "pick a frame after onset" % frozen_mcs)

    # a0(x): DECLARED, HYPOTHETICAL. Uniform over occupied interior voxels, sums to 1.
    a0 = np.zeros(dims, dtype=np.float64)
    a0[occupied] = 1.0 / n_occ

    outdir = os.path.dirname(stem) or "."
    # Refuse any existing run destination BEFORE writing a byte -- the same render-then-check
    # inversion that is a standing P2 at tools/render_label_trajectory.py:297. Checking only
    # for a non-empty directory would still let a re-run with a shorter --days leave the
    # previous run's tail frames on disk, unreferenced by the new .pvd but matching its glob.
    # The contract is the RUN DIRECTORY, not just this stem. An earlier revision refused
    # only artifacts whose name began with the stem, while its own comment claimed to refuse
    # any existing run destination -- so a directory holding a different run, or anything
    # else, was accepted and written into.
    existing = os.listdir(outdir) if os.path.isdir(outdir) else []
    if existing:
        same_stem = [f for f in existing if f.startswith(os.path.basename(stem))]
        raise SystemExit(
            "FATAL: output directory %s is not empty (%d entr%s, e.g. %s).%s\n"
            "This writer refuses any existing run destination before writing a byte. "
            "Point --out at a fresh directory."
            % (outdir, len(existing), "y" if len(existing) == 1 else "ies",
               sorted(existing)[0],
               ("\n  %d of them already carry this stem: %s"
                % (len(same_stem), sorted(same_stem)[0])) if same_stem else ""))
    os.makedirs(outdir, exist_ok=True)

    # Integer frame count. `int(days/step)+1` truncated whenever the quotient landed just
    # under an integer in binary -- int(0.3/0.1) == 2, so --days 0.3 --step 0.1 emitted
    # 0, 0.1, 0.2 and dropped the requested endpoint while printing "0 to 0.3".
    n_steps = int(round(days / step))
    if abs(n_steps * step - days) > 1e-9:
        raise SystemExit("FATAL: --days %g is not an integer multiple of --step %g "
                         "(nearest is %g); declare an endpoint the series can reach."
                         % (days, step, n_steps * step))
    times = [round(i * step, 12) for i in range(n_steps + 1)]
    assert abs(times[-1] - days) < 1e-9, "endpoint %g not included" % days
    # Unique filenames do not imply unique timesteps. An earlier revision rounded to six
    # decimals, so --days 0.000002 --step 0.0000004 wrote six distinct files advertising
    # the times [0, 0, 1e-6, 1e-6, 2e-6, 2e-6] -- a .pvd with duplicate timesteps and no
    # complaint. Validate the values that are actually emitted.
    if len(set(times)) != len(times):
        dup = sorted(t for t in set(times) if times.count(t) > 1)
        raise SystemExit("FATAL: --step %g is finer than the emitted time representation; "
                         "%d of %d timesteps collide (e.g. %g). Use a coarser step."
                         % (step, len(times) - len(set(times)), len(times), dup[0]))
    if any(times[i] >= times[i + 1] for i in range(len(times) - 1)):
        raise SystemExit("FATAL: emitted timesteps are not strictly increasing at --step %g"
                         % step)

    entries, metrics = [], []
    for idx, t in enumerate(times):
        surviving = 2.0 ** (-t / T_HALF_DAYS)
        act = a0 * surviving
        # Frame INDEX, not the formatted time. "%07.2f" cannot separate frames finer than
        # 0.01 d, so two .pvd entries named one file and the later frame overwrote the
        # earlier while the .pvd still advertised both timesteps.
        name = "%s_f%05d.vti" % (os.path.basename(stem), idx)
        open(os.path.join(os.path.dirname(stem), name), "wb").write(_vti_bytes(
            dims,
            {"species": species, "interior_mask": mask,
             "activity_fraction": act},
            {"decay_time_days": t, "frozen_mcs": float(frozen_mcs),
             "half_life_days": T_HALF_DAYS, "lambda_per_s": LAMBDA_PER_S,
             "surviving_fraction": surviving,
             "a0_occupied_voxels": float(n_occ)},
            {"what_this_is": "decay reference: a declared decay law on a frozen geometry",
             "what_this_is_not": "NO binding, NO dose, NO transport, NO source term",
             "a0_status": "HYPOTHETICAL INPUT, never a result: uniform over occupied "
                          "interior voxels of the frozen frame, normalized to sum 1",
             "activity_fraction_units": "dimensionless fraction of the initial inventory",
             "time_units": "days; the isotope clock is INDEPENDENT of MCS -- no "
                           "MCS-to-seconds conversion is established or implied",
             "geometry_units": "lattice sites; spacing 1.0/site, D-PITCH refusal intact",
             "daughter": DAUGHTER,
             "note": NOTE}))
        entries.append((t, name))
        metrics.append({"days": t, "surviving_fraction": surviving,
                        "total_activity": float(act.sum()),
                        "max_voxel_activity": float(act.max())})

    with open(stem + ".pvd", "w") as f:
        f.write('<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" '
                'byte_order="LittleEndian">\n  <Collection>\n')
        for t, name in entries:
            f.write('    <DataSet timestep="%s" group="" part="0" file="%s"/>\n' % (t, name))
        f.write("  </Collection>\n</VTKFile>\n")

    # Bind source, every output and the .pvd by CONTENT. Recording a mutable source path
    # and a git_sha string ties the receipt to a name, not to bytes.
    out_hashes = {name: _sha256(os.path.join(outdir, name)) for _, name in entries}
    out_hashes[os.path.basename(stem) + ".pvd"] = _sha256(stem + ".pvd")

    json.dump({"frozen_mcs": frozen_mcs, "source_frame": src,
               "source_frame_sha256": _sha256(src),
               "source_verified_against_manifest": manifest_ok,
               "manifest_pinned_by_receipt": receipt_path,
               "vti_origin": list(fields["_vti_origin"]),
               "vti_spacing": list(fields["_vti_spacing"]),
               "source_declared_mcs": float(embedded),
               "source_git_sha": fields.get("git_sha"),
               "source_parent_manifest_sha256": fields.get("parent_manifest_sha256"),
               "output_sha256": out_hashes,
               "days": days, "step": step, "endpoint_included": True,
               "half_life_days": T_HALF_DAYS,
               "half_life_uncertainty_days": T_HALF_UNCERTAINTY_DAYS,
               "lambda_per_s": LAMBDA_PER_S, "daughter": DAUGHTER,
               "half_life_provenance": "CONFIRM against LNHB/CEA or NNDC/ENSDF. Chosen as "
                                       "the value consistent with the pinned lambda "
                                       "1.20743e-6 1/s; 6.647 d is inconsistent with it.",
               "a0": "HYPOTHETICAL INPUT: uniform over %d occupied interior voxels, sums to 1"
                     % n_occ,
               "computes": "decay only", "does_not_compute": ["binding", "dose", "transport"],
               "time_axis": "days, independent of MCS",
               "note": NOTE, "frames": len(times), "metrics": metrics},
              open(stem + "_receipt.json", "w"), indent=2)

    print("frozen MCS          : %d" % frozen_mcs)
    print("occupied interior   : %d voxels (a0 uniform over these, sums to 1)" % n_occ)
    print("half-life           : %.4f +/- %.4f d   lambda = %.6e 1/s" %
          (T_HALF_DAYS, T_HALF_UNCERTAINTY_DAYS, LAMBDA_PER_S))
    print("frames              : %d, 0 to %g d in %g d steps" % (len(times), days, step))
    print("surviving at 1 t1/2 : %.6f   at 2: %.6f   at 3: %.6f" %
          tuple(2.0 ** (-k) for k in (1, 2, 3)))
    print("wrote               : %s.pvd + %d .vti + _receipt.json" % (stem, len(times)))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
