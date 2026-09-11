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
"""
import sys, os, json, math
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

    src = os.path.join(root, "paraview", "signal_mcs%06d.vti" % frozen_mcs)
    a, _, dims = read_vti(src)
    species, mask = a["species"], a["interior_mask"]
    occupied = (species > 0) & (mask == 1)
    n_occ = int(occupied.sum())
    if n_occ == 0:
        raise SystemExit("frozen frame MCS %d has no occupied interior voxels; "
                         "pick a frame after onset" % frozen_mcs)

    # a0(x): DECLARED, HYPOTHETICAL. Uniform over occupied interior voxels, sums to 1.
    a0 = np.zeros(dims, dtype=np.float64)
    a0[occupied] = 1.0 / n_occ

    os.makedirs(os.path.dirname(stem) or ".", exist_ok=True)
    times = [round(i * step, 6) for i in range(int(days / step) + 1)]
    entries, metrics = [], []
    for t in times:
        surviving = 2.0 ** (-t / T_HALF_DAYS)
        act = a0 * surviving
        name = "%s_d%07.2f.vti" % (os.path.basename(stem), t)
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

    json.dump({"frozen_mcs": frozen_mcs, "source_frame": src,
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
