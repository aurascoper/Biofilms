#!/usr/bin/env python3
"""Data-free tests for the exporter audit.

    python3 diagnostics/payload_census/test_census.py

Every case builds its own .vti, header and appended binary both, in a temp directory. No
trajectory, no argument, no environment variable, so these assertions cannot be skipped by
a missing data bundle.
"""
import os, shutil, struct, sys, tempfile, unittest
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from census import read_vti, audit, DT

VTK = {np.dtype(v).name: k for k, v in DT.items()}


def make_vti(path, cells_xyz, cell_arrays, field_arrays=(), point_arrays=(), compressed=False):
    """A real .vti: XML header plus raw appended binary, which is what the audit reads."""
    nx, ny, nz = cells_xyz
    ext = f"0 {nx} 0 {ny} 0 {nz}"
    blob, offs = bytearray(), {}
    for name, arr in cell_arrays.items():
        offs[name] = len(blob)
        b = np.ascontiguousarray(arr).tobytes()
        blob += struct.pack("<Q", len(b)) + b
    def decls(items):
        return "".join(f'<DataArray type="{VTK[np.asarray(a).dtype.name]}" Name="{n}" '
                       f'format="appended" offset="{offs.get(n, 0)}"/>' for n, a in items)
    comp = ' compressor="vtkZLibDataCompressor"' if compressed else ""
    pd = f"<PointData>{decls(point_arrays)}</PointData>" if point_arrays else ""
    fd = f"<FieldData>{decls(field_arrays)}</FieldData>" if field_arrays else ""
    head = (f'<?xml version="1.0"?>\n<VTKFile type="ImageData"{comp} header_type="UInt64">\n'
            f'<ImageData WholeExtent="{ext}">\n{fd}\n'
            f'<Piece Extent="{ext}">{pd}<CellData>{decls(cell_arrays.items())}</CellData></Piece>\n'
            '</ImageData>\n<AppendedData encoding="raw">_')
    with open(path, "wb") as f:
        f.write(head.encode() + bytes(blob) + b'\n</AppendedData>\n</VTKFile>\n')
    return path


class TestReader(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def f(self, name="a.vti"):
        return os.path.join(self.d, name)

    def test_cells_from_whole_extent_and_values_round_trip(self):
        # WholeExtent counts POINTS; cells are one fewer per axis.
        v = np.arange(4 * 3 * 2, dtype=np.int32)
        arrays, _, cells, _, _ = read_vti(make_vti(self.f(), (4, 3, 2), {"x": v}))
        self.assertEqual(cells, 24)
        np.testing.assert_array_equal(arrays["x"][0], v)
        self.assertEqual(arrays["x"][1], "Int32")
        self.assertEqual(arrays["x"][2], 4)

    def test_field_and_point_data_are_not_cell_arrays(self):
        # FieldData is a sibling of <Piece>; PointData sits inside it. Counting either
        # per cell inflates bytes/site by a whole array's width.
        c = {"species": np.ones(8, dtype=np.uint8)}
        arrays, _, _, _, _ = read_vti(make_vti(self.f(), (2, 2, 2), c,
                                               field_arrays=[("mcs", np.zeros(1, np.float64))],
                                               point_arrays=[("nodal", np.zeros(27, np.float64))]))
        self.assertEqual(list(arrays), ["species"])

    def test_malformed_input_is_refused_not_guessed(self):
        with open(self.f("bad.vti"), "wb") as fh:
            fh.write(b'<VTKFile type="ImageData"></VTKFile>')
        with self.assertRaises(ValueError):
            read_vti(self.f("bad.vti"))                       # no AppendedData
        with self.assertRaises(ValueError):
            read_vti(make_vti(self.f("c.vti"), (2, 2, 2),
                              {"x": np.zeros(8, np.uint8)}, compressed=True))
        with self.assertRaises(ValueError):
            read_vti(make_vti(self.f("z.vti"), (0, 2, 2), {"x": np.zeros(0, np.uint8)}))


class TestAudit(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def build(self, n_frames=3, cells=(4, 4, 4)):
        """Three frames: one array static, one constant, one varying, one a lookup."""
        N = cells[0] * cells[1] * cells[2]
        rng = np.random.default_rng(0)
        files = []
        static = rng.integers(0, 5, N).astype(np.float64)
        for t in range(n_frames):
            label = (rng.integers(1, 4, N)).astype(np.int32)     # 3 distinct labels
            files.append(make_vti(os.path.join(self.d, f"f{t}.vti"), cells, {
                "static_field": static,
                "always_zero": np.zeros(N, np.int32),
                "label": label,
                "label_squared": (label * label).astype(np.int32),  # a lookup on label
                "noise": rng.random(N),                             # near-unique
                # Two low-cardinality arrays that are INDEPENDENT by construction: every
                # one of the 3x3 value combinations occurs, so neither determines the other.
                "grid_a": (np.arange(N) % 3 + 1).astype(np.int32),
                "grid_b": (np.arange(N) // 3 % 3 + 1).astype(np.int32),
            }))
        return files

    def test_static_constant_and_distinct_counts(self):
        rows, cells, tot, gz, n = audit(self.build())
        by = {r["name"]: r for r in rows}
        self.assertEqual(n, 3)
        self.assertEqual(cells, 64)
        self.assertTrue(by["static_field"]["static_across_frames"])
        self.assertFalse(by["label"]["static_across_frames"])
        self.assertTrue(by["always_zero"]["constant_in_frame"])
        self.assertTrue(by["always_zero"]["static_across_frames"])
        self.assertFalse(by["label"]["constant_in_frame"])
        self.assertEqual(by["always_zero"]["distinct_max"], 1)
        self.assertEqual(by["label"]["distinct_max"], 3)
        self.assertGreater(gz, 0)
        self.assertGreater(tot, gz)          # this fixture is compressible

    def test_functional_dependence_is_discovered_not_declared(self):
        rows = {r["name"]: r for r in audit(self.build())[0]}
        src = [d["array"] for d in (rows["label_squared"]["derivable_from"] or [])]
        self.assertIn("label", src)
        # and the relation is not spuriously symmetric-free: label is determined by its
        # square here too, since squaring is injective on 1..3.
        self.assertIn("label_squared", [d["array"] for d in (rows["label"]["derivable_from"] or [])])

    def test_a_near_unique_source_is_marked_weak(self):
        # `noise` has ~64 distinct values over 64 sites, so everything trivially "depends"
        # on it. That must be labelled, or the byte total counts an artefact as a saving.
        rows = {r["name"]: r for r in audit(self.build())[0]}
        for d in rows["label"]["derivable_from"] or []:
            if d["array"] == "noise":
                self.assertFalse(d["informative"])
                self.assertGreater(d["distinct"], 64 // 8)
                break
        else:
            self.fail("expected `noise` to appear as a (weak) source for `label`")
        for d in rows["label_squared"]["derivable_from"] or []:
            if d["array"] == "label":
                self.assertTrue(d["informative"])

    def test_independent_arrays_are_not_reported_as_dependent(self):
        # Without this, a dependence test that says "yes" unconditionally passes every
        # other assertion in this file. A mutation check found exactly that hole.
        rows = {r["name"]: r for r in audit(self.build())[0]}
        a_src = [d["array"] for d in (rows["grid_a"]["derivable_from"] or [])]
        b_src = [d["array"] for d in (rows["grid_b"]["derivable_from"] or [])]
        self.assertNotIn("grid_b", a_src)
        self.assertNotIn("grid_a", b_src)
        # and the genuinely dependent pair is still found, so this is not vacuous
        self.assertIn("label", [d["array"] for d in rows["label_squared"]["derivable_from"]])

    def test_a_constant_array_claims_no_source(self):
        rows = {r["name"]: r for r in audit(self.build())[0]}
        self.assertIsNone(rows["always_zero"]["derivable_from"])

    def test_inventory_must_match_across_frames(self):
        files = self.build(n_frames=2)
        make_vti(files[1], (4, 4, 4), {"label": np.zeros(64, np.int32)})   # different inventory
        with self.assertRaises(ValueError):
            audit(files)


if __name__ == "__main__":
    unittest.main(verbosity=2)
