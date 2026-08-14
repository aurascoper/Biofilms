"""Shared fixtures: a synthetic Julia-layout transport snapshot, a complete
demo config (synthetic TEST values — not physical claims), and duck-typed
fake statepoints matching `dose.py`'s access pattern."""

from __future__ import annotations

import h5py
import numpy as np
import pytest

# Synthetic, clearly-demo physical values for tests ONLY. The shipped
# config/coupling_template.toml has every REQUIRED key unset on purpose.
VALID_CONFIG = """
schema_version = 1

[geometry]
voxel_pitch_cm = 0.001
origin_cm = [0.0, 0.0, 0.0]
cylinder_radius_cm = 0.004
cylinder_length_cm = 0.008
membrane_thickness_cm = 0.001

[boundaries]
axial = "vacuum"
radial_outer = "vacuum"

[time]
seconds_per_mcs = 1.0

[source]
photons_per_second = 1.0e9
spectrum_energies_eV = [1.0e6]
spectrum_probabilities = [1.0]
spatial = "line_z_axis"
angular = "isotropic"

[normalization]
hamiltonian_scale = 1.0
melanin_scale = 1.0
membrane_statistic = "mass_weighted"

[transport]
# cell-scale voxels make heating extremely sparse (~3e-4 interactions per
# 1 MeV history in an 8-um water cube) — enough histories that the smoke
# run scores with overwhelming probability
batches = 5
particles = 50000
seed = 7

[materials.medium]
density_g_cm3 = 1.0
  [materials.medium.elements]
  H = 0.111894
  O = 0.888106

[materials.baseline_biomass]
density_g_cm3 = 1.1
  [materials.baseline_biomass.elements]
  H = 0.10
  C = 0.50
  N = 0.10
  O = 0.30

[materials.membrane]
density_g_cm3 = 0.94
  [materials.membrane.elements]
  C = 0.856
  H = 0.144
"""


def make_snapshot_file(path, n=8, n_lineages=3, config_toml=""):
    """Write a synthetic transport snapshot in the JULIA storage layout:
    logical (x,y,z) arrays land in the file so h5py sees dims reversed —
    exactly what export_checkpoint.jl produces."""
    rng = np.random.default_rng(1)
    cell_id = np.zeros((n, n, n), dtype=np.int32)          # logical (x,y,z)
    # a few axis-asymmetric blobs so orientation errors are detectable
    cell_id[1:3, 2:5, 3:4] = 1
    cell_id[4:6, 1:2, 5:7] = 2
    cell_id[6:7, 5:7, 1:2] = 3
    lineage = np.where(cell_id > 0,
                       rng.integers(1, n_lineages + 1, size=cell_id.shape),
                       0).astype(np.int32)
    lineage[cell_id > 0] = ((cell_id[cell_id > 0] - 1) % n_lineages) + 1
    species = np.where(cell_id > 0, 6, 0).astype(np.int32)
    generation = np.zeros_like(cell_id)
    interior = np.ones_like(cell_id, dtype=np.uint8)

    probes = []
    for xyz in [(1, 2, 3), (4, 1, 5), (6, 5, 1), (0, 0, 0)]:
        probes.append([*xyz, int(cell_id[xyz])])
    probes = np.array(probes, dtype=np.int64)

    def put(f, name, logical):
        f[name] = np.ascontiguousarray(logical.transpose(2, 1, 0))

    with h5py.File(path, "w") as f:
        f.attrs["schema_version"] = 1
        f.attrs["logical_axis_order"] = "xyz"
        f.attrs["dataset_axis_order_h5py"] = "zyx"
        f.attrs["coordinate_index_base"] = 0
        f.attrs["cell_id_background"] = 0
        f.attrs["cell_id_wall"] = -1
        f.attrs["git_sha"] = "test"
        f.attrs["julia_version"] = "test"
        f.attrs["mcs"] = 5
        f.attrs["physical_time_s"] = float("nan")
        f.attrs["label_state_hash"] = "0" * 64
        f.attrs["material_class_source"] = "absent" if not config_toml else "config"
        f["config_toml"] = config_toml
        put(f, "lattice/cell_id", cell_id)
        put(f, "lattice/species_id", species)
        put(f, "lattice/lineage_id", lineage)
        put(f, "lattice/generation", generation)
        put(f, "lattice/interior_mask", interior)
        put(f, "fields/radiation_cpm", np.ones_like(cell_id, dtype=float))
        put(f, "fields/melanin", np.zeros_like(cell_id, dtype=float))
        put(f, "dose/accumulated_Gy", np.zeros_like(cell_id, dtype=float))
        ids = np.array([1, 2, 3], dtype=np.int64)
        f["cells/id"] = ids
        f["cells/species"] = np.full(3, 6, dtype=np.int64)
        f["cells/volume"] = np.array(
            [int((cell_id == i).sum()) for i in ids], dtype=np.int64)
        f["cells/lineage_id"] = ((ids - 1) % n_lineages) + 1
        f["cells/parent_id"] = np.zeros(3, dtype=np.int64)
        f["cells/generation"] = np.zeros(3, dtype=np.int64)
        f["cells/birth_mcs"] = np.zeros(3, dtype=np.int64)
        f["cells/lifetime_dose_Gy"] = np.zeros(3, dtype=float)
        # Julia writes (P,4); h5py view is (4,P)
        f["orientation_probes"] = probes.T
    return path


@pytest.fixture
def snapshot_path(tmp_path):
    return make_snapshot_file(tmp_path / "snap.h5")


@pytest.fixture
def config():
    from biofilm_openmc.config import load_config
    return load_config(VALID_CONFIG)


@pytest.fixture
def snapshot(snapshot_path):
    from biofilm_openmc.snapshot import load_snapshot
    return load_snapshot(snapshot_path)


class FakeTally:
    def __init__(self, mean, std_dev):
        self.mean = mean
        self.std_dev = std_dev


class FakeStatepoint:
    """Duck-type of openmc.StatePoint for dose.py: get_tally(name=...)."""

    def __init__(self, tallies):
        self._tallies = tallies

    def get_tally(self, name):
        return self._tallies[name]


def heating_statepoint_from_logical(logical_eV, sd_logical_eV=None):
    """Flatten a logical (x,y,z) heating field the way an OpenMC mesh tally
    stores bins (x fastest, z slowest) and wrap it as a FakeStatepoint."""
    flat = logical_eV.transpose(2, 1, 0).reshape(-1, 1)
    sd = (sd_logical_eV if sd_logical_eV is not None
          else np.zeros_like(logical_eV)).transpose(2, 1, 0).reshape(-1, 1)
    return FakeStatepoint({"heating": FakeTally(flat, sd)})
