"""Real-OpenMC-API contract tests (review amendment 10): a mocked module can
exercise our control flow but cannot validate the actual API. These import
the pinned openmc and round-trip the generated model WITHOUT running
transport or needing nuclear data. Skipped (never passed) when openmc is
not importable."""

import numpy as np
import pytest

openmc = pytest.importorskip("openmc")

from biofilm_openmc.materials import BIOMASS, MEDIUM, voxel_class_array
from biofilm_openmc.model import build_model

pytestmark = pytest.mark.openmc_api


def test_model_structure(snapshot, config):
    model = build_model(snapshot, config)

    assert model.settings.photon_transport is True
    assert model.settings.run_mode == "fixed source"
    assert model.settings.batches == 5
    src = model.settings.source[0]
    assert src.particle == "photon"

    heating = [t for t in model.tallies if t.name == "heating"]
    assert len(heating) == 1
    assert heating[0].scores == ["heating"]

    # material dedup: unique material objects bounded by the class table
    mats = model.geometry.get_all_materials()
    assert len(mats) <= len(config.materials)


def test_lattice_orientation_against_probes(snapshot, config):
    model = build_model(snapshot, config)
    lattices = [l for l in model.geometry.get_all_lattices().values()
                if l.name == "biofilm_voxels"]
    assert len(lattices) == 1
    lat = lattices[0]

    classes = voxel_class_array(snapshot, config)  # logical (x,y,z)
    n = classes.shape[0]
    # universes is (z, y-inverted, x): check at asymmetric probe points
    for x, y, z in [(1, 2, 3), (4, 1, 5), (6, 5, 1), (0, 0, 0)]:
        expected = classes[x, y, z]
        got = lat.universes[z][n - 1 - y][x].name
        assert got == expected, f"lattice orientation broken at {(x, y, z)}"


def test_model_xml_roundtrip(tmp_path, snapshot, config):
    model = build_model(snapshot, config)
    xml = tmp_path / "model.xml"
    model.export_to_model_xml(str(xml))
    assert xml.exists() and xml.stat().st_size > 0

    reread = openmc.Model.from_model_xml(str(xml))
    assert reread.settings.photon_transport is True
    assert len(reread.geometry.get_all_lattices()) == 1
