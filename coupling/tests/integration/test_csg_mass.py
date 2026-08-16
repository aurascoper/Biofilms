"""The mass denominator, measured against a shape whose volume is known exactly.

`voxel_mass_kg` and `phantom_mass_kg` give every tally bin the full rectangular
volume of one material. A0's report attributes its ~11% coarse-mesh aggregate
drift to that approximation and names `mesh.material_volumes()` as the upgrade
path. This is that path, checked where the right answer is pi*r^2*L rather than
an opinion.
"""

from __future__ import annotations

import math
import os
from pathlib import Path

import numpy as np
import pytest

openmc = pytest.importorskip("openmc")

_XS = os.environ.get("OPENMC_CROSS_SECTIONS", "")
pytestmark = [
    pytest.mark.openmc,
    pytest.mark.skipif(not _XS or not Path(_XS).exists(),
                       reason="OPENMC_CROSS_SECTIONS not available"),
]

from biofilm_openmc.config import WATER_PHANTOM, load_transport_config
from biofilm_openmc.materials import mesh_material_masses_kg, phantom_mass_kg
from biofilm_openmc.mesh import phantom_mesh_extent_cm
from biofilm_openmc.model import build_water_phantom_model

_A0 = Path(__file__).resolve().parents[3] / "config" / "reference_a0_water_phantom.toml"
_N_SAMPLES = 2_000_000


@pytest.fixture(scope="module")
def phantom():
    cfg = load_transport_config(_A0, kind=WATER_PHANTOM)
    model = build_water_phantom_model(cfg)
    mesh = model.tallies[0].filters[0].mesh
    return cfg, model, mesh, tuple(mesh.dimension)


def test_exact_mass_recovers_the_analytic_cylinder(phantom):
    cfg, model, mesh, dim = phantom
    exact = mesh_material_masses_kg(mesh, model, dim, n_samples=_N_SAMPLES,
                                    max_materials=2)
    analytic_kg = (math.pi * cfg.cylinder_radius_cm ** 2
                   * cfg.cylinder_length_cm * 1.0 * 1e-3)
    # Raytracing is an ESTIMATE: the tolerance is sized from n_samples, not
    # asserted as equality. 1e-3 is generous at 2e6 rays.
    assert exact.sum() == pytest.approx(analytic_kg, rel=1e-3)
    assert exact.shape == dim


def test_the_full_bin_denominator_overstates_by_exactly_four_over_pi(phantom):
    """Not an unexplained discrepancy — the full-bin mass is the CIRCUMSCRIBING
    CUBE's mass, and cube/cylinder is 4/pi. Naming the number is what makes the
    bias a known quantity rather than a mystery in an aggregate."""
    cfg, model, mesh, dim = phantom
    full = phantom_mass_kg(cfg, dim, phantom_mesh_extent_cm(cfg))
    analytic_kg = (math.pi * cfg.cylinder_radius_cm ** 2
                   * cfg.cylinder_length_cm * 1.0 * 1e-3)
    assert full.sum() / analytic_kg == pytest.approx(4.0 / math.pi, rel=1e-6)


def test_bins_are_empty_full_or_partially_cut(phantom):
    """The three populations the full-bin denominator cannot distinguish.

    Corner bins are void and it calls them water. Bins clipped by the cylinder
    wall hold a FRACTION of a bin's water and it calls them full. Only the
    interior bins are ones it gets right — which is why the error is worst
    exactly where the mesh is coarsest and the boundary layer is thickest.
    """
    cfg, model, mesh, dim = phantom
    exact = mesh_material_masses_kg(mesh, model, dim, n_samples=_N_SAMPLES,
                                    max_materials=2)
    full_bin_kg = phantom_mass_kg(cfg, dim, phantom_mesh_extent_cm(cfg)).flat[0]

    empty = exact == 0
    partial = (exact > 0) & (exact < full_bin_kg * (1 - 1e-6))
    interior = exact >= full_bin_kg * (1 - 1e-6)
    assert empty.any() and partial.any() and interior.any()
    assert empty.sum() + partial.sum() + interior.sum() == exact.size

    # The MASS fraction is pi/4 even though the fraction of bins holding any
    # mass is larger — a partially cut bin counts once but weighs less.
    assert exact.sum() / (full_bin_kg * exact.size) == pytest.approx(
        math.pi / 4, rel=1e-3)
    assert (exact > 0).mean() > math.pi / 4


def test_the_raytrace_wrapper_does_not_touch_the_transport_model(phantom):
    """The void box is a measuring instrument. If it leaked into the model the
    boundary condition would move off the cylinder and the physics would change.
    """
    cfg, model, mesh, dim = phantom
    before = len(model.geometry.get_all_cells())
    mesh_material_masses_kg(mesh, model, dim, n_samples=100_000, max_materials=2)
    assert len(model.geometry.get_all_cells()) == before
    assert not any(c.name == "raytrace_void"
                   for c in model.geometry.get_all_cells().values())
