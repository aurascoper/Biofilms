"""Reference A0 correctness, through the builder rather than around it.

`test_attenuation.py` already checks OpenMC's physics on hand-built geometry.
These check OUR chain: build_water_phantom_model -> heating tally ->
phantom_mass_kg -> per-source normalization. A phantom at dosimetry scale is
what makes that checkable at all — at cell scale a 661.7 keV photon crosses the
whole domain without interacting, and there is nothing to validate.
"""

from __future__ import annotations

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
from biofilm_openmc.dose import EV_TO_J, per_source_from_statepoint
from biofilm_openmc.materials import phantom_mass_kg
from biofilm_openmc.mesh import phantom_mesh_extent_cm, resolve_mesh_dimension
from biofilm_openmc.model import build_water_phantom_model

from conftest import WATER_PHANTOM_CONFIG

E_SRC_EV = 661655.0


def _cfg(toml: str):
    return load_transport_config(toml, kind=WATER_PHANTOM)


def _run(model, tmp_path):
    sp_path = model.run(cwd=str(tmp_path), output=False)
    return openmc.StatePoint(sp_path)


def _mass_and_dimension(cfg):
    dimension = resolve_mesh_dimension(cfg.mesh_base_dimension,
                                       cfg.mesh_coarsening_factor)
    return phantom_mass_kg(cfg, dimension, phantom_mesh_extent_cm(cfg)), dimension


def test_energy_balance_through_the_builder(tmp_path):
    """With reflective boundaries nothing escapes, so the summed mesh heating
    must recover the emitted energy per source particle.

    This validates the builder's GEOMETRY, TALLY and SOURCE: a closed geometry,
    a mesh that covers all the material, and a spectrum that emits what it
    claims. It deliberately does not validate the mass denominator — in
    `sum(field * mass)` the mass cancels exactly. The denominator is pinned
    separately below and, for the coarsening trap, by tests/test_mesh.py.
    """
    cfg = _cfg(WATER_PHANTOM_CONFIG
               .replace('axial = "vacuum"', 'axial = "reflective"')
               .replace('radial_outer = "vacuum"', 'radial_outer = "reflective"')
               .replace("batches = 5", "batches = 10")
               .replace("particles = 2000", "particles = 5000"))

    sp = _run(build_water_phantom_model(cfg), tmp_path)
    mass_kg, dimension = _mass_and_dimension(cfg)
    result = per_source_from_statepoint(sp, mass_kg)

    # specific energy (Gy/source) x mass (kg) = joules per source particle
    deposited_J = float((result.field * mass_kg).sum())
    emitted_J = E_SRC_EV * EV_TO_J
    assert abs(deposited_J - emitted_J) / emitted_J < 0.02, (
        f"deposited {deposited_J:.6e} J/source vs emitted {emitted_J:.6e}")

    # The denominator, pinned independently: the mass array must describe the
    # tallied cube at water density, so a bin volume taken from anything other
    # than extent/dimension is caught here rather than silently rescaling dose.
    ex, ey, ez = phantom_mesh_extent_cm(cfg)
    expected_total_kg = 1.0 * (ex * ey * ez) * 1e-3      # 1.0 g/cm3 -> kg
    assert np.isclose(mass_kg.sum(), expected_total_kg, rtol=1e-12)
    assert mass_kg.shape == tuple(dimension)


def test_dose_falls_off_with_radius(tmp_path):
    """Vacuum boundaries, point source on the axis: the field must be centrally
    peaked and decrease outward. A broad sanity check on placement and
    orientation, not a benchmark."""
    cfg = _cfg(WATER_PHANTOM_CONFIG
               .replace("batches = 5", "batches = 10")
               .replace("particles = 2000", "particles = 5000"))

    sp = _run(build_water_phantom_model(cfg), tmp_path)
    mass_kg, dimension = _mass_and_dimension(cfg)
    field = per_source_from_statepoint(sp, mass_kg).field

    nx, ny, nz = dimension
    extent = phantom_mesh_extent_cm(cfg)
    # radial distance of each bin centre from the cylinder axis
    xs = (np.arange(nx) + 0.5) * extent[0] / nx - extent[0] / 2.0
    ys = (np.arange(ny) + 0.5) * extent[1] / ny - extent[1] / 2.0
    r = np.sqrt(xs[:, None, None] ** 2 + ys[None, :, None] ** 2)
    r = np.broadcast_to(r, field.shape)

    edges = np.linspace(0, extent[0] / 2.0, 4)
    shells = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (r >= lo) & (r < hi) & (field > 0)
        assert sel.any(), f"no scoring bins in shell [{lo:.1f}, {hi:.1f}) cm"
        shells.append(field[sel].mean())

    assert shells[0] > shells[1] > shells[2], (
        f"dose not decreasing with radius: {shells}")


def test_coarsening_preserves_total_deposited_energy(tmp_path):
    """Physical answer invariant, resolution free: the same geometry tallied at
    two nested resolutions must deposit the same total energy."""
    fine = _cfg(WATER_PHANTOM_CONFIG
                .replace('axial = "vacuum"', 'axial = "reflective"')
                .replace('radial_outer = "vacuum"', 'radial_outer = "reflective"')
                .replace("batches = 5", "batches = 10")
                .replace("particles = 2000", "particles = 5000"))
    coarse = _cfg(WATER_PHANTOM_CONFIG
                  .replace('axial = "vacuum"', 'axial = "reflective"')
                  .replace('radial_outer = "vacuum"', 'radial_outer = "reflective"')
                  .replace("batches = 5", "batches = 10")
                  .replace("particles = 2000", "particles = 5000")
                  .replace("base_dimension = [12, 12, 12]",
                           "base_dimension = [12, 12, 12]\n  coarsening_factor = 2"))
    assert coarse.mesh_coarsening_factor == 2

    totals = []
    for i, cfg in enumerate((fine, coarse)):
        sp = _run(build_water_phantom_model(cfg), tmp_path / f"run{i}")
        mass_kg, _ = _mass_and_dimension(cfg)
        field = per_source_from_statepoint(sp, mass_kg).field
        totals.append(float((field * mass_kg).sum()))

    assert np.isclose(totals[0], totals[1], rtol=0.02), (
        f"total deposited energy moved with resolution: {totals}")
