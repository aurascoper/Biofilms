"""Transport + heating benchmarks (review amendment 11), correctly split:

- TRANSMISSION: Beer-Lambert describes uncollided flux, not heating. We tally
  the uncollided partial current (energy filter pins the un-scattered line)
  behind water slabs of several thicknesses and check ln(I) vs t is linear —
  self-consistent exponential attenuation without importing an external mu.
- HEATING ENERGY BALANCE: in a closed reflective water box every source eV
  must end up as deposited heating (photon transport + TTB), so the global
  heating tally must equal the source energy within tolerance.

Both need the executable + nuclear data; skipped (never passed) otherwise.
"""

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

E_SRC = 1.0e6  # eV


def _water():
    m = openmc.Material(name="water")
    m.set_density("g/cm3", 1.0)
    m.add_element("H", 0.111894, percent_type="wo")
    m.add_element("O", 0.888106, percent_type="wo")
    return m


def _run(model, tmp_path):
    sp_path = model.run(cwd=str(tmp_path), output=False)
    return openmc.StatePoint(sp_path)


def test_transmission_is_exponential(tmp_path):
    """Uncollided partial current through water thickness t must satisfy
    ln I(t) linear in t (single-line source; energy filter keeps only the
    unscattered line)."""
    thicknesses = [1.0, 2.0, 3.0]  # cm
    logs = []
    for i, t in enumerate(thicknesses):
        water = _water()
        slab_in = openmc.XPlane(x0=0.0, boundary_type="vacuum")
        slab_a = openmc.XPlane(x0=1.0)
        slab_b = openmc.XPlane(x0=1.0 + t)
        detector = openmc.XPlane(x0=1.0 + t + 1.0)
        world_out = openmc.XPlane(x0=1.0 + t + 2.0, boundary_type="vacuum")
        side = openmc.model.RightCircularCylinder(
            (-0.5, 0.0, 0.0), 5.0 + t, 50.0, axis="x",
            boundary_type="reflective")

        cells = [
            openmc.Cell(region=+slab_in & -slab_a & -side.cyl),
            openmc.Cell(fill=water, region=+slab_a & -slab_b & -side.cyl),
            openmc.Cell(region=+slab_b & -world_out & -side.cyl),
        ]
        geometry = openmc.Geometry(cells)

        settings = openmc.Settings()
        settings.run_mode = "fixed source"
        settings.photon_transport = True
        settings.batches = 10
        settings.particles = 20000
        settings.seed = 11 + i
        settings.source = openmc.IndependentSource(
            space=openmc.stats.Point((0.5, 0.0, 0.0)),
            angle=openmc.stats.Monodirectional((1.0, 0.0, 0.0)),
            energy=openmc.stats.Discrete([E_SRC], [1.0]),
            particle="photon")

        surf_filter = openmc.SurfaceFilter(detector)
        # keep only the unscattered line: Compton losses push scattered
        # photons out of this bin
        e_filter = openmc.EnergyFilter([0.999 * E_SRC, 1.001 * E_SRC])
        current = openmc.Tally(name="uncollided")
        current.filters = [surf_filter, e_filter]
        current.scores = ["current"]

        model = openmc.Model(geometry=geometry, settings=settings,
                             tallies=openmc.Tallies([current]))
        sp = _run(model, tmp_path / f"t{i}")
        val = float(sp.get_tally(name="uncollided").mean.ravel()[0])
        assert val > 0, "no uncollided transmission scored — beam geometry broken"
        logs.append(np.log(val))

    # ln I vs t linear: successive differences agree (exponential law)
    d1 = logs[0] - logs[1]
    d2 = logs[1] - logs[2]
    assert d1 > 0 and d2 > 0, "transmission must decrease with thickness"
    assert abs(d1 - d2) / d1 < 0.15, \
        f"ln-attenuation not linear: steps {d1:.4f} vs {d2:.4f}"


def test_heating_energy_balance(tmp_path):
    """Closed reflective water box: total heating per source particle must
    equal the source energy within tolerance."""
    water = _water()
    box = openmc.model.RectangularParallelepiped(
        -30, 30, -30, 30, -30, 30, boundary_type="reflective")
    geometry = openmc.Geometry([openmc.Cell(fill=water, region=-box)])

    settings = openmc.Settings()
    settings.run_mode = "fixed source"
    settings.photon_transport = True
    settings.batches = 10
    settings.particles = 5000
    settings.seed = 3
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Point((0.0, 0.0, 0.0)),
        angle=openmc.stats.Isotropic(),
        energy=openmc.stats.Discrete([E_SRC], [1.0]),
        particle="photon")

    heating = openmc.Tally(name="total_heating")
    heating.scores = ["heating"]

    model = openmc.Model(geometry=geometry, settings=settings,
                         tallies=openmc.Tallies([heating]))
    sp = _run(model, tmp_path)
    tally = sp.get_tally(name="total_heating")
    deposited = float(tally.mean.ravel()[0])           # eV per source particle
    rel_err = float(tally.std_dev.ravel()[0]) / deposited

    print(f"\nenergy balance: deposited {deposited:.1f} eV / source {E_SRC:.1f} eV"
          f" (rel_err {rel_err:.2e})")
    assert abs(deposited - E_SRC) / E_SRC < 0.02, \
        "closed-box heating does not balance source energy"
