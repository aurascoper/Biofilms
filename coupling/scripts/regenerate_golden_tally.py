#!/usr/bin/env python
"""Regenerate the golden-tally fixture for the gate-composition test.

    micromamba run -n openmc-biofilms python coupling/scripts/regenerate_golden_tally.py

Requires OPENMC_CROSS_SECTIONS pointing at real nuclear data (see
docs/openmc_stack.md). Runs REAL OpenMC transport -- 12 small runs (2 outer
draws x 3 replicates x {baseline, feedback}) on the same water-phantom
geometry docs/a0_transport_resolution_decision.md and test_water_phantom.py
already validate against closed-form physics -- and pins the raw per-source
heating tally each run actually produced. Nothing downstream of the raw
tally is pinned: coupling/tests/test_gate_composition.py
recomputes mass, dose, effect and verdict LIVE from these numbers using the
real repo code, in the ordinary (no-OpenMC) CI tier, the same
compare-only-in-CI split as tests/contract_csv.jl uses for the serial
fixture.

Feedback condition: density x1.35, the same DENSITY_SCALE lever
coupling/scripts/openmc_nested_pilot.py already uses -- not a new
perturbation invented for this fixture.

GUARD THE DESTINATION, matching writes_canonical_tables/refuse_partial_publish
in openmc_nested_pilot.py: refuses to overwrite the committed fixture unless
every one of the 12 runs actually completed.
"""

from __future__ import annotations

import json
import math
import os
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "coupling"))
sys.path.insert(0, str(REPO / "contract"))

FIXTURE_PATH = (REPO / "coupling" / "tests" / "fixtures"
                / "golden_tally_water_phantom.json")

N_OUTER = 2
N_REPLICATES = 3
DENSITY_BASELINE = 1.0
DENSITY_FEEDBACK = 1.0 * 1.35  # DENSITY_SCALE["lever_density_x1.35"]

CONFIG_TEMPLATE = """
schema_version = 1

[provenance]
reference_system_id = "golden_tally_water_phantom"
system_provenance = "declared"
evidence_policy = "measured_only"
execution_class = "reference_benchmark"
target_calibration = false

[geometry]
origin_cm = [0.0, 0.0, 0.0]
cylinder_radius_cm = 15.0
cylinder_length_cm = 30.0

[boundaries]
axial = "vacuum"
radial_outer = "vacuum"

[source]
photons_per_second = 1.0e9
spectrum_energies_eV = [661655.0]
spectrum_probabilities = [1.0]
spatial = "point_origin"
angular = "isotropic"

[transport]
batches = 5
particles = 2000
seed = {seed}

  [transport.mesh]
  base_dimension = [8, 8, 8]

[materials.medium]
density_g_cm3 = {density}
  [materials.medium.elements]
  H = 0.111894
  O = 0.888106
"""


# SIGNIFICANT FIGURES, NOT FULL REPR -- and the number is chosen against the
# tally's own uncertainty. `json.dumps` writes 16 digits of a heating value
# whose Monte Carlo relative standard deviation, measured on this very fixture,
# has a median of 27% per cell (min 3.8%, max 100% in the empty corners). Those
# trailing digits are not measurement, and pinning them makes the committed file
# assert a precision it does not have -- while requiring every one of them to
# reproduce on another machine for `git diff --exit-code` to pass.
#
# The precedent this fixture cites does exactly this: tests/fixtures/
# serial_seed42.csv pins 5 decimals, and its header says a different Julia minor
# version may legitimately change the bytes. ROUNDING IS THE TOLERANCE. Six
# figures sits ~5 orders of magnitude below the noise, so nothing physical is
# lost and the comparison stops depending on the last ulp.
SIGNIFICANT_FIGURES = 6


def _round_sig(x: float) -> float:
    if not math.isfinite(x) or x == 0.0:
        return x
    return round(x, SIGNIFICANT_FIGURES - 1 - int(math.floor(math.log10(abs(x)))))


def _round_nested(a):
    """Elementwise over an ndarray, returned as plain lists for JSON."""
    return np.vectorize(_round_sig)(np.asarray(a, dtype=float)).tolist()


def _nuclear_data_id() -> str:
    from openmc_nested_pilot import _nuclear_data_id as impl
    return impl()


def _run_one(density: float, seed: int, tmp_dir: Path):
    import openmc
    from biofilm_openmc.config import WATER_PHANTOM, load_transport_config
    from biofilm_openmc.dose import extract_heating
    from biofilm_openmc.materials import phantom_mass_kg
    from biofilm_openmc.mesh import phantom_mesh_extent_cm, resolve_mesh_dimension
    from biofilm_openmc.model import build_water_phantom_model

    toml = CONFIG_TEMPLATE.format(seed=seed, density=density)
    cfg = load_transport_config(toml, kind=WATER_PHANTOM)
    model = build_water_phantom_model(cfg)
    sp_path = model.run(cwd=str(tmp_dir), output=False)
    sp = openmc.StatePoint(sp_path)

    dimension = resolve_mesh_dimension(cfg.mesh_base_dimension,
                                       cfg.mesh_coarsening_factor)
    extent = phantom_mesh_extent_cm(cfg)
    mean_eV, sd_eV = extract_heating(sp, dimension)
    mass_kg = phantom_mass_kg(cfg, dimension, extent)
    return mean_eV, sd_eV, mass_kg


def main():
    xs = os.environ.get("OPENMC_CROSS_SECTIONS", "")
    if not xs or not Path(xs).exists():
        raise SystemExit(
            "OPENMC_CROSS_SECTIONS not set or does not exist -- see "
            "docs/openmc_stack.md. Refusing to write a fixture with no real "
            "nuclear data behind it.")

    runs = {"baseline": [], "feedback": []}
    completed = 0
    expected = N_OUTER * N_REPLICATES * 2
    # DIFFERENT DENSITIES HAVE DIFFERENT MASS. phantom_mass_kg depends on
    # medium.density_g_cm3, so baseline (1.0) and feedback (1.35) masses are
    # NOT the same array -- caught by inspecting the first regeneration's
    # output before trusting it (mass_kg summed to the baseline-only value
    # for both conditions). Recorded once per condition, not once overall.
    mass_kg_by_label = {}
    for outer in range(N_OUTER):
        for rep in range(N_REPLICATES):
            seed = 1000 * (outer + 1) + rep
            for label, density in (("baseline", DENSITY_BASELINE),
                                   ("feedback", DENSITY_FEEDBACK)):
                with __import__("tempfile").TemporaryDirectory() as td:
                    mean_eV, sd_eV, m = _run_one(density, seed, Path(td))
                mass_kg_by_label.setdefault(label, m)
                runs[label].append({
                    "outer": outer, "replicate": rep, "seed": seed,
                    "heating_mean_eV_per_src": _round_nested(mean_eV),
                    "heating_sd_eV_per_src": _round_nested(sd_eV),
                })
                completed += 1
                print(f"  [{completed}/{expected}] {label} outer={outer} "
                      f"rep={rep} seed={seed} done")

    if completed != expected:
        raise SystemExit(
            f"only {completed}/{expected} runs completed -- refusing to "
            "write a partial fixture over the committed one.")

    import openmc
    fixture = {
        "description": ("Golden tally: water-phantom point source, baseline "
                        "vs density x1.35 feedback (openmc_nested_pilot.py's "
                        "own DENSITY_SCALE lever), 2 outer draws x 3 "
                        "replicates. Regenerate with "
                        "coupling/scripts/regenerate_golden_tally.py under "
                        "the openmc-biofilms env."),
        "what_legitimately_changes_these_numbers": (
            "A change to OpenMC, to the nuclear data library, or to any module "
            "that produces the tally (config, model, mesh, materials, dose -- "
            "see golden-tally-verification.yml's paths filter, which the "
            "test_gate_composition.py control keeps in step with this script's "
            "real imports). NOT a change in the last digits: values are stored "
            "to 6 significant figures because the tally's own Monte Carlo "
            "relative sd has a median of 27% per cell, so digits past the "
            "sixth are not measurement. Same convention, and same reason, as "
            "tests/fixtures/serial_seed42.csv pinning 5 decimals."),
        "significant_figures": SIGNIFICANT_FIGURES,
        "openmc_version": openmc.__version__,
        "nuclear_data_id": _nuclear_data_id(),
        "n_outer": N_OUTER,
        "n_replicates": N_REPLICATES,
        "density_baseline_g_cm3": DENSITY_BASELINE,
        "density_feedback_g_cm3": DENSITY_FEEDBACK,
        # Rounded on the same convention, though these are exact CSG geometry
        # rather than sampled: one precision rule for the file, so a reader
        # does not have to know which arrays were which.
        "mass_kg_baseline": _round_nested(mass_kg_by_label["baseline"]),
        "mass_kg_feedback": _round_nested(mass_kg_by_label["feedback"]),
        "mesh_dimension": list(mass_kg_by_label["baseline"].shape),
        "runs": runs,
    }

    FIXTURE_PATH.parent.mkdir(parents=True, exist_ok=True)
    FIXTURE_PATH.write_text(json.dumps(fixture, indent=1))
    print(f"\nwrote {FIXTURE_PATH} ({FIXTURE_PATH.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
