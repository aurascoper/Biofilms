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
import os
import sys
from pathlib import Path

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
                    "heating_mean_eV_per_src": mean_eV.tolist(),
                    "heating_sd_eV_per_src": sd_eV.tolist(),
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
                        "the openmc-biofilms env. Legitimately changes only "
                        "if OpenMC or the nuclear data library changes."),
        "openmc_version": openmc.__version__,
        "nuclear_data_id": _nuclear_data_id(),
        "n_outer": N_OUTER,
        "n_replicates": N_REPLICATES,
        "density_baseline_g_cm3": DENSITY_BASELINE,
        "density_feedback_g_cm3": DENSITY_FEEDBACK,
        "mass_kg_baseline": mass_kg_by_label["baseline"].tolist(),
        "mass_kg_feedback": mass_kg_by_label["feedback"].tolist(),
        "mesh_dimension": list(mass_kg_by_label["baseline"].shape),
        "runs": runs,
    }

    FIXTURE_PATH.parent.mkdir(parents=True, exist_ok=True)
    FIXTURE_PATH.write_text(json.dumps(fixture, indent=1))
    print(f"\nwrote {FIXTURE_PATH} ({FIXTURE_PATH.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
