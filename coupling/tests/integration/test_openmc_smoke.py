"""Fixed-photon smoke run of the actual biofilm model (synthetic demo
config) + the sparsity/uncertainty feasibility report the review requires
before any lineage-resolved field is promised. Results are voxel-averaged
absorbed-dose estimates under OpenMC's charged-particle local-deposition
approximation — never single-cell microdosimetry."""

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


def test_biofilm_model_smoke_run(tmp_path, snapshot, config):
    from biofilm_openmc.dose import dose_from_statepoint, sparsity_report
    from biofilm_openmc.fingerprint import label_state_hash, transport_state_hash
    from biofilm_openmc.materials import voxel_mass_kg
    from biofilm_openmc.lineage import aggregate_by_label
    from biofilm_openmc.results import write_transport_result

    from biofilm_openmc.model import build_model
    model = build_model(snapshot, config)

    sp_path = model.run(cwd=str(tmp_path), output=False)
    sp = openmc.StatePoint(sp_path)

    mass = voxel_mass_kg(snapshot, config)
    result = dose_from_statepoint(sp, mass, config.photons_per_second)

    assert np.isfinite(result.dose_rate_mean_Gy_s).all()
    assert (result.dose_rate_mean_Gy_s >= 0).all()
    assert result.dose_rate_mean_Gy_s.max() > 0, "no heating scored anywhere"

    # feasibility numbers (review amendment 11) — printed for the report
    occupied = snapshot.cell_id > 0
    rep = sparsity_report(result, occupied)
    print(f"\nfeasibility: {rep}")
    assert 0.0 <= rep["occupied_unscored_fraction"] <= 1.0

    # mass-weighted lineage attribution runs end-to-end on real output
    agg = aggregate_by_label(result.dose_rate_mean_Gy_s,
                             snapshot.lineage_id, mass)
    assert set(agg) == set(np.unique(snapshot.lineage_id[occupied]))

    # transport result written with full run identity
    write_transport_result(
        tmp_path / "transport_result.h5", result,
        transport_hash=transport_state_hash(snapshot, config, "endfb-viii.0"),
        openmc_version=openmc.__version__,
        nuclear_data_id=os.environ["OPENMC_CROSS_SECTIONS"],
        batches=config.batches, particles=config.particles, seed=config.seed)
    assert (tmp_path / "transport_result.h5").exists()
