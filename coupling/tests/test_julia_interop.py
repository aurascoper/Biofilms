"""Cross-language leg of the orientation chain (review amendment 9): a REAL
Julia-written transport snapshot must load, pass its own asymmetric probes,
and agree with its registry. Skipped when Julia isn't available. The return
leg (Python-written dose read back in Julia) lives in the Julia test suite
with import_dose.jl."""

import shutil
import subprocess
from pathlib import Path

import pytest

julia = shutil.which("julia")
pytestmark = pytest.mark.skipif(julia is None, reason="julia not on PATH")

REPO = Path(__file__).parents[2]


@pytest.fixture(scope="module")
def julia_snapshot(tmp_path_factory):
    out = tmp_path_factory.mktemp("jl") / "snap.h5"
    subprocess.run(
        [julia, f"--project={REPO}", str(REPO / "export_checkpoint.jl"),
         "transport", str(out), "--mcs", "3", "--seed", "21"],
        check=True, capture_output=True, timeout=600)
    return out


def test_dose_written_by_python_reads_back_in_julia(tmp_path):
    """Return leg of the orientation chain: a logical (x,y,z) dose field
    written by Python must read back in Julia with identical coordinates."""
    import numpy as np
    from biofilm_openmc.snapshot import write_dose_field

    n = 6
    x, y, z = np.meshgrid(np.arange(n), np.arange(n), np.arange(n), indexing="ij")
    field = (x + 10 * y + 100 * z).astype(float)   # unique value per voxel
    path = tmp_path / "dose.h5"
    write_dose_field(path, field, field * 0.0, {"schema_version": 1})

    script = f'''
    using HDF5
    f = h5open("{path}", "r")
    d = read(f["mesh/dose_rate_mean_Gy_s"])
    # logical [x,y,z] (0-based here) must equal x + 10y + 100z
    ok = size(d) == ({n}, {n}, {n}) &&
         d[1, 1, 1] == 0.0 && d[2, 1, 1] == 1.0 &&
         d[1, 2, 1] == 10.0 && d[1, 1, 2] == 100.0 &&
         d[3, 4, 5] == 2.0 + 30.0 + 400.0
    println(ok ? "CHAIN_OK" : "CHAIN_BROKEN")
    '''
    out = subprocess.run([julia, f"--project={REPO}", "-e", script],
                         capture_output=True, text=True, timeout=300)
    assert "CHAIN_OK" in out.stdout, out.stdout + out.stderr


def test_julia_snapshot_loads_and_probes_pass(julia_snapshot):
    from biofilm_openmc.snapshot import load_snapshot
    snap = load_snapshot(julia_snapshot)  # probe verification happens inside
    assert snap.cell_id.shape == (20, 20, 20)
    assert snap.material_class_source == "absent"
    assert len(snap.label_state_hash) == 64
    # registry consistent with label arrays
    assert snap.cells["volume"].sum() == (snap.cell_id > 0).sum()
    assert set(snap.species_id[snap.cell_id > 0].tolist()) <= set(range(1, 8))
    # wall convention: -1 present outside the biological domain, never inside
    assert (snap.cell_id == -1).any()
    assert not ((snap.cell_id == -1) & snap.interior_mask).any()
