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
