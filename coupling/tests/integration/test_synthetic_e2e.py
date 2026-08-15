"""The synthetic one-way chain, run end to end.

This is the acceptance test for the whole path: snapshot -> OpenMC -> eV/source
-> exact CSG mass -> Gy/source -> Gy/s -> dose field -> cell/lineage/generation
attribution. The script records every invariant as data and returns non-zero if
any fails, so this test asserts on the verdict rather than re-deriving it.

The verdict's CLAIM fields are asserted separately from its pass/fail, because
the point of the run is not only that it worked but that it cannot be read as
having calibrated anything.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pytest

openmc = pytest.importorskip("openmc")

_XS = os.environ.get("OPENMC_CROSS_SECTIONS", "")
pytestmark = [
    pytest.mark.openmc,
    pytest.mark.skipif(not _XS or not Path(_XS).exists(),
                       reason="OPENMC_CROSS_SECTIONS not available"),
]

_REPO = Path(__file__).resolve().parents[3]
_CONFIG = _REPO / "config" / "reference_synthetic_biofilm_e2e.toml"
sys.path.insert(0, str(_REPO / "coupling" / "scripts"))

from conftest import make_snapshot_file  # noqa: E402


@pytest.fixture(scope="module")
def verdict(tmp_path_factory):
    import synthetic_e2e

    out = tmp_path_factory.mktemp("e2e")
    snap = out / "snapshot.h5"
    # n=20 to match the synthetic reference system's lattice, which is what
    # `check_lattice_congruence` and the source placement are declared against.
    make_snapshot_file(snap, n=20)
    verdict_path = out / "verdict.json"
    # Fewer rays than a real run: this checks the chain, not the convergence.
    rc = synthetic_e2e.main([
        "--snapshot", str(snap), "--config", str(_CONFIG),
        "--outdir", str(out), "--volume-samples", "500000",
        "--verdict", str(verdict_path)])
    return rc, json.loads(verdict_path.read_text())


def test_every_invariant_passed(verdict):
    rc, v = verdict
    failed = [r["invariant"] for r in v["invariants"] if not r["passed"]]
    assert not failed, f"failed invariants: {failed}"
    assert rc == 0 and v["execution_status"] == "PASSED"


def test_the_chain_covered_each_stage(verdict):
    _, v = verdict
    names = {r["invariant"] for r in v["invariants"]}
    # One per amendment the run exists to demonstrate.
    assert any(n.startswith("3_") for n in names)        # orientation / shape
    assert any(n.startswith("4_") for n in names)        # unit chain
    assert sum(n.startswith("5") for n in names) == 3    # mass denominator
    assert sum(n.startswith("6") for n in names) == 4    # attribution, incl. gen 0
    assert any(n.startswith("7_") for n in names)        # ordering
    assert any(n.startswith("8_") for n in names)        # label-only reuse
    assert any(n.startswith("9_") for n in names)        # provenance


def test_the_verdict_makes_no_biological_claim(verdict):
    """Invariant 10, asserted STRUCTURALLY. A grep over prose would let a
    wording edit flip a scientific safety check."""
    _, v = verdict
    assert v["target_calibration"] is False
    assert v["evidence_policy"] == "synthetic"
    assert v["execution_class"] == "synthetic_validation"
    assert v["claim_class"] == "integration_path_validation"
    assert v["biological_calibration"] == "NOT_EVALUATED"
    assert v["production_mesh_selected"] is False
    assert v["seconds_per_mcs_calibrated"] is False
    for applied in ("hamiltonian_applied", "melanin_applied",
                    "membrane_response_applied", "cpm_time_advanced"):
        assert v[applied] is False, f"{applied} must stay false: the one-way " \
                                    "chain stops before any response"


def test_the_result_can_be_reproduced_from_its_own_metadata(verdict):
    _, v = verdict
    prov = v["provenance"]
    assert prov["reference_system_id"] == "synthetic_biofilm_e2e"
    assert prov["target_calibration"] is False
    # Identity of the inputs, the run and the code — enough to rebuild it.
    for key in ("config_sha256", "snapshot_sha256", "transport_state_hash",
                "dose_state_hash", "label_state_hash", "git_commit"):
        assert prov[key], f"{key} is empty"
    assert prov["openmc_version"] == openmc.__version__
    assert prov["material_volume_n_samples"] == 500000
