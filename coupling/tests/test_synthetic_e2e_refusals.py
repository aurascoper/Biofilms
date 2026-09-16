"""`synthetic_e2e.py` runs at the base tally resolution and must say so.

`mesh_refinement_factor` is honoured by `build_biofilm_cylinder_model`, so a
config declaring 2 gives the model a tally twice as fine as the mass array the
script builds beside it; `extract_heating` then fails on the size mismatch,
after the histories are spent. The bundle writer, added later, resolved its
own dimension WITH the factor while the run resolved without it, so the two
disagreed about the shape of the same field. The refinement study is
`subvoxel_refinement.py`; this runner refuses the parameter by name before
transport, the way `drivers.scan` does.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_ROOT / "coupling" / "tests"))
_spec = importlib.util.spec_from_file_location(
    "synthetic_e2e", _ROOT / "coupling" / "scripts" / "synthetic_e2e.py")
e2e = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(e2e)

from conftest import VALID_CONFIG  # noqa: E402


def test_refinement_is_refused_by_name_before_transport(tmp_path):
    """THE CONTROL: openmc is not installed in this tier, so reaching the
    import at all raises ModuleNotFoundError rather than SystemExit; the
    refusal has to come first. The snapshot path does not exist either, for
    the same reason."""
    cfg = tmp_path / "refined.toml"
    cfg.write_text(VALID_CONFIG.replace(
        "[transport]", "[transport]\n  [transport.mesh]\n  refinement_factor = 2\n"))
    with pytest.raises(SystemExit, match="refinement_factor"):
        e2e.main(["--snapshot", str(tmp_path / "nope.h5"),
                  "--config", str(cfg), "--outdir", str(tmp_path / "out")])
    assert not (tmp_path / "out").exists()


def test_the_base_resolution_gets_past_the_refusal(tmp_path):
    """The control's control: factor 1 must reach the next stage. Which stage
    that is depends on the tier: in the bare tier the openmc import is
    missing; in the golden-tally job openmc is present and the next thing
    missing is the dosimetry config beside the transport one. The first
    version named only the import and went red in the tier with openmc."""
    cfg = tmp_path / "base.toml"
    cfg.write_text(VALID_CONFIG)
    with pytest.raises((ModuleNotFoundError, ImportError, FileNotFoundError)):
        e2e.main(["--snapshot", str(tmp_path / "nope.h5"),
                  "--config", str(cfg), "--outdir", str(tmp_path / "out")])
