"""The label-trajectory renderer refuses a completed build BEFORE it renders.

render() overwrites the run's verified Figure 5 files. A refusal that came after
it (inside stage_and_build) left the run's figures changed under an unchanged
manifest and receipt. This loads the tool without matplotlib, which the test
tiers do not install; render() imports it lazily for that reason.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]


def _tool():
    spec = importlib.util.spec_from_file_location(
        "render_label_trajectory", REPO / "tools" / "render_label_trajectory.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_completed_build_is_refused_before_render_runs(tmp_path, monkeypatch):
    mod = _tool()
    run = tmp_path / "run"
    (run / "manuscript").mkdir(parents=True)
    monkeypatch.setattr(mod, "verified_manifest", lambda run: {"status": "manuscript_built"})

    def render_must_not_run(run, manifest):
        raise AssertionError("render ran before the completed-build refusal")

    monkeypatch.setattr(mod, "render", render_must_not_run)
    monkeypatch.setattr(sys, "argv", ["render_label_trajectory.py", str(run)])
    with pytest.raises(FileExistsError, match="existing manuscript output"):
        mod.main()


def test_a_run_without_manuscript_output_reaches_render(tmp_path, monkeypatch):
    mod = _tool()
    run = tmp_path / "run"
    run.mkdir()
    monkeypatch.setattr(mod, "verified_manifest", lambda run: {"status": "trajectory_verified"})
    monkeypatch.setattr(mod, "render", lambda run, manifest: (_ for _ in ()).throw(RuntimeError("render reached")))
    monkeypatch.setattr(sys, "argv", ["render_label_trajectory.py", str(run)])
    with pytest.raises(RuntimeError, match="render reached"):
        mod.main()
