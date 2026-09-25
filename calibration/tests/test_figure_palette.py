"""The Python renderers draw species in the palette the model declares, read from the model.

tools/render_label_trajectory.py holds COLOURS, which tools/render_consortium.py imports for the
consortium panel and video. README says those frames are drawn "in the palette biofilms_potts.jl
declares". That was true when checked by eye and bound by nothing: COLOURS is typed, and
FIG_COLORS in biofilms_potts.jl is what the Julia figures use. This reads both and requires them
to agree, species by species. COLOURS[0] is the empty-site background and has no counterpart.

BOTH SIDES ARE READ AS SOURCE TEXT, NEITHER IS IMPORTED. The renderer imports h5py at module
level, and the calibration job installs no h5py, so importing it would fail this test in CI for
a reason that says nothing about colours. ast.literal_eval of the one assignment reads the list
the renderer would bind without running the file.

ONE COMPARISON, AND THE CONTROLS DRIVE IT. `palette_mismatch` is the whole check. The main test
calls it on the committed files and each control calls it on a copy with one colour changed, so
a control cannot pass while the check it guards has stopped comparing.
"""
from __future__ import annotations

import ast
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
MODEL = REPO / "biofilms_potts.jl"
RENDERER = REPO / "tools" / "render_label_trajectory.py"


def model_palette(source: str) -> list[str]:
    """The hex colours in biofilms_potts.jl's FIG_COLORS block, in declaration order."""
    m = re.search(r"const FIG_COLORS = \[(.*?)\]", source, re.S)
    assert m, "biofilms_potts.jl no longer declares FIG_COLORS as a literal block"
    return [h.lower() for h in re.findall(r'colorant"(#[0-9a-fA-F]{6})"', m.group(1))]


def renderer_palette(source: str) -> list[str]:
    """COLOURS[1:] from the renderer's source: the one top-level `COLOURS = [...]` assignment."""
    found = [node.value for node in ast.parse(source).body
             if isinstance(node, ast.Assign)
             and any(isinstance(t, ast.Name) and t.id == "COLOURS" for t in node.targets)]
    assert len(found) == 1, f"expected one top-level COLOURS assignment, found {len(found)}"
    return [c.lower() for c in ast.literal_eval(found[0])[1:]]


def palette_mismatch(model_source: str, renderer_source: str) -> list:
    """Every species index where the two palettes differ, plus a length problem if any."""
    model, renderer = model_palette(model_source), renderer_palette(renderer_source)
    out: list = [] if len(model) == len(renderer) == 7 else [("lengths", len(model), len(renderer))]
    out += [i for i in range(min(len(model), len(renderer))) if model[i] != renderer[i]]
    return out


def test_the_renderer_palette_is_the_model_palette():
    """THE ONE THAT MATTERS."""
    assert palette_mismatch(MODEL.read_text(encoding="utf-8"),
                            RENDERER.read_text(encoding="utf-8")) == []


def _change_one(source: str, pattern: str, template: str, victim: str) -> str:
    """`source` with the first quoted occurrence of `victim` replaced by a colour neither
    palette uses."""
    replacement = "#000001" if victim != "#000001" else "#000002"
    hits = [m for m in re.finditer(pattern, source) if m.group(1).lower() == victim]
    assert hits, f"{victim} not found to mutate"
    return source[:hits[0].start()] + template.format(replacement) + source[hits[0].end():]


def test_a_changed_model_colour_is_reported():
    """CONTROL: species index 2 changed in a copy of the model source, and only there."""
    model_src, renderer_src = MODEL.read_text(encoding="utf-8"), RENDERER.read_text(encoding="utf-8")
    victim = model_palette(model_src)[2]
    mutated = _change_one(model_src, r'colorant"(#[0-9a-fA-F]{6})"', 'colorant"{}"', victim)
    before, after = model_palette(model_src), model_palette(mutated)
    assert [i for i in range(7) if before[i] != after[i]] == [2]      # the mutation applied, once
    assert palette_mismatch(mutated, renderer_src) == [2]


def test_a_changed_renderer_colour_is_reported():
    """CONTROL: species index 2 changed in a copy of the renderer source, and only there."""
    model_src, renderer_src = MODEL.read_text(encoding="utf-8"), RENDERER.read_text(encoding="utf-8")
    victim = renderer_palette(renderer_src)[2]
    mutated = _change_one(renderer_src, r'"(#[0-9a-fA-F]{6})"', '"{}"', victim)
    before, after = renderer_palette(renderer_src), renderer_palette(mutated)
    assert [i for i in range(7) if before[i] != after[i]] == [2]
    assert palette_mismatch(model_src, mutated) == [2]
