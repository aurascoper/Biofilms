"""The vocabulary and shapes that `coupling/` and `calibration/` must agree on.

Why this package exists: `calibration/` may not import `coupling/` (undeclared
dependency, and it would drag h5py into a numpy-only package), but both sides
have to agree on what a boundary condition is called, what a closed composition
means, and how a material renders into TOML. The alternative — duplicating the
vocabularies on both sides — is a drift bug waiting for the first time someone
adds a boundary type to one and not the other.

So the contract lives here, in a package with NO dependencies at all, not even
numpy. It holds vocabulary, dataclass shapes, and pure rendering. It holds no
policy: whether a given material is ADMISSIBLE is an evidence judgement, and
that stays in `calibration/biofilm_calibration/materials/export.py`.
"""

from __future__ import annotations

import math
import subprocess
from dataclasses import dataclass, field

__all__ = [
    "ALLOWED_BOUNDARY_CONDITIONS", "ALLOWED_SOURCE_SPATIAL",
    "ALLOWED_SOURCE_ANGULAR", "EVIDENCE_POLICIES", "EXECUTION_CLASSES",
    "MaterialSpec", "closed_composition_problems", "render_material_toml",
    "git_provenance",
]

# --- vocabulary ----------------------------------------------------------
# These were `_ALLOWED_BC` / `_ALLOWED_SPATIAL` / `_ALLOWED_ANGULAR`, private to
# the coupling loader. A contract cannot be built out of another package's
# underscore-prefixed implementation details, so they are public here and the
# loader aliases them.
ALLOWED_BOUNDARY_CONDITIONS = frozenset({"vacuum", "reflective"})
ALLOWED_SOURCE_SPATIAL = frozenset({"line_z_axis", "point_origin"})
ALLOWED_SOURCE_ANGULAR = frozenset({"isotropic"})

# Whether a run may use unmeasured values. ORTHOGONAL to `system_provenance`
# (published_replica / certified_component / engineered_composite / declared),
# which says what kind of system it is, and to per-value `evidence_basis`. An
# engineered composite can hold measured, certified, derived and declared
# components at once, so "measured vs synthetic" cannot be a provenance value.
EVIDENCE_POLICIES = frozenset({"measured_only", "synthetic"})

# What kind of run this is. `target_calibration = false` alone cannot tell an A0
# benchmark from a public surrogate, a synthetic fixture, an uncalibrated
# Reference D, or an exploratory sensitivity case.
EXECUTION_CLASSES = frozenset({
    "reference_benchmark",      # A0: numerical validation against known physics
    "surrogate_validation",     # public data, not the target organism
    "synthetic_validation",     # invented values, exercises the software path
    "target_calibration",       # the real thing
    "exploratory_sensitivity",  # ranking, not calibrating
})


# --- material shape ------------------------------------------------------

@dataclass(frozen=True)
class MaterialSpec:
    """A material as the transport loader needs it: a density and closed
    elemental MASS fractions. Carries no evidence claim by construction —
    attaching one here is how a placeholder gets laundered into a measurement.
    """
    name: str
    density_g_cm3: float
    elements: dict = field(default_factory=dict)
    material_model_kind: str = "hydrated_effective_medium"


def closed_composition_problems(elements: dict) -> list[str]:
    """Why these mass fractions are not a usable composition. Empty means they
    are. The coupling loader enforces the same closure; this exists so the
    emitting side can refuse before writing a file the loader would reject.
    """
    problems: list[str] = []
    if not elements:
        return ["no elemental composition"]
    for el, frac in elements.items():
        if not isinstance(frac, (int, float)) or isinstance(frac, bool):
            problems.append(f"mass fraction for {el} is not a number: {frac!r}")
        elif frac < 0:
            problems.append(f"negative mass fraction for {el}: {frac}")
    total = sum(v for v in elements.values() if isinstance(v, (int, float)))
    if not math.isclose(total, 1.0, rel_tol=1e-6):
        problems.append(f"elemental mass fractions sum to {total}, not 1 — the "
                        "coupling loader requires a closed composition")
    return problems


def render_material_toml(spec: MaterialSpec, class_name: str,
                         header: str = "") -> str:
    """Render `[materials.<class_name>]`, or raise on an unusable composition.

    Rendering only. It does NOT ask where the numbers came from — the caller
    decides whether this material is allowed to be written at all.
    """
    problems = closed_composition_problems(spec.elements)
    if spec.density_g_cm3 is None or spec.density_g_cm3 <= 0:
        problems.append(f"density {spec.density_g_cm3} is not positive")
    if problems:
        raise ValueError(f"refusing to render material {spec.name!r}:\n  - "
                         + "\n  - ".join(problems))

    lines = []
    if header:
        lines += [f"# {line}" for line in header.splitlines()]
    lines += [f"[materials.{class_name}]",
              f"density_g_cm3 = {spec.density_g_cm3!r}",
              f"  [materials.{class_name}.elements]"]
    lines += [f"  {el} = {frac!r}" for el, frac in sorted(spec.elements.items())]
    return "\n".join(lines) + "\n"


# --- run provenance ------------------------------------------------------

def git_provenance(root: str | None = None) -> dict:
    """The commit this ran at, and whether the tree was dirty.

    A dirty tree recorded as a clean commit is a false provenance claim, so the
    marker is part of the commit string rather than a separate field somebody
    can drop. `root` defaults to the enclosing repository, discovered rather
    than assumed from this file's depth — this module is imported from two
    packages at different depths.
    """
    def git(*args) -> str:
        cmd = ["git"] + (["-C", root] if root else []) + list(args)
        return subprocess.run(cmd, capture_output=True, text=True,
                              check=True).stdout.strip()

    try:
        head = git("rev-parse", "HEAD")
        # Tracked modifications only. Untracked files are not a difference
        # between the code and the commit — and the scripts that call this
        # WRITE some of them, so counting them self-reports every run as dirty.
        dirty = bool(git("status", "--porcelain", "--untracked-files=no"))
    except (OSError, subprocess.CalledProcessError):
        return {"git_commit": None, "git_dirty": None}
    return {"git_commit": head + ("-dirty" if dirty else ""), "git_dirty": dirty}
