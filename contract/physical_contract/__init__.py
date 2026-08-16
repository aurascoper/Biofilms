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
    "SYSTEM_PROVENANCE",
    "MaterialSpec", "closed_composition_problems", "render_material_toml",
    "source_placement_problems",
    "git_provenance",
    "UNCERTAINTY_TYPES", "SAMPLING_ROLES", "EFFECT_DIRECTIONS",
    "FEEDBACK_GATE_VERDICTS", "distribution_row_problems",
    "DEFAULT_EXCEEDANCE_PROBABILITY", "DEFAULT_NUMERICAL_BUDGET_FRACTION",
    "DEFAULT_VOLUME_BUDGET_FRACTION", "DEFAULT_SURROGATE_BUDGET_FRACTION",
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

# What kind of system this is. The ledger has enforced this vocabulary since
# before the contract existed; it lives here now so the loader, the emitter and
# the ledger test all read one definition.
SYSTEM_PROVENANCE = frozenset({"published_replica", "certified_component",
                               "engineered_composite", "declared"})

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


# --- source placement ----------------------------------------------------

def source_placement_problems(position, origin, pitch_cm: float, n: int,
                              cylinder_radius_cm: float,
                              cylinder_length_cm: float) -> list[str]:
    """Why this source position is degenerate or out of bounds. Empty is fine.

    Lives here because the transport builder must REFUSE a bad placement and
    the emitting side must not WRITE one, and two implementations of "is this
    point on a lattice plane" are two chances to disagree about it.

    A point exactly on a lattice plane is a point on a surface: which cell the
    particle starts in is then settled by floating-point tie-breaking rather
    than by the geometry.
    """
    problems: list[str] = []
    x0, y0, z0 = origin
    side = n * pitch_cm
    cx, cy = x0 + side / 2.0, y0 + side / 2.0

    for axis, (p, lo) in enumerate(zip(position, (x0, y0, z0))):
        offset = (p - lo) % pitch_cm
        if min(offset, pitch_cm - offset) <= 1e-9 * pitch_cm:
            problems.append(
                f"axis {'xyz'[axis]}: source at {p} lies on a lattice plane "
                f"(origin {lo}, pitch {pitch_cm})")

    r = math.hypot(position[0] - cx, position[1] - cy)
    if r > cylinder_radius_cm * (1 - 1e-9):
        problems.append(
            f"source is {r} cm from the axis, outside the biological domain "
            f"(radius {cylinder_radius_cm})")
    z_lo, z_hi = z0, z0 + cylinder_length_cm
    if not (z_lo < position[2] < z_hi):
        problems.append(
            f"source z = {position[2]} is outside the domain [{z_lo}, {z_hi}]")
    return problems


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


# --- feedback gates ------------------------------------------------------
#
# Two gates, answering different questions. Conflating them would make
# "feedback matters" indistinguishable from "the software can perform
# feedback", which is the single most important distinction this repository
# has left to protect.
#
#   OFFLINE  a counterfactual transport experiment on immutable snapshots.
#            Does changing a biologically controlled material state produce a
#            transport effect distinguishable from Monte Carlo noise, numerical
#            discretization, calibration uncertainty and biological-state
#            uncertainty? It never advances the CPM and never feeds dose back
#            into biology.
#   ONLINE   an AUTHORIZATION condition on the real two-way loop. It may open
#            only after the offline gate passes, the calibration prerequisites
#            are ready, the current state is inside the validated envelope, and
#            the effect still exceeds a predeclared threshold.

# What KIND of uncertainty a quantity carries. The separation is load-bearing:
# putting a probability distribution on a numerical control or a model choice
# mixes epistemic uncertainty with engineering decisions and produces a
# probability nobody can interpret.
UNCERTAINTY_TYPES = frozenset({
    "physical_measurement",    # a posterior over a measured quantity
    "biological_posterior",    # a fitted response parameter
    "source_metrology",        # certificate/assay, propagated to date
    "monte_carlo_estimator",   # transport stochasticity; nested replicates
    "numerical_convergence",   # mesh, histories, ray counts; bounded, not sampled
    "discrete_model_choice",   # occupancy map, mass denominator; branched
    "nuclear_data_model",      # library version; sensitivity, not a Gaussian
    "declared_exact",
    "unsupported",
})

# HOW a quantity enters the experiment. `outer_random` is the only role that
# draws from a probability distribution, and it must be refused for numerical
# controls and model choices unless a separately reviewed probability model
# genuinely exists.
SAMPLING_ROLES = frozenset({
    "outer_random",              # drawn from its distribution per outer point
    "inner_transport_replicate", # an independent OpenMC seed
    "convergence_axis",          # swept to bound a residual, never sampled
    "scenario_branch",           # an explicit named alternative
    "fixed",
    "excluded",
})

# Roles that may never be assigned to these uncertainty types, because doing so
# would launder an engineering decision into an epistemic probability.
_NEVER_RANDOM = frozenset({"numerical_convergence", "discrete_model_choice",
                           "nuclear_data_model"})

# The direction a scientifically meaningful effect is expected to move in.
# Declared with the threshold, before results are inspected: a two-sided test
# would pass on an effect of the wrong sign.
EFFECT_DIRECTIONS = frozenset({"increase", "decrease"})

# One state machine, not a scatter of if-statements. Every terminal state says
# WHY, so a blocked gate is diagnostic rather than merely closed.
FEEDBACK_GATE_VERDICTS = (
    "NOT_EVALUATED",
    "BLOCKED_ON_PHYSICAL_CALIBRATION",
    "BLOCKED_ON_NUMERICAL_RESOLUTION",
    "BLOCKED_ON_BIOLOGICAL_CALIBRATION",
    "OFFLINE_EFFECT_BELOW_THRESHOLD",
    "OFFLINE_UNCERTAINTY_TOO_LARGE",
    "OFFLINE_PASS",
    "ONLINE_OUT_OF_DOMAIN",
    "ONLINE_UNCERTAINTY_TOO_LARGE",
    "ONLINE_ENABLED",
    "UNSUPPORTED_BY_CURRENT_MODEL",
)

# Conservative project defaults, versioned in an acceptance policy rather than
# presented as metrological constants. The scientific effect threshold itself
# gets NO default: until someone declares what magnitude of feedback would
# matter, the gate stays NOT_EVALUATED, which is the honest answer.
DEFAULT_EXCEEDANCE_PROBABILITY = 0.99
DEFAULT_NUMERICAL_BUDGET_FRACTION = 0.25   # U99(transport+numerical) <= 0.25*delta
DEFAULT_VOLUME_BUDGET_FRACTION = 0.10      # CSG mass is denominator-critical
DEFAULT_SURROGATE_BUDGET_FRACTION = 0.10   # emulator error << decision margin


def distribution_row_problems(row: dict) -> list[str]:
    """Why this uncertainty-ledger row may not be sampled as written.

    The refusal that matters: a mesh factor or a mass-denominator method is an
    engineering choice, not a random variable. Sampling one produces a
    probability that looks like evidence and is not.
    """
    problems: list[str] = []
    kind = (row.get("uncertainty_type") or "").strip()
    role = (row.get("sampling_role") or "").strip()
    if kind not in UNCERTAINTY_TYPES:
        problems.append(f"uncertainty_type {kind!r} not in {sorted(UNCERTAINTY_TYPES)}")
    if role not in SAMPLING_ROLES:
        problems.append(f"sampling_role {role!r} not in {sorted(SAMPLING_ROLES)}")
    if kind in _NEVER_RANDOM and role == "outer_random":
        problems.append(
            f"uncertainty_type={kind} may not take sampling_role=outer_random — "
            "assigning a probability distribution to a numerical control or a "
            "model choice mixes engineering decisions with epistemic "
            "uncertainty and makes the resulting probability uninterpretable. "
            "Use convergence_axis or scenario_branch")
    if kind == "unsupported" and role not in ("excluded", "fixed"):
        problems.append("an unsupported quantity may only be excluded or fixed")
    return problems
