#!/usr/bin/env python3
"""Verify the slab (biofilm-depth) stencil of ``biofilms_radiodialysis.R``.

WHAT THIS COVERS, AND WHAT IT DOES NOT
--------------------------------------
This script rebuilds the slab face weights in numpy and checks them against a
closed-form solution.  ``biofilms_radiodialysis.R`` is what ships, and it is
covered separately by ``verify_biofilm_depth_profile.R`` (see the last section).
So this file establishes ONE of two independent codings of one scheme; the
shipping implementation is verified by its companion, not here.

Do not report "the geometry change is verified against the analytic solution"
without saying which of the two codings was verified.

THE CONTRACT BEING CHECKED
--------------------------
Setting ds/dt = 0 in the immobile-phase equation gives

    s = U c / (k_des + k_loss),      U = X_total * (k_ads + k_red*f_red_active)

and substituting it back shows the sorption/desorption pair is a *reversible
buffer* that cancels.  The only true steady sink is precipitation:

    k_eff = U * k_loss / (k_des + k_loss)

so the steady mobile equation is D c'' = k_eff c, whose solution on 0 <= z <= L
with no-flux at the substratum and a Dirichlet bulk face is

    c(z) = c_b * cosh(z / lam) / cosh(L / lam),     lam = sqrt(D / k_eff)

Four assertions:

  1. accuracy  -- max relative error < 1e-5 on a 40-node grid
  2. order     -- error falls ~4x per grid halving (second order)
  3. NEGATIVE CONTROL -- the *cylindrical* weights, run against this same slab
     solution, must FAIL assertion 1.  Without this the suite would pass on a
     stencil that ignored geometry entirely (AGENTS.md rule 1).
  4. non-flat regime -- 1 and 2 repeated at phi = L/lam ~ 3.7, so the operator
     is exercised where a gradient actually exists rather than only in the
     nearly-flat default regime.

WHAT THIS DOES *NOT* ESTABLISH
------------------------------
Whether the model predicts a *measurable* gradient across a real film.  That is
downstream of ``RADIODIALYSIS: BLOCKED`` (README.md:350-353): U is proportional
to X_total and lam ~ 1/sqrt(U), so a 20-100x correction to the biomass basis
grows lam by 4.5-10x and flattens the profile.  The X_total values used here
are TEST VALUES chosen to exercise the numerics.  They are not claims about a
biofilm and no profile number from this script may be quoted as one.

Nor does it cover the transient, which is the regime an SECM profile would
actually be compared against; separation of variables would give an analytic
form for it.  Assertion 4 covers the spatial operator in the non-flat regime.
A transient check would additionally cover the time integration.

THE COMPANION CHECK ON THE R ITSELF
-----------------------------------
``analysis/verify_biofilm_depth_profile.R`` now covers the shipping R code --
its assembled operator, its steady profile, its lsoda time integration, and a
regression proving the cylindrical path is unchanged -- with mutation-based
negative controls.  Run both::

    Rscript analysis/verify_biofilm_depth_profile.R
    coupling/.venv/bin/python analysis/verify_biofilm_depth_profile.py

They agree to the digit on the shared quantity (max rel err 2.31e-07 at N=40),
which is worth more than either alone: this file is an independent re-coding,
so agreement rules out a shared transcription error, and it still runs where R
is absent.  Keep both.  Neither supersedes the other.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np

# Tolerances.
#
# A fixed absolute tolerance cannot serve here: the truncation error scales as
# (dz/lam)^2 = (phi/(N-1))^2, so any single number either fails a strong
# gradient or is vacuous for a weak one.  The bound is therefore on the
# NORMALISED constant  C = err * ((N-1)/phi)^2,  which is N-independent (that
# is what "second order" means) and depends only weakly on phi.  Measured:
#
#     phi   0.30    1.00    2.00    3.70    6.00
#     C     0.0036  0.0317  0.0803  0.1539  0.2501     -> saturating near 1/4
#
# 0.3 bounds all of them.  It is deliberately not tighter: the sharp check is
# the refinement ratio, which sits at 4.05-4.21 and breaks on any stencil error.
ERROR_CONSTANT_MAX = 0.3
ORDER_TOL = 3.6  # second order is 4.0; allow for the coarsest-grid pre-asymptote
GRID_N = 40


@dataclass(frozen=True)
class Kinetics:
    """Rate constants as written in biofilms_radiodialysis.R:141-148."""

    k_ads: float = 0.05
    k_red: float = 0.02
    k_des: float = 0.005
    k_loss: float = 0.001
    f_red_active: float = 0.3  # ACTIVE-reducer fraction OF X_total, in [0,1]
    X_total: float = 1.0  # TEST VALUE -- gated, see module docstring

    def validate(self) -> None:
        for name, value in asdict(self).items():
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and non-negative, got {value!r}")
        if self.f_red_active > 1.0:
            raise ValueError(
                "f_red_active is a fraction OF X_total and must lie in [0, 1], "
                f"got {self.f_red_active!r}. It is not a second biomass density; "
                "the previous fixed X_red let the reducing mass exceed X_total.")
        if self.k_des + self.k_loss <= 0.0:
            raise ValueError("k_des + k_loss must be positive for a steady state to exist")

    @property
    def uptake(self) -> float:
        """U, the transient sink (s^-1)."""
        self.validate()
        return self.X_total * (self.k_ads + self.k_red * self.f_red_active)

    @property
    def k_eff(self) -> float:
        """The only true steady sink: sorption/desorption cancels."""
        return self.uptake * self.k_loss / (self.k_des + self.k_loss)


def face_weights(grid: np.ndarray, geom: str) -> tuple[np.ndarray, np.ndarray]:
    """Port of face_weights() in biofilms_radiodialysis.R.

    Node 0 is NaN for the same reason it is NA in the R source: it uses the
    reflected-ghost form and reads no metric weight.  NaN rather than a
    plausible number so that indexing it propagates loudly.
    """
    if geom not in ("cylindrical", "slab"):
        raise ValueError(f"geom must be 'cylindrical' or 'slab', got {geom!r}")
    dz = grid[1] - grid[0]
    if geom == "slab":
        w_plus = np.ones_like(grid)
        w_minus = np.ones_like(grid)
    else:
        with np.errstate(divide="ignore", invalid="ignore"):
            w_plus = (grid + 0.5 * dz) / grid
            w_minus = (grid - 0.5 * dz) / grid
    w_plus[0] = np.nan
    w_minus[0] = np.nan
    return w_plus, w_minus


def steady_profile(n: int, length: float, diffusivity: float, sink: float,
                   c_bulk: float, geom: str) -> tuple[np.ndarray, np.ndarray]:
    """Solve D*L(c) = sink*c on the given grid, Dirichlet at the outer face.

    Uses the same operator as radiodialysis_rhs: reflected ghost at node 0
    (factor 2, identical for the cylindrical axis and the slab substratum),
    face-weighted interior.  The Dirichlet outer face is the k_L -> infinity
    limit of the code's Robin condition, which is what has a closed form.
    """
    if n < 4:
        raise ValueError(f"need at least 4 nodes, got {n}")
    grid = np.linspace(0.0, length, n)
    dz = grid[1] - grid[0]
    w_plus, w_minus = face_weights(grid, geom)
    operator = np.zeros((n, n))
    rhs = np.zeros(n)

    # node 0: reflected ghost c[-1] = c[1]  ->  2*(c[1] - c[0])/dz^2
    operator[0, 0] = -(2.0 * diffusivity / dz**2 + sink)
    operator[0, 1] = 2.0 * diffusivity / dz**2

    for i in range(1, n - 1):
        operator[i, i - 1] = diffusivity * w_minus[i] / dz**2
        operator[i, i] = -(diffusivity * (w_plus[i] + w_minus[i]) / dz**2 + sink)
        operator[i, i + 1] = diffusivity * w_plus[i] / dz**2

    operator[-1, -1] = 1.0
    rhs[-1] = c_bulk
    return np.linalg.solve(operator, rhs), grid


def analytic(grid: np.ndarray, length: float, lam: float, c_bulk: float) -> np.ndarray:
    return c_bulk * np.cosh(grid / lam) / np.cosh(length / lam)


def max_rel_error(n: int, length: float, diffusivity: float, sink: float,
                  c_bulk: float, geom: str) -> float:
    numeric, grid = steady_profile(n, length, diffusivity, sink, c_bulk, geom)
    exact = analytic(grid, length, math.sqrt(diffusivity / sink), c_bulk)
    return float(np.max(np.abs(numeric - exact) / exact))


def check_cylindrical_refactor_equivalence(n: int = 40, radius: float = 1.0,
                                           diffusivity: float = 1e-3) -> dict[str, Any]:
    """Regression guard: the weight form must equal the original inline stencil.

    The R change hoisted r_plus/r_i and r_minus/r_i out of radiodialysis_rhs
    into face_weights().  That touches the CYLINDRICAL path, which ships and is
    cross-validated against biofilms_potts.jl:1153-1163.  The slab assertions
    cannot protect it: there w_plus == w_minus == 1, so a weight-ordering error
    is invisible.  A mutation test confirmed exactly that -- swapping the two
    weights survived every other check in this file.

    This compares the refactored form against the pre-change arithmetic

        D * (r_plus*(c[i+1]-c[i]) - r_minus*(c[i]-c[i-1])) / (r_i * dr^2)

    node by node, and is asymmetric under a swap.
    """
    grid = np.linspace(0.0, radius, n)
    dr = grid[1] - grid[0]
    w_plus, w_minus = face_weights(grid, "cylindrical")
    # A profile with structure at every node, so no term cancels by accident.
    c = np.cos(3.0 * grid / radius) + 0.5 * grid / radius

    original = np.empty(n - 2)
    refactored = np.empty(n - 2)
    for i in range(1, n - 1):
        r_i = grid[i]
        r_plus, r_minus = r_i + 0.5 * dr, r_i - 0.5 * dr
        original[i - 1] = diffusivity * (
            r_plus * (c[i + 1] - c[i]) - r_minus * (c[i] - c[i - 1])
        ) / (r_i * dr**2)
        refactored[i - 1] = diffusivity * (
            w_plus[i] * (c[i + 1] - c[i]) - w_minus[i] * (c[i] - c[i - 1])
        ) / dr**2

    scale = np.max(np.abs(original))
    residual = float(np.max(np.abs(original - refactored)) / scale)

    # Negative control: the same comparison with the weights swapped must NOT
    # come out equal, or this guard proves nothing.
    swapped = np.empty(n - 2)
    for i in range(1, n - 1):
        swapped[i - 1] = diffusivity * (
            w_minus[i] * (c[i + 1] - c[i]) - w_plus[i] * (c[i] - c[i - 1])
        ) / dr**2
    swap_residual = float(np.max(np.abs(original - swapped)) / scale)

    # The baseline verdict GATES the control verdict.  A negative control can
    # fail to measure anything in two opposite ways: never firing (obviously
    # broken) or ALWAYS firing, which looks like perfect detection and is the
    # dangerous one.  A wrong reference makes every mutant read as caught -- it
    # happened in the R harness, where an omitted Robin term put the baseline
    # at 9.8e-01 and three mutants read as caught over a harness measuring
    # nothing.  So if the baseline does not match, the swap control is reported
    # as INVALID rather than as a pass.
    baseline_ok = residual < 1e-13
    checks = {"matches_pre_change_stencil": baseline_ok}
    if baseline_ok:
        checks["swap_control_differs"] = swap_residual > 1e-6
    else:
        checks["swap_control_INVALID_baseline_failed"] = False

    return {
        "case": "cylindrical_refactor_equivalence",
        "note": "protects the shipping cylindrical path and its Julia cross-validation",
        "residual": residual,
        "swap_control_residual": swap_residual,
        "baseline_gated": True,
        "checks": checks,
        "passed": all(checks.values()),
    }


def run_case(name: str, length: float, diffusivity: float, sink: float,
             note: str, c_bulk: float = 1.0) -> dict[str, Any]:
    lam = math.sqrt(diffusivity / sink)
    phi = length / lam
    errors = {n: max_rel_error(n, length, diffusivity, sink, c_bulk, "slab")
              for n in (GRID_N // 2, GRID_N, GRID_N * 2, GRID_N * 4)}
    ratios = [errors[n] / errors[2 * n] for n in sorted(errors)[:-1]]
    control = max_rel_error(GRID_N, length, diffusivity, sink, c_bulk, "cylindrical")
    profile, _ = steady_profile(GRID_N, length, diffusivity, sink, c_bulk, "slab")

    def constant(err: float, n: int) -> float:
        return err * ((n - 1) / phi) ** 2

    bound = ERROR_CONSTANT_MAX * (phi / (GRID_N - 1)) ** 2
    checks = {
        "accuracy": constant(errors[GRID_N], GRID_N) < ERROR_CONSTANT_MAX,
        "second_order": all(r > ORDER_TOL for r in ratios),
        "negative_control_fails_slab_bound": control > bound,
    }
    return {
        "case": name,
        "note": note,
        "thiele_modulus": phi,
        "lambda_um": lam * 1e4,
        "profile_ratio_c0_over_cL": float(profile[0] / profile[-1]),
        "errors": {str(k): v for k, v in errors.items()},
        "error_constant": constant(errors[GRID_N], GRID_N),
        "refinement_ratios": ratios,
        "slab_bound_at_N": bound,
        "negative_control_error": control,
        "checks": checks,
        "passed": all(checks.values()),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--report", type=Path, default=None,
                        help="write the machine-readable report here")
    args = parser.parse_args()

    length = 0.01  # 100 um film
    diffusivity = 1e-5
    kin = Kinetics()
    cases = [
        run_case(
            "default_regime", length, diffusivity, kin.k_eff,
            note="k_eff from the file's own rate constants at X_total=1.0 "
                 "(a TEST VALUE: the biomass basis is gated). Nearly flat.",
        ),
        # Assertion 4.  The sink is set directly from a target phi rather than
        # by inventing a biomass number: this case exists to exercise the
        # operator where a gradient exists, and corresponds to no parameter set.
        run_case(
            "strong_gradient", length, diffusivity,
            sink=diffusivity / (length / 3.7) ** 2,
            note="phi = 3.7 imposed directly on the operator. NOT a physical "
                 "parameter set and not a prediction about any biofilm.",
        ),
        check_cylindrical_refactor_equivalence(),
    ]

    report = {
        "tolerances": {"error_constant_max": ERROR_CONSTANT_MAX,
                       "order_ratio": ORDER_TOL, "N": GRID_N},
        "k_eff_over_U": Kinetics().k_eff / Kinetics().uptake,
        "covers": "the numpy specification of the slab stencil",
        "does_not_cover": [
            "biofilms_radiodialysis.R itself (covered by verify_biofilm_depth_profile.R)",
            "the transient regime",
            "whether any gradient is measurable (downstream of RADIODIALYSIS: BLOCKED)",
        ],
        "cases": cases,
        "passed": all(c["passed"] for c in cases),
    }

    for case in cases:
        print(f"[{'PASS' if case['passed'] else 'FAIL'}] {case['case']}")
        if "thiele_modulus" in case:
            print(f"       phi={case['thiele_modulus']:.3f} "
                  f"lambda={case['lambda_um']:.1f}um "
                  f"c(0)/c(L)={case['profile_ratio_c0_over_cL']:.4f}")
            print(f"       err(N={GRID_N})={case['errors'][str(GRID_N)]:.2e} "
                  f"C={case['error_constant']:.4f}/{ERROR_CONSTANT_MAX} "
                  f"ratios={['%.2f' % r for r in case['refinement_ratios']]} "
                  f"control={case['negative_control_error']:.2e} "
                  f"(bound {case['slab_bound_at_N']:.2e})")
        else:
            print(f"       residual={case['residual']:.2e} "
                  f"swap_control={case['swap_control_residual']:.2e}")
        for check, ok in case["checks"].items():
            if not ok:
                print(f"       FAILED: {check}")

    print(f"\nk_eff/U = {report['k_eff_over_U']:.6f} = 1/{1/report['k_eff_over_U']:.0f} "
          f"(structural, independent of X_total)")
    print("COVERS: the numpy specification, NOT biofilms_radiodialysis.R.")

    if args.report:
        args.report.write_text(json.dumps(report, indent=2) + "\n")

    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
