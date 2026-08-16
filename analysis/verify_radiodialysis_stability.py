#!/usr/bin/env python3
"""Verify the radiodialysis linear-stability and time-step contracts.

This audit deliberately separates three statements that are easy to conflate:

1. The linearized reaction-diffusion *continuous/semi-discrete operator* can be
   asymptotically stable.
2. The loss-free spatially uniform mode can be only semistable (one zero
   eigenvalue), so "strictly stable for all non-negative parameters" is false.
3. A stable continuous operator does not make an explicit-Euler update stable;
   the time step must satisfy its own amplification condition.

The finite cylindrical domain is represented by a non-negative Laplacian mode
parameter for the analytical 2x2 calculation and by the repository's stated
radial finite-volume stencil for the numerical audit.  Plane waves are not
claimed to be the eigenfunctions of a finite cylinder with Robin boundaries.

The script writes a machine-readable report and exits non-zero only when the
analytical/semidiscrete contracts themselves are violated.  Pass
--require-explicit-euler-stability to promote the explicit time-step diagnostic
to a blocking check after the production solver owns a declared substep policy.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np


@dataclass(frozen=True)
class Rates:
    diffusion: float
    forward: float
    desorption: float
    loss: float

    def validate(self) -> None:
        for name, value in asdict(self).items():
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and non-negative, got {value!r}")


def mode_matrix(rates: Rates, laplacian_eigenvalue: float) -> np.ndarray:
    """2x2 perturbation matrix for a mode with -Laplacian eigenvalue mu >= 0."""
    rates.validate()
    mu = float(laplacian_eigenvalue)
    if not math.isfinite(mu) or mu < 0.0:
        raise ValueError(f"laplacian_eigenvalue must be finite and >= 0, got {mu!r}")
    return np.array(
        [
            [-(rates.diffusion * mu + rates.forward), rates.desorption],
            [rates.forward, -(rates.desorption + rates.loss)],
        ],
        dtype=float,
    )


def analytic_trace_determinant(rates: Rates, mu: float) -> tuple[float, float]:
    trace = -(
        rates.diffusion * mu
        + rates.forward
        + rates.desorption
        + rates.loss
    )
    determinant = (
        rates.diffusion * mu * (rates.desorption + rates.loss)
        + rates.forward * rates.loss
    )
    return trace, determinant


def classify_eigenvalues(eigenvalues: np.ndarray, tol: float = 1.0e-11) -> str:
    max_real = float(np.max(np.real(eigenvalues)))
    if max_real < -tol:
        return "strictly_asymptotically_stable"
    if max_real <= tol:
        return "semistable_or_neutral"
    return "unstable"


def build_radial_diffusion_matrix(
    *, nr: int, radius: float, diffusion: float, permeability: float
) -> tuple[np.ndarray, float]:
    """Repository finite-volume/ghost-point diffusion operator for perturbations.

    At r=R the non-homogeneous external concentration disappears after
    subtracting the steady state, leaving the homogeneous Robin perturbation
    boundary used here.
    """
    if nr < 3:
        raise ValueError("nr must be >= 3")
    if radius <= 0.0 or diffusion <= 0.0 or permeability < 0.0:
        raise ValueError("radius/diffusion must be positive; permeability >= 0")

    r = np.linspace(0.0, radius, nr)
    dr = float(r[1] - r[0])
    op = np.zeros((nr, nr), dtype=float)

    # r=0: L'Hopital limit used by the manuscript implementation.
    op[0, 0] = -2.0 * diffusion / dr**2
    op[0, 1] = 2.0 * diffusion / dr**2

    for i in range(1, nr - 1):
        ri = float(r[i])
        r_plus = ri + 0.5 * dr
        r_minus = ri - 0.5 * dr
        op[i, i + 1] = diffusion * r_plus / (ri * dr**2)
        op[i, i] = -diffusion * (r_plus + r_minus) / (ri * dr**2)
        op[i, i - 1] = diffusion * r_minus / (ri * dr**2)

    # Homogeneous perturbation ghost point:
    # u_ghost = u[N-1] - 2 dr P u[N] / D.
    i = nr - 1
    ri = float(r[i])
    r_plus = ri + 0.5 * dr
    r_minus = ri - 0.5 * dr
    op[i, i - 1] = diffusion * (r_plus + r_minus) / (ri * dr**2)
    op[i, i] = diffusion * (
        -r_plus * (1.0 + 2.0 * dr * permeability / diffusion) - r_minus
    ) / (ri * dr**2)
    return op, dr


def build_semidiscrete_operator(
    *, nr: int, radius: float, rates: Rates, permeability: float
) -> tuple[np.ndarray, float]:
    diffusion, dr = build_radial_diffusion_matrix(
        nr=nr,
        radius=radius,
        diffusion=rates.diffusion,
        permeability=permeability,
    )
    eye = np.eye(nr)
    operator = np.block(
        [
            [diffusion - rates.forward * eye, rates.desorption * eye],
            [rates.forward * eye, -(rates.desorption + rates.loss) * eye],
        ]
    )
    return operator, dr


def explicit_euler_diagnostic(eigenvalues: np.ndarray, dt: float) -> dict[str, Any]:
    amplification = np.abs(1.0 + dt * eigenvalues)
    max_amplification = float(np.max(amplification))
    negative_real = -np.real(eigenvalues[np.real(eigenvalues) < 0.0])
    conservative_real_axis_limit = (
        float(2.0 / np.max(negative_real)) if negative_real.size else math.inf
    )
    return {
        "dt": float(dt),
        "max_amplification": max_amplification,
        "stable": bool(max_amplification <= 1.0 + 1.0e-10),
        "conservative_real_axis_dt_limit": conservative_real_axis_limit,
    }


def run_audit() -> dict[str, Any]:
    # Manuscript/default radiodialysis values.
    rates = Rates(
        diffusion=1.0e-3,
        forward=0.05 * 1.0 + 0.02 * 0.3,
        desorption=0.005,
        loss=0.001,
    )

    mode_rows: list[dict[str, Any]] = []
    analytic_ok = True
    for mu in (0.0, 1.0e-6, 1.0, 100.0):
        matrix = mode_matrix(rates, mu)
        eig = np.linalg.eigvals(matrix)
        trace, det = analytic_trace_determinant(rates, mu)
        row = {
            "laplacian_eigenvalue": mu,
            "trace_analytic": trace,
            "trace_numeric": float(np.trace(matrix)),
            "determinant_analytic": det,
            "determinant_numeric": float(np.linalg.det(matrix)),
            "eigenvalues": [
                {"real": float(v.real), "imag": float(v.imag)} for v in eig
            ],
            "classification": classify_eigenvalues(eig),
        }
        row_ok = (
            math.isclose(trace, float(np.trace(matrix)), rel_tol=1e-12, abs_tol=1e-12)
            and math.isclose(det, float(np.linalg.det(matrix)), rel_tol=1e-10, abs_tol=1e-12)
            and row["classification"] == "strictly_asymptotically_stable"
        )
        row["contract_passed"] = row_ok
        analytic_ok = analytic_ok and row_ok
        mode_rows.append(row)

    # Edge case disproving an unconditional *strict* stability statement:
    # without irreversible loss, the spatially uniform mode conserves total
    # mobile+immobile material and has one zero eigenvalue.
    no_loss = Rates(diffusion=1.0e-3, forward=0.056, desorption=0.005, loss=0.0)
    no_loss_eig = np.linalg.eigvals(mode_matrix(no_loss, 0.0))
    no_loss_class = classify_eigenvalues(no_loss_eig)
    edge_case_ok = no_loss_class == "semistable_or_neutral" and np.any(
        np.isclose(no_loss_eig, 0.0, atol=1.0e-11)
    )

    # TWO GEOMETRIES, because the struct default and the as-run value differ and
    # only one of them is the grid an explicit step is ever taken on.
    #
    #   radius = 1.0   RadiolysisParams' own default. NO CALL SITE USES IT. It is
    #                  also default_parms()' grid in biofilms_radiodialysis.R,
    #                  which integrates with deSolve::ode(method="lsoda") —
    #                  adaptive and stiff-capable, taking no explicit step at all.
    #   radius = 20.0  what every Julia call site passes: R = N/2 in LATTICE
    #                  units at N = 40 (biofilms_potts.jl:1441, :1807;
    #                  biofilms_potts_jacc.jl:387). This is the grid the explicit
    #                  stepper actually runs, and the grid Figures 3 and 4 ran on.
    #
    # Reporting only the first would say dt = 0.5 is unstable for a solver that
    # never sees that grid. Reporting only the second would hide the trap: the
    # explicit path is stable ONLY BECAUSE r is in lattice units while D_eff is
    # in cm^2/s — the documented units error behind RADIODIALYSIS: BLOCKED.
    # Correct the units so R = 1.0 cm and dt = 0.5 becomes the unstable case.
    geometries = []
    for label, radius, in_use in (
        ("struct_default_unused", 1.0, False),
        ("as_run_lattice_units", 20.0, True),
    ):
        op, dr_g = build_semidiscrete_operator(
            nr=40, radius=radius, rates=rates, permeability=0.01)
        eig = np.linalg.eigvals(op)
        cls = classify_eigenvalues(eig)
        geometries.append({
            "label": label,
            "radius": radius,
            "dr": dr_g,
            "explicit_euler_path_runs_here": in_use,
            "max_real_eigenvalue": float(np.max(np.real(eig))),
            "min_real_eigenvalue": float(np.min(np.real(eig))),
            "classification": cls,
            "contract_passed": cls == "strictly_asymptotically_stable",
            "explicit_euler": explicit_euler_diagnostic(eig, dt=0.5),
        })

    # The semidiscrete CONTRACT must hold on both; only the as-run geometry
    # licenses a statement about the time step actually taken.
    semi_ok = all(g["contract_passed"] for g in geometries)
    as_run = next(g for g in geometries if g["explicit_euler_path_runs_here"])
    semi, dr = build_semidiscrete_operator(
        nr=40, radius=1.0, rates=rates, permeability=0.01)
    semi_eig = np.linalg.eigvals(semi)
    semi_class = classify_eigenvalues(semi_eig)
    euler = as_run["explicit_euler"]

    return {
        "schema_version": 1,
        "scope": "radiodialysis_linear_stability_and_time_step_contract",
        "interpretation": {
            "base_state": "perturbations about the boundary-forced steady state",
            "spatial_modes": (
                "non-negative eigenvalues of -Laplacian; finite-cylinder modes "
                "are not asserted to be plane waves"
            ),
            "continuous_vs_numerical": (
                "continuous/semidiscrete stability does not imply explicit-Euler stability"
            ),
        },
        "parameters": {
            **asdict(rates),
            "nr": 40,
            "radius": 1.0,
            "permeability": 0.01,
            "explicit_euler_dt": 0.5,
            "explicit_euler_dt_source": "biofilms_potts.jl:1208 dt_rd, s per CPM MCS",
            "production_substep_rule": "0.4 * dr^2 / (2 * D_eff), biofilms_potts.jl:1245",
        },
        "geometries": geometries,
        "analytical_modes": mode_rows,
        "loss_free_zero_mode": {
            "eigenvalues": [
                {"real": float(v.real), "imag": float(v.imag)} for v in no_loss_eig
            ],
            "classification": no_loss_class,
            "expected": "semistable_or_neutral",
            "contract_passed": bool(edge_case_ok),
        },
        "finite_volume_semidiscrete": {
            "dr": dr,
            "max_real_eigenvalue": float(np.max(np.real(semi_eig))),
            "min_real_eigenvalue": float(np.min(np.real(semi_eig))),
            "classification": semi_class,
            "contract_passed": bool(semi_ok),
        },
        "explicit_euler": euler,
        "overall_continuous_contract_passed": bool(analytic_ok and edge_case_ok and semi_ok),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("radiodialysis-stability.json"))
    parser.add_argument(
        "--require-explicit-euler-stability",
        action="store_true",
        help="Fail when the declared dt is outside the Euler stability region.",
    )
    args = parser.parse_args()

    report = run_audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    print(json.dumps(report, indent=2, sort_keys=True))
    if not report["overall_continuous_contract_passed"]:
        return 1
    if args.require_explicit_euler_stability and not report["explicit_euler"]["stable"]:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
