#!/usr/bin/env python3
"""Reproduce §3.4's overdamped-regime numbers.

    python3 analysis/overdamped_regime.py
    python3 analysis/overdamped_regime.py --report out.json

WHY THIS FILE EXISTS. Section 3.4 states that no independently-integrated momentum
state exists, and now states the stronger claim that one would be physically
inappropriate: at cell scale the regime is overdamped by roughly ten orders of
magnitude. Those are published numbers, and a published number nobody can re-run is
the defect PP-62-11 records and the P1 Codex raised on pull request #23.

WHAT IS COMPUTED, AND ITS INPUTS ARE ASSUMPTIONS RATHER THAN MEASUREMENTS.

    Re    = rho * U * L / eta          inertial vs viscous forces
    tau_p = 2 * rho * a^2 / (9 * eta)  inertial (momentum) relaxation time

rho and eta are water at room temperature, carried from the literature. L, U and a
are ASSUMED scales for a bacterial cell and a biofilm feature; none is measured in
this repository, and nothing here should be read as a measurement of this system.
Re is a property of the flow and not of the organism, so every figure is reported
with the length and speed that produced it -- a Reynolds number quoted without them
is not reproducible, which is the pinned-inputs discipline 0a exists to enforce.

WHAT THE CLAIM ACTUALLY IS, BECAUSE THE OBVIOUS CONTROL TESTS THE WRONG QUANTITY.
The claim is not "Re < 1". It is that the inertial memory is negligible against the
timescales the model addresses -- that T_bio / tau_p is enormous. Those come apart:
a millimetre object at cm/s gives Re = 10, which looks like a failing input, while
its T_bio/tau_p is 1.9e4 -- four hundred times larger than the genuine control's 49.
So crossing Re = 1 says little about where the ratio sits, and a control built on Re
would be testing a different quantity from the one claimed. The control below
straddles the RATIO, by taking an object large enough that tau_p approaches a
biological timescale.
"""
from __future__ import annotations

import json
import sys

RHO = 1.0e3      # kg m^-3, water
ETA = 1.0e-3     # Pa s,    water at room temperature
T_BIO = 1080.0   # s, the 18-minute biofilm elastic relaxation time (literature)

# The ratio below which the overdamped claim would not hold as stated. Ten orders is
# what §3.4 says; 1e8 is the floor asserted, so the claim has margin rather than
# sitting on its own threshold.
RATIO_FLOOR = 1.0e8

checks_run = 0
failures = 0


def report(ok: bool, msg: str, detail: str = "") -> bool:
    global checks_run, failures
    checks_run += 1
    if not ok:
        failures += 1
    print(f"[{'PASS' if ok else 'FAIL'}] {msg}" + (f"  ({detail})" if detail else ""))
    return ok


def reynolds(L: float, U: float) -> float:
    return RHO * U * L / ETA


def tau_p(a: float) -> float:
    return 2.0 * RHO * a * a / (9.0 * ETA)


CASES = [  # (label, L m, U m/s)
    ("cell, 1 um at 20 um/s",            1e-6,  20e-6),
    ("cell, 10 um at 20 um/s",           10e-6, 20e-6),
    ("biofilm feature, 100 um at 20 um/s", 100e-6, 20e-6),
]
RADII = [("0.5 um", 0.5e-6), ("1 um", 1e-6), ("100 um", 100e-6)]

print("=== inputs (assumed scales, not measurements) ===")
print(f"  rho = {RHO:.3g} kg/m3   eta = {ETA:.3g} Pa s   T_bio = {T_BIO:.0f} s (18 min)\n")

print("=== Reynolds number, each with the L and U that produced it ===")
for label, L, U in CASES:
    print(f"  {label:36} Re = {reynolds(L, U):.2e}")

print("\n=== inertial relaxation and the ratio the claim is about ===")
for label, a in RADII:
    t = tau_p(a)
    print(f"  a = {label:8} tau_p = {t:.2e} s   T_bio/tau_p = {T_BIO / t:.2e}")

# --- the assertions -------------------------------------------------------
report(all(reynolds(L, U) < 1.0 for _, L, U in CASES),
       "every assumed scale is below Re = 1",
       f"largest {max(reynolds(L, U) for _, L, U in CASES):.2e}")

cell_ratios = [T_BIO / tau_p(a) for label, a in RADII if label != "100 um"]
report(min(cell_ratios) >= RATIO_FLOOR,
       f"cell-scale inertial memory is at least {RATIO_FLOOR:.0e} below T_bio",
       f"smallest ratio {min(cell_ratios):.2e}")

# --- CONTROL: straddle the RATIO, not the Reynolds number ------------------
print("\n=== control: an input where tau_p approaches a biological timescale ===")
a_ctrl = 1e-2                      # 1 cm
t_ctrl = tau_p(a_ctrl)
r_ctrl = T_BIO / t_ctrl
print(f"  a = 1 cm    tau_p = {t_ctrl:.2f} s   T_bio/tau_p = {r_ctrl:.1f}")

# and the input that LOOKS like a control and is not, kept because the distinction
# is the point: it fails Re < 1 while its ratio stays far above the real control's.
L_bad, U_bad, a_bad = 1e-3, 1e-2, 0.5e-3
print(f"  the tempting control (1 mm at 1 cm/s): Re = {reynolds(L_bad, U_bad):.1f} "
      f"but T_bio/tau_p = {T_BIO / tau_p(a_bad):.1e}")

report(r_ctrl < RATIO_FLOOR,
       "the control input fails the ratio floor, so the assertion is not vacuous",
       f"1 cm gives {r_ctrl:.1f}, far below {RATIO_FLOOR:.0e}")
# CORRECTED WHILE WRITING THIS: the first version asserted that the Re>1 input
# still passes RATIO_FLOOR, and it does not -- 1.9e4 is well under 1e8, and the
# assertion failed on its first run. The overstatement is instructive and the
# weaker true claim is the one that makes the point: crossing Re = 1 says little
# about where the RATIO sits, because this input's ratio is still hundreds of
# times larger than the genuine control's. Re and the ratio are different
# quantities; that is why the control straddles the ratio.
ratio_bad = T_BIO / tau_p(a_bad)
report(ratio_bad > 100.0 * r_ctrl,
       "a Re>1 input's ratio is still orders above the control's, so Re is the wrong control",
       f"Re = {reynolds(L_bad, U_bad):.0f}, ratio {ratio_bad:.1e} vs control {r_ctrl:.1f} "
       f"({ratio_bad / r_ctrl:.0f}x)")

verdict = "ALL PASS" if failures == 0 else f"FAILURES: {failures}"
print(f"\n{verdict} ({checks_run} checks run, {failures} failure"
      f"{'' if failures == 1 else 's'}, 0 skipped)")

if len(sys.argv) >= 3 and sys.argv[1] == "--report":
    with open(sys.argv[2], "w", encoding="utf-8") as fh:
        # THE RECEIPT CARRIES EVERY COMPUTED VALUE, NOT A HAND-PICKED SUBSET.
        # It used to emit three scalars -- Re_cell_1um, tau_p_0p5um, ratio_0p5um --
        # while section 3.4 states FIVE numbers, so 2e-4, 2e-3 and 2.2e-7 existed in
        # stdout only. A guard written from that receipt would have covered two of
        # five and looked complete. That is the same defect as asserting on a union
        # where categories were computed, moved to the SERIALISATION boundary: the
        # discriminating data was computed and thrown away on the way out. Serialise
        # the sets CASES and RADII already define, so a sixth stated number cannot be
        # missed by a receipt that never mentioned it.
        json.dump({"checks_run": checks_run, "failures": failures, "skipped": 0,
                   "skips": [],
                   "inputs": {"rho": RHO, "eta": ETA, "T_bio": T_BIO},
                   "reynolds": {label: reynolds(L, U) for label, L, U in CASES},
                   "lengths": {label: L for label, L, _ in CASES},
                   "speeds": {label: U for label, _, U in CASES},
                   "radii": dict(RADII),
                   "tau_p": {label: tau_p(a) for label, a in RADII},
                   "ratio": {label: T_BIO / tau_p(a) for label, a in RADII},
                   "complete": failures == 0}, fh)

sys.exit(0 if failures == 0 else 1)
