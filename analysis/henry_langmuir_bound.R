#!/usr/bin/env Rscript
# ============================================================================
#  Where the shipped Henry sorption term stops agreeing with Langmuir.
#
#    Rscript analysis/henry_langmuir_bound.R
#    Rscript analysis/henry_langmuir_bound.R --report out.json
#
#  WHY THIS FILE EXISTS.  Section 3.12 states that three sorption forms are in
#  play -- irreversible sink, Henry, Langmuir -- and that "they are not
#  interchangeable".  That claim is structural and correct.  Nowhere in this
#  repository was it ever said BY HOW MUCH, over the concentration range the
#  solver actually integrates.  This produces that number, so the ledger row
#  quoting it can be re-run rather than believed.  Codex raised a P1 on pull
#  request #23 -- "Ship the producer for the restated proposal fraction" -- for
#  a published number nobody could reproduce; this is written to not be that.
#
#  WHAT IS COMPUTED
#    The shipped immobile-phase equation (biofilms_radiodialysis.R:144) is
#        ds/dt = U c - (k_des + k_loss) s,     U = X_total (k_ads + k_red f_red)
#    first-order in c with no site limitation: a Henry isotherm.  The kinetic
#    Langmuir form of the same reaction is
#        ds/dt = U c (1 - s/X_max) - (k_des + k_loss) s
#    whose equilibrium satisfies  s_L = s_H / (1 + s_H / X_max).  So the
#    relative error of the shipped form is  s_H / (X_max + s_H), and the two
#    diverge by more than eps once  c_ref > X_max (1 - eps) / (eps s_H_norm).
#
#  WHAT IS *NOT* CONCLUDED, AND THE DISTINCTION IS THE POINT.  This does NOT
#  say the Henry form is adequate.  s is normalised -- c_ext = 1.0, and section
#  3.12 says the sorbate has "no chemical identity" -- so converting s to a
#  capacity needs a reference concentration that, per
#  data/calibration/suspended_isotherm_proposal.csv, exists nowhere in this
#  repository.  Nothing establishes that a trace-contaminant c_ref is the
#  relevant one.  The result is therefore a THRESHOLD and an undetermined side:
#  the forms diverge above c_ref = X_max / K, and whether this solver's regime
#  sits above or below that is undetermined.  Written this way so it cannot be
#  quoted as "the Henry form is fine".
#
#  X_max IS A PRIOR ON A PRIOR.  X_max = q_max * rho_dry multiplies a
#  suspended-MEASURED capacity by a biofilm-density PRIOR, so the product
#  inherits the weaker class (suspended_isotherm_proposal.csv:42-44).  Every
#  crossover printed below is therefore a prior, not a measurement.
#
#  WHAT THE EQUILIBRIUM CANCELLATION DOES AND DOES NOT BUY.  s_H is
#  proportional to X_total and X_max = q_max * rho_dry, so at equilibrium
#  X_total cancels against rho_dry and the ratio reduces to
#  (k_ads + k_red f_red) c c_ref / (L q_max).  That makes the EQUILIBRIUM
#  crossover immune to PP-66-08, the open needs_calibration row saying X_total
#  is a site-occupancy fraction in a g/cm^3 slot so the rate is "wrong by an
#  unstated rho".  It cancels X_total against rho_dry and does NOTHING for the
#  q_max uncertainty.  The transient crossover cancels only partly, so both are
#  reported and the difference between them is the size of that partiality.
#
#  Loads only the model expressions, skipping library() calls and the Shiny app
#  at biofilms_radiodialysis.R:611, the same idiom as
#  analysis/verify_biofilm_depth_profile.R.  deSolve is required; without it
#  every check reports SKIPPED rather than being represented as passed
#  (AGENTS.md rule 2).
# ============================================================================

MODEL <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  p <- if (length(a)) file.path(dirname(dirname(normalizePath(
         sub("^--file=", "", a[1])))), "biofilms_radiodialysis.R") else ""
  if (nzchar(p) && file.exists(p)) p else "biofilms_radiodialysis.R"
})

WANTED <- c("radiodialysis_rhs", "face_weights", "default_parms",
            "slab_parms", "penetration_depth", "run_radiodialysis",
            "uptake_rate_of")
for (e in parse(MODEL)) {
  if (is.call(e) && identical(as.character(e[[1]]), "<-") &&
      as.character(e[[2]]) %in% WANTED) eval(e, envir = globalenv())
}
stopifnot(all(vapply(WANTED, exists, logical(1))))

failures   <- 0L
skipped    <- character(0)
checks_run <- 0L
report <- function(ok, msg, detail = "") {
  checks_run <<- checks_run + 1L
  if (!ok) failures <<- failures + 1L
  cat(sprintf("[%s] %s%s\n", if (ok) "PASS" else "FAIL", msg,
              if (nzchar(detail)) paste0("  (", detail, ")") else ""))
  invisible(ok)
}

# ATTACHED, not merely available: run_radiodialysis() calls `ode` unqualified,
# so requireNamespace() alone leaves it unresolvable and the model errors out
# mid-run rather than reporting a clean SKIP.
HAVE_DESOLVE <- suppressWarnings(require(deSolve, quietly = TRUE,
                                         warn.conflicts = FALSE))

# --- the shipped constants, read from the model rather than restated ---------
p  <- default_parms()
U  <- p$X_total * (p$k_ads + p$k_red * p$f_red_active)
L  <- p$k_des + p$k_loss
RATIO_EQ <- U / L                      # Henry equilibrium s/c

cat("=== shipped constants (read from default_parms()) ===\n")
cat(sprintf("  k_ads %.4g  k_red %.4g  k_des %.4g  k_loss %.4g\n",
            p$k_ads, p$k_red, p$k_des, p$k_loss))
cat(sprintf("  X_total %.4g  f_red_active %.4g  c_ext %.4g\n",
            p$X_total, p$f_red_active, p$c_ext))
cat(sprintf("  U = X_total*(k_ads + k_red*f_red) = %.4f 1/s\n", U))
cat(sprintf("  k_des + k_loss                    = %.4f 1/s\n", L))
cat(sprintf("  Henry equilibrium ratio s_eq/c    = %.4f\n\n", RATIO_EQ))

# --- the visited range ------------------------------------------------------
T_END <- 100
if (!HAVE_DESOLVE) {
  skipped <- c(skipped, "every check: deSolve absent, no trajectory to measure")
  s_max <- NA_real_; c_max <- NA_real_
} else {
  o <- run_radiodialysis(p, t_end = T_END, n_out = 200)
  c_max <- max(o$c_mat); s_max <- max(o$s_mat)
  cat("=== range the solver actually visits ===\n")
  cat(sprintf("  c: min %.6f  max %.6f   (c_ext = %.1f)\n",
              min(o$c_mat), c_max, p$c_ext))
  cat(sprintf("  s: min %.6f  max %.6f\n", min(o$s_mat), s_max))
  cat(sprintf("  fraction of Henry equilibrium reached: %.3f",
              s_max / (RATIO_EQ * c_max)))
  cat(sprintf("   (relaxation 1/L = %.1f s vs t_end = %d s)\n\n", 1 / L, T_END))
}

# --- the bound --------------------------------------------------------------
# eps-divergence crossover.  With s_H(physical) = s_norm * c_ref,
#   rel_err = s_norm c_ref / (X_max + s_norm c_ref) > eps
#     <=>  s_norm c_ref (1 - eps) > eps X_max
#     <=>  c_ref > X_max / K,     K = (1 - eps) s_norm / eps
# K is what the ledger row quotes, so it is derived here once and control 2
# pins it from both sides rather than trusting this line.
crossover_K <- function(s_norm, eps) ((1 - eps) * s_norm) / eps
EPS <- 0.10

if (HAVE_DESOLVE) {
  K_transient   <- crossover_K(s_max, EPS)
  K_equilibrium <- crossover_K(RATIO_EQ * c_max, EPS)
  cat("=== the bound: forms diverge by more than 10% above c_ref = X_max / K ===\n")
  cat(sprintf("  K (transient, t_end = %d s)   = %.2f\n", T_END, K_transient))
  cat(sprintf("  K (equilibrium limit)         = %.2f", K_equilibrium))
  cat("   <- immune to PP-66-08; q_max uncertainty unaffected\n\n")
  cat(sprintf("  %-14s %-12s %s\n", "c_ref (mg/L)", "X_max = 50", "X_max = 5"))
  for (cr_L in c(1, 10, 100, 1000, 10000)) {
    cr <- cr_L / 1000                                   # mg/L -> mg/cm^3
    e <- function(xm) 100 * s_max * cr / (xm + s_max * cr)
    cat(sprintf("  %-14g %-12s %s\n", cr_L,
                sprintf("%.3f%%", e(50)), sprintf("%.3f%%", e(5))))
  }
  cat(sprintf("\n  crossover: X_max = 50 -> c_ref = %.1f mg/L; X_max = 5 -> %.1f mg/L\n\n",
              1000 * 50 / K_transient, 1000 * 5 / K_transient))
}

# ============================================================================
#  CONTROL 1 -- the closed form against a NUMERICAL SOLVE, not against algebra.
#
#  s_L = s_H/(1 + s_H/X_max) is a derivation.  Checking it by re-deriving it
#  proves nothing: a control that never meets the pipeline confirms a guard
#  against an input the pipeline cannot produce, which is the defect AGENTS.md
#  rule 1 now records from the figure-sidecar case.  So the Langmuir RHS is
#  integrated on the same grid with the same constants and its equilibrium is
#  compared against the closed form.
# ============================================================================
langmuir_equilibrium_numeric <- function(parms, X_max, c_fixed, t_end = 5000) {
  # Well-mixed single node at fixed c: ds/dt = U c (1 - s/X_max) - L s.
  rhs <- function(t, y, pr) list(U * c_fixed * (1 - y[1] / X_max) - L * y[1])
  out <- deSolve::ode(y = c(s = 0), times = c(0, t_end), func = rhs,
                      parms = NULL, method = "lsoda",
                      rtol = 1e-10, atol = 1e-12)
  as.numeric(out[nrow(out), 2])
}

if (HAVE_DESOLVE) {
  cat("=== control 1: closed form vs numerical Langmuir integration ===\n")
  worst <- 0
  for (X_max in c(0.5, 5, 50)) {
    for (cc in c(0.1, 0.5, 0.925)) {
      sH    <- RATIO_EQ * cc                       # Henry equilibrium at this c
      closed <- sH / (1 + sH / X_max)
      numer  <- langmuir_equilibrium_numeric(p, X_max, cc)
      rel    <- abs(closed - numer) / numer
      worst  <- max(worst, rel)
      cat(sprintf("    X_max %5.1f  c %5.3f : closed %.6f  numeric %.6f  rel %.2e\n",
                  X_max, cc, closed, numer, rel))
    }
  }
  report(worst < 1e-6,
         "closed form s_L = s_H/(1 + s_H/X_max) matches numerical integration",
         sprintf("worst relative disagreement %.2e over 9 (X_max, c) pairs", worst))
} else {
  report(FALSE, "control 1 could not run")   # unreachable: skip path set above
}

# ============================================================================
#  CONTROL 2 -- TWO-SIDED on the crossover.
#
#  A one-sided assertion ("error exceeds 10% somewhere above") passes on any
#  monotone curve regardless of where the crossover sits.  Both sides are
#  asserted so the check fails if the threshold moves in EITHER direction.
# ============================================================================
if (HAVE_DESOLVE) {
  cat("\n=== control 2: the crossover is pinned from both sides ===\n")
  X_max <- 5
  c_star <- X_max / K_transient
  rel <- function(cr) s_max * cr / (X_max + s_max * cr)
  below <- rel(c_star * 0.5)
  above <- rel(c_star * 2.0)
  at    <- rel(c_star)
  cat(sprintf("    c_ref = 0.5 c*  -> %.4f   (must be < 0.10)\n", below))
  cat(sprintf("    c_ref =     c*  -> %.4f   (must be ~ 0.10)\n", at))
  cat(sprintf("    c_ref = 2.0 c*  -> %.4f   (must be > 0.10)\n", above))
  report(below < EPS && above > EPS && abs(at - EPS) < 1e-9,
         "crossover c* = X_max/K separates <10% from >10% in both directions",
         sprintf("below %.4f, at %.4f, above %.4f", below, at, above))
}

# ============================================================================
#  CONTROL 3 -- the conclusion is CONDITIONAL, asserted rather than asserted-of.
#
#  The ledger row says the bound depends on c_ref and that the solver's regime
#  is undetermined.  If the divergence were large for every plausible c_ref
#  that would be a fact about the FORMS and would belong in the manuscript
#  instead.  This fails if the span ever stops straddling the threshold, which
#  is the condition under which the row's framing would have to change.
# ============================================================================
if (HAVE_DESOLVE) {
  cat("\n=== control 3: the result is conditional, not unconditional ===\n")
  lo <- s_max * (1 / 1000)  / (50 + s_max * (1 / 1000))     # 1 mg/L,  X_max 50
  hi <- s_max * (10000/1000) / (5 + s_max * (10000/1000))   # 10 g/L,  X_max 5
  cat(sprintf("    smallest case (1 mg/L, X_max 50): %.5f%%\n", 100 * lo))
  cat(sprintf("    largest  case (10 g/L, X_max  5): %.2f%%\n",  100 * hi))
  report(lo < 0.01 && hi > 0.5,
         "divergence straddles the threshold, so the bound is conditional on c_ref",
         sprintf("%.4f%% to %.1f%% -- a fact about a missing parameter, not the forms",
                 100 * lo, 100 * hi))
}

# --- verdict ----------------------------------------------------------------
# A skip is uncovered surface, not a neutral fact (AGENTS.md rule 2), so a run
# that skipped is never reported as a clean pass.
verdict <- if (failures > 0L) {
  sprintf("FAILURES: %d", failures)
} else if (length(skipped)) {
  "PASSED WHAT RAN -- NOT A CLEAN RUN"
} else {
  "ALL PASS"
}
cat(sprintf("\n%s (%d check%s run, %d failure%s, %d skipped)\n", verdict,
            checks_run, if (checks_run == 1L) "" else "s", failures,
            if (failures == 1L) "" else "s", length(skipped)))
for (s in skipped) cat("  UNCOVERED:", s, "\n")

args <- commandArgs(TRUE)
if (length(args) >= 2L && args[1] == "--report") {
  # `skips` is emitted even when empty, so the CI receipt gate can name what was
  # uncovered rather than reporting a count with nothing behind it.
  writeLines(sprintf(
    '{"checks_run": %d, "failures": %d, "skipped": %d, "skips": [%s], "K_transient": %s, "K_equilibrium": %s, "complete": %s}',
    checks_run, failures, length(skipped),
    paste(sprintf('"%s"', skipped), collapse = ", "),
    if (HAVE_DESOLVE) sprintf("%.4f", K_transient) else "null",
    if (HAVE_DESOLVE) sprintf("%.4f", K_equilibrium) else "null",
    tolower(as.character(failures == 0L && !length(skipped)))), args[2])
}
quit(status = if (failures == 0L && !length(skipped)) 0L else 1L)
