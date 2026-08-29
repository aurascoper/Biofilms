#!/usr/bin/env Rscript
# ============================================================================
#  Verify the slab (biofilm-depth) geometry of biofilms_radiodialysis.R.
#
#  This checks the R code that SHIPS.  Its numpy sibling,
#  analysis/verify_biofilm_depth_profile.py, checks an independent re-coding of
#  the same scheme; run both -- agreement between them is worth more than
#  either alone, and the numpy one runs where R is absent.
#
#  Loads only the model expressions, skipping the library() calls and the Shiny
#  app, so the stencil checks need no third-party packages at all.  The
#  time-integration check needs deSolve and reports SKIPPED without it rather
#  than being represented as passed (AGENTS.md rule 2).
#
#  WHAT IS CHECKED
#    1. slab_parms() refuses a defaulted X_total and names the gate.
#    2. The operator radiodialysis_rhs actually assembles equals an independent
#       assembly, node by node, including the Robin row.
#    3. The steady slab profile matches the closed-form Robin solution.
#    4. penetration_depth() reproduces lambda_steady/lambda_transient = sqrt(6).
#    5. REGRESSION: the cylindrical path -- which ships and is cross-validated
#       against biofilms_potts.jl:1153-1163 -- is byte-equivalent to the
#       pre-change inline arithmetic across every row the refactor touched.
#    6. NEGATIVE CONTROLS: check 5 is re-run under three mutations of
#       face_weights(), each of which must be caught.  Without this the suite
#       could pass on a stencil that ignored geometry (AGENTS.md rule 1) --
#       and an earlier draft of check 5 did exactly that, because in slab
#       geometry w_plus == w_minus == 1 makes a swap invisible.
#
#  WHAT IS *NOT* CHECKED
#    Whether the model predicts a MEASURABLE gradient in a real biofilm.  That
#    is downstream of RADIODIALYSIS: BLOCKED (README.md:350-353): U scales with
#    X_total and lambda ~ 1/sqrt(U), so a 20-100x correction to the biomass
#    basis grows lambda by 4.5-10x and flattens the profile.  Every X_total
#    below is a TEST VALUE chosen to exercise the numerics.  No profile number
#    produced here may be quoted as a claim about a biofilm.
# ============================================================================

MODEL <- file.path(dirname(dirname(normalizePath(
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))),
  "biofilms_radiodialysis.R")
if (!file.exists(MODEL)) MODEL <- "biofilms_radiodialysis.R"

WANTED <- c("radiodialysis_rhs", "face_weights", "default_parms",
            "slab_parms", "penetration_depth", "run_radiodialysis")
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

# ---------------------------------------------------------------------------
#  Mutation campaign, with the baseline guard built in.
#
#  A mutation harness has TWO ways to measure nothing, and they look opposite:
#
#    always-passes : no mutant is ever caught -- obviously broken, you notice.
#    always-FAILS  : every mutant is "caught" -- looks like perfect detection.
#
#  The second is the dangerous one, because universal failure is the shape you
#  are hoping for.  It happened here: an earlier reference matrix omitted the
#  D(w+ + w-)/dg^2 term in the Robin row, the baseline sat at 9.8e-01, and all
#  three mutants read as caught while the harness measured nothing at all.
#
#  So the baseline verdict GATES the mutant verdicts.  If an unmutated run does
#  not pass, no mutant result is reported as meaningful -- the campaign is
#  declared invalid instead.  Do not rely on remembering this.
# ---------------------------------------------------------------------------
mutation_campaign <- function(label, probe, tol, restore, mutants) {
  base <- probe()
  if (!(base < tol)) {
    report(FALSE, sprintf("%s: MUTATION HARNESS INVALID", label),
           sprintf("baseline residual %.2e is not < %.0e, so every mutant would
             read as caught regardless of the code. No mutant verdict reported.",
                   base, tol))
    return(invisible(FALSE))
  }
  cat(sprintf("\n  %s -- baseline clean at %.2e; each mutation must be caught:\n",
              label, base))
  for (nm in names(mutants)) {
    on.exit(restore(), add = TRUE)
    mutants[[nm]]()
    r <- tryCatch(probe(), error = function(e) Inf)
    restore()
    report(r > tol, sprintf("  caught: %s", nm), sprintf("%.2e", r))
  }
  report(probe() < tol, "  baseline restored after mutation")
}

# --- the operator radiodialysis_rhs assembles: rhs(c) = A c + b -------------
operator_of <- function(p) {
  Nr <- length(p$r_grid)
  zero <- c(rep(0, Nr), rep(0, Nr), 1.0)
  b <- radiodialysis_rhs(0, zero, p)[[1]][seq_len(Nr)]
  A <- sapply(seq_len(Nr), function(i) {
    y <- zero; y[i] <- 1.0
    radiodialysis_rhs(0, y, p)[[1]][seq_len(Nr)] - b
  })
  list(A = A, b = b)
}

# --- independent assembly from the pre-change inline arithmetic -------------
# Derived from the ghost substitution, not guessed: substituting
#   c_ghost = c[N-1] - 2 dr P (c[N] - c_ext)/D
# into the outer stencil gives the A[N,N] and b[N] below.  An earlier draft
# omitted the D(w+ + w-)/dr^2 term here; the baseline then failed at 9.8e-01
# and every mutant looked "caught" for the wrong reason.
reference <- function(p, geom) {
  g <- p$r_grid; N <- length(g); dg <- g[2] - g[1]; D <- p$D_eff
  U <- p$k_ads * p$X_total + p$k_red * p$X_red
  bc <- if (identical(geom, "slab")) p$k_L else p$P0
  wp <- if (identical(geom, "slab")) rep(1, N) else (g + 0.5 * dg) / g
  wm <- if (identical(geom, "slab")) rep(1, N) else (g - 0.5 * dg) / g
  A <- matrix(0, N, N); b <- numeric(N)
  for (i in 2:(N - 1)) {
    A[i, i - 1] <- D * wm[i] / dg^2
    A[i, i]     <- -(D * (wp[i] + wm[i]) / dg^2 + U)
    A[i, i + 1] <- D * wp[i] / dg^2
  }
  A[N, N - 1] <- D * (wp[N] + wm[N]) / dg^2
  A[N, N]     <- -(D * (wp[N] + wm[N]) / dg^2 + 2 * wp[N] * bc / dg + U)
  b[N]        <- 2 * wp[N] * bc * p$c_ext / dg
  list(A = A, b = b)
}

residual <- function(p, geom) {
  o <- operator_of(p); r <- reference(p, geom); rows <- 2:nrow(o$A)
  max(abs(o$A[rows, ] - r$A[rows, ])) / max(abs(r$A[rows, ]))
}

# --- 1. the refuse-to-default guard ----------------------------------------
report(tryCatch({ slab_parms(); FALSE },
                error = function(e) grepl("RADIODIALYSIS: BLOCKED",
                                          conditionMessage(e), fixed = TRUE)),
       "slab_parms() refuses a defaulted X_total and names the gate")

# --- 1b. the geometry dispatch refuses what it does not recognise ----------
# Codex P2 on #19.  The old `if (geom == "slab") k_L else P_eff` read every
# unrecognised value as cylindrical.  Both halves are asserted: an unknown
# geom must ERROR, and BOTH supported values must still run -- a dispatch that
# refuses everything would satisfy the first half alone.
local({
  p_bad <- slab_parms(X_total = 1.0); p_bad$geom <- "slabb"
  y <- c(rep(0, length(p_bad$r_grid)), rep(0, length(p_bad$r_grid)), 1)
  refused <- tryCatch({ radiodialysis_rhs(0, y, p_bad); FALSE },
                      error = function(e) grepl("unsupported geom",
                                                conditionMessage(e), fixed = TRUE))
  accepts <- vapply(list(slab_parms(X_total = 1.0), default_parms()), function(p) {
    yy <- c(rep(0, length(p$r_grid)), rep(0, length(p$r_grid)), 1)
    tryCatch({ radiodialysis_rhs(0, yy, p); TRUE }, error = function(e) FALSE)
  }, logical(1))
  report(refused && all(accepts),
         "geom dispatch refuses an unknown value and still accepts both known ones",
         sprintf("refused=%s slab=%s cylindrical=%s",
                 refused, accepts[1], accepts[2]))
})

# --- 2-3. slab operator and steady profile ---------------------------------
L_f <- 0.01; D <- 1e-5; N <- 40; kL <- 1e2
k_eff <- 0.056 * 0.001 / 0.006
p <- slab_parms(Nr = N, L_f = L_f, X_total = 1.0, D_eff = D, k_L = kL)
p$k_ads <- k_eff; p$k_red <- 0; p$k_des <- 0; p$k_loss <- 0; p$X_red <- 0
res <- residual(p, "slab")
report(res < 1e-13, "slab operator matches independent assembly",
       sprintf("%.2e", res))

op <- operator_of(p)
c_num <- solve(op$A, -op$b)
lam <- sqrt(D / k_eff)
c_exact <- p$c_ext * cosh(p$r_grid / lam) /
  (cosh(L_f / lam) + (D / (kL * lam)) * sinh(L_f / lam))
err <- max(abs(c_num - c_exact) / c_exact)
report(err < 1e-5, "steady slab profile matches the analytic Robin solution",
       sprintf("max rel err %.2e, lambda %.1f um, phi %.3f, c(0)/c(L) %.4f",
               err, lam * 1e4, L_f / lam, c_num[1] / c_num[N]))

# --- 4. the structural sqrt(6) ---------------------------------------------
p2 <- slab_parms(Nr = N, L_f = L_f, X_total = 1.0, D_eff = D)
ratio <- penetration_depth(p2, "steady") / penetration_depth(p2, "transient")
report(abs(ratio - sqrt(6)) < 1e-12,
       "lambda_steady/lambda_transient = sqrt(6), independent of X_total",
       sprintf("%.9f", ratio))

# --- 5. cylindrical regression ---------------------------------------------
cyl <- function() default_parms(Nr = 40, R = 1.0)
res_c <- residual(cyl(), "cylindrical")
report(res_c < 1e-13,
       "cylindrical path byte-equivalent to the pre-change stencil",
       sprintf("%.2e", res_c))

# --- 6. negative controls ---------------------------------------------------
orig <- face_weights
mutation_campaign(
  label   = "negative controls",
  probe   = function() residual(cyl(), "cylindrical"),
  tol     = 1e-13,
  restore = function() face_weights <<- orig,
  mutants = list(
    "face_weights swapped" = function()
      face_weights <<- function(g, gm) {
        w <- orig(g, gm); list(w_plus = w$w_minus, w_minus = w$w_plus) },
    "w_plus[5] perturbed by 0.1%" = function()
      face_weights <<- function(g, gm) {
        w <- orig(g, gm); w$w_plus[5] <- w$w_plus[5] * 1.001; w },
    "weights forced to slab (= 1)" = function()
      face_weights <<- function(g, gm)
        list(w_plus = rep(1, length(g)), w_minus = rep(1, length(g)))
  )
)

# The gate itself needs a negative control, or it is just another check that
# cannot fail.  Deliberately break the reference the campaign measures against
# and confirm the campaign REFUSES rather than reporting three caught mutants.
local({
  good <- reference
  reference <<- function(p, geom) {          # drop the Robin diffusion term:
    r <- good(p, geom); N <- nrow(r$A)       # exactly the bug that occurred
    g <- p$r_grid; dg <- g[2] - g[1]
    wp <- (g[N] + 0.5 * dg) / g[N]; wm <- (g[N] - 0.5 * dg) / g[N]
    r$A[N, N] <- r$A[N, N] + p$D_eff * (wp + wm) / dg^2
    r
  }
  before <- failures
  invisible(capture.output(mutation_campaign(
    "self-test", function() residual(cyl(), "cylindrical"), 1e-13,
    function() NULL, list("noop" = function() NULL))))
  reference <<- good
  fired <- failures > before
  failures <<- before                        # do not count the deliberate break
  report(fired,
         "  the baseline gate itself fires on a broken reference")
})

# --- 7. time integration (needs deSolve) ------------------------------------
cat("\n")
if (!requireNamespace("deSolve", quietly = TRUE)) {
  skipped <- c(skipped, "time integration (deSolve not installed)")
  cat("[SKIP] time integration: deSolve not installed.",
      "This is uncovered surface, not a pass.\n")
} else {
  library(deSolve)
  ps <- slab_parms(Nr = N, L_f = L_f, X_total = 1.0, D_eff = D, k_L = kL)
  out <- run_radiodialysis(ps, t_end = 2e5, n_out = 50)
  final <- out$c_mat[nrow(out$c_mat), ]
  # full system: sorption chain on, so the steady sink is k_eff
  ke <- (ps$k_ads * ps$X_total + ps$k_red * ps$X_red) *
        ps$k_loss / (ps$k_des + ps$k_loss)
  lm2 <- sqrt(ps$D_eff / ke)
  ex <- ps$c_ext * cosh(ps$r_grid / lm2) /
    (cosh(L_f / lm2) + (ps$D_eff / (kL * lm2)) * sinh(L_f / lm2))
  e2 <- max(abs(final - ex) / ex)
  report(e2 < 1e-3,
         "lsoda steady state matches the analytic Robin solution",
         sprintf("max rel err %.2e at t=2e5 s", e2))
  report(abs(out$m_vec[length(out$m_vec)] - 1.0) < 1e-12,
         "membrane integrity m stays 1 in slab geometry (radiation off)")
}

# A skip is a question, not a pass (AGENTS.md rule 2).  Never print a bare
# "ALL PASS" over a suite that did not run everything: the summary line is what
# gets quoted, and "215 passed, 8 skipped" is how a dead function shipped here.
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

# Machine-readable receipt.  CI must assert on THIS, not on the exit status:
# a job whose R setup step failed or was omitted can exit 0 having run nothing,
# which is the ALL-PASS-over-skip bug relocated from R into YAML.
args <- commandArgs(TRUE)
if (length(args) >= 2L && args[1] == "--report") {
  writeLines(sprintf(
    '{"checks_run": %d, "failures": %d, "skipped": %d, "skips": [%s], "complete": %s}',
    checks_run, failures, length(skipped),
    paste(sprintf('"%s"', skipped), collapse = ", "),
    tolower(as.character(failures == 0L && !length(skipped)))), args[2])
}
quit(status = if (failures == 0L && !length(skipped)) 0L else 1L)
