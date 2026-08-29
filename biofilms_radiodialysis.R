# ============================================================
#  Radiodialysis Membrane Transport — Cylindrical Bioreactor
#  "Modeling Radioresistance and Radiotropic Fitness" (Kinder & Faulkner, 2026)
#
#  Implements the minimal 3-equation PDE system derived from
#  the Deep Research synthesis (April 2026):
#
#  (1) Mobile species: cylindrical reaction-diffusion PDE
#        ∂c/∂t = (1/r) ∂/∂r (r D_eff ∂c/∂r)
#                - X (k_ads + k_red f_red) c + k_des s
#
#  (2) Immobile phase ODE (biosorption / bioreduction):
#        ∂s/∂t = X (k_ads + k_red f_red) c - (k_des + k_loss) s
#
#  (3) Membrane damage ODE (radiation-driven permeability):
#        dm/dt = -k_dam Ḋ(R) m
#        P_eff(t) = P₀ exp(α D_cum(t)),  D_cum(t) = Ḋ(R) · t
#
#  Robin BC at r = R (outer membrane):
#        -D_eff ∂c/∂r|_{r=R} = P_eff(t) · (c(R,t) − c_ext)
#
#  Zero-flux (symmetry) at r = 0.
#
#  TWO GEOMETRIES, one stencil.  default_parms() gives the cylindrical
#  reactor above (biofilm lumped into a uniform scalar sink over a 1 cm
#  radius).  slab_parms() gives depth across a single biofilm — z = 0
#  substratum, z = L_f bulk-liquid face — which is the only one of the two
#  that can be compared with micron-scale in-film measurements.  Geometry
#  enters solely through the face weights built by face_weights(); the
#  axis/substratum node is identical in both and is unchanged.
#
#  Spatial discretisation: finite-volume method of lines (Nr cells).
#  Time integration: deSolve::ode (LSODA adaptive stiff solver).
#
#  Key references:
#   Fox et al. (2009) SRNL/Nafion Donnan dialysis
#   Lara et al. (2023) Membranes — radiation-evolving permeability
#   Renslow et al. (2017) Shewanella uranium biofilm
#   Aydogan Gokturk et al. (2022) Nat. Commun. Donnan dialysis
# ============================================================

library(deSolve)
library(shiny)
library(plotly)

# ------------------------------------------------------------
#  Right-hand side: method-of-lines ODE for the coupled system
#
#  State vector y (length 2*Nr + 1):
#    y[1 .. Nr]        = c_i  (mobile species concentration)
#    y[Nr+1 .. 2*Nr]   = s_i  (immobile phase concentration)
#    y[2*Nr+1]         = m    (membrane integrity, 1 = intact)
#
#  parms: named list (see defaults below)
# ------------------------------------------------------------
radiodialysis_rhs <- function(t, y, parms) {
  with(parms, {
    Nr    <- length(r_grid)
    dr    <- r_grid[2] - r_grid[1]

    c_vec <- y[seq_len(Nr)]
    s_vec <- y[Nr + seq_len(Nr)]
    m_val <- y[2 * Nr + 1]

    # Cumulative absorbed dose at the membrane (Gy), linear in t
    D_cum <- Ddot_R * t

    # Outer-face transfer coefficient (cm s⁻¹).  The producer declares which
    # physical quantity this is; the stencil below must not assume one.
    #   cylindrical : radiation-evolving Nafion permeability P_eff(t)
    #                 (Lara et al. 2023, Eq. 6)
    #   slab        : external liquid-film mass-transfer coefficient k_L,
    #                 time-invariant.  This is NOT P_eff and must not be named
    #                 or reported as one — commit e4021c0 withdrew P_eff from
    #                 every producer for carrying a fabricated unit.
    # switch() with an unnamed default, NOT `if slab else cylindrical`: the
    # else form treats EVERY unrecognised value as cylindrical, so a misspelt
    # or future geom silently acquires a scientific interpretation instead of
    # refusing (AGENTS.md rule 3).  face_weights() validates, but the weights
    # reach this function precomputed in parms, so it never runs here.
    bc_coef <- switch(
      geom,
      slab        = k_L,
      cylindrical = P0 * exp(alpha_P * D_cum),
      stop("radiodialysis_rhs(): unsupported geom ", sQuote(geom),
           ". Supported: 'slab', 'cylindrical'.  An unrecognised geometry ",
           "must refuse, not fall through to one of them.", call. = FALSE)
    )

    # --------------------------------------------------------
    # (3) Membrane damage ODE
    # --------------------------------------------------------
    dm_dt <- -k_dam * Ddot_R * m_val

    # --------------------------------------------------------
    # Source/sink term (identical for mobile and immobile)
    # --------------------------------------------------------
    uptake_rate <- uptake_rate_of(X_total, f_red_active, k_ads, k_red)

    # --------------------------------------------------------
    # (1) Mobile species — diffusion + reaction (geometry via w_plus/w_minus)
    # --------------------------------------------------------
    dc_dt <- numeric(Nr)

    # i = 1: identical in both geometries, and unchanged by the slab work.
    #   cylindrical  r = 0, symmetry, L'Hôpital limit → 2 D_eff ∂²c/∂r²
    #   slab         z = 0, substratum, no-flux ghost c₀ = c₂ → 2 D_eff ∂²c/∂z²
    # The reflected ghost is c₀ = c₂, NOT c₀ = c₁; the latter drops the whole
    # scheme from second order to first.  Reads no metric weight (both are NA).
    dc_dt[1] <- D_eff * 2.0 * (c_vec[2] - c_vec[1]) / dr^2 +
                (-uptake_rate * c_vec[1] + k_des * s_vec[1])

    # i = 2 .. Nr-1: interior finite-volume stencil.
    # Geometry enters only through the face weights w_plus / w_minus, built
    # once in the parameter constructor:
    #   cylindrical  w± = (r ± dr/2) / r     slab  w± = 1
    if (Nr > 2) {
      for (i in 2:(Nr - 1)) {
        diff_op <- D_eff *
          (w_plus[i]  * (c_vec[i + 1] - c_vec[i]) -
           w_minus[i] * (c_vec[i]     - c_vec[i - 1])) / dr^2
        dc_dt[i] <- diff_op +
          (-uptake_rate * c_vec[i] + k_des * s_vec[i])
      }
    }

    # i = Nr: Robin BC via ghost-point.
    #   cylindrical  r = R, the Donnan / Nafion membrane
    #   slab         z = L_f, the biofilm / bulk-liquid interface
    #   ghost: c[Nr+1] = c[Nr-1] - 2 dr bc_coef (c[Nr] - c_ext) / D_eff
    {
      c_ghost <- c_vec[Nr - 1] -
        2.0 * dr * bc_coef * (c_vec[Nr] - c_ext) / D_eff
      diff_op <- D_eff *
        (w_plus[Nr]  * (c_ghost   - c_vec[Nr]) -
         w_minus[Nr] * (c_vec[Nr] - c_vec[Nr - 1])) / dr^2
      dc_dt[Nr] <- diff_op +
        (-uptake_rate * c_vec[Nr] + k_des * s_vec[Nr])
    }

    # --------------------------------------------------------
    # (2) Immobile phase ODE (no spatial flux)
    # --------------------------------------------------------
    ds_dt <- uptake_rate * c_vec - (k_des + k_loss) * s_vec

    list(c(dc_dt, ds_dt, dm_dt))
  })
}

# ------------------------------------------------------------
#  Default parameters
#  Values drawn from Table 2 of the paper and Deep Research
#  synthesis (April 2026) where available; otherwise reasonable
#  literature-consistent estimates for low-level nuclear waste
#  remediation context.
# ------------------------------------------------------------

# ------------------------------------------------------------
#  Uptake rate U (s⁻¹).
#
#  f_red_active is a FRACTION of X_total, not a second density:
#
#      X_red = f_red_active * X_total,  f_red_active in [0, 1]
#
#  which is the contract already written in
#  docs/research/radiotrophic_calibration_map.md:395 as
#  X_red = f_red,dry * X_total.  Two consequences, both structural:
#
#    * X_red <= X_total holds by construction rather than by an assertion
#      nobody wrote.  The previous form fixed X_red at 0.3 while X_total was
#      caller-declared, so a corrected low-biomass basis produced a reducing
#      mass EXCEEDING total biomass, and U did not scale with the declared
#      basis at all (Codex P1 on #19).
#    * U = X_total * (k_ads + k_red * f_red_active) is proportional to
#      X_total, so the whole uptake term inherits the X_total gate.  There is
#      no second gated quantity to forget.
#
#  WHAT f_red_active IS NOT.  It is an ACTIVE-reducer dry-mass fraction, not a
#  taxonomic one.  The default 0.3 is inherited, labelled "(Shewanella proxy)"
#  since 45de4ba, and is a TAXONOMIC proxy standing in an activity slot -- the
#  very substitution active_from_taxonomic() refuses without a measured
#  activity fraction ("taxonomic abundance is not functional activity").  It is
#  not 2/7 of seven species: the coupled path counts one species, S. oneidensis
#  (biofilms_potts.jl:1377).  Treat 0.3 as an unvalidated placeholder gated by
#  RADIODIALYSIS: BLOCKED, not as a composition.
# ------------------------------------------------------------
uptake_rate_of <- function(X_total, f_red_active, k_ads, k_red) {
  # X_total is validated HERE too (Codex on #19): the structural bound
  # f*X_total <= X_total assumes a finite non-negative scalar, and a hand-built
  # parms list carrying X_total = -1 reached both callers, returned a negative
  # uptake, and broke the very bound this helper is cited for.  A guard that
  # checks one of its two operands proves nothing about their product.
  finite_scalar <- function(x, nm, hi = Inf) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0 || x > hi)
      stop(nm, " must be a single finite number in [0, ",
           if (is.finite(hi)) hi else "Inf", "], got ",
           paste(format(x), collapse = ", "), call. = FALSE)
  }
  finite_scalar(X_total, "X_total")
  finite_scalar(f_red_active,
                "f_red_active (a fraction OF X_total, not a second density)", 1)
  X_total * (k_ads + k_red * f_red_active)
}

# ------------------------------------------------------------
#  Diffusion face weights.  The ONLY thing geometry changes.
#
#    cylindrical (1/r) ∂/∂r (r D ∂c/∂r)  ->  w± = (r ± dr/2) / r
#    slab (Cartesian)  ∂/∂z (D ∂c/∂z)    ->  w± = 1
#
#  Node 1 is deliberately NA: it uses the reflected-ghost form
#  (c₀ = c₂ -> factor 2), which is identical for the cylindrical
#  axis (L'Hôpital) and the slab substratum (no-flux), and reads
#  no metric weight.  NA so that a future reader who *does* index
#  it fails loudly rather than picking up a silently wrong number.
# ------------------------------------------------------------
face_weights <- function(r_grid, geom) {
  stopifnot(geom %in% c("cylindrical", "slab"))
  dr <- r_grid[2] - r_grid[1]
  if (identical(geom, "slab")) {
    w_plus <- rep(1.0, length(r_grid)); w_minus <- rep(1.0, length(r_grid))
  } else {
    w_plus  <- (r_grid + 0.5 * dr) / r_grid
    w_minus <- (r_grid - 0.5 * dr) / r_grid
  }
  w_plus[1] <- NA_real_; w_minus[1] <- NA_real_
  list(w_plus = w_plus, w_minus = w_minus)
}

default_parms <- function(Nr = 40, R = 1.0) {
  r_grid <- seq(0, R, length.out = Nr)
  w      <- face_weights(r_grid, "cylindrical")
  list(
    # --- Spatial grid ---
    r_grid  = r_grid,
    R       = R,
    geom    = "cylindrical",
    w_plus  = w$w_plus,
    w_minus = w$w_minus,

    # --- Transport ---
    D_eff   = 1e-3,    # effective diffusivity (cm² s⁻¹), Table 2 range 1e-4..1e-2

    # --- Biosorption / bioreduction (Renslow et al. 2017) ---
    k_ads   = 0.05,    # adsorption rate constant (cm³ g⁻¹ s⁻¹)
    k_red   = 0.02,    # bioreduction rate constant (cm³ g⁻¹ s⁻¹)
    k_des   = 0.005,   # desorption rate constant (s⁻¹)
    k_loss  = 0.001,   # immobile-phase loss / precipitation (s⁻¹)

    # --- Biomass ---
    X_total = 1.0,     # total biofilm dry mass density (g cm⁻³)
    f_red_active = 0.3, # ACTIVE-reducer fraction OF X_total, in [0,1].
                       # Unvalidated placeholder; see uptake_rate_of() above.

    # --- Membrane (Nafion / Donnan) ---
    P0      = 0.01,    # baseline permeability (cm s⁻¹), Fox et al. 2009
    alpha_P = 0.02,    # radiation-damage permeability coefficient (Gy⁻¹)
                       # Lara et al. 2023: permeability ×2-4 under 10-50 Gy
    k_dam   = 0.005,   # membrane structural damage rate (Gy⁻¹)

    # --- Radiation at membrane r = R ---
    Ddot_R  = 1.0,     # dose rate at membrane (Gy s⁻¹)
                       # scales with gamma_intensity / exp(-kappa * R)

    # --- Boundary ---
    c_ext   = 1.0      # external contaminant concentration (normalised)
  )
}

# ------------------------------------------------------------
#  Slab preset: depth across a biofilm, NOT radius across a reactor.
#
#    z = 0    substratum          (no-flux)
#    z = L_f  bulk-liquid face    (Robin, external mass transfer k_L)
#
#  Same three equations, same stencil, Cartesian face weights.  Written to
#  be comparable with micron-scale in-film measurements (e.g. SECM metal
#  profiles), which the cylindrical preset cannot produce: there the biofilm
#  is a spatially uniform scalar sink over a 1 cm reactor radius.
#
#  Three semantics are declared here rather than inherited (AGENTS.md rule 4):
#    * the Robin coefficient is k_L, a liquid-film mass-transfer coefficient,
#      NOT the membrane permeability P_eff;
#    * the radiation terms are switched off — there is no membrane at a
#      biofilm/liquid interface, so k_dam = Ddot_R = 0 and m(t) stays 1;
#    * D_eff is a MOLECULAR diffusivity.  The cylindrical preset's 1e-3 cm² s⁻¹
#      is ~100x the self-diffusivity of water: it is a reactor-scale dispersion
#      coefficient and is meaningless on a 100 µm domain.
#
#  X_total has NO default on purpose.  See the stop() below.
# ------------------------------------------------------------
slab_parms <- function(Nr = 40, L_f = 0.01, X_total,
                       D_eff = 1e-5, k_L = 1e-3, f_red_active = 0.3) {
  if (missing(X_total)) {
    stop(
      "slab_parms(): X_total must be stated explicitly, not defaulted.\n",
      "  The biofilm dry-mass basis is exactly what RADIODIALYSIS: BLOCKED\n",
      "  names (README.md:350-353).  It sets the uptake rate\n",
      "    U = X_total * (k_ads + k_red*f_red_active),\n",
      "  and the penetration depth lambda = sqrt(D_eff/k_eff) scales as\n",
      "  1/sqrt(U), so whether ANY gradient is resolvable across L_f is\n",
      "  downstream of that gate.  Pass a value and label it a test value.",
      call. = FALSE
    )
  }
  z_grid <- seq(0, L_f, length.out = Nr)
  w      <- face_weights(z_grid, "slab")
  list(
    # --- Spatial grid: depth, not radius ---
    r_grid  = z_grid,   # name kept so run_radiodialysis() is reused unchanged
    R       = L_f,      # film thickness (cm)
    geom    = "slab",
    w_plus  = w$w_plus,
    w_minus = w$w_minus,

    # --- Transport ---
    D_eff   = D_eff,   # molecular diffusivity in the film (cm² s⁻¹)

    # --- Biosorption / bioreduction (unchanged from the cylindrical preset) ---
    k_ads   = 0.05,
    k_red   = 0.02,
    k_des   = 0.005,
    k_loss  = 0.001,   # the ONLY true sink at steady state; see note below

    # --- Biomass ---
    X_total = X_total, # caller-declared; gated, see stop() above
    f_red_active = f_red_active,  # fraction OF X_total; U scales with the
                       # declared basis, so it inherits that same gate.

    # --- Outer face: liquid-film mass transfer, NOT membrane permeability ---
    k_L     = k_L,     # (cm s⁻¹).  P0 / alpha_P are deliberately absent.

    # --- Radiation: off.  No membrane exists at a biofilm/liquid interface ---
    k_dam   = 0.0,
    Ddot_R  = 0.0,

    # --- Boundary ---
    c_ext   = 1.0      # bulk-liquid concentration (normalised)
  )
}

# ------------------------------------------------------------
#  Steady-state penetration depth.
#
#  Setting ds/dt = 0 gives s = U c / (k_des + k_loss), and substituting into
#  the mobile equation shows the sorption/desorption pair is a REVERSIBLE
#  BUFFER that cancels.  The only true steady sink is precipitation:
#
#      k_eff = U * k_loss / (k_des + k_loss)
#
#  With the file's constants k_eff/U = 1/6 exactly, so lambda_steady is
#  sqrt(6) ~ 2.45x the transient value.  This ratio is structural and does
#  not depend on X_total; the absolute lambda does.
# ------------------------------------------------------------
penetration_depth <- function(parms, phase = c("steady", "transient")) {
  phase <- match.arg(phase)
  with(parms, {
    U <- uptake_rate_of(X_total, f_red_active, k_ads, k_red)
    k <- if (phase == "steady") U * k_loss / (k_des + k_loss) else U
    sqrt(D_eff / k)
  })
}

# ------------------------------------------------------------
#  Run solver: returns list with times, c-matrix, s-matrix, m-vec
# ------------------------------------------------------------
run_radiodialysis <- function(parms, t_end = 100, n_out = 200,
                               c0 = NULL, s0 = NULL) {
  Nr <- length(parms$r_grid)

  # Initial conditions: contaminant absent inside, full membrane integrity
  if (is.null(c0)) c0 <- rep(0.0, Nr)
  if (is.null(s0)) s0 <- rep(0.0, Nr)
  m0 <- 1.0

  y0    <- c(c0, s0, m0)
  times <- seq(0, t_end, length.out = n_out)

  out <- ode(y = y0, times = times, func = radiodialysis_rhs,
             parms = parms, method = "lsoda",
             rtol = 1e-6, atol = 1e-8)

  list(
    times  = out[, "time"],
    c_mat  = out[, 2:(Nr + 1)],             # Nr columns
    s_mat  = out[, (Nr + 2):(2 * Nr + 1)],  # Nr columns
    m_vec  = out[, 2 * Nr + 2],
    parms  = parms
  )
}

# ------------------------------------------------------------
#  Shiny UI
# ------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Radiodialysis Membrane Transport — Cylindrical Bioreactor"),

  sidebarLayout(
    sidebarPanel(width = 3,
      h4("Transport"),
      sliderInput("D_eff", "D_eff (cm² s⁻¹, ×10⁻³)",
                  min = 0.1, max = 20, value = 1.0, step = 0.1),

      h4("Biosorption / Bioreduction"),
      sliderInput("k_ads",   "k_ads (cm³ g⁻¹ s⁻¹)",
                  min = 0, max = 0.2, value = 0.05, step = 0.005),
      sliderInput("k_red",   "k_red (cm³ g⁻¹ s⁻¹)",
                  min = 0, max = 0.1, value = 0.02, step = 0.002),
      sliderInput("k_des",   "k_des (s⁻¹)",
                  min = 0, max = 0.05, value = 0.005, step = 0.001),
      sliderInput("f_red_active", "Active metal-reducing fraction of X_total",
                  min = 0, max = 1, value = 0.3, step = 0.05),

      h4("Membrane"),
      sliderInput("P0",      "P₀ (cm s⁻¹, baseline)",
                  min = 0.001, max = 0.05, value = 0.01, step = 0.001),
      sliderInput("alpha_P", "α (Gy⁻¹, permeability gain)",
                  min = 0, max = 0.1, value = 0.02, step = 0.002),
      sliderInput("k_dam",   "k_dam (Gy⁻¹, structural damage)",
                  min = 0, max = 0.05, value = 0.005, step = 0.001),

      h4("Radiation"),
      sliderInput("Ddot_R",  "Ḋ(R) at membrane (Gy s⁻¹)",
                  min = 0.1, max = 5, value = 1.0, step = 0.1),
      sliderInput("c_ext",   "c_ext (external conc., normalised)",
                  min = 0.1, max = 2, value = 1.0, step = 0.1),

      h4("Simulation"),
      sliderInput("t_end", "Duration (s)", min = 10, max = 500, value = 100, step = 10),
      sliderInput("Nr",    "Radial grid points", min = 10, max = 80, value = 40, step = 5),

      br(),
      actionButton("run", "Run Simulation", class = "btn-primary"),
      br(), br(),
      verbatimTextOutput("status")
    ),

    mainPanel(width = 9,
      tabsetPanel(
        tabPanel("Mobile concentration c(r,t)",
                 plotlyOutput("plot_c", height = "450px")),
        tabPanel("Immobile phase s(r,t)",
                 plotlyOutput("plot_s", height = "450px")),
        tabPanel("Membrane / P_eff(t)",
                 plotlyOutput("plot_m", height = "450px")),
        tabPanel("Radial profiles",
                 plotlyOutput("plot_profiles", height = "450px"))
      )
    )
  )
)

# ------------------------------------------------------------
#  Shiny server
# ------------------------------------------------------------
server <- function(input, output, session) {

  result <- reactiveVal(NULL)

  observeEvent(input$run, {
    parms <- default_parms(Nr = input$Nr, R = 1.0)
    parms$D_eff   <- input$D_eff  * 1e-3
    parms$k_ads   <- input$k_ads
    parms$k_red   <- input$k_red
    parms$k_des   <- input$k_des
    parms$f_red_active <- input$f_red_active
    parms$P0      <- input$P0
    parms$alpha_P <- input$alpha_P
    parms$k_dam   <- input$k_dam
    parms$Ddot_R  <- input$Ddot_R
    parms$c_ext   <- input$c_ext

    res <- run_radiodialysis(parms, t_end = input$t_end, n_out = 150)
    result(res)

    output$status <- renderText({
      m_final  <- tail(res$m_vec, 1)
      P_final  <- parms$P0 * exp(parms$alpha_P * parms$Ddot_R * input$t_end)
      c_centre <- res$c_mat[nrow(res$c_mat), 1]
      c_wall   <- res$c_mat[nrow(res$c_mat), input$Nr]
      sprintf(
        # P_eff/P0, dimensionless. The absolute value was printed with
        # `cm/s` attached, which is a fabricated unit on an uncalibrated
        # placeholder: P0 is a literature-anchored default, not a measurement
        # of this system, so P0 * exp(...) is not a permeability in cm/s.
        "t_end=%.0f s  |  m(t_end)=%.4f  |  P_eff(t_end)/P0=%.4f (dimensionless)\nc(0,t_end)=%.4f  |  c(R,t_end)=%.4f  |  Nr=%d",
        input$t_end, m_final, P_final / parms$P0, c_centre, c_wall, input$Nr
      )
    })
  })

  # ---- Heatmap: c(r,t) ----
  output$plot_c <- renderPlotly({
    res <- result()
    req(res)
    p <- res$parms
    plot_ly(
      x = p$r_grid,
      y = res$times,
      z = res$c_mat,
      type   = "heatmap",
      colorscale = "Viridis",
      colorbar = list(title = "c (norm.)")
    ) %>% layout(
      title  = "Mobile contaminant concentration c(r,t)",
      xaxis  = list(title = "r (cm)"),
      yaxis  = list(title = "t (s)")
    )
  })

  # ---- Heatmap: s(r,t) ----
  output$plot_s <- renderPlotly({
    res <- result()
    req(res)
    p <- res$parms
    plot_ly(
      x = p$r_grid,
      y = res$times,
      z = res$s_mat,
      type   = "heatmap",
      colorscale = "Hot",
      colorbar = list(title = "s (norm.)")
    ) %>% layout(
      title  = "Immobile phase s(r,t) — sorbed + reduced contaminant",
      xaxis  = list(title = "r (cm)"),
      yaxis  = list(title = "t (s)")
    )
  })

  # ---- Time series: m(t) and P_eff(t) ----
  output$plot_m <- renderPlotly({
    res <- result()
    req(res)
    p       <- res$parms
    P_eff_t <- p$P0 * exp(p$alpha_P * p$Ddot_R * res$times)
    D_cum_t <- p$Ddot_R * res$times

    plot_ly() %>%
      add_trace(x = res$times, y = res$m_vec,
                name = "m(t) membrane integrity",
                type = "scatter", mode = "lines",
                line = list(color = "steelblue", width = 2.5)) %>%
      add_trace(x = res$times, y = P_eff_t / p$P0,
                name = "P_eff(t) / P₀ (normalised)",
                type = "scatter", mode = "lines",
                line = list(color = "tomato", width = 2.5, dash = "dash"),
                yaxis = "y2") %>%
      add_trace(x = res$times, y = D_cum_t,
                name = "D_cum(t) (Gy)",
                type = "scatter", mode = "lines",
                line = list(color = "goldenrod", width = 1.5, dash = "dot"),
                yaxis = "y3") %>%
      layout(
        title  = "Membrane damage m(t), permeability P_eff(t), and cumulative dose D_cum(t)",
        xaxis  = list(title = "t (s)"),
        yaxis  = list(title = "m(t) (integrity)", range = c(0, 1.05)),
        yaxis2 = list(title = "P_eff / P₀", overlaying = "y", side = "right",
                      showgrid = FALSE),
        yaxis3 = list(title = "D_cum (Gy)", overlaying = "y", side = "right",
                      anchor = "free", position = 1.0, showgrid = FALSE),
        legend = list(x = 0.05, y = 0.95)
      )
  })

  # ---- Radial profiles at selected times ----
  output$plot_profiles <- renderPlotly({
    res <- result()
    req(res)
    p    <- res$parms
    nout <- length(res$times)
    t_idx <- round(seq(1, nout, length.out = 5))

    pal <- c("navy", "steelblue", "seagreen", "orange", "firebrick")

    fig <- plot_ly()
    for (k in seq_along(t_idx)) {
      idx <- t_idx[k]
      fig <- fig %>%
        add_trace(
          x    = p$r_grid,
          y    = res$c_mat[idx, ],
          name = sprintf("c  t=%.1f s", res$times[idx]),
          type = "scatter", mode = "lines",
          line = list(color = pal[k], width = 2.5)
        ) %>%
        add_trace(
          x    = p$r_grid,
          y    = res$s_mat[idx, ],
          name = sprintf("s  t=%.1f s", res$times[idx]),
          type = "scatter", mode = "lines",
          line = list(color = pal[k], width = 1.5, dash = "dot")
        )
    }
    fig %>% layout(
      title  = "Radial profiles: solid=c(r,t) (mobile), dotted=s(r,t) (immobile)",
      xaxis  = list(title = "r (cm)"),
      yaxis  = list(title = "Concentration (normalised)"),
      legend = list(x = 0.01, y = 0.99)
    )
  })
}

app <- shinyApp(ui = ui, server = server)

# Run interactively (VS Code R extension / RStudio) or as Rscript
if (interactive()) {
  app
} else {
  port <- as.integer(Sys.getenv("SHINY_PORT", "7799"))
  message(sprintf("Radiodialysis app → http://127.0.0.1:%d", port))
  shiny::runApp(app, host = "127.0.0.1", port = port, launch.browser = FALSE)
}
