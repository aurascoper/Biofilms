# ============================================================
#  Radiodialysis Membrane Transport — Cylindrical Bioreactor
#  "Modeling Radiotrophic Fitness" (Kinder & Faulkner, 2026)
#
#  Implements the minimal 3-equation PDE system derived from
#  the Deep Research synthesis (April 2026):
#
#  (1) Mobile species: cylindrical reaction-diffusion PDE
#        ∂c/∂t = (1/r) ∂/∂r (r D_eff ∂c/∂r)
#                - (k_ads X + k_red X_red) c + k_des s
#
#  (2) Immobile phase ODE (biosorption / bioreduction):
#        ∂s/∂t = (k_ads X + k_red X_red) c - (k_des + k_loss) s
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

    # Radiation-driven effective permeability (Lara et al. 2023, Eq. 6)
    P_eff <- P0 * exp(alpha_P * D_cum)

    # --------------------------------------------------------
    # (3) Membrane damage ODE
    # --------------------------------------------------------
    dm_dt <- -k_dam * Ddot_R * m_val

    # --------------------------------------------------------
    # Source/sink term (identical for mobile and immobile)
    # --------------------------------------------------------
    uptake_rate <- k_ads * X_total + k_red * X_red  # s⁻¹

    # --------------------------------------------------------
    # (1) Mobile species — cylindrical diffusion + reaction
    # --------------------------------------------------------
    dc_dt <- numeric(Nr)

    # i = 1 (r = 0): L'Hôpital limit → ∂c/∂t = 2 D_eff ∂²c/∂r²
    dc_dt[1] <- D_eff * 2.0 * (c_vec[2] - c_vec[1]) / dr^2 +
                (-uptake_rate * c_vec[1] + k_des * s_vec[1])

    # i = 2 .. Nr-1: interior finite-volume stencil
    if (Nr > 2) {
      for (i in 2:(Nr - 1)) {
        r_i      <- r_grid[i]
        r_plus   <- r_i + 0.5 * dr
        r_minus  <- r_i - 0.5 * dr
        diff_cyl <- D_eff *
          (r_plus  * (c_vec[i + 1] - c_vec[i]) -
           r_minus * (c_vec[i]     - c_vec[i - 1])) /
          (r_i * dr^2)
        dc_dt[i] <- diff_cyl +
          (-uptake_rate * c_vec[i] + k_des * s_vec[i])
      }
    }

    # i = Nr (r = R): Robin BC via ghost-point (Donnan / Nafion membrane)
    #   ghost: c[Nr+1] = c[Nr-1] - 2 dr P_eff (c[Nr] - c_ext) / D_eff
    {
      i       <- Nr
      r_i     <- r_grid[i]
      r_plus  <- r_i + 0.5 * dr
      r_minus <- r_i - 0.5 * dr
      c_ghost <- c_vec[Nr - 1] -
        2.0 * dr * P_eff * (c_vec[Nr] - c_ext) / D_eff
      diff_cyl <- D_eff *
        (r_plus  * (c_ghost   - c_vec[Nr]) -
         r_minus * (c_vec[Nr] - c_vec[Nr - 1])) /
        (r_i * dr^2)
      dc_dt[Nr] <- diff_cyl +
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
default_parms <- function(Nr = 40, R = 1.0) {
  r_grid <- seq(0, R, length.out = Nr)
  list(
    # --- Spatial grid ---
    r_grid  = r_grid,
    R       = R,

    # --- Transport ---
    D_eff   = 1e-3,    # effective diffusivity (cm² s⁻¹), Table 2 range 1e-4..1e-2

    # --- Biosorption / bioreduction (Renslow et al. 2017) ---
    k_ads   = 0.05,    # adsorption rate constant (cm³ g⁻¹ s⁻¹)
    k_red   = 0.02,    # bioreduction rate constant (cm³ g⁻¹ s⁻¹)
    k_des   = 0.005,   # desorption rate constant (s⁻¹)
    k_loss  = 0.001,   # immobile-phase loss / precipitation (s⁻¹)

    # --- Biomass ---
    X_total = 1.0,     # total biofilm dry mass density (g cm⁻³)
    X_red   = 0.3,     # metal-reducing fraction (Shewanella proxy)

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
      sliderInput("X_red",   "Metal-reducing biomass fraction",
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
    parms$X_red   <- input$X_red
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
        "t_end=%.0f s  |  m(t_end)=%.4f  |  P_eff(t_end)=%.4f cm/s\nc(0,t_end)=%.4f  |  c(R,t_end)=%.4f  |  Nr=%d",
        input$t_end, m_final, P_final, c_centre, c_wall, input$Nr
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
