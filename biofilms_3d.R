# Load necessary libraries
library(deSolve)
library(stats)
library(plotly)
library(shiny)

# ---------------------------------------------------------------------------
# Species parameters
# 7-species list from the paper (Table 1).
# radiotrophic = TRUE  -> species moves TOWARD the radiation axis (positive
#                         gamma sensitivity, melanin-mediated radiotropism)
# radiotrophic = FALSE -> species moves AWAY from the axis (radiosensitive)
# ---------------------------------------------------------------------------
species_params <- list(
  list(
    name            = "Pseudoalteromonas",
    mu              = 0.09,
    D               = 0.005,
    rad_sensitivity = 0.06,
    quad_coeff      = 0.02,
    radiotrophic    = FALSE
  ),
  list(
    name            = "Shewanella oneidensis",
    mu              = 0.11,
    D               = 0.010,
    rad_sensitivity = 0.04,
    quad_coeff      = 0.02,
    radiotrophic    = FALSE
  ),
  list(
    name            = "Deinococcus radiodurans",
    mu              = 0.08,
    D               = 0.008,
    rad_sensitivity = 0.05,
    quad_coeff      = 0.015,
    radiotrophic    = FALSE
  ),
  list(
    name            = "Bacillus subtilis",
    mu              = 0.10,
    D               = 0.009,
    rad_sensitivity = 0.03,
    quad_coeff      = 0.018,
    radiotrophic    = FALSE
  ),
  list(
    name            = "Aspergillus niger",
    mu              = 0.07,
    D               = 0.007,
    rad_sensitivity = 0.025,
    quad_coeff      = 0.012,
    radiotrophic    = FALSE
  ),
  list(
    name            = "Cryptococcus neoformans",
    mu              = 0.09,
    D               = 0.008,
    rad_sensitivity = 0.055,
    quad_coeff      = 0.020,
    radiotrophic    = TRUE   # melanin-mediated radiotropism
  ),
  list(
    name            = "Cladosporium sphaerospermum",
    mu              = 0.085,
    D               = 0.007,
    rad_sensitivity = 0.050,
    quad_coeff      = 0.018,
    radiotrophic    = TRUE   # melanin-mediated radiotropism
  )
)

# ---------------------------------------------------------------------------
# Helper: cylindrical elastic wall reflection
#
# Given a proposed new position (pos_new) and velocity (vel), enforce:
#   radial:  r <= R   — reflect radial velocity component
#   axial:   0 <= z <= L — reflect z velocity
#
# Returns a list(pos = c(x,y,z), vel = c(vx,vy,vz))
# ---------------------------------------------------------------------------
reflect_cylinder <- function(pos_new, vel, R, L) {
  x  <- pos_new[1]
  y  <- pos_new[2]
  z  <- pos_new[3]
  vx <- vel[1]
  vy <- vel[2]
  vz <- vel[3]

  # --- radial reflection ---
  r_new <- sqrt(x^2 + y^2)
  if (r_new > R && r_new > 0) {
    # unit radial vector
    nx <- x / r_new
    ny <- y / r_new
    # reflect velocity: v_r <- -(v . n_hat) * n_hat  (flip radial component)
    v_dot_n <- vx * nx + vy * ny
    vx <- vx - 2 * v_dot_n * nx
    vy <- vy - 2 * v_dot_n * ny
    # project position back onto the cylinder wall
    x <- R * nx
    y <- R * ny
  }

  # --- axial reflection ---
  if (z < 0) {
    z  <- -z
    vz <- -vz
  }
  if (z > L) {
    z  <- 2 * L - z
    vz <- -vz
  }

  list(pos = c(x, y, z), vel = c(vx, vy, vz))
}

# ---------------------------------------------------------------------------
# Radiation intensity along the central axis (Beer–Lambert, radial decay)
#   I_gamma(r) = gamma_intensity * exp(-kappa * r)
# ---------------------------------------------------------------------------
I_gamma <- function(r, gamma_intensity, kappa) {
  gamma_intensity * exp(-kappa * r)
}

# ---------------------------------------------------------------------------
# Radial nutrient concentration — nutrient enters from the outer membrane,
# depleted near the source axis:
#   nutrient_conc(r) = C0 * (r / R)
# ---------------------------------------------------------------------------
nutrient_conc <- function(r, C0, R) {
  C0 * (r / R)
}

# ---------------------------------------------------------------------------
# Langevin dynamics (cylindrical geometry)
#
# Forces acting on species i at position state = c(x, y, z):
#   F_rad   — radial radiation avoidance / attraction force
#   F_nutr  — radial nutrient chemotaxis (toward outer wall)
#   F_heat  — isotropic diffusion boost from heat
#   F_quad  — pairwise quadratic interaction with other species
#   F_tho   — thorium decay reduces effective nutrient contribution
#   noise   — Brownian (Wiener) noise scaled by sqrt(2D)
# ---------------------------------------------------------------------------
langevin_dynamics_cyl <- function(state, parameters, positions,
                                   gamma_intensity, kappa,
                                   thorium_intensity,
                                   C0, R, L,
                                   heat_intensity) {
  x <- state[1]; y <- state[2]; z <- state[3]
  r <- sqrt(x^2 + y^2)

  # Brownian noise
  noise_x <- rnorm(1, mean = 0, sd = sqrt(2 * parameters$D))
  noise_y <- rnorm(1, mean = 0, sd = sqrt(2 * parameters$D))
  noise_z <- rnorm(1, mean = 0, sd = sqrt(2 * parameters$D))

  # --- Radiation force (radial vector) ---
  # Radiotrophic: move TOWARD axis (sign = +1 inward = -radial direction
  #               expressed as negative of outward unit vector)
  # Radiosensitive: move AWAY from axis (sign = -1, i.e. outward)
  I_r   <- I_gamma(r, gamma_intensity, kappa)
  F_rad <- c(0, 0, 0)
  if (r > 1e-10) {
    radial_unit <- c(x / r, y / r, 0)
    if (parameters$radiotrophic) {
      # radiotrophic: attracted toward source (axis), so force points inward
      F_rad <- parameters$rad_sensitivity * I_r * (-radial_unit)
    } else {
      # radiosensitive: repelled from source, force points outward
      F_rad <- -parameters$rad_sensitivity * I_r * radial_unit
    }
  }

  # --- Nutrient chemotaxis (toward outer wall = radially outward) ---
  nutr_grad_magnitude <- C0 / R   # d/dr [C0 * r/R] = C0/R (constant gradient)
  F_nutr <- c(0, 0, 0)
  if (r > 1e-10) {
    radial_unit <- c(x / r, y / r, 0)
    F_nutr <- nutr_grad_magnitude * radial_unit
  }

  # --- Heat: isotropic diffusion boost ---
  F_heat_scalar <- heat_intensity * parameters$D

  # --- Thorium: suppresses nutrient uptake (scalar penalty) ---
  nutrient_effect <- -thorium_intensity * nutrient_conc(r, C0, R)

  # --- Pairwise quadratic interaction ---
  F_quad <- c(0, 0, 0)
  num_sp <- nrow(positions)
  for (j in seq_len(num_sp)) {
    if (j != parameters$idx) {
      diff_vec <- positions[j, ] - state
      dist_ij  <- sqrt(sum(diff_vec^2))
      if (dist_ij > 1e-10) {
        interaction <- parameters$quad_coeff / dist_ij
        F_quad      <- F_quad + interaction * diff_vec
      }
    }
  }

  # --- Assemble total force and compute velocity increments ---
  # F_heat and nutrient_effect are scalars added uniformly to all components
  F_total_x <- F_rad[1] + F_nutr[1] + F_quad[1] + F_heat_scalar + nutrient_effect
  F_total_y <- F_rad[2] + F_nutr[2] + F_quad[2] + F_heat_scalar + nutrient_effect
  F_total_z <- F_rad[3] + F_nutr[3] + F_quad[3] + F_heat_scalar + nutrient_effect

  dx_dt <- parameters$mu * F_total_x + noise_x
  dy_dt <- parameters$mu * F_total_y + noise_y
  dz_dt <- parameters$mu * F_total_z + noise_z

  return(c(dx_dt, dy_dt, dz_dt))
}

# ---------------------------------------------------------------------------
# Simulation: integrate Langevin steps inside the cylindrical bioreactor
# ---------------------------------------------------------------------------
simulate_species_motility <- function(num_species, num_steps, dt, species_params,
                                       gamma_intensity, kappa,
                                       thorium_intensity,
                                       C0, R, L,
                                       heat_intensity) {

  # Initial positions: uniform in cylinder using rejection sampling
  positions <- matrix(NA_real_, nrow = num_species, ncol = 3)
  for (i in seq_len(num_species)) {
    repeat {
      x_try <- runif(1, -R, R)
      y_try <- runif(1, -R, R)
      if (sqrt(x_try^2 + y_try^2) <= R) {
        positions[i, ] <- c(x_try, y_try, runif(1, 0, L))
        break
      }
    }
  }

  # Velocity state (initialised to zero; updated each step for reflection)
  velocities <- matrix(0, nrow = num_species, ncol = 3)

  trajectory_data <- array(NA_real_, dim = c(num_steps, num_species, 3))

  for (step in seq_len(num_steps)) {
    for (i in seq_len(num_species)) {
      state  <- positions[i, ]
      params <- species_params[[i]]
      params$idx <- i

      # Compute Langevin velocity increment
      delta_v <- langevin_dynamics_cyl(
        state          = state,
        parameters     = params,
        positions      = positions,
        gamma_intensity = gamma_intensity,
        kappa          = kappa,
        thorium_intensity = thorium_intensity,
        C0             = C0,
        R              = R,
        L              = L,
        heat_intensity = heat_intensity
      )

      # Proposed new position
      pos_new <- state + delta_v * dt

      # Apply cylindrical reflection
      reflected       <- reflect_cylinder(pos_new, delta_v, R, L)
      positions[i, ]  <- reflected$pos
      velocities[i, ] <- reflected$vel

      trajectory_data[step, i, ] <- positions[i, ]
    }
  }

  return(trajectory_data)
}

# ---------------------------------------------------------------------------
# Helper: build cylinder wireframe traces for plotly
# Returns a list of add_trace calls (top circle, bottom circle, vertical lines)
# ---------------------------------------------------------------------------
cylinder_wireframe_traces <- function(R, L, n_circle = 60, n_vert = 12) {
  theta <- seq(0, 2 * pi, length.out = n_circle)

  # Bottom circle (z = 0)
  bot <- list(
    x    = R * cos(theta),
    y    = R * sin(theta),
    z    = rep(0, n_circle),
    type = "scatter3d",
    mode = "lines",
    line = list(color = "rgba(180,180,255,0.35)", width = 2),
    showlegend = FALSE,
    name = "Cylinder"
  )

  # Top circle (z = L)
  top <- list(
    x    = R * cos(theta),
    y    = R * sin(theta),
    z    = rep(L, n_circle),
    type = "scatter3d",
    mode = "lines",
    line = list(color = "rgba(180,180,255,0.35)", width = 2),
    showlegend = FALSE,
    name = "Cylinder"
  )

  # Vertical lines at n_vert equally-spaced angles
  phi_v <- seq(0, 2 * pi, length.out = n_vert + 1)[seq_len(n_vert)]
  vert_traces <- lapply(phi_v, function(ph) {
    list(
      x    = c(R * cos(ph), R * cos(ph), NA),
      y    = c(R * sin(ph), R * sin(ph), NA),
      z    = c(0, L, NA),
      type = "scatter3d",
      mode = "lines",
      line = list(color = "rgba(180,180,255,0.25)", width = 1),
      showlegend = FALSE,
      name = "Cylinder"
    )
  })

  c(list(bot, top), vert_traces)
}

# ---------------------------------------------------------------------------
# Shiny UI
# ---------------------------------------------------------------------------
SPECIES_COLORS <- c("red", "dodgerblue", "green", "darkorchid",
                    "orange", "deeppink", "saddlebrown")

ui <- fluidPage(
  titlePanel("3D Radiotrophic Biofilm — Cylindrical Bioreactor"),

  sidebarLayout(
    sidebarPanel(
      h4("Radiation"),
      sliderInput("gamma_intensity",
                  "Gamma Radiation Intensity (I\u2080)",
                  min = 0, max = 1, value = 0.05, step = 0.05),
      sliderInput("kappa",
                  "Attenuation Coefficient (\u03ba)",
                  min = 0.5, max = 5.0, value = 2.0, step = 0.1),
      sliderInput("thorium_intensity",
                  "Thorium Intensity",
                  min = 0, max = 1, value = 0.05, step = 0.05),

      h4("Environment"),
      sliderInput("heat_intensity",
                  "Heat Intensity",
                  min = 0, max = 1, value = 0.1, step = 0.05),
      sliderInput("C0",
                  "Outer-wall Nutrient Concentration (C\u2080)",
                  min = 0, max = 1, value = 0.5, step = 0.05),

      br(),
      actionButton("run",   "Run 45-Step Simulation",
                   class = "btn-primary"),
      actionButton("reset", "Reset Simulation")
    ),

    mainPanel(
      plotlyOutput("plot", height = "600px"),
      verbatimTextOutput("status")
    )
  )
)

# ---------------------------------------------------------------------------
# Shiny server
# ---------------------------------------------------------------------------
server <- function(input, output, session) {

  # Fixed bioreactor geometry (normalised units)
  R_bio <- 1.0
  L_bio <- 2.0

  num_species <- length(species_params)
  num_steps   <- 45
  dt          <- 0.1

  # Reactive: run simulation on button press
  observeEvent(input$run, {
    gamma_intensity   <- input$gamma_intensity
    kappa_val         <- input$kappa
    thorium_intensity <- input$thorium_intensity
    heat_intensity    <- input$heat_intensity
    C0_val            <- input$C0

    traj <- simulate_species_motility(
      num_species       = num_species,
      num_steps         = num_steps,
      dt                = dt,
      species_params    = species_params,
      gamma_intensity   = gamma_intensity,
      kappa             = kappa_val,
      thorium_intensity = thorium_intensity,
      C0                = C0_val,
      R                 = R_bio,
      L                 = L_bio,
      heat_intensity    = heat_intensity
    )

    # Build plotly figure
    p <- plot_ly()

    # Add cylinder wireframe
    wf_traces <- cylinder_wireframe_traces(R_bio, L_bio)
    for (tr in wf_traces) {
      p <- p %>% add_trace(
        x          = tr$x,
        y          = tr$y,
        z          = tr$z,
        type       = tr$type,
        mode       = tr$mode,
        line       = tr$line,
        showlegend = tr$showlegend,
        name       = tr$name
      )
    }

    # Add central axis (radiation source)
    p <- p %>% add_trace(
      x          = c(0, 0),
      y          = c(0, 0),
      z          = c(0, L_bio),
      type       = "scatter3d",
      mode       = "lines",
      line       = list(color = "rgba(255,220,0,0.7)", width = 4,
                        dash = "dash"),
      name       = "Radiation axis",
      showlegend = TRUE
    )

    # Add species trajectories
    for (i in seq_len(num_species)) {
      label_suffix <- if (species_params[[i]]$radiotrophic) " (radiotrophic)" else ""
      p <- p %>% add_trace(
        x    = traj[, i, 1],
        y    = traj[, i, 2],
        z    = traj[, i, 3],
        type = "scatter3d",
        mode = "lines+markers",
        line = list(width = 3, color = SPECIES_COLORS[i]),
        marker = list(
          size   = 3,
          color  = SPECIES_COLORS[i],
          symbol = if (species_params[[i]]$radiotrophic) "diamond" else "circle"
        ),
        name = paste0(species_params[[i]]$name, label_suffix)
      )
    }

    plot_title <- sprintf(
      "Cylindrical Bioreactor  |  I\u2080=%.2f  \u03ba=%.1f  C\u2080=%.2f  Thorium=%.2f  Heat=%.2f",
      gamma_intensity, kappa_val, C0_val, thorium_intensity, heat_intensity
    )

    p <- p %>% layout(
      title = plot_title,
      scene = list(
        xaxis = list(title = "x", range = c(-R_bio * 1.1, R_bio * 1.1)),
        yaxis = list(title = "y", range = c(-R_bio * 1.1, R_bio * 1.1)),
        zaxis = list(title = "z (axial)", range = c(-0.1, L_bio + 0.1)),
        aspectmode = "manual",
        aspectratio = list(x = 1, y = 1, z = L_bio / R_bio)
      ),
      legend = list(x = 0.01, y = 0.99)
    )

    output$plot   <- renderPlotly({ p })
    output$status <- renderText({
      sprintf("Simulation complete: %d species, %d steps, dt=%.2f\nGeometry: cylinder R=%.1f L=%.1f",
              num_species, num_steps, dt, R_bio, L_bio)
    })
  })

  # Reset: clear the plot
  observeEvent(input$reset, {
    output$plot   <- renderPlotly({ plot_ly() })
    output$status <- renderText({ "Simulation reset." })
  })
}

# ---------------------------------------------------------------------------
# Launch
# ---------------------------------------------------------------------------
shinyApp(ui = ui, server = server)
