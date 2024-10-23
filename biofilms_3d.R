# Load necessary libraries
library(deSolve)
library(stats)
library(plotly)
library(shiny)

# Define species-specific motility, radiation parameters, and quadratic interaction coefficients
species_params <- list(
  list(name = "Pseudoalteromonas", mu = 0.09, D = 0.005, rad_sensitivity = 0.06, quad_coeff = 0.02),
  list(name = "Shewanella Oneidensis", mu = 0.11, D = 0.01, rad_sensitivity = 0.04, quad_coeff = 0.02),
  list(name = "Polaribacter", mu = 0.1, D = 0.01, rad_sensitivity = 0.03, quad_coeff = 0.01),
  list(name = "Flavobacterium", mu = 0.07, D = 0.01, rad_sensitivity = 0.02, quad_coeff = 0.02),
  list(name = "Ochrobactrum intermedium AM7", mu = 0.06, D = 0.01, rad_sensitivity = 0.07, quad_coeff = 0.015)
)

# Langevin dynamics function with heat and pressure effects
langevin_dynamics <- function(state, parameters, positions, gamma_intensity, thorium_intensity, nutrient_grid, heat_intensity, pressure_field) {
  noise_x <- rnorm(1, mean = 0, sd = sqrt(2 * parameters$D))
  noise_y <- rnorm(1, mean = 0, sd = sqrt(2 * parameters$D))
  noise_z <- rnorm(1, mean = 0, sd = sqrt(2 * parameters$D))
  
  # Radiation sensitivity force
  F_rad_avoidance <- -parameters$rad_sensitivity * gamma_intensity
  
  # Heat influences diffusion
  F_heat <- heat_intensity * parameters$D
  
  # Interaction forces with other species (quadratic interaction)
  F_quad <- c(0, 0, 0)
  for (j in 1:length(species_params)) {
    if (j != parameters$idx) {
      distance <- sqrt(sum((positions[j, ] - state)^2))
      if (distance > 0) {
        interaction <- parameters$quad_coeff / distance
        F_quad <- F_quad + interaction * (positions[j, ] - state)
      }
    }
  }
  
  # Thorium decay affecting nutrient uptake
  nutrient_effect <- -thorium_intensity * nutrient_grid[parameters$idx]
  
  # Pressure field creates additional movement
  F_pressure <- pressure_field[parameters$idx, ]
  
  # Update velocities with noise, forces, heat, and pressure
  dx_dt <- parameters$mu * (F_rad_avoidance + F_heat + F_quad[1] + nutrient_effect + F_pressure[1]) + noise_x
  dy_dt <- parameters$mu * (F_rad_avoidance + F_heat + F_quad[2] + nutrient_effect + F_pressure[2]) + noise_y
  dz_dt <- parameters$mu * (F_rad_avoidance + F_heat + F_quad[3] + nutrient_effect + F_pressure[3]) + noise_z
  
  return(c(dx_dt, dy_dt, dz_dt))
}

# Simulation function for species motility with heat and pressure dynamics
simulate_species_motility <- function(num_species, num_steps, dt, species_params, gamma_intensity, thorium_intensity, nutrient_grid, heat_intensity, pressure_field) {
  positions <- matrix(runif(num_species * 3), ncol = 3)  # Initial (x, y, z) positions
  
  # Store trajectory data for each species
  trajectory_data <- array(NA, dim = c(num_steps, num_species, 3))
  
  for (step in 1:num_steps) {
    for (i in 1:num_species) {
      state <- positions[i, ]
      params <- species_params[[i]]
      params$gamma_intensity <- gamma_intensity
      params$idx <- i  # Current species index
      
      # Perform Langevin dynamics to get the rate of change in positions
      result <- langevin_dynamics(state, params, positions, gamma_intensity, thorium_intensity, nutrient_grid, heat_intensity, pressure_field)
      
      # Update positions with the results
      positions[i, 1] <- max(0, min(1, state[1] + result[1] * dt))
      positions[i, 2] <- max(0, min(1, state[2] + result[2] * dt))
      positions[i, 3] <- max(0, min(1, state[3] + result[3] * dt))
      
      # Store trajectory for each species
      trajectory_data[step, i, ] <- positions[i, ]
    }
  }
  
  return(trajectory_data)
}

# Define UI for the Shiny app
ui <- fluidPage(
  titlePanel("3D Biofilm Species Dynamics with Heat, Pressure, and Radiation Control"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("gamma_intensity", "Gamma Radiation Intensity", min = 0, max = 1, value = 0.05, step = 0.05),
      sliderInput("thorium_intensity", "Thorium Intensity", min = 0, max = 1, value = 0.05, step = 0.05),
      sliderInput("heat_intensity", "Heat Intensity", min = 0, max = 1, value = 0.1, step = 0.05),
      actionButton("run", "Run 45-Step Simulation"),
      actionButton("reset", "Reset Simulation")
    ),
    
    mainPanel(
      plotlyOutput("plot"),
      verbatimTextOutput("status")
    )
  )
)

# Define server logic for the Shiny app
server <- function(input, output, session) {
  # Initialize simulation parameters
  num_species <- length(species_params)
  num_steps <- 45
  dt <- 0.1
  
  # Create nutrient grid and pressure field
  nutrient_grid <- runif(num_species)
  pressure_field <- matrix(runif(num_species * 3), ncol = 3)
  
  # Reactively run the simulation when the button is clicked
  observeEvent(input$run, {
    gamma_intensity <- input$gamma_intensity
    thorium_intensity <- input$thorium_intensity
    heat_intensity <- input$heat_intensity
    
    trajectory_data <- simulate_species_motility(num_species, num_steps, dt, species_params, gamma_intensity, thorium_intensity, nutrient_grid, heat_intensity, pressure_field)
    
    # Generate a plotly plot with 3D scatter for the trajectories
    p <- plot_ly() %>%
      add_trace(
        x = trajectory_data[, 1, 1], y = trajectory_data[, 1, 2], z = trajectory_data[, 1, 3],
        type = 'scatter3d', mode = 'lines+markers', line = list(width = 4),
        name = species_params[[1]]$name, marker = list(color = 'red')
      ) %>%
      add_trace(
        x = trajectory_data[, 2, 1], y = trajectory_data[, 2, 2], z = trajectory_data[, 2, 3],
        type = 'scatter3d', mode = 'lines+markers', line = list(width = 4),
        name = species_params[[2]]$name, marker = list(color = 'blue')
      ) %>%
      add_trace(
        x = trajectory_data[, 3, 1], y = trajectory_data[, 3, 2], z = trajectory_data[, 3, 3],
        type = 'scatter3d', mode = 'lines+markers', line = list(width = 4),
        name = species_params[[3]]$name, marker = list(color = 'green')
      ) %>%
      add_trace(
        x = trajectory_data[, 4, 1], y = trajectory_data[, 4, 2], z = trajectory_data[, 4, 3],
        type = 'scatter3d', mode = 'lines+markers', line = list(width = 4),
        name = species_params[[4]]$name, marker = list(color = 'purple')
      ) %>%
      add_trace(
        x = trajectory_data[, 5, 1], y = trajectory_data[, 5, 2], z = trajectory_data[, 5, 3],
        type = 'scatter3d', mode = 'lines+markers', line = list(width = 4),
        name = species_params[[5]]$name, marker = list(color = 'orange')
      )
    
    output$plot <- renderPlotly({ p })
  })
  
  # Reset button clears the plot
  observeEvent(input$reset, {
    output$plot <- renderPlotly({
      plot_ly()
    })
  })
}

# Run the application
shinyApp(ui = ui, server = server)