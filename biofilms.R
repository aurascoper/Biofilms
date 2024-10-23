# Load necessary libraries
library(deSolve)
library(stats)

# Define species-specific motility, radiation parameters, and quadratic interaction coefficients
species_params <- list(
  list(name = "New Extremophile", mu = 0.09, D = 0.005, rad_aversion = 0.06, quad_coeff = 0.02),  # Species 1
  list(name = "Cryptococcus neoformans", mu = 0.06, D = 0.01, rad_sensitivity = 0.07, quad_coeff = 0.015),  # Species 2
  list(name = "Deinococcus radiodurans", mu = 0.12, D = 0.01, rad_sensitivity = 0.0, quad_coeff = 0.02),  # Species 3
  list(name = "Bacillus subtilis", mu = 0.1, D = 0.01, rad_sensitivity = 0.03, quad_coeff = 0.01),  # Species 4
  list(name = "Cladosporium sphaerospermum", mu = 0.05, D = 0.01, rad_sensitivity = 0.05, quad_coeff = 0.015),  # Species 5
  list(name = "Aspergillus niger", mu = 0.07, D = 0.01, rad_sensitivity = 0.02, quad_coeff = 0.02),  # Species 6
  list(name = "Shewanella oneidensis", mu = 0.11, D = 0.01, rad_sensitivity = 0.04, quad_coeff = 0.02)  # Species 7
)

# Define the Langevin dynamics function for each species with quadratic terms
langevin_dynamics <- function(t, state, parameters, positions) {
  with(as.list(c(state, parameters)), {
    # Define the random noise term (Gaussian noise) for x and y directions
    noise_x <- rnorm(1, mean = 0, sd = sqrt(2 * D))
    noise_y <- rnorm(1, mean = 0, sd = sqrt(2 * D))
    
    # Handle radiation aversion or sensitivity
    F_rad_avoidance <- 0
    if (!is.null(parameters$rad_aversion)) {
      F_rad_avoidance <- -parameters$rad_aversion * parameters$gamma_intensity
    } else if (!is.null(parameters$rad_sensitivity)) {
      F_rad_avoidance <- -parameters$rad_sensitivity * parameters$gamma_intensity
    }
    
    # Calculate the quadratic interaction forces with all other species
    F_quad <- c(0, 0)
    for (j in 1:nrow(positions)) {
      if (j != parameters$idx) {  # Do not calculate self-interaction
        distance <- sqrt(sum((positions[j, ] - state)^2))
        interaction <- parameters$quad_coeff * distance
        F_quad <- F_quad + interaction * (positions[j, ] - state)
      }
    }
    
    # Calculate the rate of change of position
    dx_dt <- mu * (F_rad_avoidance + F_quad[1]) + noise_x
    dy_dt <- mu * (F_rad_avoidance + F_quad[2]) + noise_y
    
    # Return the rate of change in both x and y directions, limit to within [0, 1]
    new_x <- max(0, min(1, state[1] + dx_dt * dt))  # Limit x to [0, 1]
    new_y <- max(0, min(1, state[2] + dy_dt * dt))  # Limit y to [0, 1]
    return(list(c(new_x - state[1], new_y - state[2])))
  })
}

# Define a function to simulate the motility of all species, including radiation effects and clustering analysis
simulate_species_motility_with_clustering <- function(num_species, num_steps, dt, species_params, gamma_intensity) {
  # Initialize positions of species randomly within a unit square
  positions <- matrix(runif(num_species * 2), ncol = 2)  # Initial (x, y) positions
  
  # Store trajectory data for visualization
  trajectory_data <- list()
  for (i in 1:num_species) {
    trajectory_data[[i]] <- matrix(NA, nrow = num_steps, ncol = 2)
  }
  
  # Iterate through each time step to update positions and perform k-means clustering
  for (step in 1:num_steps) {
    for (i in 1:num_species) {
      # Extract the current position of species i
      state <- positions[i, ]
      
      # Set the parameters for the current species
      params <- species_params[[i]]
      params$gamma_intensity <- gamma_intensity  # Constant gamma radiation intensity
      params$idx <- i  # Track current species index
      
      # Solve the Langevin dynamics for species i using Euler's method
      result <- langevin_dynamics(t = step * dt, state = state, parameters = params, positions = positions)
      
      # Update the position of species i
      positions[i, 1] <- positions[i, 1] + result[[1]][1]  # Update x-coordinate
      positions[i, 2] <- positions[i, 2] + result[[1]][2]  # Update y-coordinate
      
      # Store trajectory for each species
      trajectory_data[[i]][step, ] <- positions[i, ]
    }
    
    # Perform k-means clustering every 10 steps for efficiency in visualization
    if (step %% 10 == 0) {
      kmeans_result <- kmeans(positions, centers = 4)  # Choose 4 clusters to observe group behavior
      
      # Visualize the positions and clusters at the current time step
      plot(positions, col = kmeans_result$cluster, pch = 16, xlim = c(0, 1), ylim = c(0, 1),
           xlab = "X Position", ylab = "Y Position", main = paste("Step:", step, "- K-Means Clustering"))
      points(kmeans_result$centers, col = 1:4, pch = 8, cex = 2)  # Mark cluster centers
      Sys.sleep(0.02)  # Pause for visualization (shorter pause for efficiency)
    }
  }
  
  return(trajectory_data)
}

# Set simulation parameters
num_species <- 7      # Number of species in the hypothetical SCOBY
num_steps <- 1000     # Number of time steps
dt <- 0.1             # Time step size
gamma_intensity <- 0.2  # Constant intensity of gamma radiation due to 14.05 billion year half-life

# Run the simulation of species motility with quadratic interactions and clustering analysis
trajectory_data <- simulate_species_motility_with_clustering(num_species, num_steps, dt, species_params, gamma_intensity)

# Final visualization: Plot the overall trajectories
colors <- c("red", "blue", "purple", "green", "orange", "brown", "pink")
par(mar = c(5, 5, 4, 8))  # Adjust margins to leave space for the legend

# Initialize the plot before adding lines
plot(NULL, xlim = c(0, 1), ylim = c(0, 1), xlab = "X Position", ylab = "Y Position", main = "Species Trajectories with K-Means Clustering Analysis")

for (i in 1:num_species) {
  lines(trajectory_data[[i]][, 1], trajectory_data[[i]][, 2], col = colors[i], lwd = 2)
}

# Add the legend after plotting the lines
legend("topright", inset = c(-0.2, 0), legend = sapply(species_params, function(x) x$name), 
       col = colors, lty = 1, lwd = 2, xpd = TRUE, bty = "n")  # Legend outside the plot
