# Load required libraries
library(ggplot2)
library(dplyr)
library(MASS)
library(class)  # For k-NN algorithm

# Step 1: Data setup
data <- data.frame(
  Element = c("Thorium", "Americium", "Curium", "Neptunium"),
  Cost = c(7500, 1500000, 3000000, 10000),                # Cost in USD/kg
  EnergyDensity = c(79400, 110000, 3000000, 110000),       # Energy Density in MJ/kg
  RTG = c(0, 1, 1, 1),                                    # Recommended for RTG (1/0)
  Reactor = c(1, 0, 0, 1)                                 # Recommended for Reactor (1/0)
)

# View the dataset
print(data)

# Logistic regression for RTG use case
rtg_model <- glm(RTG ~ Cost + EnergyDensity, data = data, family = binomial)

# Logistic regression for Reactor use case
reactor_model <- glm(Reactor ~ Cost + EnergyDensity, data = data, family = binomial)

# Stepwise regression for both models
stepwise_rtg <- stepAIC(rtg_model, direction = "both")
stepwise_reactor <- stepAIC(reactor_model, direction = "both")

# Step 2: Scatter plot of samples
ggplot(data, aes(x = Cost, y = EnergyDensity, color = as.factor(RTG))) +
  geom_point(size = 4) +
  labs(title = "Scatter Plot: Cost vs Energy Density for RTG Use Case",
       x = "Cost (USD/kg)",
       y = "Energy Density (MJ/kg)",
       color = "RTG Recommendation") +
  theme_minimal()

# Step 3: k-NN classification on the dataset

# Define feature matrix and target vector for RTG
features <- data[, c("Cost", "EnergyDensity")]
target_rtg <- data$RTG

# Normalize the data
features_scaled <- as.data.frame(scale(features))

# Create grid for decision boundary
grid <- expand.grid(Cost = seq(min(features_scaled$Cost), max(features_scaled$Cost), length = 100),
                    EnergyDensity = seq(min(features_scaled$EnergyDensity), max(features_scaled$EnergyDensity), length = 100))

# Apply k-NN with k = 3
k <- 3
knn_predictions <- knn(train = features_scaled, test = grid, cl = target_rtg, k = k)

# Convert back to original scale for plotting
grid$Cost <- grid$Cost * sd(features$Cost) + mean(features$Cost)
grid$EnergyDensity <- grid$EnergyDensity * sd(features$EnergyDensity) + mean(features$EnergyDensity)

# Plot k-NN decision boundary
ggplot(grid, aes(x = Cost, y = EnergyDensity)) +
  geom_tile(aes(fill = knn_predictions), alpha = 0.3) +
  geom_point(data = data, aes(x = Cost, y = EnergyDensity, color = as.factor(RTG)), size = 4) +
  labs(title = "k-NN Decision Boundary with k = 3 for RTG Use Case",
       x = "Cost (USD/kg)",
       y = "Energy Density (MJ/kg)",
       fill = "k-NN Prediction",
       color = "RTG Recommendation") +
  theme_minimal()

# Optionally, repeat the k-NN classification for Reactor use case by adjusting target
target_reactor <- data$Reactor
knn_predictions_reactor <- knn(train = features_scaled, test = grid, cl = target_reactor, k = k)

# Plot k-NN decision boundary for Reactor
ggplot(grid, aes(x = Cost, y = EnergyDensity)) +
  geom_tile(aes(fill = knn_predictions_reactor), alpha = 0.3) +
  geom_point(data = data, aes(x = Cost, y = EnergyDensity, color = as.factor(Reactor)), size = 4) +
  labs(title = "k-NN Decision Boundary with k = 3 for Reactor Use Case",
       x = "Cost (USD/kg)",
       y = "Energy Density (MJ/kg)",
       fill = "k-NN Prediction",
       color = "Reactor Recommendation") +
  theme_minimal()
