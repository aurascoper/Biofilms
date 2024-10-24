## Biofilm Dynamics Simulation under Radiation Stress

## Introduction

This repository contains scripts for simulating biofilm dynamics under radiation stress, evaluating radiotrophic microbial communities, and analyzing microbial interactions under varying environmental conditions. These simulations model nutrient uptake, microbial fitness, and motility in 2D and 3D environments under both radiation exposure and nutrient gradients and heat/pressure gradients.

### Mid-Level Overview

At the core of this model, we simulate **multi-species biofilm dynamics** under **thorium decay** and **gamma radiation**. The system models cooperative growth using a **Langevin dynamics** framework, with species-specific motility, sensitivity to radiation, and interspecies interactions. The Hamiltonian formalism is introduced to capture **phase-locking kernels** for radiation-driven microbial adaptations.

The simulations rely on **partial stochastic differential equations (SPDEs)** for microbial fitness, incorporating **diffusion coefficients**, **mutual interactions**, and **nutrient uptake efficiencies**. The **subcellular localization** data is used to map metabolic functions under these stress conditions.

## Equations

The general Hamiltonian-based system is written with the factors of: 

Species motility, radiation sensitivity, and nutrient uptake are supervised-dynamically and updated through these Hamiltonian and Langevin, and Euler-Langrangian interactions.

## Scripts

### reactor_decision_tree.R

This R script can simulate synthetic radiation noise, which can be used to select a reactor model for energy or bioremediation. It uses **radiated environments** and **decision tree analysis**. It can leverage the **thorium decay constant** and subcellular localization RNA data to model microbial fitness and metabolic activities under radiation stress.

**Usage**:
1. Load the dataset (e.g., `subcellular_locations.tsv`).
2. Run the decision tree to analyze microbial fitness under **gamma radiation** and **thorium interactions**.

The decision tree helps identify optimal survival strategies of **extremophilic microbes** by analyzing nutrient uptake, radiation resistance, and species-specific growth models.

### biofilms.R

This R script models **biofilm growth dynamics** of multiple species exposed to radiation. It incorporates:
- **Species-specific motility**
- **Radiation sensitivity**
- **Nutrient uptake models**

**Usage**:
1. Define the number of species and set parameters such as **growth rates**, **radiation sensitivities**, and **nutrient uptake efficiencies**.
2. Visualize biofilm growth using **ggplot2** for 2D plots or **Plotly** for interactive 3D plots.

This model explores how different species **cooperate** or **compete** for resources, factoring in radiation-induced stress.

### biofilms_3d.R

<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>KMeans Species Trajectory</title>
</head>
<body>
    <h1>KMeans Species Trajectory</h1>
    <p>This animation demonstrates the species' motility and clustering behavior over time using KMeans clustering.</p>
    <img src="https://drive.google.com/uc?export=view&id=1IgOHbf3GupcqlWmakA0wooLYh-orulZ9" alt="KMeans Species Trajectory">
</body>
</html>


This R script enables **3D visualization** of biofilm growth using **PlotlyJS**. It can also work with julia preprocessing or ML and Python as well. It simulates species interactions in a **structured 3D grid**, modeling their growth and nutrient uptake under **radiation and nutrient gradients**.

**Usage**:
1. Define parameters such as species-specific **motility** and **diffusion** rates.
2. Use **PlotlyJS** for 3D visualization of the simulation.

It can incorporate **k-means clustering** and **decision tree analysis**, with the aim to model and help predict biofilm structures and community interactions over time.
