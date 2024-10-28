## Biofilm Dynamics Simulation under Radiation Stress

## Introduction

This repository contains scripts for simulating biofilm dynamics under radiation stress, evaluating radiotrophic microbial communities, and analyzing microbial interactions under varying environmental conditions. These simulations model nutrient uptake, microbial fitness, and motility in 2D and 3D environments under both radiation exposure and nutrient gradients and heat/pressure gradients.

### Mid-Level Overview

At the core of this model, we simulate **multi-species biofilm dynamics** under **thorium decay** and **gamma radiation**. The system models cooperative growth using a **Langevin dynamics** framework, with species-specific motility, sensitivity to radiation, and interspecies interactions. The Hamiltonian formalism is introduced to capture **phase-locking kernels** for radiation-driven microbial adaptations.

The simulations rely on **partial stochastic differential equations (PSDEs)** for microbial fitness, incorporating **diffusion coefficients**, **mutual interactions**, and **nutrient uptake efficiencies**. The **subcellular localization** data is used to map metabolic functions under these stress conditions.

We can include Subcellular data, for potential thermodynamically stable molecular pathways using RadioDyalysis with mutation rates. Included is the Human Subcellular RNA Sequence Numbers per Location in our Cells.

![Subcellular RNA Locations.png](https://github.com/aurascoper/Biofilms/blob/8f6bcbbdea8bd81300404d61aa0f5749b5334f92/Subcellular%20RNA%20Locations.png)

## Equations

The general Hamiltonian-based system is written with the factors of: 

Species motility, radiation sensitivity, and nutrient uptake are supervised-dynamically and updated through these Hamiltonian and Langevin, and Euler-Langrangian interactions.

$$
\frac{dq}{dt} = \frac{\partial H}{\partial p}, \quad \frac{dp}{dt} = -\frac{\partial H}{\partial q} + \eta(t)
$$

## Scripts

### reactor_decision_tree.R

This R script can synthesize radiosympletic noise, which can be used to select a reactor model for bioenergy, particular biological product supplies, or bioremediation. It uses **radiated environments** and **decision tree analysis**. It can leverage the **thorium decay constant** and subcellular localization RNA data to model microbial fitness and metabolic activities under radiation stress.

**Usage Optional**:
1. Load the dataset (e.g., `subcellular_locations.tsv`).
2. Run a simple decision tree analysis to analyze cellular automata that could sync or **lock phases** with microbial fitness factors under **gamma radiation** and **thorium or your own radioactive or radiodialytic interactions**.

The decision tree helps identify optimal survival strategies of **extremophilic microbes** by analyzing nutrient uptake, radiation resistance, and species-specific growth models.

### biofilms.R

![kmeans_species_trajectory (1).gif](https://github.com/aurascoper/Biofilms/blob/b7a111904fc0d8f70b4df84e1f13eb9728e00ce5/kmeans_species_trajectory%20(1).gif)
![biofilm_dynamics_7_species.gif](https://github.com/aurascoper/Biofilms/blob/91ded6274b16aa950569f49d9ec51f23d4f729e1/biofilm_dynamics_7_species.gif)
![biofilm_dynamics.gif](https://github.com/aurascoper/Biofilms/blob/91ded6274b16aa950569f49d9ec51f23d4f729e1/biofilm_dynamics.gif)
![Functor.jl Coordinate .gif](https://github.com/aurascoper/Biofilms/blob/b87615d25b630a393ea00211e73d69dfe3d196b5/12_step_24_second_animation_sped_up.gif)

This R script models **biofilm growth dynamics** of multiple species exposed to radiation. It incorporates:
- **Species-specific motility**
- **Radiation sensitivity**
- **Nutrient uptake models**

The same can be said for the Optimization Problems that can be done with Bioreactor.jl

This Growth curve was calculated with JuMP and Cbc in Julia and is currently operating in JupyterLab. It can be interpreted as a scalar graph for scaling molecular growth rates of interes (Billions of CFUs (colony forming units)).

![Bioreactor Growth Curves.png.png](https://github.com/aurascoper/Biofilms/blob/abcaeb99cc8d98cbae08995b925836c86cc0993c/Bioreactor%20Growth%20Curves.png.png)

**Usage**:
1. Define the number of species and set parameters such as **growth rates**, **radiation sensitivities**, and **nutrient uptake efficiencies**.
2. Visualize biofilm growth using **ggplot2** for 2D plots or **PlotlyJS** for interactive 3D plots.

This model explores how different species **cooperate** or **compete** for resources, factoring in radiation-induced stress.

### biofilms_3d.R

This R script enables **3D visualization** of biofilm growth using **PlotlyJS**. It can also work with julia preprocessing or ML and Python as well. It simulates species interactions in a **structured 3D grid**, modeling their growth and nutrient uptake under **radiation and nutrient gradients**.

**Usage**:
1. Define parameters such as species-specific **motility** and **diffusion** rates.
2. Use **PlotlyJS** for 3D visualization of the simulation.

It can incorporate **k-means clustering** and **decision tree analysis** based on a simple Hamiltonian = k-nn network, with the aim to model and help predict biofilm structures and community interactions over time.
