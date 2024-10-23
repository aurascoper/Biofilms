## Introduction

This repository contains scripts for simulating biofilm dynamics under radiation stress, evaluating radiotrophic microbial communities, and analyzing microbial interactions under varying environmental conditions. These simulations can model nutrient uptake, microbial fitness, and motility in 2D and 3D environments, under both radiation exposure and nutrient gradients.

### Mid-Level Technical Overview

At the core of this model, we simulate **multi-species biofilm dynamics** under **thorium decay** and **gamma radiation**. The system models cooperative growth using a **Langevin dynamics** framework with species-specific motility, sensitivity to radiation, and interspecies interactions. The Hamiltonian formalism is introduced to capture **phase-locking kernels** for radiation-driven microbial adaptations.

The simulations rely on **stochastic partial differential equations** (SPDEs) for microbial fitness, incorporating **diffusion coefficients**, **mutual interactions**, and **nutrient uptake efficiencies**. The **subcellular localization** data is used to map metabolic functions under these stress conditions.
`;

`
The general Hamiltonian-based system is written as:

$$H_{k-NN} = \sum_{j=1}^{n} P_{s_j}(t) F_j(t, x) \Gamma_s(t, x) $$

Where:
- $$P_{s_j}(t) $$: transition probabilities between species
- $$F_j(t, x) $$: fitness function of species at position x
- $$\Gamma_s(t, x)$$: phase-locked kernel affecting transitions

Time evolution of the system follows:

$$\frac{dq}{dt}$$ $$=$$ $$\frac{\partial H}{\partial p}$$ $$ $$-$$ $$\Gamma_s(t, x)$$ $$F_s(t, x)$$

$$\frac{dp}{dt}$$ = $$-\frac{\partial H}{\partial q}$$ + $$\Gamma_s(t, x)$$ $$F_s(t, x)$$

Species motility, radiation sensitivity, and nutrient uptake are dynamically updated through these Hamiltonian interactions.

`;
`
### reactor_decision_tree.R

This R script simulates microbial interactions in **radiotrophic environments** using **decision tree analysis**. It leverages **thorium decay** and subcellular localization data to model microbial fitness and metabolic activities under radiation stress.

**Usage**:
1. Load the dataset (e.g., \`subcellular_locations.tsv\`).
2. Run the decision tree to analyze microbial fitness under **gamma radiation** and **thorium interactions**.

The decision tree helps identify optimal survival strategies of **extremophilic microbes** by analyzing nutrient uptake, radiation resistance, and species-specific growth models.
`;

`
### biofilms.R

This R script models **biofilm growth dynamics** of multiple species exposed to radiation. It incorporates:
- **Species-specific motility**
- **Radiation sensitivity**
- **Nutrient uptake models**

**Usage**:
1. Define the number of species and set parameters such as **growth rates**, **radiation sensitivities**, and **nutrient uptake efficiencies**.
2. Visualize biofilm growth using **ggplot2** for 2D plots or **Plotly** for interactive 3D plots.

This model explores how different species **cooperate** or **compete** for resources, factoring in radiation-induced stress.
`;

// 3D Biofilm Simulations Overview
`
### biofilms_3d.R

This Julia script enables **3D visualization** of biofilm growth using **PlotlyJS**. It simulates species interactions in a **structured 3D grid**, modeling their growth and nutrient uptake under **radiation and nutrient gradients**.

**Usage**:
1. Define parameters such as species-specific **motility** and **diffusion** rates.
2. Use **PlotlyJS** for 3D visualization of the simulation.

Incorporating **k-means clustering** and **decision tree analysis**, this model helps predict biofilm structures and community interactions over time.
`;

`
### Subcellular Location Data

The subcellular location data provided in \`subcellular_locations.tsv\` and \`subcellular_location_data.tsv\` is essential for mapping **microbial metabolic processes** under **radiation stress**. It includes key biological processes localized in various subcellular compartments and their fitness implications in extreme environments.

**Applications**:
- Understanding **radiotrophic behavior** at the subcellular level.
- Mapping **nutrient uptake efficiency** and **DNA repair** in stressed conditions.
`;

`
### License

This repository is licensed under the Creative Commons. See the \`LICENSE.txt\` file for more information.

---

For any questions or collaboration inquiries, feel free to contact \`hkinder@stlteach.org\`.
`;
