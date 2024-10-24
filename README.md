
Species motility, radiation sensitivity, and nutrient uptake are dynamically updated through these Hamiltonian interactions.

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

This Julia script enables **3D visualization** of biofilm growth using **PlotlyJS**. It simulates species interactions in a **structured 3D grid**, modeling their growth and nutrient uptake under **radiation and nutrient gradients**.

**Usage**:
1. Define parameters such as species-specific **motility** and **diffusion** rates.
2. Use **PlotlyJS** for 3D visualization of the simulation.

Incorporating **k-means clustering** and **decision tree analysis**, this model helps predict biofilm structures and community interactions over time.

### TypeScript Plot Integration (For Web)

To render interactive plots for biofilm dynamics using Plotly.js, you can use the following TypeScript code.

For questions, collaboration interest, or help, please email hkinder@stlteach.org.