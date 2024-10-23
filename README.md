To run the R scripts, make sure you have the following packages installed:

R
Copy code
install.packages(c("deSolve", "plotly", "kmeans", "data.table", "ggplot2"))
Julia Environment
For the 3D simulations, install the necessary Julia packages by running:

julia
Copy code
using Pkg
Pkg.add(["DifferentialEquations", "PlotlyJS", "MLJ", "Clustering"])

The R script reactor_decision_tree.R generates the invariance in the decisions on reactor design and location which can be used to generate noise interactions in geospatial and radioactive and radiotrophic environments using decision tree analysis. It utilizes thorium decay data and CAN use subcellular localization to assess microbial fitness and metabolic activity under radiation stress.

Usage:

Load the dataset for subcellular locations (e.g., subcellular_locations.tsv).
Run the decision tree to analyze fitness functions of metagenomic analyses for screening of extremophilic microbials under gamma radiation and thorium interactions.

biofilms.R Simulates the biofilm growth dynamics of multiple species exposed to radiation. Includes nutrient uptake and species-specific motility models.

Usage:

Set the number of species, general sensitivity and motility gradients other parameters in the fitness function like specific growth rates and radiation sensitivities in the model.
Visualize biofilm growth using ggplot2 for 2D plots or plotly for interactive 3D plots.
biofilms_3d.R
This script uses plotly for 3D biofilm growth simulations. It visualizes the biofilm as species grow and interact under radiation and nutrient gradients.

Usage:

Specify parameters like species motility and diffusion.
Use PlotlyJS to visualize the 3D simulation.
Subcellular Location Data
The subcellular location data is provided in subcellular_locations.tsv and subcellular_location_data.tsv. These datasets contain key biological processes mapped to their respective subcellular locations and can be used to map metabolic functions under radiation stress.

License
This repository is licensed under the Creative Commons. See the LICENSE.txt file for more information.

If there are questions or collaborator interest, please email hkinder@stlteach.org 
