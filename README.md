# Biofilm Dynamics Simulation under Radiation Stress

## Introduction

This repository contains scripts for simulating biofilm dynamics under radiation stress, evaluating radiotrophic microbial communities, and analyzing microbial interactions under varying environmental conditions. These simulations model nutrient uptake, microbial fitness, and motility in 2D and 3D environments under both radiation exposure and nutrient gradients and heat/pressure gradients.

### Mid-Level Overview

At the core of this model, we simulate **multi-species biofilm dynamics** under **thorium decay** and **gamma radiation**. The system models cooperative growth using a **Langevin dynamics** framework, with species-specific motility, sensitivity to radiation, and interspecies interactions. The Hamiltonian formalism is introduced to capture **phase-locking kernels** for radiation-driven microbial adaptations.

The simulations rely on **partial stochastic differential equations (PSDEs)** for microbial fitness, incorporating **diffusion coefficients**, **mutual interactions**, and **nutrient uptake efficiencies**. The **subcellular localization** data is used to map metabolic functions under these stress conditions.

## Equations

The general Hamiltonian-based system is written as:

$$
H_{k-NN} = \sum_{j=1}^{n} P_{s_j}(t) F_j(t, x) \Gamma_s(t, x)
$$

Where:
- $$\(P_{s_j}(t)\)$$ : transition probabilities between species
- $$\(F_j(t, x)\)$$ : fitness function of species at position \(x\)
- $$\(\Gamma_s(t, x)\)$$: phase-locked kernel affecting transitions

The time evolution of the system follows:

$$
\frac{dq}{dt} = \frac{\partial H}{\partial p} - \Gamma_s(t, x) F_s(t, x)
$$

$$
\frac{dp}{dt} = -\frac{\partial H}{\partial q} + \Gamma_s(t, x) F_s(t, x)
$$

// Import Three.js modules
import * as THREE from 'three';

// Setup the scene
const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
const renderer = new THREE.WebGLRenderer();

// Set the size of the renderer
renderer.setSize(window.innerWidth, window.innerHeight);
document.body.appendChild(renderer.domElement);

// Add lighting
const ambientLight = new THREE.AmbientLight(0x404040);
scene.add(ambientLight);

// Example: Add a sphere to represent biofilm behavior
const geometry = new THREE.SphereGeometry(1, 32, 32);
const material = new THREE.MeshStandardMaterial({ color: 0x7f7fff });
const sphere = new THREE.Mesh(geometry, material);
scene.add(sphere);

// Adjust camera position
camera.position.z = 5;

// Animation loop to represent biofilm growth or motion
function animate() {
    requestAnimationFrame(animate);

    // Example of dynamic biofilm behavior (rotation)
    sphere.rotation.x += 0.01;
    sphere.rotation.y += 0.01;

    renderer.render(scene, camera);
}

// Start the animation
animate();

Species motility, radiation sensitivity, and nutrient uptake are dynamically updated through these Hamiltonian interactions.

## Scripts

### reactor_decision_tree.R

This R script can simulate synthetic radiation noise which can be used to select reactor model path for energy or bioremediation, it uses **radiated environments** using **decision tree analysis**. It can leverage the **various nuclear decay constants** and the subcellular localization RNA data to model microbial fitness and metabolic activities under radiation stress.

**Usage**:
1. Load your own dataset or use extremophiles as is and select sequence from datasets (e.g., `subcellular_locations.tsv, subcellular_locations_data`, etc). This model uses various radioprotective functions like Cryptococcus Neoformans and Deinonoccus Radiodurans, and Shewanella algae as well as other extremophiles.
2. Fill in the decision tree from scratch for bioreactor engineering using energy/fitness landscapes to analyze bioreactor design strategiess that use **gamma radiation** and **thorium interactions** or others.

The decision tree helps identify optimal survival strategies of **extremophilic microbes** by analyzing nutrient uptake, radiation resistance, and species-specific growth models.

### biofilms.R

This R script models **biofilm growth dynamics** of multiple species exposed to radiation. It incorporates:
- **Species-specific motility**
- **Radiation sensitivity**
- **Nutrient uptake models**

**Usage**:
1. Define the number of species and set parameters such as **growth rates**, **radiation sensitivities**, and **nutrient uptake efficiencies**.
2. Visualize biofilm growth using **ggplot2** for 2D plots or **Plotly** for interactive 3D plots.

<div>
  <h2>KMeans Species Trajectory</h2>
  <p>This GIF shows the simulation of species motility and clustering behavior over time using KMeans clustering:</p>
  <img src="https://i.giphy.com/media/v1.Y2lkPTc5MGI3NjExeHpqNmNvYTFjMHIyeGNpa3RtOHE3cmI1ajZhbHJib2hjODZ2bWdpeSZlcD12MV9pbnRlcm5hbF9naWZfYnlfaWQmY3Q9Zw/oKfaA8IEcAoHQMlZG3/giphy.gif" alt="KMeans Species Trajectory Animation">
</div>

This model explores how different species **cooperate** or **compete** for resources, factoring in radiation-induced stress.

### biofilms_3d.R

## Project Demo

[![Watch the video](https://https://drive.google.com/file/d/1huR7BzRkdLneXVtaaHBpmnxO7F_77v_a/view?usp=drive_link/)

This Julia script enables **3D visualization** of biofilm growth using **PlotlyJS**. It simulates species interactions in a **structured 3D grid**, modeling their growth and nutrient uptake under **radiation and nutrient gradients**.

**Usage**:
1. Define parameters such as species-specific **motility** and **diffusion** rates.
2. Use **PlotlyJS** for 3D visualization of the simulation.

Incorporating **k-means clustering** and **decision tree analysis**, this model helps predict biofilm structures and community interactions over time.

For questions, collaboration interest of help, please email hkinder@stlteach.org
