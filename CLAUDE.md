# Biofilms — Modeling Radiotrophic Fitness

Simulation suite supporting **Kinder & Faulkner (2026)**, a bioRxiv preprint on radiotrophic biofilm communities in cylindrical bioreactors for nuclear bioremediation. Multi-scale stack: Langevin PSDE dynamics in R, a 3D Cellular Potts Model in Julia coupled to a melanin/radiation field, and a radiodialysis membrane transport PDE system. Includes a Hamiltonian-kNN reactor decision tree and a Shiny interactive bioreactor visualizer.

## Stack
- **R** (Shiny, deSolve, stats) — original Langevin PSDE simulator and the two interactive apps
- **Julia** — 3D Cellular Potts Model with coupled fields (CairoMakie for figures)
- LaTeX (preprint manuscript)
- Tableau (`.twb`) — primary fuel source geolocator viz
- No build system, no test suite — these are research scripts producing figures and animations

## Run / test / build
```sh
# R scripts — each is self-contained, run interactively or via Rscript
Rscript biofilms.R                          # original flat-domain Langevin PSDE
Rscript biofilms_3d.R                       # cylindrical bioreactor (launches Shiny app)
Rscript biofilms_radiodialysis.R            # radiodialysis PDE (Shiny app)
Rscript reactor_decision_tree.R             # Hamiltonian kNN reactor selection

# Julia simulations
julia biofilms_potts.jl                     # 3D Cellular Potts Model + radiodialysis coupling

# Preprint
cd preprint && pdflatex modeling_radiotrophic_fitness.tex
```

## Layout
- `biofilms.R` — flat-domain Langevin PSDE with 7 species + k-means clustering (original)
- `biofilms_3d.R` — cylindrical bioreactor Shiny app
- `biofilms_potts.jl` — 3D CPM (`H_adhesion + H_volume + H_radiation + H_pairwise + H_melanin`) with coupled melanin RD and radiation fields per preprint §3.3, §3.5, §3.7
- `biofilms_radiodialysis.R` — Shiny app for the membrane transport PDE system (§3.9)
- `reactor_decision_tree.R` — Hamiltonian kNN over reactor design space
- `scoby_3d.jl` — placeholder (empty)
- `preprint/` — LaTeX source, compiled PDF, submission zip, `figures/` (PDF + PNG outputs from the simulations)
- `assets/` — README preview images
- `subcellular_location*.tsv` — supporting dataset (subcellular RNA locations)
- `Citations.md`, `nihms-983981.pdf`, `Wolfram's Scoby.pdf`, `file.pdf` — reference material
- Pre-rendered animations: `biofilm_dynamics*.gif`, `biofilms_3d.mp4.mp4` (yes, double `.mp4`), `12_step_24_second_animation_sped_up.gif`, `kmeans_species_trajectory (1).gif`
- `BiofilmMining.mp4.mp4`, `Primary Fuel Source Geolocator.mp4` + `.twb` — Tableau workbook + recording

## Conventions
- **Species params live as named lists at the top of each R file** — `species_params <- list(list(name=..., mu=..., D=..., rad_sensitivity=..., quad_coeff=...))`. Add new species here; downstream code iterates the list.
- Julia CPM uses **MCS (Monte Carlo Sweeps)** as the time unit; coupled fields (melanin, nutrient, radiation) update once per MCS per Eq. 7 (melanin RD) and Eq. 5 (radiation field) of the preprint
- Figures land in `preprint/figures/` as both PDF (for LaTeX `\includegraphics`) and PNG (for README previews)
- The preprint and the simulations are versioned together — when math changes in code, update the corresponding section in `preprint/modeling_radiotrophic_fitness.tex`

## Gotchas
- **Two languages, two ecosystems.** R scripts are independent of the Julia CPM — they model different scales (PSDE/macroscopic vs. lattice/cellular). Don't try to unify them.
- **Empty file:** `scoby_3d.jl` is a 0-byte placeholder. The active Julia simulation is `biofilms_potts.jl`.
- **Duplicate-extension files (`*.mp4.mp4`)** are intentional artifacts from the recording tool — don't rename, the `README.md` and the preprint LaTeX may reference them by that exact name.
- **Tableau workbook lockfile** (`~Primary Fuel Source Geolocator__17960.twbr`) is a session marker — safe to ignore, don't commit your own.
- **Subcellular location data is duplicated** (`subcellular_location.tsv`, `subcellular_location_data.tsv`, `subcellular_locations.tsv`) — check which one a given script reads before "cleaning up". The first two are identical (1.9 MB); the third is a small summary (1.2 KB).
- No CI, no tests, no lint config — this is a paper-supporting research repo. Changes should reproduce the paper's figures byte-for-byte unless you're explicitly revising the manuscript.
- Large binary assets (videos, PDFs, animations) live in the repo without Git LFS — `git clone` is heavy.
