# Modeling Radiotrophic Fitness

**Hunter Kinder**, B.A., M.A.T.L. — Independent Researcher, Missouri, USA  
**Brett Faulkner**, B.Sc, BioPhys — Independent Researcher

*Preprint — submitted to bioRxiv (Systems Biology)*  
*Version 1.0 — April 2026*

---

## Abstract

This work presents a mathematical framework to describe the radiotrophic and adaptive dynamics of microbial communities exposed to gamma radiation. The developed model integrates Langevin dynamics with reaction-diffusion equations to simulate the spatial motility, growth, and nonlinear interactions of radiotrophic fungi, extremophiles, and radiation-sensitive averse species. In this framework, radiation gradients and chemical factors, such as melanin synthesis, are introduced as modifiable parameters to analyze their influence on biofilm growth and adaptive responses. The model leverages stochastic differential equations to capture complex dependencies between environmental conditions and microbial interactions. Numerical simulations of seven interacting species — including *Cryptococcus neoformans*, *Deinococcus radiodurans*, *Cladosporium sphaerospermum*, *Bacillus subtilis*, *Aspergillus niger*, *Shewanella oneidensis*, and *Ochrobactrum intermedium* AM7 — demonstrate that radiation gradients drive spatial niche partitioning, with melanized radiotrophic fungi concentrating in high-radiation zones while radiosensitive metal-reducing bacteria occupy attenuated peripheries. This work contributes a predictive modeling tool capable of elucidating the nonlinear adaptive behaviors of radiotrophic communities in extreme environments, with applications to nuclear bioremediation.

---

## 1. Introduction

Microbial communities possess an extraordinary ability to adapt and survive under extreme environmental stressors, such as elevated radiation levels. Within such communities, radiotrophic fungi and other extremophilic microorganisms have evolved specialized mechanisms that allow them to thrive in high-radiation environments. This study aims to develop a comprehensive mathematical framework that captures the complex interplay between these organisms' nonlinear interactions, radiation sensitivities, and chemotactic behaviors within biofilm communities. By leveraging a combined approach that integrates Langevin dynamics, reaction-diffusion modeling, and a Hamiltonian-based formulation, we simulate and explore the spatiotemporal evolution of these communities under radiation stress and other external influences.

The practical motivation for this work lies in nuclear bioremediation. Contaminated sites such as the Chernobyl Exclusion Zone and the Hanford Nuclear Reservation host complex microbial communities in which radiotrophic species coexist with metal-reducing bacteria capable of immobilizing radionuclides. Mathematical models that predict community dynamics under spatially heterogeneous radiation fields are essential for designing and optimizing bioremediation strategies. The framework developed here provides the theoretical foundation for such predictions.

---

## 2. Related Work

### 2.1 Radiotrophic Organisms and Melanin-Mediated Radiotrophy

The discovery that melanized fungi thrive in high-radiation environments has fundamentally reshaped our understanding of biological energy transduction. Zhdanova et al. [1] first documented the prevalence of darkly pigmented fungi within the damaged reactor at the Chernobyl Nuclear Power Plant, noting that species such as *Cladosporium sphaerospermum* dominated the mycobiota of the most contaminated structures. Subsequent work by Dadachova et al. [2] demonstrated that ionizing radiation alters the electronic properties of melanin, enhancing electron-transfer activity in melanized cells of *Cryptococcus neoformans* and *Wangiella dermatitidis*. Melanized cells exposed to radiation levels approximately 500 times above background grew significantly faster, accumulated more biomass, and incorporated three-fold more ¹⁴C-acetate than non-irradiated melanized cells or irradiated albino mutants. Dadachova and Casadevall [3] subsequently proposed the term "radiosynthesis" to describe this phenomenon, drawing an analogy to photosynthesis in which melanin serves a role loosely comparable to chlorophyll. More recently, Shunk et al. [4] cultivated *C. sphaerospermum* aboard the International Space Station for 26 days, demonstrating a growth advantage of approximately 21% under space radiation conditions and measurable attenuation of ionizing radiation beneath the fungal biomass.

The radiation resistance of *Deinococcus radiodurans* operates through fundamentally different mechanisms. Rather than harnessing radiation for metabolic benefit, *D. radiodurans* survives doses exceeding 12 kGy through a combination of efficient DNA double-strand break repair via extended synthesis-dependent strand annealing (ESDSA), multiple genome copies enabling homologous recombination, and a potent manganese-based antioxidant system that protects proteins from oxidative damage [5, 6]. Daly et al. [7] showed that small-molecule Mn²⁺-metabolite complexes specifically protect the proteome against radiation-induced carbonylation, establishing that protein protection, rather than DNA repair capacity alone, governs extreme radioresistance. *Aspergillus niger*, while not radiotrophic in the strict sense, produces melanin that confers substantial radioprotection, with melanized fungi generally exhibiting LD₁₀ values approaching or exceeding 1 kGy [3].

### 2.2 Mathematical Models of Biofilm Communities

The mathematical modeling of biofilm communities has a rich history beginning with the foundational one-dimensional multispecies model of Wanner and Gujer [8], which coupled reaction-diffusion equations for substrate transport with biomass conservation laws to predict biofilm thickness and spatial species distributions. Xavier et al. [9] extended this framework to multiple spatial dimensions, providing a deterministic continuum approach to heterogeneous biofilm development. Individual-based modeling approaches, exemplified by the iDynoMiCS platform of Lardon et al. [10], complement continuum models by resolving cell-level stochastic interactions. Cross-diffusion biofilm models have been analyzed by Sonner et al. [11], who established mathematical properties of volume-filling multispecies systems that exhibit porous-medium-type degeneracy. Our PSDE framework extends these precedents by incorporating radiation-dependent fitness terms and phase-locking dynamics that couple species interactions to external radiation fields, moving beyond the substrate-limited growth paradigms that dominate existing biofilm models.

### 2.3 Hamiltonian and Stochastic Approaches in Theoretical Ecology

The application of Hamiltonian mechanics to ecological systems remains comparatively underexplored, despite the natural conservation laws that govern energy flow in ecosystems. Symplectic integration methods, which preserve the geometric structure of Hamiltonian flows, have been advocated for ecological modeling by Diele et al. [12] to ensure that numerical solutions maintain essential qualitative dynamics. The Kuramoto model of coupled phase oscillators [13] provides a well-established framework for synchronization phenomena in biological systems, from neural networks to circadian rhythms. Acebrón et al. [14] reviewed the model as a paradigm for synchronization, noting its generalization to systems with time delays, frequency-weighted coupling, and heterogeneous interaction topologies. Our phase-locking kernel Γ_s(t,x) draws directly on this tradition, extending it to multispecies microbial communities where species-level oscillatory dynamics in growth and resource acquisition couple through shared radiation and nutrient fields. The use of Langevin stochastic dynamics to describe microbial population fluctuations follows naturally from the Fokker-Planck formalism and has been applied to competing microbial populations in chemostat systems [15], where stochastic differential equations revealed long-term competitive outcomes differing from deterministic predictions.

### 2.4 Bioremediation Context

The practical motivation for modeling radiotrophic communities lies in nuclear bioremediation. *Shewanella oneidensis* MR-1 is a model dissimilatory metal-reducing bacterium capable of reducing soluble U(VI) to insoluble U(IV) via extracellular electron transfer through c-type cytochromes [16, 17]. *Ochrobactrum intermedium* AM7, isolated from soil near the Kakrapar Atomic Power Station in India, tolerates high concentrations of thorium and produces exopolysaccharides (EPS) that complex with Th(IV), offering a pathway toward bioremediation of actinide-contaminated environments [18]. The combination of radiotrophic fungi that actively benefit from radiation with metal-reducing bacteria that immobilize radionuclides suggests the possibility of engineered biofilm consortia for in situ remediation. The fitness landscape framework introduced by Wright [19] and formalized through the NK model of Kauffman [20] provides a conceptual foundation for understanding how epistatic interactions among species traits shape adaptation under radiation stress, a perspective that our Hamiltonian kNN decision tree operationalizes quantitatively.

---

## 3. Mathematical Framework

### 3.1 Symbols and Notation

| Symbol | Description |
|--------|-------------|
| F_s(t,x) | Fitness function of species s, evolving over time and space |
| H(p,q) | Hamiltonian function defining total system energy |
| p, q | Phase-space variables (momentum and position) |
| η_s(t,x) | Random environmental noise term |
| D_s, μ_s | Species-specific diffusion and motility coefficients |
| P_sj(t) | Permutation matrix describing species transitions |
| γ(t,x) | Gamma irradiation field |
| k_B | Boltzmann constant |
| S[q(t)] | Action functional for species trajectories |
| ∇ | Gradient operator |
| Γ_s(t,x) | Phase-locked kernel for species s |
| β_s,ion | Ionizing radiation sensitivity coefficient |
| α_s,nir | Non-ionizing radiation coupling coefficient |
| θ_s | Adaptive phase-locking coefficient |

### 3.2 Partial Stochastic Differential Equation (PSDE) for Species Fitness

The fitness of each species is governed by the following PSDE:

$$\partial_t F_s(t,x) = \nabla \cdot (D_s \nabla F_s) - \nabla \cdot \left(\mu_s \sum_{j=1}^n P_{sj}(t) F_j\right) + R_s + \eta_s - \beta_{s,\text{ion}} I F_s + \gamma_s \Delta_s - \alpha_{s,\text{nir}} N F_s + \theta_s H_s + C_s$$

The terms encode:
- **Diffusion**: D_s ∇F_s — passive spatial spreading
- **Directed motility**: μ_s Σ_j P_sj F_j — chemotaxis and radiation-gradient-directed movement
- **Deterministic reaction**: R_s(t,x) — nutrient-dependent growth
- **Stochastic noise**: η_s(t,x) = σ_s ξ(t,x) — white noise perturbation
- **Ionizing damage**: −β_s,ion I(t,x) F_s — radiation-induced fitness reduction
- **Phase-locking adjustment**: γ_s Δ_s — adaptive synchronization gain
- **Non-ionizing effect**: −α_s,nir N(t,x) F_s — UV/heat-driven fitness coupling
- **Hamiltonian interaction**: θ_s H_s(t,x) — energy-conserving inter-species force
- **External forcing**: C_s(t,x) — environmental control inputs

### 3.3 Hamiltonian Framework

The total Hamiltonian for the multi-species biofilm is:

$$H = \sum_{i=1}^N \left[\frac{1}{2} m_i v_i^2 + U_i(x_i)\right] + \sum_{i \neq j} V_{ij}(r_{ij}) + \sum_{k=1}^K W_k(t,x)$$

where m_i is effective species biomass density, U_i(x_i) is the potential energy in the nutrient concentration field, V_ij(r_ij) is the pairwise interaction energy (mutualistic: V_mutual = −γ exp(−r²/σ²)), and W_k(t,x) captures external energy contributions from radiation exposure and stochastic effects.

### 3.4 Symplectic Integration

Symplectic integrators maintain the structure of the Hamiltonian system, ensuring energy and momentum conservation over extended simulation windows. Phase-space trajectories evolve via the phase-locked canonical equations:

$$\frac{dq}{dt} = \frac{\partial H}{\partial p} - \Gamma_s(t,x) F_s(t,x), \quad \frac{dp}{dt} = -\frac{\partial H}{\partial q} + \Gamma_s(t,x) F_s(t,x)$$

The leapfrog/Verlet scheme is employed numerically: (1) evaluate H_s based on current state; (2) adjust phase-locking coefficient Δ_s; (3) integrate via symplectic step.

### 3.5 Radiation Field Models

**Ionizing radiation (gamma):**
$$I(t,x) = I_\gamma \exp(-\kappa x)$$

where I_γ is the initial intensity and κ is the attenuation coefficient scaled by biofilm density.

**Non-ionizing radiation (UV):**
$$N(t,x) = N_{UV} \cos(\omega t)$$

### 3.6 Melanin Reaction-Diffusion

The diffusion-driven melanin production by radiotrophic fungi is modeled as:

$$\frac{\partial M}{\partial t} = D_M \nabla^2 M + \alpha_M \cdot N_{\text{RadioF}} \cdot R(t,x)$$

where D_M is the melanin diffusion coefficient, α_M the radiation-driven production rate, and N_RadioF the local radiotrophic cell density.

### 3.7 Nutrient Uptake with Adaptive Feedback

$$R_s(t,x) = \alpha_s C_s(t,x)\left(1 - \frac{F_s(t,x)}{K_s}\right) + A_s(t) C_s(t,x)$$

$$A_s(t) = \gamma_s \phi_s(t) - \beta_{s,\text{ion}} I(t)$$

where φ_s(t) = cos(ω_s t − θ_s) is the phase alignment function and A_s(t) is the adaptive feedback driven by phase-locking and radiation stress.

### 3.8 Hamiltonian kNN Decision Tree

To model species transition probabilities modulated by neighborhood phase states:

$$H_{\text{k-NN}} = \frac{1}{\sigma_n} \sum_{j=1}^n P_{sj}(t)\, F_j(t,x)\, \Gamma_s(t,x)$$

### 3.9 Radiation Damage and Resilience

$$D_{\text{rad},s}(t,x) = -\beta_{s,\text{rad}} I(t,x) F_s(t,x)$$

Species with higher melanin production or stronger DNA repair exhibit lower β_s,rad values.

### 3.10 Biofilm Mechanical Properties

The EPS matrix is modeled as a Finite Extensible Nonlinear Elastic (FENE) entropic spring:

$$U(r) = -\frac{1}{2} k R^2 \ln\left(1 - \frac{r^2}{R^2}\right)$$

Viscoelastic creep follows the Kelvin-Voigt model:

$$\sigma(t) = E\epsilon(t) + \eta \frac{d\epsilon}{dt}$$

with stress relaxation time τ = η/E governing biofilm restructuring dynamics.

---

## 4. Parameter Estimation

Biologically justified parameter ranges for the six primary modeled species, derived from published experimental data:

| Parameter | Symbol | Species | Value Range | Units | Biological Basis | Ref |
|-----------|--------|---------|-------------|-------|-----------------|-----|
| Diffusion coefficient | D_s | *C. neoformans* | 0.01–0.10 | μm²/s | Passive Brownian diffusion, ~5 μm diameter yeast | [21] |
| Diffusion coefficient | D_s | *D. radiodurans* | 0.05–0.50 | μm²/s | Small coccoid (~1.5 μm); higher D from smaller size | [5] |
| Diffusion coefficient | D_s | *B. subtilis* | 0.10–1.00 | μm²/s | Rod-shaped motile cells; active motility augments D | [22] |
| Diffusion coefficient | D_s | *C. sphaerospermum* | 0.005–0.05 | μm²/s | Filamentous hyphae; low effective diffusion | [4] |
| Diffusion coefficient | D_s | *A. niger* | 0.005–0.05 | μm²/s | Filamentous; comparable to *C. sphaerospermum* | [3] |
| Diffusion coefficient | D_s | *S. oneidensis* | 0.10–0.80 | μm²/s | Facultative anaerobe, rod-shaped, flagellated | [16] |
| Motility coefficient | μ_s | *C. neoformans* | 0.00–0.05 | μm/s | Non-motile yeast; residual drift from EPS flow | [2] |
| Motility coefficient | μ_s | *D. radiodurans* | 0.00–0.01 | μm/s | Non-motile coccus | [5] |
| Motility coefficient | μ_s | *B. subtilis* | 15–45 | μm/s | Flagellum-driven; ~25 μm/s at 30°C | [22] |
| Motility coefficient | μ_s | *C. sphaerospermum* | 0.00–0.005 | μm/s | Hyphal extension only; ~0.1–5 μm/min tip growth | [4] |
| Motility coefficient | μ_s | *S. oneidensis* | 10–40 | μm/s | Polar flagellum; comparable to *E. coli* | [16] |
| Ionizing rad. sensitivity | β_s,ion | *C. neoformans* | 1×10⁻⁵–1×10⁻⁴ | Gy⁻¹ | Melanized; LD₁₀ ~2–5 kGy | [2, 3] |
| Ionizing rad. sensitivity | β_s,ion | *D. radiodurans* | 5×10⁻⁶–5×10⁻⁵ | Gy⁻¹ | D₁₀ ~12 kGy; most radioresistant known organism | [5, 6] |
| Ionizing rad. sensitivity | β_s,ion | *B. subtilis* | 1×10⁻³–5×10⁻³ | Gy⁻¹ | Spore D₁₀ ~1–2 kGy; vegetative D₁₀ ~0.2–0.6 kGy | [23] |
| Ionizing rad. sensitivity | β_s,ion | *C. sphaerospermum* | 1×10⁻⁵–1×10⁻⁴ | Gy⁻¹ | Melanized Chernobyl isolate; comparable to *C. neoformans* | [1, 4] |
| Ionizing rad. sensitivity | β_s,ion | *A. niger* | 5×10⁻⁵–5×10⁻⁴ | Gy⁻¹ | Melanized but less resistant than Chernobyl isolates | [3] |
| Ionizing rad. sensitivity | β_s,ion | *S. oneidensis* | 5×10⁻²–1×10⁻¹ | Gy⁻¹ | D₁₀ ~0.07 kGy; extremely radiosensitive | [6] |
| Melanin production rate | α_M | *C. neoformans* | 0.05–0.15 | μg/cell/Gy | Substrate-dependent; enhanced by radiation | [2] |
| Melanin production rate | α_M | *C. sphaerospermum* | 0.08–0.20 | μg/cell/Gy | Constitutively melanized; high production | [1, 4] |
| Melanin production rate | α_M | *A. niger* | 0.03–0.10 | μg/cell/Gy | DHN-melanin pathway | [3] |
| Carrying capacity | K_s | *C. neoformans* | 10⁴–10⁶ | cells/mm³ | Yeast biofilm density under nutrient-rich conditions | [2] |
| Carrying capacity | K_s | *D. radiodurans* | 10⁵–10⁷ | cells/mm³ | Dense tetrads; high packing efficiency | [5] |
| Carrying capacity | K_s | *B. subtilis* | 10⁵–10⁷ | cells/mm³ | Dense biofilm with EPS matrix | [22] |
| Carrying capacity | K_s | *S. oneidensis* | 10⁵–10⁷ | cells/mm³ | Planktonic and biofilm modes | [16] |
| Phase-locking frequency | ω_s | All species | 0.01–1.0 | rad/hr | Circadian/ultradian metabolic oscillation; species-specific | [13, 14] |
| Noise intensity | σ_s | All species | 0.001–0.05 | — | Thermal and demographic stochasticity | [15] |

---

## 5. Computational Methods

Simulations were implemented in Julia using `DifferentialEquations.jl` for numerical integration of the PSDE system. Spatial species distributions were initialized with small random perturbations around a central point and evolved on a [0,1]² domain. The Euler-Maruyama scheme was applied for stochastic integration; Runge-Kutta (RK4) was used for deterministic validation runs. Symplectic leapfrog integration was employed for the Hamiltonian phase-space trajectories.

Optimization of species-specific parameters was performed with `JuMP.jl` and `Ipopt`. Sensitivity analysis employed Sobol indices (global variance decomposition) and Morris screening to identify the highest-influence parameters (μ_s, D_s, β_s,ion, θ_s). Visualization used `PlotlyJS.jl` for interactive 3D plots, `Agents.jl` for agent-based validation runs, and `Plots.jl` with `GR` backend for static trajectory outputs. The interactive simulation includes sliders for gamma radiation intensity, thorium intensity, and heat intensity at 45-step resolution.

---

## 6. Results

### 6.1 Species Trajectory Clustering Under Radiation Gradients

The Langevin diffusion simulations on the [0,1]² domain reveal pronounced spatial phenotype separation under the imposed radiation gradient (Figure 1). k-means clustering of species trajectories identifies three distinct spatial niches. The first cluster, concentrated in the high-radiation zone, is dominated by the melanized species *C. neoformans* and *C. sphaerospermum*, whose positive gamma sensitivity coefficients (γ_s ≈ 0.07 and 0.06 respectively) enable them to exploit ionizing radiation as a metabolic stimulus. The second cluster occupies the intermediate-radiation zone and contains *D. radiodurans* and *A. niger* — species that tolerate radiation without deriving direct metabolic benefit. The third cluster, in the low-radiation periphery, comprises *B. subtilis*, *S. oneidensis*, and the marine taxa (*Pseudoalteromonas* sp., *Polaribacter* sp., *Flavobacterium* sp.), whose higher radiation sensitivities confine them to regions where the ionizing flux has attenuated below their damage thresholds. This spatial segregation is consistent with radiation-mediated niche partitioning, in which ionizing radiation simultaneously serves as stressor and resource — a dual role with no direct analogue in classical substrate-competition biofilm models.

### 6.2 Thorium-232 Decay Integration

The Th-232 decay model (Figure 2) compares actual versus predicted radioactive decay on a logarithmic scale spanning billions of years, incorporating stochastic perturbation around the deterministic exponential decay law. The 14.05-billion-year half-life produces a quasi-static radiation source on biological timescales (hours to years), but the stochastic perturbation term captures quantum-mechanical uncertainty in individual decay events, which manifests as Poisson-distributed dose-rate fluctuations. For biofilm modeling purposes, Th-232 provides an effectively constant ionizing background, while stochastic perturbation introduces dose-rate heterogeneity at the spatial scale of microcolonies (tens of micrometers), driving the phase-space fluctuations captured by the Langevin dynamics.

### 6.3 Motility-Diffusion-Gamma Sensitivity Relationships

The three-dimensional scatter relating motility (μ_s), diffusion coefficient (D_s), and gamma sensitivity across species (Figure 3) reveals a striking inverse relationship between radiation sensitivity and motility. *C. neoformans*, at the highest gamma sensitivity (≈0.07), exhibits effectively zero motility, consistent with the hypothesis that radiotrophic organisms invest resources in melanin production rather than active locomotion. *D. radiodurans*, with the lowest gamma sensitivity, is likewise non-motile but for different mechanistic reasons: radioresistance derives from molecular protection rather than radiation exploitation. The motile species (*B. subtilis*, *S. oneidensis*) occupy the high-motility, low-gamma-sensitivity region, suggesting a trade-off between radiation exploitation and behavioral avoidance. As radiation intensity rises, the competitive advantage shifts from motile radiation-avoiders to sessile radiotrophic beneficiaries.

### 6.4 Phase-Locking and Hamiltonian kNN Dynamics

The Hamiltonian kNN decision tree (Figure 4) reveals that synchronized metabolic oscillations emerge among species occupying the same spatial niche. Within the high-radiation cluster, *C. neoformans* and *C. sphaerospermum* exhibit strong phase coherence in their melanin production cycles, suggesting cooperative dynamics in which shared melanin-derived radioprotection benefits both. Species in different radiation niches show weak phase coupling, indicating competitive exclusion at niche boundaries. The kNN decision tree correctly predicts species dominance transitions: as simulated radiation intensity increases beyond 1 kGy, *B. subtilis* populations collapse, *D. radiodurans* maintains steady-state abundance, and *C. neoformans* increases in fitness — reproducing the biologically expected ordering of radiation tolerance.

---

## 7. Discussion

### 7.1 The Hamiltonian Framework Versus Classical Biofilm Models

The phase-locked Hamiltonian framework differs fundamentally from prior reaction-diffusion biofilm models in its treatment of inter-species interactions and energy conservation. Classical models following the Wanner-Gujer tradition [8] describe biofilm dynamics through coupled PDEs governing substrate consumption and biomass growth, with species interactions mediated exclusively through competition for shared substrates. Our Hamiltonian H(p,q) explicitly conserves total system energy across kinetic, potential, and interaction terms, ensuring that numerical integration via the symplectic leapfrog scheme does not introduce artificial energy drift. The stochastic Langevin layer adds thermal and demographic noise while preserving the underlying Hamiltonian structure through the fluctuation-dissipation theorem. This dual structure captures phenomena that ODE/PDE models inherently miss: synchronized metabolic oscillations through phase-locking, energy-conserving transitions between fitness states, and the distinction between radiation as a damaging perturbation versus radiation as an energy source positively coupled to the Hamiltonian.

### 7.2 Biological Interpretation of *C. neoformans* and *D. radiodurans* Dynamics

The simulation trajectories for *C. neoformans* and *D. radiodurans* align with known biology while illuminating an underappreciated distinction between radiotrophy and radioresistance. *C. neoformans* in the model increases its fitness monotonically with radiation intensity up to a species-specific threshold, consistent with the experimental findings of Dadachova et al. [2] that melanized cells show enhanced metabolic activity and faster growth under ionizing radiation. The model captures the mechanism through the positive gamma sensitivity term γ_s, which converts absorbed dose into fitness gain via the melanin-mediated electron-transfer pathway [24, 38]. *D. radiodurans*, by contrast, maintains approximately constant fitness across a wide dose range, declining only at extreme doses. This flat fitness profile emerges from the balance between the radiation damage term β_s,ion and the organism's molecular protection coefficient, reflecting the Mn²⁺-dependent proteome shielding identified by Daly et al. [7]. The model thus correctly predicts that *C. neoformans* outcompetes *D. radiodurans* in moderately irradiated environments where radiation provides metabolic benefit, while *D. radiodurans* prevails at extreme doses where even melanin-mediated protection is overwhelmed.

### 7.3 Implications for Nuclear Bioremediation

The model's predictions for a mixed *S. oneidensis* and *O. intermedium* AM7 biofilm in a Th-contaminated environment suggest a synergistic remediation strategy. *S. oneidensis*, with its capacity for extracellular electron transfer to metal oxides [16, 17], can reduce soluble actinide species and immobilize them as insoluble precipitates. However, its extreme radiosensitivity (D₁₀ ~0.07 kGy) limits viability in high-radiation zones. The model predicts that *O. intermedium* AM7, with its thorium-tolerant EPS production [18], establishes a protective biofilm matrix in the high-radiation zone that attenuates local dose rates sufficiently to permit *S. oneidensis* colonization in the intermediate zone. This spatial partitioning — emergent from the Langevin dynamics rather than imposed as a boundary condition — suggests that self-organizing biofilm architectures may naturally optimize remediation performance. The potential for engineering *D. radiodurans* strains for combined radioresistance and metal reduction [28] offers an additional path toward single-organism bioremediation in extreme environments.

### 7.4 Mechanical Properties and Radiation Survival

The inclusion of FENE entropic spring and Kelvin-Voigt viscoelastic components addresses a dimension of biofilm physics largely absent from ecological radiation models. Biofilm EPS matrices exhibit viscoelastic behavior that mediates mechanical stress transmission, nutrient diffusion resistance, and protection from external perturbations. The FENE potential provides finite extensibility that prevents unphysical deformation under radiation-induced swelling or thermal expansion, while the Maxwell element captures stress relaxation on timescales relevant to biofilm restructuring. Species with higher EPS production rates (*O. intermedium* AM7, *B. subtilis*) generate mechanically stiffer local environments that resist radiation-induced structural degradation, creating protected microniches for radiosensitive community members.

### 7.5 Limitations and Future Directions

Several simplifying assumptions limit the current model's applicability to real contaminated sites. The radiation field is treated as spatially homogeneous with exponential attenuation, whereas actual environments at Chernobyl and Hanford exhibit complex three-dimensional dose distributions governed by heterogeneous source geometries and shielding materials. The UV oscillatory term N(t,x) = N_UV cos(ωt) assumes a purely sinusoidal diurnal cycle, neglecting atmospheric absorption spectra and seasonal variation. The model does not incorporate horizontal gene transfer, which is known to spread radiation resistance determinants among biofilm community members. Future work should couple the fitness model to experimentally measured dose-rate maps from contaminated facilities, incorporate spatially resolved nutrient transport using Darcy flow through the EPS matrix, and validate predicted species clustering patterns against metagenomic surveys of radiotrophic biofilm communities at nuclear sites.

---

## 8. Conclusion

This work introduces a Hamiltonian-Langevin framework for modeling radiotrophic fitness in multispecies biofilm communities exposed to ionizing radiation. The PSDE fitness equation integrates radiation-dependent growth, melanin-mediated energy transduction, stochastic demographic noise, and viscoelastic biofilm mechanics into a unified formalism that conserves total system energy through symplectic integration. The central result is that species niche partitioning under radiation gradients emerges naturally from the phase-locked Hamiltonian dynamics, with radiotrophic fungi (*C. neoformans*, *C. sphaerospermum*) clustering in high-radiation zones, radioresistant bacteria (*D. radiodurans*) occupying intermediate regions, and radiosensitive metal-reducing species (*S. oneidensis*) confined to attenuated peripheries. These predictions are consistent with known radiobiology and suggest design principles for engineered biofilm consortia targeting nuclear bioremediation. Experimental validation through controlled irradiation of mixed-species biofilms, combined with spatially resolved transcriptomics and metabolomics, represents the natural next step toward translating this theoretical framework into practical bioremediation strategies.

---

## References

1. Zhdanova, N.N., Zakharchenko, V.A., Vember, V.V., Nakonechnaya, L.T. (2000). Fungi from Chernobyl: mycobiota of the inner regions of the containment structures of the damaged nuclear reactor. *Mycological Research*, 104(12), 1421–1426. https://doi.org/10.1017/S0953756200002756

2. Dadachova, E., Bryan, R.A., Huang, X., Moadel, T., Schweitzer, A.D., Aisen, P., Nosanchuk, J.D., Casadevall, A. (2007). Ionizing radiation changes the electronic properties of melanin and enhances the growth of melanized fungi. *PLoS ONE*, 2(5), e457. https://doi.org/10.1371/journal.pone.0000457

3. Dadachova, E., Casadevall, A. (2008). Ionizing radiation: how fungi cope, adapt, and exploit with the help of melanin. *Current Opinion in Microbiology*, 11(6), 525–531. https://doi.org/10.1016/j.mib.2008.09.013

4. Shunk, G.K., Gomez, X.R., Kern, C., Averesch, N.J.H. (2022). Cultivation of the dematiaceous fungus *Cladosporium sphaerospermum* aboard the International Space Station and effects of ionizing radiation. *Frontiers in Microbiology*, 13, 877625. https://doi.org/10.3389/fmicb.2022.877625

5. Slade, D., Radman, M. (2011). Oxidative stress resistance in *Deinococcus radiodurans*. *Microbiology and Molecular Biology Reviews*, 75(1), 133–191. https://doi.org/10.1128/MMBR.00015-10

6. Daly, M.J. (2009). A new perspective on radiation resistance based on *Deinococcus radiodurans*. *Nature Reviews Microbiology*, 7(3), 237–245. https://doi.org/10.1038/nrmicro2073

7. Daly, M.J., et al. (2010). Small-molecule antioxidant proteome-shields in *Deinococcus radiodurans*. *PLoS ONE*, 5(9), e12570. https://doi.org/10.1371/journal.pone.0012570

8. Wanner, O., Gujer, W. (1986). A multispecies biofilm model. *Biotechnology and Bioengineering*, 28(3), 314–328. https://doi.org/10.1002/bit.260280304

9. Xavier, J.B., Picioreanu, C., van Loosdrecht, M.C.M. (2005). A framework for multidimensional modelling of activity and structure of multispecies biofilms. *Environmental Microbiology*, 7(8), 1085–1103. https://doi.org/10.1111/j.1462-2920.2005.00787.x

10. Lardon, L.A., et al. (2011). iDynoMiCS: next-generation individual-based modelling of biofilms. *Environmental Microbiology*, 13(9), 2416–2434. https://doi.org/10.1111/j.1462-2920.2011.02414.x

11. Sonner, S., Efendiev, M.A., Eberl, H.J. (2015). On the well-posedness of a mathematical model of quorum-sensing in patchy biofilm communities. *Mathematical Methods in the Applied Sciences*, 38(3), 3037–3042. https://doi.org/10.1002/mma.3237

12. Diele, F., Marangi, C., Ragni, S. (2015). Geometric numerical integration in ecological modelling. *Mathematics and Computers in Simulation*, 110, 40–52. https://doi.org/10.1016/j.matcom.2014.02.006

13. Kuramoto, Y. (1984). *Chemical Oscillations, Waves, and Turbulence*. Springer-Verlag, Berlin. https://doi.org/10.1007/978-3-642-69689-3

14. Acebrón, J.A., Bonilla, L.L., Pérez Vicente, C.J., Ritort, F., Spigler, R. (2005). The Kuramoto model: a simple paradigm for synchronization phenomena. *Reviews of Modern Physics*, 77(1), 137–185. https://doi.org/10.1103/RevModPhys.77.137

15. Campillo, F., Joannides, M., Larramendy-Valverde, I. (2017). Stochastic analysis of a full system of two competing populations in a chemostat. *Chemical Engineering Science*, 175, 424–440. https://doi.org/10.1016/j.ces.2017.10.052

16. Heidelberg, J.F., et al. (2002). Genome sequence of the dissimilatory metal ion-reducing bacterium *Shewanella oneidensis*. *Nature Biotechnology*, 20(11), 1118–1123. https://doi.org/10.1038/nbt749

17. Veeramani, H., et al. (2011). Products of abiotic U(VI) reduction by biogenic magnetite and vivianite. *Geochimica et Cosmochimica Acta*, 75(9), 2512–2528. https://doi.org/10.1016/j.gca.2011.02.024

18. Shukla, A., Parmar, P., Goswami, D., Patel, B., Saraf, M. (2020). Characterization of novel thorium tolerant *Ochrobactrum intermedium* AM7. *Journal of Hazardous Materials*, 388, 122047. https://doi.org/10.1016/j.jhazmat.2020.122047

19. Wright, S. (1932). The roles of mutation, inbreeding, crossbreeding and selection in evolution. *Proceedings of the Sixth International Congress of Genetics*, 1, 356–366.

20. Kauffman, S.A., Weinberger, E.D. (1989). The NK model of rugged fitness landscapes. *Journal of Theoretical Biology*, 141(2), 211–245. https://doi.org/10.1016/S0022-5193(89)80019-0

21. Berg, H.C. (1993). *Random Walks in Biology*. Princeton University Press.

22. Guttenplan, S.B., Shaw, S., Kearns, D.B. (2013). The cell biology of peritrichous flagella in *Bacillus subtilis*. *Molecular Microbiology*, 87(1), 211–229. https://doi.org/10.1111/mmi.12103

23. Nicholson, W.L., et al. (2000). Resistance of *Bacillus* endospores to extreme terrestrial and extraterrestrial environments. *Microbiology and Molecular Biology Reviews*, 64(3), 548–572. https://doi.org/10.1128/MMBR.64.3.548-572.2000

24. Turick, C.E., et al. (2011). Gamma radiation interacts with melanin to alter its oxidation-reduction potential and results in electric current production. *Bioelectrochemistry*, 82(1), 69–73. https://doi.org/10.1016/j.bioelechem.2011.04.009

25. Robertson, K.L., et al. (2012). Adaptation of the black yeast *Wangiella dermatitidis* to ionizing radiation. *PLoS ONE*, 7(11), e48674. https://doi.org/10.1371/journal.pone.0048674

26. Malo, M.E., et al. (2018). Morphological changes in melanized and non-melanized *Cryptococcus neoformans* cells post exposure to ionizing radiation. *Fungal Biology*, 122(6), 449–456. https://doi.org/10.1016/j.funbio.2017.08.012

27. Newsome, L., Morris, K., Lloyd, J.R. (2014). The biogeochemistry and bioremediation of uranium and other priority radionuclides. *Chemical Geology*, 363, 164–184. https://doi.org/10.1016/j.chemgeo.2013.10.034

28. Brim, H., et al. (2000). Engineering *Deinococcus radiodurans* for metal remediation in radioactive mixed waste environments. *Nature Biotechnology*, 18(1), 85–90. https://doi.org/10.1038/71986

29. Kazy, S.K., D'Souza, S.F., Sar, P. (2009). Uranium and thorium sequestration by a *Pseudomonas* sp. *Journal of Hazardous Materials*, 163(1), 65–72. https://doi.org/10.1016/j.jhazmat.2008.06.076

30. Blasius, B., Huppert, A., Stone, L. (1999). Complex dynamics and phase synchronization in spatially extended ecological systems. *Nature*, 399, 354–359. https://doi.org/10.1038/20676

31. Alpkvist, E., Klapper, I. (2007). A multidimensional multispecies continuum model for heterogeneous biofilm development. *Bulletin of Mathematical Biology*, 69(2), 765–789. https://doi.org/10.1007/s11538-006-9168-7

32. Eberl, H.J., Parker, D.F., van Loosdrecht, M.C.M. (2001). A new deterministic spatio-temporal continuum model for biofilm development. *Journal of Theoretical Medicine*, 3(3), 161–175. https://doi.org/10.1080/10273660108833072

33. Hairer, E., Lubich, C., Wanner, G. (2006). *Geometric Numerical Integration*, 2nd ed. Springer-Verlag. https://doi.org/10.1007/3-540-30666-8

34. Eisenman, H.C., Casadevall, A. (2012). Synthesis and assembly of fungal melanin. *Applied Microbiology and Biotechnology*, 93(3), 931–940. https://doi.org/10.1007/s00253-011-3777-2

35. Lloyd, J.R., Renshaw, J.C. (2005). Bioremediation of radioactive waste. *Current Opinion in Biotechnology*, 16(3), 254–260. https://doi.org/10.1016/j.copbio.2005.04.012

36. Khajo, A., et al. (2011). Protection of melanized *Cryptococcus neoformans* from lethal dose gamma irradiation. *PLoS ONE*, 6(9), e25092. https://doi.org/10.1371/journal.pone.0025092

37. Battista, J.R. (1997). Against all odds: the survival strategies of *Deinococcus radiodurans*. *Annual Review of Microbiology*, 51, 203–224. https://doi.org/10.1146/annurev.micro.51.1.203

38. Casadevall, A., et al. (2017). Melanin, radiation, and energy transduction in fungi. *Microbiology Spectrum*, 5(2), FUNK-0037-2016. https://doi.org/10.1128/microbiolspec.FUNK-0037-2016

---

## Notation Fixes Required Before Submission

1. **Species index**: Standardize to subscript *s* throughout; use *j* only for neighbor/interaction species in summations
2. **Radiation intensity**: Rename I(t,x) → I_γ(t,x) everywhere to prevent confusion with the identity matrix
3. **Noise term**: Explicitly state η_s(t,x) = σ_s ξ(t,x) where ξ is spatiotemporal white noise
4. **Hamiltonian mass**: Define m_i explicitly as effective biomass density [cells/μm³ × effective mass analog]
5. **Melanin units**: Verify dimensional consistency in ∂M/∂t = D_M ∇²M + α_M · N_RadioF · R; each term should carry [μg/μm³/s]
6. **Adaptive feedback**: Provide explicit form of A_s(t) = γ_s φ_s(t) − β_s,ion I(t) in the nutrient uptake section
7. **Section numbering**: Renumber sequentially; remove duplicate §2, §3, §8 headers from original draft
