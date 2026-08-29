#!/usr/bin/env julia
# ============================================================
#  Cellular Potts Model (CPM) for a Radiotropic Biofilm System
#  Based on: "Modeling Radioresistance and Radiotropic Fitness"
#  (Kinder & Faulkner, 2026)
#
#  PHENOTYPE, PRECISELY. This model represents RADIOTROPISM: a directional
#  spatial preference in a dose gradient, expressed as the sign and magnitude
#  of beta_ion. It does NOT represent radiotrophy (radiation measurably
#  ENHANCING growth or metabolism) and it does NOT represent radioresistance
#  (elevated SURVIVAL), because it has neither growth nor death: divide_cell!
#  has no trigger and no cell is ever removed. Radioresistance evidence
#  motivates the per-species magnitudes below; it is not what the term computes.
#
#  Maps the paper's Hamiltonian (Eq. 2, Section 3.3) to CPM energy:
#    H_CPM = H_adhesion + H_volume + H_radiation + H_pairwise + H_melanin
#
#  Coupled fields (melanin, nutrient, radiation) updated each MCS
#  per Eq. 7 (melanin RD) and Eq. 5 (radiation field).
#
#  Pure Julia — requires only: LinearAlgebra, Statistics, Random, Printf
# ============================================================

using LinearAlgebra, Statistics, Random, Printf

# ============================================================
#  1. Species definitions — from Table 2 of the paper
# ============================================================

const SPECIES_NAMES = [
    "C. neoformans",       # 1 — melanized yeast, radiotropic in this model
    "D. radiodurans",      # 2 — extreme radioresistant coccus
    "C. sphaerospermum",   # 3 — melanized filamentous fungus, radiotropic here
    "B. subtilis",         # 4 — motile rod, moderate radiosensitivity
    "A. niger",            # 5 — melanized filamentous fungus
    "S. oneidensis",       # 6 — metal-reducing bacterium, radiosensitive
    "O. intermedium AM7",  # 7 — thorium-tolerant, EPS producer
]

const N_SPECIES = 7

# Species indices for convenience
const CN = 1; const DR = 2; const CS = 3; const BS = 4
const AN = 5; const SO = 6; const OI = 7

# Species drawn TOWARD the source (negative beta_ion). A tropism, not a
# metabolism: occupying a high-dose site lowers this cell's Hamiltonian, which
# is a Metropolis acceptance functional in arbitrary units, not an energy budget.
const RADIOTROPIC = Set([CN, CS])

# Melanin-producing species (a superset: producing melanin does not imply
# being drawn toward the source)
const MELANIN_PRODUCERS = Set([CN, CS, AN])

"""
Parameters struct — values from Table 2 where available.
β_ion: TROPISM coefficient, per species — sign sets direction, magnitude sets
  strength. Negative beta means high-dose sites lower the Hamiltonian, so the
  cell drifts up the dose gradient; positive beta means it drifts down.
  The magnitudes are anchored to published radioresistance (D10) evidence, but
  the quantity the model computes is a DIRECTIONAL PREFERENCE, not survival and
  not growth. Units are nominally Gy⁻¹ from Table 2, but no dose calibration
  exists, so beta is dimensionless in practice — see the calibration contracts.
α_M: melanin production rate (μg cell⁻¹ Gy⁻¹) — Table 2
D_s: diffusion coefficient (μm² s⁻¹) → mapped to CPM motility temperature
"""
Base.@kwdef struct CPMParams
    # --- Lattice ---
    N::Int = 60                          # grid dimension (N×N×N)

    # --- CPM dynamics ---
    T_cpm::Float64 = 5.0                 # CPM Boltzmann temperature (acceptance)
    λ_V::Float64 = 10.0                  # volume constraint strength (H_volume)
    V_target::Int = 120                  # target volume per cell (sites)

    # --- Radiation field (Eq. 5: I_γ = I₀ exp(-κ‖x‖)) ---
    I0::Float64 = 1.0                    # source intensity at axis
    κ::Float64 = 2.0                     # attenuation coefficient

    # --- Melanin diffusion (Eq. 7) ---
    D_M::Float64 = 0.1                   # melanin diffusion coefficient
    dt_field::Float64 = 0.5              # time step for field updates

    # --- Nutrient field ---
    D_C::Float64 = 0.2                   # nutrient diffusion coefficient
    C_wall::Float64 = 1.0                # nutrient concentration at wall boundary

    # --- Species-specific β_ion (Table 2, midpoint of range) ---
    # Sign convention: positive = drifts DOWN the dose gradient (away from the
    # source), negative = drifts UP it (toward the source). Not damage, not
    # benefit: the model has no death and no growth to express either in.
    β_ion::Vector{Float64} = [
        -5e-5,   # CN: drawn toward the source. CONTESTED — see
                 #     docs/research/radiotrophic_compatibility_audit.md;
                 #     radiotrophy is not established for this species.
         2.5e-5, # DR: radioresistant evidence -> weak avoidance here
        -5e-5,   # CS: drawn toward the source. CONTESTED on the same grounds;
                 #     Zhdanova bears on radioTROPISM, which is this term.
         3e-3,   # BS: moderate radiosensitivity (Nicholson 2000)
         2.5e-4, # AN: melanized, drifts away from the source
         7.5e-2, # SO: extremely radiosensitive (Daly 2009)
         1e-2,   # OI: moderate (estimated, not in Table 2)
    ]

    # --- Species-specific melanin rate α_M (Table 2, midpoint) ---
    α_M_species::Vector{Float64} = [
        0.10,  # CN
        0.0,   # DR (no melanin)
        0.14,  # CS
        0.0,   # BS
        0.065, # AN
        0.0,   # SO
        0.0,   # OI
    ]

    # --- Species diffusion D_s (Table 2, midpoint) → used for motility bias ---
    D_s::Vector{Float64} = [
        0.05,  # CN
        0.25,  # DR
        0.025, # CS
        0.50,  # BS
        0.025, # AN
        0.40,  # SO
        0.20,  # OI (estimated)
    ]

    # --- Nutrient uptake rates per species ---
    uptake::Vector{Float64} = [
        0.01, 0.02, 0.01, 0.03, 0.01, 0.03, 0.02
    ]

    # --- Pairwise mutualistic interaction (Eq. 3: V_ij = -γ exp(-r²/σ²)) ---
    γ_mutual::Float64 = 2.0              # attraction strength
    σ_mutual::Float64 = 10.0             # interaction range (lattice units)

    # --- Number of cells per species at initialization ---
    n_cells_per_species::Int = 8

    # --- Physical time contract (REQUIRED for any physical coupling).
    #     NaN = undefined: dose accrual and physical-time tracking refuse to
    #     run. Never default this to a number — see audit §5. ---
    seconds_per_mcs::Float64 = NaN

    # --- Snapshot interval ---
    snapshot_interval::Int = 10
end

# ============================================================
#  2. Adhesion energy matrix J — derived from paper's V_ij
#     Lower J = more adhesive (mutualistic); higher J = repulsive
#     Index 0 → medium (row/col 1), species 1..7 → rows/cols 2..8
# ============================================================

"""
Build the (N_SPECIES+1) × (N_SPECIES+1) adhesion matrix.
Index 1 = medium, indices 2..8 = species 1..7.

Mutualistic pairs (lower J): CN↔DR, CN↔CS, CS↔AN, SO↔OI
Antagonistic/neutral pairs get higher J.
Species–medium contact has moderate J (cells prefer to stick together).
"""
function build_J_matrix()
    # Default: moderate adhesion
    J = fill(12.0, N_SPECIES + 1, N_SPECIES + 1)

    # Self-adhesion (same species) — low energy, cells stick to own type
    for i in 1:(N_SPECIES + 1)
        J[i, i] = 0.0
    end

    # Species–medium (index 1 = medium)
    for s in 2:(N_SPECIES + 1)
        J[1, s] = 16.0
        J[s, 1] = 16.0
    end

    # Mutualistic pairs — lower J (attractive, from paper's V_ij Eq. 3)
    mutualistic = [(CN, DR), (CN, CS), (CS, AN), (SO, OI), (DR, BS)]
    for (a, b) in mutualistic
        J[a + 1, b + 1] = 4.0
        J[b + 1, a + 1] = 4.0
    end

    # Antagonistic pairs — higher J
    antagonistic = [(CN, SO), (CS, SO), (BS, AN)]
    for (a, b) in antagonistic
        J[a + 1, b + 1] = 20.0
        J[b + 1, a + 1] = 20.0
    end

    return J
end

# ============================================================
#  3. Data structures
# ============================================================

mutable struct CellInfo
    species::Int          # 1..7
    volume::Int           # number of occupied lattice sites
    com::Vector{Float64}  # center of mass [x, y, z]
    # --- genealogy labels (biological identity; NEVER a material identity) ---
    lineage_id::Int       # founder's cell id, inherited through divisions
    parent_id::Int        # 0 for founders
    generation::Int       # 0 for founders
    birth_mcs::Int
    lifetime_dose_Gy::Float64   # accrued while alive; finalized at division/death
end

# Legacy 3-arg form (founder-agnostic; genealogy fields zeroed)
CellInfo(species::Int, volume::Int, com::Vector{Float64}) =
    CellInfo(species, volume, com, 0, 0, 0, 0, 0.0)

# ---- Lifecycle event records (typed, archival) ----------------
# Division is binary fission: the parent identity ENDS and two daughters
# with fresh IDs are created. No automatic scheduler exists — see §12b.

struct DivisionEvent
    mcs::Int
    physical_time_s::Float64      # NaN until seconds_per_mcs is configured
    parent_id::Int
    daughter_a_id::Int
    daughter_b_id::Int
    lineage_id::Int
    parent_generation::Int
    daughter_generation::Int
    parent_volume::Int
    daughter_a_volume::Int
    daughter_b_volume::Int
end

struct DeathEvent
    mcs::Int
    physical_time_s::Float64
    cell_id::Int
    species::Int
    lineage_id::Int
    generation::Int
    lifetime_dose_Gy::Float64
end

struct ArchivedCell
    id::Int
    species::Int
    lineage_id::Int
    parent_id::Int
    generation::Int
    birth_mcs::Int
    end_mcs::Int
    fate::Symbol                  # :divided | :died
    lifetime_dose_Gy::Float64
end

mutable struct CPMState
    # Lattice: 0 = medium, -1 = wall, >0 = cell ID
    lattice::Array{Int32, 3}

    # Cell registry: cell_id → CellInfo
    cells::Dict{Int, CellInfo}

    # Coupled fields (same dimensions as lattice)
    radiation::Array{Float64, 3}  # I_γ(x,y,z) — Eq. 5
    melanin::Array{Float64, 3}    # M(x,y,z)   — Eq. 7
    nutrient::Array{Float64, 3}   # C(x,y,z)

    # Cylinder mask: true = interior (updatable), false = wall
    interior::BitArray{3}

    # Next cell ID counter
    next_id::Int

    # Parameters
    params::CPMParams

    # Adhesion matrix
    J::Matrix{Float64}

    # --- Simulation clock (single time authority; see audit §5) ---
    current_mcs::Int              # incremented by mcs_step!
    physical_time_s::Float64      # NaN until seconds_per_mcs is configured

    # --- Lifecycle archive ---
    division_log::Vector{DivisionEvent}
    death_log::Vector{DeathEvent}
    archive::Vector{ArchivedCell}

    # --- Physical dose contract (Gy); inactive until a dose field is imported.
    #     `radiation` above remains the dimensionless Hamiltonian signal;
    #     physical dose is NEVER stored into it implicitly. ---
    dose_active::Bool
    dose_rate_mean_Gy_s::Array{Float64, 3}
    dose_rate_sd_Gy_s::Array{Float64, 3}
    dose_increment_Gy::Array{Float64, 3}
    accumulated_dose_Gy::Array{Float64, 3}
    membrane_dose_rate_Gy_s::Float64   # NaN until the m-vs-P_eff coupling is
                                       # declared (audit §6 STOP CONDITION)

    # --- Per-consumer signal split: melanin production reads its own drive
    #     field (initialized identical to `radiation`, so legacy behavior is
    #     bit-for-bit unchanged until a transform is explicitly imported) ---
    melanin_drive::Array{Float64, 3}
end

# ============================================================
#  4. Initialization
# ============================================================

"""
Compute cylindrical distance from the z-axis (center of the lattice).
"""
@inline function radial_dist(x, y, N)
    cx = N / 2.0
    cy = N / 2.0
    return sqrt((x - cx)^2 + (y - cy)^2)
end

"""
Initialize the radiation field: I_γ(r) = I₀ exp(-κ r/N).
Radiation source along z-axis (center), attenuated radially (Eq. 5).
"""
function init_radiation!(state::CPMState)
    p = state.params
    N = p.N
    @inbounds for z in 1:N, y in 1:N, x in 1:N
        r = radial_dist(x, y, N)
        state.radiation[x, y, z] = p.I0 * exp(-p.κ * r / N)
    end
end

"""
Initialize nutrient field: gradient from cylinder wall (high) to center (low).
"""
function init_nutrient!(state::CPMState)
    p = state.params
    N = p.N
    R = N / 2.0
    @inbounds for z in 1:N, y in 1:N, x in 1:N
        r = radial_dist(x, y, N)
        # Nutrient supplied from wall, linear gradient
        state.nutrient[x, y, z] = p.C_wall * clamp(r / R, 0.0, 1.0)
    end
end

"""
Place a single cell as a small sphere of random radius 2–3 at a given center.
Returns number of sites placed.
"""
function place_cell!(lattice::Array{Int32, 3}, interior::BitArray{3},
                     cell_id::Int32, cx::Int, cy::Int, cz::Int,
                     radius::Int, N::Int)
    count = 0
    for dz in -radius:radius, dy in -radius:radius, dx in -radius:radius
        if dx^2 + dy^2 + dz^2 <= radius^2
            x = cx + dx; y = cy + dy; z = cz + dz
            if 1 <= x <= N && 1 <= y <= N && 1 <= z <= N
                if interior[x, y, z] && lattice[x, y, z] == 0
                    lattice[x, y, z] = cell_id
                    count += 1
                end
            end
        end
    end
    return count
end

"""
Initialize the full CPM state: lattice, cells, fields.
"""
function init_state(params::CPMParams; seed::Int = 42)
    rng = MersenneTwister(seed)
    N = params.N

    # Create lattice (all medium = 0)
    lattice = zeros(Int32, N, N, N)

    # Create interior mask (cylindrical geometry)
    interior = falses(N, N, N)
    R = N / 2.0
    @inbounds for z in 1:N, y in 1:N, x in 1:N
        if radial_dist(x, y, N) <= R
            interior[x, y, z] = true
        else
            lattice[x, y, z] = Int32(-1)  # wall
        end
    end

    # Place cells
    cells = Dict{Int, CellInfo}()
    next_id = 1

    for s in 1:N_SPECIES
        for _ in 1:params.n_cells_per_species
            # Random position inside cylinder
            placed = false
            attempts = 0
            while !placed && attempts < 100
                cx = rand(rng, 5:(N - 4))
                cy = rand(rng, 5:(N - 4))
                cz = rand(rng, 5:(N - 4))
                if radial_dist(cx, cy, N) < R - 4
                    radius = rand(rng, 2:3)
                    cid = Int32(next_id)
                    vol = place_cell!(lattice, interior, cid, cx, cy, cz, radius, N)
                    if vol > 0
                        # founder: lineage = own id, generation 0
                        cells[next_id] = CellInfo(s, vol,
                            Float64[cx, cy, cz], next_id, 0, 0, 0, 0.0)
                        next_id += 1
                        placed = true
                    end
                end
                attempts += 1
            end
        end
    end

    # Fields
    radiation = zeros(Float64, N, N, N)
    melanin = zeros(Float64, N, N, N)
    nutrient = zeros(Float64, N, N, N)

    J = build_J_matrix()

    state = CPMState(lattice, cells, radiation, melanin, nutrient,
                     interior, next_id, params, J,
                     0, NaN,                                 # clock
                     DivisionEvent[], DeathEvent[], ArchivedCell[],
                     false,                                  # dose_active
                     zeros(Float64, N, N, N), zeros(Float64, N, N, N),
                     zeros(Float64, N, N, N), zeros(Float64, N, N, N),
                     NaN,                                    # membrane dose rate
                     zeros(Float64, N, N, N))                # melanin_drive

    init_radiation!(state)
    init_nutrient!(state)
    copyto!(state.melanin_drive, state.radiation)  # legacy: identical signals

    return state
end

# ============================================================
#  5. Hamiltonian terms — energy computation
# ============================================================

# 26-connected neighborhood offsets in 3D
const NEIGHBORS_26 = Tuple{Int,Int,Int}[
    (dx, dy, dz)
    for dx in -1:1 for dy in -1:1 for dz in -1:1
    if !(dx == 0 && dy == 0 && dz == 0)
]

# 6-connected (face) neighbors for field diffusion
const NEIGHBORS_6 = [
    (1,0,0), (-1,0,0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1)
]

"""
Get the species index (1..7) for a lattice value, or 0 for medium.
Wall (-1) is treated as medium for adhesion purposes.
"""
@inline function species_of(σ::Int32, cells::Dict{Int, CellInfo})
    σ <= 0 && return 0
    c = get(cells, Int(σ), nothing)
    c === nothing && return 0
    return c.species
end

"""
J matrix index: 1 = medium, 2..8 = species 1..7.
"""
@inline j_idx(sp::Int) = sp + 1

"""
Compute the adhesion energy contribution for a single site.
H_adhesion = Σ_{neighbors} J[σ(x), σ(x')] * (1 - δ(σ(x), σ(x')))
(Section 3.3 → mapped to CPM contact energies)
"""
function site_adhesion(lattice, cells, J, x, y, z, N)
    σ_center = lattice[x, y, z]
    sp_center = species_of(σ_center, cells)
    E = 0.0
    @inbounds for (dx, dy, dz) in NEIGHBORS_26
        nx, ny, nz = x + dx, y + dy, z + dz
        if 1 <= nx <= N && 1 <= ny <= N && 1 <= nz <= N
            σ_nb = lattice[nx, ny, nz]
            if σ_nb != σ_center  # (1 - δ) factor
                sp_nb = species_of(σ_nb, cells)
                E += J[j_idx(sp_center), j_idx(sp_nb)]
            end
        else
            # boundary treated as medium
            E += J[j_idx(sp_center), 1]
        end
    end
    return E
end

"""
    compute_delta_H_terms(state, sx, sy, sz, tx, ty, tz)

The four Hamiltonian contributions to a copy attempt, separately.

THIS FUNCTION IS THE OLD BODY OF `compute_delta_H`, WHICH ALREADY COMPUTED ALL
FOUR AND THEN DISCARDED THEM on its return line. Nothing new is calculated here.
`compute_delta_H` now sums this in the same order it used to add them, so the
scalar it returns is bit-identical by construction rather than by luck — same
operands, same associativity, same result.

It is split rather than retyped because three callers outside this file want the
scalar, and one of them is `biofilms_potts_jacc.jl`'s cross-implementation check
against the parallel port. Changing what that reads to prove a rendering point
would be trading a real guard for a picture.
"""
function compute_delta_H_terms(state::CPMState, sx::Int, sy::Int, sz::Int,
                               tx::Int, ty::Int, tz::Int)
    lat = state.lattice
    cells = state.cells
    J = state.J
    p = state.params
    N = p.N

    σ_source = lat[sx, sy, sz]  # cell ID being copied TO target
    σ_target = lat[tx, ty, tz]  # cell ID being replaced AT target

    # ---- H_adhesion: ΔH_adh ----
    # Before: target site has σ_target; After: target site has σ_source
    H_adh_before = site_adhesion(lat, cells, J, tx, ty, tz, N)
    # Temporarily set target to source value
    lat[tx, ty, tz] = σ_source
    H_adh_after = site_adhesion(lat, cells, J, tx, ty, tz, N)
    lat[tx, ty, tz] = σ_target  # restore
    ΔH_adh = H_adh_after - H_adh_before

    # ---- H_volume: λ_V * (V - V_target)² ----
    # Volume constraint (Section 3.3 kinetic term analog)
    ΔH_vol = 0.0
    if σ_source > 0
        c = cells[Int(σ_source)]
        V = c.volume
        # After copy: source cell gains 1 site
        ΔH_vol += p.λ_V * ((V + 1 - p.V_target)^2 - (V - p.V_target)^2)
    end
    if σ_target > 0
        c = cells[Int(σ_target)]
        V = c.volume
        # After copy: target cell loses 1 site
        ΔH_vol += p.λ_V * ((V - 1 - p.V_target)^2 - (V - p.V_target)^2)
    end

    # ---- H_radiation: β_{s,ion} * I_γ(r) per site ----
    # (Eq. 5 + Table 2: negative β = drawn up the dose gradient)
    ΔH_rad = 0.0
    I_local = state.radiation[tx, ty, tz]
    if σ_source > 0
        sp = cells[Int(σ_source)].species
        ΔH_rad += p.β_ion[sp] * I_local  # source gains a site in radiation
    end
    if σ_target > 0
        sp = cells[Int(σ_target)].species
        ΔH_rad -= p.β_ion[sp] * I_local  # target loses a site in radiation
    end

    # ---- H_melanin: α_M * M(x) modulates local energy ----
    # (Eq. 7: melanin-rich sites are preferred by the radiotropic species —
    #  again a spatial preference, not energy transduction)
    ΔH_mel = 0.0
    M_local = state.melanin[tx, ty, tz]
    # THE 0.5 BELOW IS THE COEFFICIENT THAT ACTUALLY MOVES THIS MODEL, and it is
    # hard-coded here rather than tabulated. At the shipped I0 = 1.0 and
    # T_cpm = 5.0, the radiation term for a radiotropic species is
    # β_ion·I = -5e-5, an acceptance bias of 1.000010 — one part in 1e5. This
    # term at the reported M = 1.44 is -0.72, a bias of 1.155. Four orders of
    # magnitude. The radial stratification is therefore melanin-mediated, not
    # β_ion-mediated. Radiation still drives it, but only indirectly:
    # melanin_drive is copied from the radiation field, so
    #   radiation -> production (α_M, tabulated) -> M -> here -> tropism.
    # Ledgered as cpm.melanin_coupling in data/parameter_provenance.csv.
    if σ_source > 0 && cells[Int(σ_source)].species in RADIOTROPIC
        ΔH_mel -= 0.5 * M_local  # melanin-rich site preferred
    end
    if σ_target > 0 && cells[Int(σ_target)].species in RADIOTROPIC
        ΔH_mel += 0.5 * M_local  # losing a melanin-rich site costs
    end

    return (adh = ΔH_adh, vol = ΔH_vol, rad = ΔH_rad, mel = ΔH_mel)
end

"""
Compute ΔH for a proposed copy: source site (sx,sy,sz) copies into target (tx,ty,tz).
Only recomputes terms affected by the single-site change.
Returns the total energy change ΔH.

Unchanged in value, type and summation order; it now adds up
`compute_delta_H_terms` instead of four locals with the same names.
"""
function compute_delta_H(state::CPMState, sx::Int, sy::Int, sz::Int,
                         tx::Int, ty::Int, tz::Int)
    t = compute_delta_H_terms(state, sx, sy, sz, tx, ty, tz)
    return t.adh + t.vol + t.rad + t.mel
end

# ============================================================
#  6. Monte Carlo Step (MCS)
# ============================================================

"""
Perform one Monte Carlo Step = N³ attempted copy operations.
Standard CPM Metropolis dynamics (Section 3.4 analog → stochastic accept/reject).
"""
function mcs_step!(state::CPMState, rng::AbstractRNG)
    p = state.params
    N = p.N
    lat = state.lattice
    n_attempts = N^3
    state.current_mcs += 1

    for _ in 1:n_attempts
        # 1. Pick random source site inside cylinder
        sx = rand(rng, 1:N)
        sy = rand(rng, 1:N)
        sz = rand(rng, 1:N)
        if !state.interior[sx, sy, sz]
            continue
        end

        # 2. Pick random neighbor as target
        dx, dy, dz = NEIGHBORS_26[rand(rng, 1:26)]
        tx = sx + dx; ty = sy + dy; tz = sz + dz
        if !(1 <= tx <= N && 1 <= ty <= N && 1 <= tz <= N)
            continue
        end
        if !state.interior[tx, ty, tz]
            continue
        end

        # 3. Skip if same cell ID
        σ_s = lat[sx, sy, sz]
        σ_t = lat[tx, ty, tz]
        if σ_s == σ_t
            continue
        end
        # Don't copy medium into medium
        if σ_s <= 0 && σ_t <= 0
            continue
        end

        # 4. Compute ΔH
        ΔH = compute_delta_H(state, sx, sy, sz, tx, ty, tz)

        # 5. Metropolis acceptance
        accept = if ΔH <= 0
            true
        else
            rand(rng) < exp(-ΔH / p.T_cpm)
        end

        if accept
            # Execute copy: target site gets source cell ID
            # Update volumes
            if σ_t > 0 && haskey(state.cells, Int(σ_t))
                state.cells[Int(σ_t)].volume -= 1
            end
            if σ_s > 0 && haskey(state.cells, Int(σ_s))
                state.cells[Int(σ_s)].volume += 1
            end
            lat[tx, ty, tz] = σ_s
        end
    end

    # Clean up dead cells (volume = 0) — archived, never silently dropped
    for (id, cell) in state.cells
        if cell.volume <= 0
            push!(state.death_log,
                  DeathEvent(state.current_mcs, state.physical_time_s, id,
                             cell.species, cell.lineage_id, cell.generation,
                             cell.lifetime_dose_Gy))
            push!(state.archive,
                  ArchivedCell(id, cell.species, cell.lineage_id,
                               cell.parent_id, cell.generation, cell.birth_mcs,
                               state.current_mcs, :died, cell.lifetime_dose_Gy))
            delete!(state.cells, id)
        end
    end
end

# ============================================================
#  7. Coupled field updates (each MCS)
# ============================================================

"""
Discrete Laplacian (6-connected) at site (x,y,z).
"""
function laplacian_3d(field::Array{Float64, 3}, x, y, z, N)
    val = field[x, y, z]
    L = 0.0
    @inbounds for (dx, dy, dz) in NEIGHBORS_6
        nx, ny, nz = x + dx, y + dy, z + dz
        if 1 <= nx <= N && 1 <= ny <= N && 1 <= nz <= N
            L += field[nx, ny, nz] - val
        end
        # Neumann BC: if outside, contributes 0
    end
    return L
end

"""
Count cells of radiotropic species at site (x,y,z).
Returns 1.0 if a radiotropic cell occupies the site, else 0.0.
"""
@inline function rf_density(lattice, cells, x, y, z)
    σ = lattice[x, y, z]
    σ <= 0 && return 0.0
    c = get(cells, Int(σ), nothing)
    c === nothing && return 0.0
    return c.species in MELANIN_PRODUCERS ? 1.0 : 0.0
end

"""
Count cells of any species at site for nutrient consumption.
Returns uptake rate if cell present, else 0.
"""
@inline function local_uptake(lattice, cells, uptake, x, y, z)
    σ = lattice[x, y, z]
    σ <= 0 && return 0.0
    c = get(cells, Int(σ), nothing)
    c === nothing && return 0.0
    return uptake[c.species]
end

"""
Update melanin field M(x,y,z) via reaction-diffusion (Eq. 7):
  ∂M/∂t = D_M ∇²M + α_M · n_RF(x) · I_γ(x)
"""
function update_melanin!(state::CPMState)
    p = state.params
    N = p.N
    dt = p.dt_field
    M = state.melanin
    M_new = copy(M)

    @inbounds for z in 1:N, y in 1:N, x in 1:N
        if !state.interior[x, y, z]
            continue
        end
        # Diffusion
        diff = p.D_M * laplacian_3d(M, x, y, z, N)
        # Production: melanin-producing cells in radiation field
        n_rf = rf_density(state.lattice, state.cells, x, y, z)
        σ = state.lattice[x, y, z]
        α = 0.0
        if σ > 0 && haskey(state.cells, Int(σ))
            sp = state.cells[Int(σ)].species
            α = p.α_M_species[sp]
        end
        prod = α * n_rf * state.melanin_drive[x, y, z]
        M_new[x, y, z] = max(0.0, M[x, y, z] + dt * (diff + prod))
    end

    copyto!(state.melanin, M_new)
end

"""
Update nutrient field C(x,y,z):
  ∂C/∂t = D_C ∇²C - Σ_s uptake_s · n_s(x)
Boundary: C = C_wall at wall sites.
"""
function update_nutrient!(state::CPMState)
    p = state.params
    N = p.N
    dt = p.dt_field
    C = state.nutrient
    C_new = copy(C)

    @inbounds for z in 1:N, y in 1:N, x in 1:N
        if !state.interior[x, y, z]
            # Wall sites maintain nutrient supply
            C_new[x, y, z] = p.C_wall
            continue
        end
        diff = p.D_C * laplacian_3d(C, x, y, z, N)
        consumption = local_uptake(state.lattice, state.cells, p.uptake, x, y, z)
        C_new[x, y, z] = max(0.0, C[x, y, z] + dt * (diff - consumption))
    end

    copyto!(state.nutrient, C_new)
end

"""
Recompute center of mass for all living cells.
"""
function update_centers_of_mass!(state::CPMState)
    N = state.params.N
    lat = state.lattice

    # Reset accumulators
    sums = Dict{Int, Vector{Float64}}()
    counts = Dict{Int, Int}()
    for id in keys(state.cells)
        sums[id] = [0.0, 0.0, 0.0]
        counts[id] = 0
    end

    @inbounds for z in 1:N, y in 1:N, x in 1:N
        σ = lat[x, y, z]
        if σ > 0 && haskey(sums, Int(σ))
            id = Int(σ)
            sums[id][1] += x
            sums[id][2] += y
            sums[id][3] += z
            counts[id] += 1
        end
    end

    for (id, cell) in state.cells
        n = counts[id]
        if n > 0
            cell.com .= sums[id] ./ n
            cell.volume = n
        end
    end
end

# ============================================================
#  8. Pairwise interaction energy (Eq. 3)
#     V_ij = -γ exp(-r_ij²/σ²) for mutualistic pairs
#     Computed between cell centers, not per-site.
# ============================================================

# Mutualistic species pairs (from biological context in paper)
const MUTUALISTIC_PAIRS = Set([
    (CN, DR), (DR, CN),
    (CN, CS), (CS, CN),
    (CS, AN), (AN, CS),
    (SO, OI), (OI, SO),
    (DR, BS), (BS, DR),
])

"""
Total pairwise interaction energy between all cell pairs.
V_ij^mutual = -γ exp(-r_ij²/σ²)  (Eq. 3)
"""
function total_pairwise_energy(state::CPMState)
    p = state.params
    E = 0.0
    cell_list = collect(state.cells)
    n = length(cell_list)
    for i in 1:n
        id_i, ci = cell_list[i]
        for j in (i + 1):n
            id_j, cj = cell_list[j]
            if (ci.species, cj.species) in MUTUALISTIC_PAIRS
                r2 = sum((ci.com .- cj.com) .^ 2)
                E += -p.γ_mutual * exp(-r2 / p.σ_mutual^2)
            end
        end
    end
    return E
end

# ============================================================
#  9. Snapshot and trajectory recording
# ============================================================

struct SpeciesSnapshot
    species::Int
    name::String
    total_volume::Int
    n_cells::Int
    mean_r::Float64      # mean radial distance from z-axis
    mean_melanin::Float64 # mean melanin at cell sites
end

struct Snapshot
    mcs::Int
    species_data::Vector{SpeciesSnapshot}
    total_energy_pairwise::Float64
end

"""
Collect a snapshot of the current state.
"""
function take_snapshot(state::CPMState, mcs::Int)
    p = state.params
    N = p.N

    # Per-species accumulators
    vol = zeros(Int, N_SPECIES)
    ncells = zeros(Int, N_SPECIES)
    r_sum = zeros(Float64, N_SPECIES)
    r_count = zeros(Int, N_SPECIES)
    mel_sum = zeros(Float64, N_SPECIES)

    for (id, cell) in state.cells
        s = cell.species
        ncells[s] += 1
        vol[s] += cell.volume
        r = radial_dist(cell.com[1], cell.com[2], N)
        r_sum[s] += r
        r_count[s] += 1
    end

    # Mean melanin per species from lattice
    lat = state.lattice
    species_site_count = zeros(Int, N_SPECIES)
    @inbounds for z in 1:N, y in 1:N, x in 1:N
        σ = lat[x, y, z]
        if σ > 0 && haskey(state.cells, Int(σ))
            s = state.cells[Int(σ)].species
            mel_sum[s] += state.melanin[x, y, z]
            species_site_count[s] += 1
        end
    end

    species_data = SpeciesSnapshot[]
    for s in 1:N_SPECIES
        mean_r = r_count[s] > 0 ? r_sum[s] / r_count[s] : 0.0
        mean_mel = species_site_count[s] > 0 ? mel_sum[s] / species_site_count[s] : 0.0
        push!(species_data, SpeciesSnapshot(s, SPECIES_NAMES[s], vol[s], ncells[s],
                                             mean_r, mean_mel))
    end

    E_pair = total_pairwise_energy(state)
    return Snapshot(mcs, species_data, E_pair)
end

"""
Print a snapshot as a formatted table.
"""
function print_snapshot(snap::Snapshot)
    @printf("╔══════════════════════════════════════════════════════════════════════════╗\n")
    @printf("║  MCS = %-6d                      Pairwise E = %10.3f              ║\n",
            snap.mcs, snap.total_energy_pairwise)
    @printf("╠══════════════════════════════════════════════════════════════════════════╣\n")
    @printf("║  %-22s  %6s  %5s  %8s  %10s  ║\n",
            "Species", "Volume", "Cells", "Mean r", "Melanin")
    @printf("╠══════════════════════════════════════════════════════════════════════════╣\n")
    for sd in snap.species_data
        @printf("║  %-22s  %6d  %5d  %8.2f  %10.4f  ║\n",
                sd.name, sd.total_volume, sd.n_cells, sd.mean_r, sd.mean_melanin)
    end
    @printf("╚══════════════════════════════════════════════════════════════════════════╝\n")
end

# ============================================================
#  10. K-means clustering of cell centers (k=3)
#      Identifies spatial niches in (x,y,z) space
# ============================================================

"""
Simple k-means clustering of cell centers into k clusters.
Returns (assignments, centroids) where assignments maps cell_id → cluster.
"""
function kmeans_cells(state::CPMState; k::Int = 3, max_iter::Int = 50,
                       seed::Int = 123)
    rng = MersenneTwister(seed)
    cell_ids = collect(keys(state.cells))
    n = length(cell_ids)

    if n == 0
        return Dict{Int,Int}(), Matrix{Float64}(undef, 0, 0)
    end

    # Collect cell centers as matrix (n × 3)
    points = zeros(Float64, n, 3)
    for (i, id) in enumerate(cell_ids)
        points[i, :] .= state.cells[id].com
    end

    # Initialize centroids randomly from data points
    idx = randperm(rng, n)[1:min(k, n)]
    centroids = points[idx, :]
    actual_k = size(centroids, 1)

    assignments = zeros(Int, n)

    for _ in 1:max_iter
        # Assign each point to nearest centroid
        changed = false
        for i in 1:n
            best_c = 1
            best_d = Inf
            for c in 1:actual_k
                d = sum((points[i, :] .- centroids[c, :]) .^ 2)
                if d < best_d
                    best_d = d
                    best_c = c
                end
            end
            if assignments[i] != best_c
                assignments[i] = best_c
                changed = true
            end
        end

        !changed && break

        # Update centroids
        for c in 1:actual_k
            members = findall(==(c), assignments)
            if !isempty(members)
                centroids[c, :] .= mean(points[members, :], dims=1)[1, :]
            end
        end
    end

    # Build result dict
    result = Dict{Int, Int}()
    for (i, id) in enumerate(cell_ids)
        result[id] = assignments[i]
    end

    return result, centroids
end

"""
Print cluster analysis with species composition per cluster.
"""
function print_cluster_analysis(state::CPMState, assignments::Dict{Int, Int},
                                 centroids::Matrix{Float64})
    k = size(centroids, 1)
    N = state.params.N

    println("\n" * "="^72)
    println("  SPATIAL NICHE ANALYSIS (k-means, k=$k)")
    println("="^72)

    for c in 1:k
        members = [id for (id, cl) in assignments if cl == c]
        if isempty(members)
            continue
        end

        # Cluster centroid
        cx, cy, cz = centroids[c, :]
        r = radial_dist(cx, cy, N)

        @printf("\n  Cluster %d: centroid=(%.1f, %.1f, %.1f), radial=%.1f\n",
                c, cx, cy, cz, r)
        @printf("  %-22s  %6s  %8s\n", "Species", "Cells", "Volume")
        println("  " * "-"^40)

        # Species breakdown
        sp_count = zeros(Int, N_SPECIES)
        sp_vol = zeros(Int, N_SPECIES)
        for id in members
            cell = state.cells[id]
            sp_count[cell.species] += 1
            sp_vol[cell.species] += cell.volume
        end
        for s in 1:N_SPECIES
            if sp_count[s] > 0
                @printf("  %-22s  %6d  %8d\n",
                        SPECIES_NAMES[s], sp_count[s], sp_vol[s])
            end
        end
    end
    println("\n" * "="^72)
end

# ============================================================
#  11. Main simulation driver
# ============================================================

"""
    run_simulation(params, n_mcs; seed=42)

Run the CPM simulation for `n_mcs` Monte Carlo Steps.
Returns (state, trajectory) where trajectory is a Vector{Snapshot}.
"""
function run_simulation(params::CPMParams, n_mcs::Int; seed::Int = 42)
    rng = MersenneTwister(seed)
    state = init_state(params; seed=seed)

    trajectory = Snapshot[]

    # Initial snapshot
    update_centers_of_mass!(state)
    snap = take_snapshot(state, 0)
    push!(trajectory, snap)
    print_snapshot(snap)

    for mcs in 1:n_mcs
        # --- CPM Monte Carlo step ---
        mcs_step!(state, rng)

        # --- Update coupled fields ---
        update_melanin!(state)
        update_nutrient!(state)

        # --- Update cell statistics ---
        update_centers_of_mass!(state)

        # --- Snapshot ---
        if mcs % params.snapshot_interval == 0 || mcs == n_mcs
            snap = take_snapshot(state, mcs)
            push!(trajectory, snap)
            print_snapshot(snap)
        end
    end

    return state, trajectory
end

"""
    main()

Run 200 MCS with default parameters and print cluster analysis.
"""
function main()
    println("="^72)
    println("  Cellular Potts Model — Radiotropic Biofilm System")
    println("  Based on Kinder & Faulkner (2026)")
    println("  Hamiltonian: H = H_adh + H_vol + H_rad + H_mel  (H_pair is a diagnostic, not in the acceptance path)")
    println("="^72)
    println()

    params = CPMParams()

    @printf("  Lattice:  %d x %d x %d (cylindrical)\n", params.N, params.N, params.N)
    @printf("  Species:  %d\n", N_SPECIES)
    @printf("  Cells:    %d per species (%d total)\n",
            params.n_cells_per_species, params.n_cells_per_species * N_SPECIES)
    @printf("  T_cpm:    %.1f\n", params.T_cpm)
    @printf("  lambda_V: %.1f  (V_target = %d)\n", params.λ_V, params.V_target)
    @printf("  kappa:    %.1f  (radiation attenuation)\n", params.κ)
    @printf("  D_M:      %.2f  (melanin diffusion)\n", params.D_M)
    @printf("  D_C:      %.2f  (nutrient diffusion)\n", params.D_C)
    println()

    n_mcs = 200
    @printf("  Running %d Monte Carlo Steps...\n\n", n_mcs)

    t_start = time()
    state, trajectory = run_simulation(params, n_mcs)
    elapsed = time() - t_start

    @printf("\n  Simulation completed in %.1f seconds.\n", elapsed)
    @printf("  Surviving cells: %d\n", length(state.cells))

    # --- Cluster analysis ---
    assignments, centroids = kmeans_cells(state; k=3)
    if !isempty(assignments)
        print_cluster_analysis(state, assignments, centroids)
    else
        println("\n  No cells survived — cannot perform cluster analysis.")
    end

    # --- Summary: radial stratification ---
    println("\n  RADIAL STRATIFICATION SUMMARY")
    println("  " * "-"^50)
    for s in 1:N_SPECIES
        cells_s = [(id, c) for (id, c) in state.cells if c.species == s]
        if !isempty(cells_s)
            rs = [radial_dist(c.com[1], c.com[2], params.N) for (_, c) in cells_s]
            @printf("  %-22s  mean_r=%6.2f  (n=%d)\n",
                    SPECIES_NAMES[s], mean(rs), length(cells_s))
        else
            @printf("  %-22s  EXTINCT\n", SPECIES_NAMES[s])
        end
    end
    println()

    # --- Export figures (Figs 1–2 only; no radiodialysis data) ---
    println("  Exporting figures (plain CPM — Figs 1 & 2 only)...")
    export_figures(trajectory, nothing, params)

    return state, trajectory
end

# ============================================================
#  12. Radiodialysis Membrane Transport Coupling
#
#  Implements the minimal 3-equation PDE system from the Deep
#  Research synthesis (April 2026):
#
#  (1) Mobile contaminant: cylindrical reaction-diffusion PDE
#        ∂c/∂t = (1/r) ∂/∂r (r D_eff ∂c/∂r)
#                - (k_ads X + k_red X_red) c + k_des s
#
#  (2) Immobile phase ODE (biosorption + bioreduction):
#        ∂s/∂t = (k_ads X + k_red X_red) c - (k_des + k_loss) s
#
#  (3) Membrane damage ODE:
#        dm/dt = -k_dam Ḋ(R) m
#        P_eff(t) = P₀ exp(α D_cum),  D_cum = Ḋ(R)·t
#
#  Robin BC at r = R (membrane wall):
#        -D_eff ∂c/∂r|_{R} = P_eff(t) · (c(R,t) - c_ext)
#
#  CPM Coupling:
#    • c(r,t) projected onto 3D lattice as `contaminant` field
#    • s(r,t) tracked as radial-band averages
#    • P_eff(t) scales the wall boundary in update_nutrient_coupled!
#      — membrane damage reduces effective nutrient import
#
#  Key refs:
#    Fox et al. 2009 (SRNL/Nafion Donnan dialysis)
#    Lara et al. 2023 (Membranes, radiation-evolving permeability)
#    Renslow et al. 2017 (Shewanella uranium biofilm)
#    Aydogan Gokturk et al. 2022 (Nat. Commun., Donnan)
# ============================================================

# ---- Radiodialysis parameters --------------------------------

Base.@kwdef struct RadiolysisParams
    Nr::Int = 40               # radial grid points (0..R)

    # Transport
    D_eff::Float64 = 1e-3      # effective diffusivity (cm² s⁻¹)

    # Biosorption / bioreduction (Renslow et al. 2017)
    k_ads::Float64 = 0.05      # adsorption rate (cm³ g⁻¹ s⁻¹)
    k_red::Float64 = 0.02      # bioreduction rate (cm³ g⁻¹ s⁻¹)
    k_des::Float64 = 0.005     # desorption rate (s⁻¹)
    k_loss::Float64 = 0.001    # immobile-phase loss (s⁻¹)

    # Biomass
    X_total::Float64 = 1.0     # total dry-mass density (g cm⁻³)
    X_red::Float64 = 0.3       # metal-reducing fraction (Shewanella proxy)

    # Basis gate acknowledgement.  false refuses any X_total != 1.0.
    # Bool rather than a Symbol because this field is HDF5-serialised by
    # export_checkpoint.jl and must survive a restart; the exemption is binary
    # anyway, and WHICH sites hold it is pinned by the census test rather than
    # by a symbol name.  See _assert_basis_gate.
    basis_gate_ack::Bool = false

    # Membrane (Nafion / Donnan — Fox et al. 2009, Lara et al. 2023)
    P0::Float64 = 0.01         # baseline permeability (cm s⁻¹)
    alpha_P::Float64 = 0.02    # radiation-damage coefficient (Gy⁻¹)
    k_dam::Float64 = 0.005     # structural damage rate (Gy⁻¹)

    # Radiation dose rate at membrane
    Ddot_R::Float64 = 1.0      # Gy s⁻¹ (scaled with CPM radiation field)

    # External contaminant concentration (normalised)
    c_ext::Float64 = 1.0

    # ODE time step
    dt_rd::Float64 = 0.5       # s per CPM MCS (tune to match CPM time scale)
end

# ---- Radiodialysis state ------------------------------------

mutable struct RadiolysisState
    r_grid::Vector{Float64}   # Nr radial points  [0, R]
    c::Vector{Float64}        # mobile contaminant c(r)
    s::Vector{Float64}        # immobile phase     s(r)
    m::Float64                # membrane integrity  m ∈ [0,1]
    t::Float64                # simulation time
    params::RadiolysisParams
    # Sticky basis provenance.  NOT inferred from X_total: a coupled state can
    # legitimately produce mean(X_tot) == 1.0 (every sampled interior site
    # occupied), which would collide with the standalone default and slip the
    # gate; and c/s are PATH-dependent, so once a step has run on an occupancy
    # basis they stay gated even if X_total later returns to 1.0.  Rule 4: the
    # producer that installs the basis declares it, and the flag never clears.
    basis_from_occupancy::Bool
end

# Six-argument form: provenance defaults to false, so a standalone state is
# unmarked and only the coupled installers set it.
RadiolysisState(r_grid, c, s, m, t, params) =
    RadiolysisState(r_grid, c, s, m, t, params, false)

"""
Initialise RadiolysisState: clean interior, intact membrane.
"""
function init_radiolysis(rp::RadiolysisParams; R::Float64 = 1.0)
    r_grid = range(0.0, R, length = rp.Nr) |> collect
    RadiolysisState(r_grid,
                    zeros(Float64, rp.Nr),  # c: initially zero inside
                    zeros(Float64, rp.Nr),  # s: initially zero
                    1.0,                    # m: intact
                    0.0,                    # t
                    rp)
end

# ---- Finite-volume method-of-lines step ---------------------

"""
Advance the radiodialysis PDE system by one time step dt.
Substeps automatically if dt exceeds the explicit-Euler diffusion stability
bound (0.4·dr²/2D_eff).

n_sub = 1 AT THE GEOMETRY EVERY CALL SITE USES, not at this struct's own
default. Every caller passes R = N/2 in LATTICE units (`:1441`, `:1807`), giving
dr = 0.513 and dt_stable = 52.6 at N = 40, so legacy trajectories are bit-for-bit
unchanged. At the RadiolysisParams default R = 1.0 — which no call site uses, and
which is also `biofilms_radiodialysis.R`'s grid, integrated there by LSODA rather
than by this stepper — dt_stable = 0.132 and n_sub = 4.

So the safety margin here is an artifact of the documented r-in-lattice-units
error (RADIODIALYSIS: BLOCKED), not of dt_rd being small. Correcting the units to
R = 1.0 cm makes dt_rd = 0.5 genuinely unstable, and this guard is what absorbs
it. `biofilms_potts_jacc.jl` carries the identical wrapper for the same reason.
"""
function step_radiolysis!(rd::RadiolysisState, dt::Float64)
    _assert_basis_gate(rd.params.X_total, rd.params.basis_gate_ack,
                       rd.basis_from_occupancy)
    dr = rd.r_grid[2] - rd.r_grid[1]
    dt_stable = 0.4 * dr^2 / (2.0 * rd.params.D_eff)
    n_sub = max(1, ceil(Int, dt / dt_stable))
    dt_sub = dt / n_sub
    for _ in 1:n_sub
        _step_radiolysis_euler!(rd, dt_sub)
    end
end

"""
    _assert_basis_gate(X_total)

Refuse to integrate on a coupled biomass basis, EXPLICITLY.

RADIODIALYSIS: BLOCKED gates the biomass basis fed into this coupling, and
until now nothing here enforced it. The block was being done by an accident:
`uptake = k_ads*X_total + k_red*X_red` is the pre-fraction additive form, which
happens to misbehave at non-unit `X_total`, and that side effect was standing in
for a guard. An accidental tripwire is not a gate -- it can be removed by
someone tidying the arithmetic, leaving nothing to say the basis was blocked.

So the refusal is stated. `X_total == 1.0` is the standalone default and stays
allowed; anything else means a real basis was supplied, which is what the gate
covers.

WHY THE ARITHMETIC IS NOT ALSO FIXED. `biofilms_radiodialysis.R` now derives
`X_red = f_red_active * X_total` (`uptake_rate_of()`), and mirroring that here
would make this path *conformant*. It would not make it *correct*: the `X_red`
reaching it is `red_cells[i] / counts[i]` (`compute_radial_biomass`), one
species' occupied sites over ALL interior sites, which README.md:344 records as
"neither a biomass fraction nor a reducer fraction". Conformant arithmetic over
a quantity the repository has already refused is worse than visibly
non-conformant arithmetic over the same one: the defect would stop being visible
while staying just as gated. The arithmetic stays as a marker that this
reconciliation is unfinished.

Nor can the fraction be derived from parcel counts. Counts give a TAXONOMIC
fraction, and `active_from_taxonomic()` refuses converting one to an
active-reducer fraction without a measured activity fraction; `D-XRED` in
`data/calibration/reference_d_requirements.csv` records that as blocked by this
units error rather than by missing data.
"""
function _assert_basis_gate(X_total::Float64, ack::Bool = false,
                            from_occupancy::Bool = false)
    (X_total == 1.0 && !from_occupancy) && return nothing
    # THE ONE EXEMPTION.  validate_serial.jl steps this path only to reproduce a
    # bit-for-bit CPM trajectory, and records no radiodialysis quantity: its CSV
    # carries CPM columns plus rd.m, whose ODE (dm/dt = -k_dam*Ddot_R*m) has no
    # X_total or X_red in it.  That independence is not taken on trust -- it is
    # asserted in tests/radiodialysis_basis_gate.jl by running the harness at
    # two different gated bases and requiring byte-identical CSV.  Widen this
    # exemption and that test is what should stop you.
    ack && return nothing
    error("""
        RADIODIALYSIS: BLOCKED -- refusing to integrate at X_total = $X_total.

        Only the standalone default X_total == 1.0 is allowed here.

        A coupled basis reaches this path as mean(compute_radial_biomass(...)):
        one species' occupied sites over all interior sites, which is
        neither a biomass fraction nor a reducer fraction
        (README.md:344). The uptake arithmetic here is also still the
        pre-fraction additive form, unlike biofilms_radiodialysis.R's
        uptake_rate_of().

        This refusal is deliberate and is asserted by
        tests/radiodialysis_basis_gate.jl. Removing it to let a coupled run
        proceed re-opens the defect the gate names; repair the quantity first.
        """)
end

"""
One forward-Euler step of equations (1)–(3) with ghost-point Robin BC at r = R.
"""
function _step_radiolysis_euler!(rd::RadiolysisState, dt::Float64)
    rp  = rd.params
    Nr  = rp.Nr
    dr  = rd.r_grid[2] - rd.r_grid[1]

    c   = rd.c
    s   = rd.s
    m   = rd.m
    t   = rd.t

    # Current radiation-driven quantities
    D_cum  = rp.Ddot_R * t
    P_eff  = rp.P0 * exp(rp.alpha_P * D_cum)

    # (3) Membrane damage ODE
    dm_dt = -rp.k_dam * rp.Ddot_R * m

    uptake = rp.k_ads * rp.X_total + rp.k_red * rp.X_red

    dc = zeros(Float64, Nr)
    ds = zeros(Float64, Nr)

    # i = 1: r = 0 — L'Hôpital limit for (1/r)∂/∂r(r ∂/∂r) → 2 ∂²/∂r²
    dc[1] = rp.D_eff * 2.0 * (c[2] - c[1]) / dr^2 +
            (-uptake * c[1] + rp.k_des * s[1])

    # i = 2..Nr-1: interior finite-volume
    for i in 2:(Nr - 1)
        r_i    = rd.r_grid[i]
        r_plus  = r_i + 0.5 * dr
        r_minus = r_i - 0.5 * dr
        diff_cyl = rp.D_eff *
            (r_plus  * (c[i + 1] - c[i]) -
             r_minus * (c[i]     - c[i - 1])) /
            (r_i * dr^2)
        dc[i] = diff_cyl + (-uptake * c[i] + rp.k_des * s[i])
    end

    # i = Nr: r = R — Robin BC via ghost point
    #   c_ghost = c[Nr-1] - 2 dr P_eff (c[Nr] - c_ext) / D_eff
    let i = Nr
        r_i     = rd.r_grid[i]
        r_plus  = r_i + 0.5 * dr
        r_minus = r_i - 0.5 * dr
        c_ghost = c[Nr - 1] -
            2.0 * dr * P_eff * (c[Nr] - rp.c_ext) / rp.D_eff
        diff_cyl = rp.D_eff *
            (r_plus  * (c_ghost - c[Nr]) -
             r_minus * (c[Nr]   - c[Nr - 1])) /
            (r_i * dr^2)
        dc[i] = diff_cyl + (-uptake * c[Nr] + rp.k_des * s[Nr])
    end

    # Immobile phase: no spatial flux
    @. ds = uptake * c - (rp.k_des + rp.k_loss) * s

    # Forward Euler update
    @. rd.c = max(0.0, c + dt * dc)
    @. rd.s = max(0.0, s + dt * ds)
    rd.m = max(0.0, m + dt * dm_dt)
    rd.t = t + dt
end

# ---- 3D projection: radial → lattice ------------------------

"""
Project the 1D radial contaminant field c(r) onto the 3D lattice.
Each lattice site gets the c value of the nearest radial bin.
"""
function radial_to_3d!(field3d::Array{Float64, 3}, c_radial::Vector{Float64},
                        r_grid::Vector{Float64}, interior::BitArray{3},
                        N::Int)
    R_total = r_grid[end]
    dr = r_grid[2] - r_grid[1]
    Nr = length(r_grid)
    @inbounds for z in 1:N, y in 1:N, x in 1:N
        if interior[x, y, z]
            r = radial_dist(x, y, N) / (N / 2.0) * R_total
            idx = clamp(round(Int, r / dr) + 1, 1, Nr)
            field3d[x, y, z] = c_radial[idx]
        end
    end
end

"""
Compute radial-band-averaged biomass density from the CPM lattice.
Used to feed X_total / X_red into the radiolysis PDE.
"""
function compute_radial_biomass(state::CPMState, Nr::Int)
    p = state.params
    N = p.N
    R_total = N / 2.0
    dr = R_total / (Nr - 1)

    total_cells = zeros(Float64, Nr)
    red_cells   = zeros(Float64, Nr)
    counts      = zeros(Int, Nr)

    lat = state.lattice
    @inbounds for z in 1:N, y in 1:N, x in 1:N
        if !state.interior[x, y, z]
            continue
        end
        σ = lat[x, y, z]
        r = radial_dist(x, y, N)
        idx = clamp(round(Int, r / dr) + 1, 1, Nr)
        counts[idx] += 1
        if σ > 0
            total_cells[idx] += 1.0
            c = get(state.cells, Int(σ), nothing)
            if c !== nothing && c.species == SO  # Shewanella = metal-reducer
                red_cells[idx] += 1.0
            end
        end
    end

    # Normalise to density (fraction of sites occupied)
    X_total = [counts[i] > 0 ? total_cells[i] / counts[i] : 0.0 for i in 1:Nr]
    X_red   = [counts[i] > 0 ? red_cells[i]   / counts[i] : 0.0 for i in 1:Nr]
    return X_total, X_red
end

# ---- Coupled nutrient update (P_eff scales wall supply) -----

"""
Like update_nutrient! but the wall boundary is scaled by membrane integrity m(t).
Membrane damage reduces the effective nutrient import from outside.
"""
function update_nutrient_coupled!(state::CPMState, m::Float64)
    p = state.params
    N = p.N
    dt = p.dt_field
    C = state.nutrient
    C_new = copy(C)

    # Effective wall concentration: membrane damage reduces import
    C_wall_eff = p.C_wall * m

    @inbounds for z in 1:N, y in 1:N, x in 1:N
        if !state.interior[x, y, z]
            C_new[x, y, z] = C_wall_eff
            continue
        end
        diff = p.D_C * laplacian_3d(C, x, y, z, N)
        consumption = local_uptake(state.lattice, state.cells, p.uptake, x, y, z)
        C_new[x, y, z] = max(0.0, C[x, y, z] + dt * (diff - consumption))
    end

    copyto!(state.nutrient, C_new)
end

# ---- Extended snapshot with radiodialysis state -------------

struct CoupledSnapshot
    cpm_snap::Snapshot
    mcs::Int
    m::Float64                  # membrane integrity
    P_eff::Float64              # effective permeability
    c_mean::Float64             # mean mobile contaminant
    s_mean::Float64             # mean sorbed contaminant
    c_wall::Float64             # c at r = R (membrane face)
end

# ---- Coupled simulation driver ------------------------------

"""
    run_simulation_coupled(params, rp, n_mcs; seed=42)

Run the CPM simulation coupled with the radiodialysis membrane
transport PDE system. Returns `(state, rd, trajectory)`.

Coupling loop per MCS:
  1. CPM Metropolis step (mcs_step!)
  2. Advance radiodialysis PDE (step_radiolysis!) — X_total/X_red
     updated from lattice biomass every 10 MCS for efficiency
  3. Project c(r) onto 3D contaminant field (radial_to_3d!)
  4. Update melanin field (update_melanin!)
  5. Update nutrient field with P_eff scaling (update_nutrient_coupled!)
  6. Update cell centers of mass
  7. Snapshot every params.snapshot_interval MCS
"""
function run_simulation_coupled(params::CPMParams, rp::RadiolysisParams,
                                 n_mcs::Int; seed::Int = 42)
    rng = MersenneTwister(seed)
    state = init_state(params; seed = seed)
    rd    = init_radiolysis(rp; R = Float64(params.N) / 2.0)

    # Add contaminant field to the state — store in nutrient-slot shadow
    # (we'll track it separately as a local 3D array since CPMState is immutable)
    N = params.N
    contaminant_3d = zeros(Float64, N, N, N)

    trajectory = CoupledSnapshot[]

    update_centers_of_mass!(state)
    snap0 = take_snapshot(state, 0)
    cs0 = CoupledSnapshot(snap0, 0, rd.m,
                           rp.P0 * exp(rp.alpha_P * rp.Ddot_R * 0.0),
                           mean(rd.c), mean(rd.s),
                           rd.c[end])
    push!(trajectory, cs0)
    print_snapshot(snap0)
    # RATIO, NOT AN ABSOLUTE. P0 is an uncalibrated placeholder, so
    # `P0 * exp(...)` carries its arbitrariness into a number that looks
    # measured. P_eff/P0 is dimensionless and is what the manuscript reports.
    @printf("  [RD] t=%.1f  m=%.4f  P_eff/P0=%.5f  c_wall=%.4f  c_mean=%.4f\n\n",
            rd.t, rd.m, exp(rp.alpha_P * rd.t * rp.Ddot_R),
            rd.c[end], mean(rd.c))

    for mcs in 1:n_mcs
        # 1. CPM step
        mcs_step!(state, rng)

        # 2. Radiodialysis PDE step (biomass updated every 10 MCS)
        if mcs % 10 == 1
            X_tot, X_rd = compute_radial_biomass(state, rp.Nr)
            # Update spatially averaged biomass in params for this step
            rd.params = RadiolysisParams(
                Nr      = rp.Nr,
                D_eff   = rp.D_eff,
                k_ads   = rp.k_ads,
                k_red   = rp.k_red,
                k_des   = rp.k_des,
                k_loss  = rp.k_loss,
                X_total = mean(X_tot),
                X_red   = mean(X_rd),
                basis_gate_ack = rp.basis_gate_ack,
                P0      = rp.P0,
                alpha_P = rp.alpha_P,
                k_dam   = rp.k_dam,
                Ddot_R  = rp.Ddot_R,
                c_ext   = rp.c_ext,
                dt_rd   = rp.dt_rd
            )
            # Declared where the occupancy basis is INSTALLED, and never
            # cleared: c and s are path-dependent from here on, so provenance
            # cannot be re-derived from the current X_total (rule 4).
            rd.basis_from_occupancy = true
        end
        step_radiolysis!(rd, rp.dt_rd)

        # 3. Project c(r) → 3D contaminant field
        radial_to_3d!(contaminant_3d, rd.c, rd.r_grid, state.interior, N)

        # 4. Melanin reaction-diffusion
        update_melanin!(state)

        # 5. Nutrient with membrane integrity coupling
        update_nutrient_coupled!(state, rd.m)

        # 6. Update cell CoMs
        update_centers_of_mass!(state)

        # 7. Snapshot
        if mcs % params.snapshot_interval == 0 || mcs == n_mcs
            snap = take_snapshot(state, mcs)
            P_eff_now = rp.P0 * exp(rp.alpha_P * rp.Ddot_R * rd.t)
            cs = CoupledSnapshot(snap, mcs, rd.m, P_eff_now,
                                  mean(rd.c), mean(rd.s), rd.c[end])
            push!(trajectory, cs)
            print_snapshot(snap)
            @printf("  [RD] t=%.1f  m=%.4f  P_eff/P0=%.5f  c_wall=%.4f  c_mean=%.4f  s_mean=%.4f\n\n",
                    rd.t, rd.m, P_eff_now / rp.P0, rd.c[end], mean(rd.c), mean(rd.s))
        end
    end

    return state, rd, contaminant_3d, trajectory
end

# ---- Radiodialysis-coupled main -----------------------------

"""
    print_membrane_report(rp, rd, n_mcs)

Print the membrane-integrity/permeability/dose summary. Split out of
`main_coupled()` so it's testable without a full simulation run, and so a
negative control can capture exactly what it prints: no fabricated physical
units (Gy, cm/s) attached to a quantity this repository cannot calibrate —
see docs/research/session_claims_2026-08-24_redteam.md and PR #12.
"""
function print_membrane_report(rp::RadiolysisParams, rd::RadiolysisState, n_mcs::Int)
    @printf("  Final membrane integrity: m = %.4f\n", rd.m)
    @printf("  Final P_eff / P0 = %.4f  (dimensionless ratio, exact by construction from alpha_P·Ddot_R·dt_rd·n_MCS — not a measurement)\n",
            exp(rp.alpha_P * rp.Ddot_R * rd.t))
    @printf("  Final P_eff, absolute: not computed in this work (P0 = %.3g cm/s is a Nafion-117 literature prior, not a calibration)\n", rp.P0)
    @printf("  D_cum (dimensionless placeholder) after %d MCS: %.1f  (no seconds_per_mcs conversion exists — not Gy)\n", n_mcs, rp.Ddot_R * rd.t)
end

"""
    main_coupled()

Run the fully coupled CPM + radiodialysis simulation and print
the radial stratification alongside contaminant uptake summary.
"""
function main_coupled()
    println("="^72)
    println("  CPM + Radiodialysis Membrane Transport (Coupled)")
    println("  Kinder & Faulkner (2026) — Equations (1)–(3)")
    println("  H = H_adh + H_vol + H_rad + H_mel  (H_pair is a diagnostic, not in the acceptance path)")
    println("  PDE: ∂c/∂t = ∇·(D ∇c) - uptake·c + k_des·s (cylindrical)")
    println("="^72)
    println()

    params = CPMParams(N = 40, n_cells_per_species = 6, snapshot_interval = 20)
    rp     = RadiolysisParams(Nr = 40, Ddot_R = 1.0, c_ext = 1.0)

    @printf("  CPM lattice:  %d³  |  Species: %d  |  Cells: %d each\n",
            params.N, N_SPECIES, params.n_cells_per_species)
    @printf("  Radiolysis:   Nr=%d  D_eff=%.3g  P₀=%.3g  α=%.3g\n",
            rp.Nr, rp.D_eff, rp.P0, rp.alpha_P)
    @printf("  Membrane:     Ḋ(R)=%.1f (placeholder, no seconds_per_mcs conversion — not Gy/s)  c_ext=%.2f  k_dam=%.3g\n\n",
            rp.Ddot_R, rp.c_ext, rp.k_dam)

    n_mcs = 100
    @printf("  Running %d Monte Carlo Steps (coupled)...\n\n", n_mcs)

    t_start = time()
    state, rd, contaminant_3d, traj = run_simulation_coupled(params, rp, n_mcs)
    elapsed = time() - t_start

    @printf("\n  Simulation completed in %.1f seconds.\n", elapsed)
    @printf("  Surviving cells: %d\n", length(state.cells))
    print_membrane_report(rp, rd, n_mcs)

    # Contaminant uptake summary
    println("\n  CONTAMINANT UPTAKE SUMMARY")
    println("  " * "-"^50)
    total_sorbed = sum(rd.s) * (rd.r_grid[2] - rd.r_grid[1]) * params.N^2
    @printf("  Total sorbed (∫s dr × area):  %.4f (normalised units)\n",
            total_sorbed)
    @printf("  c at axis (r=0):  %.4f\n", rd.c[1])
    @printf("  c at membrane (r=R):  %.4f\n", rd.c[end])

    # Cluster analysis
    assignments, centroids = kmeans_cells(state; k = 3)
    if !isempty(assignments)
        print_cluster_analysis(state, assignments, centroids)
    end

    # Radial stratification
    println("\n  RADIAL STRATIFICATION (CPM + radiodialysis)")
    println("  " * "-"^60)
    for s in 1:N_SPECIES
        cells_s = [(id, c) for (id, c) in state.cells if c.species == s]
        if !isempty(cells_s)
            rs = [radial_dist(c.com[1], c.com[2], params.N) for (_, c) in cells_s]
            mr = mean(rs)
            R_total = params.N / 2.0
            r_norm = mr / R_total
            Nr = rp.Nr
            idx = clamp(round(Int, r_norm * (Nr - 1)) + 1, 1, Nr)
            c_local = rd.c[idx]
            s_local = rd.s[idx]
            @printf("  %-22s  mean_r=%5.1f  c(r)=%.4f  s(r)=%.4f  (n=%d)\n",
                    SPECIES_NAMES[s], mr, c_local, s_local, length(cells_s))
        else
            @printf("  %-22s  EXTINCT\n", SPECIES_NAMES[s])
        end
    end
    println()

    # --- Export figures ---
    println("  Exporting figures...")
    export_figures(traj, traj, params)

    return state, rd, traj
end

# ============================================================
#  12b. Coupling foundation: dose contract, lifecycle primitive,
#       windowed simulation API   (feat/openmc-dose-coupling)
#
#  - Physical dose (Gy) is stored in dedicated CPMState fields and is
#    NEVER written into the dimensionless `radiation` Hamiltonian signal
#    or the melanin drive except through an explicitly supplied
#    transform (import_dose_field!). See docs/audit_biofilms_potts.md §5.
#  - divide_cell! is a MANUALLY INVOKED demonstration primitive. There
#    is deliberately no automatic division scheduler: the CPM has no
#    growth dynamics from which biological generations could be claimed
#    (audit §2).
# ============================================================

"""
    import_dose_field!(state, dose_mean_Gy_s, dose_sd_Gy_s;
                       hamiltonian_transform = nothing,
                       melanin_transform = nothing)

Install a physical dose-rate field (Gy/s, CPM lattice shape). The
Hamiltonian signal (`state.radiation`) and melanin drive are updated
ONLY if the corresponding transform is supplied — passing `nothing`
keeps the legacy dimensionless signals untouched, so importing dose
never silently changes the dynamics.
"""
function import_dose_field!(state::CPMState,
                            dose_mean_Gy_s::Array{Float64, 3},
                            dose_sd_Gy_s::Array{Float64, 3};
                            hamiltonian_transform = nothing,
                            melanin_transform = nothing)
    size(dose_mean_Gy_s) == size(state.lattice) ||
        error("import_dose_field!: dose field shape $(size(dose_mean_Gy_s)) ≠ lattice $(size(state.lattice))")
    copyto!(state.dose_rate_mean_Gy_s, dose_mean_Gy_s)
    copyto!(state.dose_rate_sd_Gy_s, dose_sd_Gy_s)
    if hamiltonian_transform !== nothing
        state.radiation .= hamiltonian_transform.(dose_mean_Gy_s)
    end
    if melanin_transform !== nothing
        state.melanin_drive .= melanin_transform.(dose_mean_Gy_s)
    end
    state.dose_active = true
    return state
end

"""
    accrue_dose!(state)

Per-MCS dose attribution: ΔD_v = Ḋ_v · seconds_per_mcs accrued to the
CURRENT voxel labels (cells move between transport solves; attributing
only at transport checkpoints would misassign the interval's exposure to
endpoint geometry). Cell lifetime dose accrues the unweighted site-mean
increment.
"""
function accrue_dose!(state::CPMState)
    state.dose_active || return state
    dt = state.params.seconds_per_mcs
    isnan(dt) &&
        error("accrue_dose!: seconds_per_mcs is not configured (physical unit contract, audit §5)")
    N = state.params.N
    # ponytail: unweighted site-mean cell dose — equal voxel masses assumed
    # until the material config supplies densities; replace with mass weights then.
    @inbounds for z in 1:N, y in 1:N, x in 1:N
        ΔD = state.dose_rate_mean_Gy_s[x, y, z] * dt
        state.dose_increment_Gy[x, y, z] = ΔD
        state.accumulated_dose_Gy[x, y, z] += ΔD
        σ = state.lattice[x, y, z]
        if σ > 0
            c = get(state.cells, Int(σ), nothing)
            if c !== nothing
                c.lifetime_dose_Gy += ΔD / max(c.volume, 1)
            end
        end
    end
    return state
end

# ---- 26-connectivity check for daughter integrity -------------

function _is_connected_26(sites::Vector{NTuple{3, Int}})
    isempty(sites) && return false
    site_set = Set(sites)
    seen = Set{NTuple{3, Int}}()
    stack = [sites[1]]
    push!(seen, sites[1])
    while !isempty(stack)
        (x, y, z) = pop!(stack)
        for (dx, dy, dz) in NEIGHBORS_26
            nb = (x + dx, y + dy, z + dz)
            if nb in site_set && !(nb in seen)
                push!(seen, nb)
                push!(stack, nb)
            end
        end
    end
    return length(seen) == length(sites)
end

"""
    divide_cell!(state, id; min_daughter_volume = 8)
        -> DivisionEvent | Symbol

Binary fission of cell `id`: the parent identity ENDS (archived with
fate `:divided`, lifetime dose finalized) and two daughters with fresh
IDs are created (generation = parent + 1, lineage inherited, exposure
records start at zero).

Deterministic and RNG-free: the center of mass and second-moment matrix
are computed FRESH from the current lattice sites (stored COM may be
stale mid-MCS), and the split plane passes through the COM perpendicular
to the principal axis (largest eigenvalue of the symmetric second-moment
matrix). Sites exactly on the plane go to daughter A.

Fails WITHOUT mutating state, returning a reason Symbol:
`:no_such_cell`, `:below_min_volume`, `:empty_daughter`,
`:disconnected_daughter` (26-connectivity on each daughter).
"""
function divide_cell!(state::CPMState, id::Int; min_daughter_volume::Int = 8)
    haskey(state.cells, id) || return :no_such_cell
    parent = state.cells[id]
    N = state.params.N

    sites = NTuple{3, Int}[]
    for z in 1:N, y in 1:N, x in 1:N
        state.lattice[x, y, z] == Int32(id) && push!(sites, (x, y, z))
    end
    isempty(sites) && return :no_such_cell

    com = [mean(s[1] for s in sites),
           mean(s[2] for s in sites),
           mean(s[3] for s in sites)]
    S = zeros(3, 3)
    for s in sites
        d = [s[1] - com[1], s[2] - com[2], s[3] - com[3]]
        S .+= d * d'
    end
    # eigen(Symmetric) returns ascending eigenvalues; principal axis is last.
    axis = eigen(Symmetric(S)).vectors[:, 3]

    a = NTuple{3, Int}[]
    b = NTuple{3, Int}[]
    for s in sites
        proj = (s[1] - com[1]) * axis[1] + (s[2] - com[2]) * axis[2] +
               (s[3] - com[3]) * axis[3]
        proj >= 0 ? push!(a, s) : push!(b, s)
    end

    (isempty(a) || isempty(b)) && return :empty_daughter
    (length(a) < min_daughter_volume || length(b) < min_daughter_volume) &&
        return :below_min_volume
    (_is_connected_26(a) && _is_connected_26(b)) || return :disconnected_daughter

    # ---- Commit point: mutate state ----
    da = state.next_id
    db = state.next_id + 1
    state.next_id += 2

    for s in a
        state.lattice[s...] = Int32(da)
    end
    for s in b
        state.lattice[s...] = Int32(db)
    end

    gen = parent.generation + 1
    mcs = state.current_mcs
    state.cells[da] = CellInfo(parent.species, length(a),
        [mean(s[1] for s in a), mean(s[2] for s in a), mean(s[3] for s in a)],
        parent.lineage_id, id, gen, mcs, 0.0)
    state.cells[db] = CellInfo(parent.species, length(b),
        [mean(s[1] for s in b), mean(s[2] for s in b), mean(s[3] for s in b)],
        parent.lineage_id, id, gen, mcs, 0.0)

    push!(state.archive,
          ArchivedCell(id, parent.species, parent.lineage_id, parent.parent_id,
                       parent.generation, parent.birth_mcs, mcs, :divided,
                       parent.lifetime_dose_Gy))
    delete!(state.cells, id)

    event = DivisionEvent(mcs, state.physical_time_s, id, da, db,
                          parent.lineage_id, parent.generation, gen,
                          length(sites), length(a), length(b))
    push!(state.division_log, event)
    return event
end

# ---- Windowed simulation API (RNG ownership; restart foundation) ----
# The simulation struct owns the exact stochastic state. Mirrors
# run_simulation_coupled's initialization and per-MCS call order exactly,
# so an unbroken legacy run and a windowed run produce identical
# trajectories (tested in tests/genealogy_tests.jl).

mutable struct CoupledSimulation
    state::CPMState
    rd::RadiolysisState
    contaminant::Array{Float64, 3}
    rng::MersenneTwister
    mcs::Int
    physical_time_s::Float64
end

function init_coupled_simulation(params::CPMParams, rp::RadiolysisParams;
                                 seed::Int = 42)
    rng = MersenneTwister(seed)
    state = init_state(params; seed = seed)
    rd = init_radiolysis(rp; R = Float64(params.N) / 2.0)
    contaminant = zeros(Float64, params.N, params.N, params.N)
    update_centers_of_mass!(state)
    return CoupledSimulation(state, rd, contaminant, rng, 0, 0.0)
end

function advance_window!(sim::CoupledSimulation, n_mcs::Int)
    state = sim.state
    rd = sim.rd
    N = state.params.N
    for _ in 1:n_mcs
        sim.mcs += 1
        mcs = sim.mcs

        mcs_step!(state, sim.rng)

        if mcs % 10 == 1
            X_tot, X_rd = compute_radial_biomass(state, rd.params.Nr)
            rp = rd.params
            rd.params = RadiolysisParams(
                Nr = rp.Nr, D_eff = rp.D_eff, k_ads = rp.k_ads,
                k_red = rp.k_red, k_des = rp.k_des, k_loss = rp.k_loss,
                X_total = mean(X_tot), X_red = mean(X_rd),
                basis_gate_ack = rp.basis_gate_ack,
                P0 = rp.P0, alpha_P = rp.alpha_P, k_dam = rp.k_dam,
                Ddot_R = rp.Ddot_R, c_ext = rp.c_ext, dt_rd = rp.dt_rd)
            # Declared where the occupancy basis is INSTALLED, and never
            # cleared: c and s are path-dependent from here on, so provenance
            # cannot be re-derived from the current X_total (rule 4).
            rd.basis_from_occupancy = true
        end
        step_radiolysis!(rd, rd.params.dt_rd)
        radial_to_3d!(sim.contaminant, rd.c, rd.r_grid, state.interior, N)
        update_melanin!(state)
        update_nutrient_coupled!(state, rd.m)
        update_centers_of_mass!(state)

        state.dose_active && accrue_dose!(state)
        spm = state.params.seconds_per_mcs
        if !isnan(spm)
            sim.physical_time_s += spm
            state.physical_time_s = sim.physical_time_s
        end
    end
    return sim
end

# ============================================================
#  13. Figure export (CairoMakie — headless, publication-ready)
#
#  Four figures produced by export_figures():
#
#  Fig 1 — Radial stratification over MCS
#           Line plot: mean radial distance per species vs. MCS.
#           Shows radiotropic species drifting toward the axis,
#           radiosensitive species drifting toward the outer wall.
#
#  Fig 2 — Melanin accumulation
#           Line plot: mean melanin field value at occupied sites
#           for melanin-producing species (CN, CS, AN) vs. MCS.
#
#  Fig 3 — Membrane integrity and permeability
#           Dual-axis: m(t) integrity [left] and P_eff(t)/P₀ [right]
#           vs. simulation time. Derived from CoupledSnapshot log.
#
#  Fig 4 — Contaminant penetration
#           Dual-axis: c_wall and c_mean vs. time. Quantifies how
#           much contaminant entered and how far it penetrated.
#
#  All outputs saved under outdir (default: preprint/figures/).
# ============================================================

using CairoMakie

# Publication theme: clean, LaTeX-compatible
function _set_theme!()
    set_theme!(
        Theme(
            fontsize  = 13,
            fonts     = (; regular = "Arial", bold = "Arial Bold"),
            Axis      = (
                spinewidth       = 1.2,
                xgridvisible     = false,
                ygridvisible     = false,
                xminorticksvisible = true,
                yminorticksvisible = true,
                xtickalign       = 1,
                ytickalign       = 1,
            ),
            Lines     = (linewidth = 2.2,),
            Legend    = (framevisible = false, padding = (4,4,4,4)),
        )
    )
end

# Palette: one colour per species (matches SPECIES_NAMES order)
const FIG_COLORS = [
    colorant"#e6194b",   # CN  — red
    colorant"#3cb44b",   # DR  — green
    colorant"#4363d8",   # CS  — blue
    colorant"#f58231",   # BS  — orange
    colorant"#911eb4",   # AN  — purple
    colorant"#42d4f4",   # SO  — cyan
    colorant"#f032e6",   # OI  — magenta
]

# Short labels for legend
const FIG_LABELS = [
    "C. neoformans",
    "D. radiodurans",
    "C. sphaerospermum",
    "B. subtilis",
    "A. niger",
    "S. oneidensis",
    "O. intermedium",
]

"""
    export_figures(trajectory, coupled_traj; outdir, params)

Generate and save the four publication figures from simulation output.
`trajectory` can be Vector{Snapshot} or Vector{CoupledSnapshot}.
`coupled_traj` is Vector{CoupledSnapshot} (pass `nothing` to skip Figs 3–4).
"""
function export_figures(trajectory::Vector{<:Any},
                         coupled_traj,
                         params::CPMParams;
                         outdir::String = joinpath(@__DIR__, "preprint", "figures"))

    mkpath(outdir)
    _set_theme!()

    # Unpack snapshots (handle both Snapshot and CoupledSnapshot)
    snaps = if eltype(trajectory) == CoupledSnapshot
        [cs.cpm_snap for cs in trajectory]
    else
        trajectory
    end

    mcs_vec = [s.mcs for s in snaps]
    N = params.N
    R_wall = N / 2.0   # wall in lattice units

    # ----------------------------------------------------------------
    # Figure 1: Radial stratification — normalised r/R, coloured bands
    # ----------------------------------------------------------------
    fig1 = Figure(size = (860, 500), figure_padding = (12, 20, 10, 10))
    ax1  = Axis(fig1[1,1],
                xlabel      = "Monte Carlo Steps",
                ylabel      = "Mean radial position  r / R",
                title       = "Mean radial position — cylindrical CPM bioreactor",
                titlesize   = 15,
                xlabelsize  = 13,
                ylabelsize  = 13,
                yticks      = LinearTicks(6),
                limits      = (nothing, (0.0, 1.15)))

    # Shaded zones labeled by the radiation field only (I(r) is maximal on the
    # axis and decays outward — see init_radiation!). Biological occupancy of
    # each zone is an outcome of the run, not asserted here.
    inner_zone = (0.0, 0.45)    # inner high-radiation region (r/R)
    outer_zone = (0.65, 1.0)    # outer lower-radiation region (r/R)
    poly!(ax1,
          Point2f[(mcs_vec[1], inner_zone[1]), (mcs_vec[end], inner_zone[1]),
                  (mcs_vec[end], inner_zone[2]), (mcs_vec[1], inner_zone[2])];
          color = (:firebrick, 0.08), strokewidth = 0)
    poly!(ax1,
          Point2f[(mcs_vec[1], outer_zone[1]), (mcs_vec[end], outer_zone[1]),
                  (mcs_vec[end], outer_zone[2]), (mcs_vec[1], outer_zone[2])];
          color = (:royalblue, 0.08), strokewidth = 0)

    # Dashed wall line at r/R = 1.0
    hlines!(ax1, [1.0]; color = (:black, 0.40), linestyle = :dash, linewidth = 1.5)

    # Zone annotation text (neutral physical labels)
    text!(ax1, mcs_vec[end] * 0.02, 0.22;
          text = "inner\nhigh-radiation region", fontsize = 10,
          color = (:firebrick, 0.70), align = (:left, :center))
    text!(ax1, mcs_vec[end] * 0.02, 0.82;
          text = "outer\nlower-radiation region", fontsize = 10,
          color = (:royalblue, 0.70), align = (:left, :center))

    line_handles = []
    for s in 1:N_SPECIES
        mean_r_vec = Float64[]
        for snap in snaps
            sd = snap.species_data[s]
            r_norm = sd.n_cells > 0 ? sd.mean_r / R_wall : NaN
            push!(mean_r_vec, r_norm)
        end
        lw = s in [CN, CS] ? 3.2 : 2.2
        ls = s in [CN, CS] ? :solid : (s in [BS, SO] ? :dash : :dot)
        h = lines!(ax1, mcs_vec, mean_r_vec;
                   color     = FIG_COLORS[s],
                   label     = FIG_LABELS[s],
                   linewidth = lw,
                   linestyle = ls)
        push!(line_handles, h)
    end

    # Annotate final r/R values for key species
    for s in [CN, CS, BS]
        last_r = NaN
        for snap in reverse(snaps)
            sd = snap.species_data[s]
            if sd.n_cells > 0
                last_r = sd.mean_r / R_wall
                break
            end
        end
        isnan(last_r) && continue
        text!(ax1, Float64(mcs_vec[end]) + Float64(mcs_vec[end]) * 0.01, last_r;
              text      = @sprintf("%.2f", last_r),
              fontsize  = 9,
              color     = FIG_COLORS[s],
              align     = (:left, :center))
    end

    Legend(fig1[1,2], ax1;
           labelsize   = 11,
           rowgap      = 4,
           framevisible = false,
           title       = "Species",
           titlesize   = 12)
    save(joinpath(outdir, "fig1_radial_stratification.pdf"), fig1)
    save(joinpath(outdir, "fig1_radial_stratification.png"), fig1; px_per_unit = 2)
    @printf("  Saved fig1_radial_stratification (.pdf + .png)\n")

    # ----------------------------------------------------------------
    # Figure 2: Melanin accumulation — filled area + final annotations
    # ----------------------------------------------------------------
    mel_producers = [CN, CS, AN]

    fig2 = Figure(size = (780, 460), figure_padding = (12, 20, 10, 10))
    ax2  = Axis(fig2[1,1],
                xlabel    = "Monte Carlo Steps",
                ylabel    = "Mean melanin field value at occupied sites",
                title     = "Melanin field at occupied sites — dimensionless model units",
                titlesize = 15,
                xlabelsize = 13,
                ylabelsize = 13)

    for s in mel_producers
        mel_vec = Float64[snap.species_data[s].mean_melanin for snap in snaps]
        # Filled band from 0 to line
        band!(ax2, mcs_vec, zeros(length(mcs_vec)), mel_vec;
              color = (FIG_COLORS[s], 0.18))
        lines!(ax2, mcs_vec, mel_vec;
               color     = FIG_COLORS[s],
               linewidth = 2.8,
               label     = FIG_LABELS[s])
        # Annotate final value
        final_val = mel_vec[end]
        text!(ax2, Float64(mcs_vec[end]) + Float64(mcs_vec[end]) * 0.01, final_val;
              text     = @sprintf("%.2f", final_val),
              fontsize = 10,
              color    = FIG_COLORS[s],
              align    = (:left, :center))
    end

    # Mark radiotropic species with a label in the legend title
    Legend(fig2[1,2], ax2;
           labelsize    = 11,
           rowgap       = 4,
           framevisible = false,
           title        = "Melanin producers",
           titlesize    = 11)
    # Override first two labels to add star
    text!(ax2, mcs_vec[1], -0.05;
          text = "Ordering follows the hand-specified alpha_M_species scales; not a measured quantity",
          fontsize = 8.5, color = (:gray40, 1.0), align = (:left, :top))

    save(joinpath(outdir, "fig2_melanin_accumulation.pdf"), fig2)
    save(joinpath(outdir, "fig2_melanin_accumulation.png"), fig2; px_per_unit = 2)
    @printf("  Saved fig2_melanin_accumulation (.pdf + .png)\n")

    # ----------------------------------------------------------------
    # Figures 3 & 4 require CoupledSnapshot data
    # ----------------------------------------------------------------
    if coupled_traj === nothing
        @printf("  Skipping Figs 3–4 (no coupled_traj provided)\n")
        return
    end

    ct = coupled_traj   # Vector{CoupledSnapshot}
    t_vec      = [cs.cpm_snap.mcs for cs in ct]   # MCS axis
    m_vec      = [cs.m        for cs in ct]
    Peff_vec   = [cs.P_eff    for cs in ct]
    P0_val     = Peff_vec[1]   # baseline (t=0)
    Peff_norm  = Peff_vec ./ P0_val
    c_wall_vec = [cs.c_wall   for cs in ct]
    c_mean_vec = [cs.c_mean   for cs in ct]
    s_mean_vec = [cs.s_mean   for cs in ct]

    # ----------------------------------------------------------------
    # Figure 3: Membrane integrity and effective permeability — dual axis
    # ----------------------------------------------------------------
    fig3  = Figure(size = (820, 500), figure_padding = (12, 24, 10, 10))

    ax3l  = Axis(fig3[1,1],
                 xlabel          = "Monte Carlo Steps",
                 ylabel          = "Membrane integrity  m(t)",
                 title           = "Membrane damage and radiation-driven permeability",
                 titlesize       = 15,
                 xlabelsize      = 13,
                 ylabelsize      = 13,
                 yticklabelcolor = colorant"#1f77b4",
                 ylabelcolor     = colorant"#1f77b4",
                 limits          = (nothing, (0.0, 1.1)))

    ax3r  = Axis(fig3[1,1],
                 ylabel          = "P_eff / P₀",
                 yaxisposition   = :right,
                 yticklabelcolor = colorant"#d62728",
                 ylabelcolor     = colorant"#d62728")
    hidespines!(ax3r)
    hidexdecorations!(ax3r)

    # Shaded region showing intact vs damaged membrane
    band!(ax3l, t_vec, m_vec, ones(length(t_vec));
          color = (colorant"#1f77b4", 0.12))

    l1 = lines!(ax3l, t_vec, m_vec;
                color     = colorant"#1f77b4",
                linewidth = 3.0,
                label     = "m(t)  membrane integrity")
    l2 = lines!(ax3r, t_vec, Peff_norm;
                color     = colorant"#d62728",
                linewidth = 2.8,
                linestyle = :dash,
                label     = "P_eff / P₀  permeability ratio")

    # Key annotations
    m_final   = m_vec[end]
    Peff_final = Peff_norm[end]
    text!(ax3l, Float64(t_vec[end]) * 0.55, m_final + 0.05;
          text    = @sprintf("m = %.3f", m_final),
          fontsize = 10, color = colorant"#1f77b4", align = (:center, :bottom))
    text!(ax3r, Float64(t_vec[end]) * 0.55, Peff_final * 0.88;
          text    = @sprintf("%.1f× baseline", Peff_final),
          fontsize = 10, color = colorant"#d62728", align = (:center, :top))

    Legend(fig3[1,2], [l1, l2],
           ["m(t)  integrity", "P_eff / P₀  permeability"];
           labelsize    = 11,
           rowgap       = 6,
           framevisible = false)
    save(joinpath(outdir, "fig3_membrane_transport.pdf"), fig3)
    save(joinpath(outdir, "fig3_membrane_transport.png"), fig3; px_per_unit = 2)
    @printf("  Saved fig3_membrane_transport (.pdf + .png)\n")

    # ----------------------------------------------------------------
    # Figure 4: Contaminant penetration — c_wall, c_mean, s_mean
    #           s_mean is on right axis; scale matched to left for visibility
    # ----------------------------------------------------------------
    fig4  = Figure(size = (820, 500), figure_padding = (12, 24, 10, 10))

    ax4l  = Axis(fig4[1,1],
                 xlabel      = "Monte Carlo Steps",
                 ylabel      = "Mobile contaminant  c / c_ext",
                 title       = "Contaminant penetration and biosorption",
                 titlesize   = 15,
                 xlabelsize  = 13,
                 ylabelsize  = 13,
                 limits      = (nothing, (-0.02, 1.05)))

    # Scale s_mean to [0,1] for dual-axis legibility; annotate true scale on right
    s_max   = max(maximum(s_mean_vec), 1e-12)
    s_norm  = s_mean_vec ./ s_max

    ax4r  = Axis(fig4[1,1],
                 ylabel          = @sprintf("Sorbed phase  s_mean  (×%.4f c_ext)", s_max),
                 yaxisposition   = :right,
                 yticklabelcolor = colorant"#2ca02c",
                 ylabelcolor     = colorant"#2ca02c",
                 limits          = (nothing, (-0.02, 1.05)))
    hidespines!(ax4r)
    hidexdecorations!(ax4r)

    # Filled bands
    band!(ax4l, t_vec, zeros(length(t_vec)), c_wall_vec;
          color = (colorant"#d62728", 0.12))
    band!(ax4l, t_vec, zeros(length(t_vec)), c_mean_vec;
          color = (colorant"#1f77b4", 0.12))

    l3 = lines!(ax4l, t_vec, c_wall_vec;
                color     = colorant"#d62728",
                linewidth = 3.0,
                label     = "c(R,t)  membrane wall")
    l4 = lines!(ax4l, t_vec, c_mean_vec;
                color     = colorant"#1f77b4",
                linewidth = 2.5,
                linestyle = :dash,
                label     = "c_mean  interior average")
    l5 = lines!(ax4r, t_vec, s_norm;
                color     = colorant"#2ca02c",
                linewidth = 2.5,
                linestyle = :dashdot,
                label     = "s_mean  sorbed (immobile)")

    # Annotations
    c_wall_final = c_wall_vec[end]
    c_mean_final = c_mean_vec[end]
    text!(ax4l, Float64(t_vec[end]) * 0.55, c_wall_final + 0.04;
          text    = @sprintf("c(R) = %.0f%% c_ext", c_wall_final * 100),
          fontsize = 10, color = colorant"#d62728", align = (:center, :bottom))
    text!(ax4l, Float64(t_vec[end]) * 0.55, c_mean_final + 0.04;
          text    = @sprintf("c_mean = %.1f%% c_ext  (%.0f%% non-penetration)",
                             c_mean_final * 100, (1.0 - c_mean_final) * 100),
          fontsize = 10, color = colorant"#1f77b4", align = (:center, :bottom))

    Legend(fig4[1,2], [l3, l4, l5],
           ["c(R,t)  at wall", "c_mean  interior", "s_mean  sorbed"];
           labelsize    = 11,
           rowgap       = 6,
           framevisible = false)
    save(joinpath(outdir, "fig4_contaminant_penetration.pdf"), fig4)
    save(joinpath(outdir, "fig4_contaminant_penetration.png"), fig4; px_per_unit = 2)
    @printf("  Saved fig4_contaminant_penetration (.pdf + .png)\n")

    @printf("\n  All figures written to: %s\n", outdir)
end

# ============================================================
#  Run if executed directly
# ============================================================
if abspath(PROGRAM_FILE) == @__FILE__
    # Run coupled simulation by default; fall back to plain CPM
    # if --no-radiolysis is passed as first argument.
    if "--no-radiolysis" in ARGS
        main()
    else
        main_coupled()
    end
end
