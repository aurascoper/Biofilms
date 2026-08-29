#!/usr/bin/env julia
# ============================================================
#  JACC.jl port of biofilms_potts.jl — GPU-parallel CPM
#
#  8-color (2×2×2) checkerboard Metropolis: Moore-26 adhesion means
#  no two same-color sites share a face/edge/corner, so each color
#  pass is race-free on the lattice. Each thread owns one TARGET
#  site and pulls a random Moore neighbor as source (role inversion
#  of the serial copy-attempt; identical ordered-pair distribution).
#
#  Backend is a JACC preference (restart to switch):
#    julia --project=. -e 'using JACC; JACC.set_backend("threads")'
#    julia --project=. -e 'using JACC; JACC.set_backend("amdgpu"; storage=:host)'
#
#  Usage:
#    julia --project=. [-t auto] biofilms_potts_jacc.jl [seeds...]
#    julia --project=. biofilms_potts_jacc.jl --selftest
#
#  Validation is statistical (checkerboard ≠ random-sequential
#  dynamics): compare the CSV lines against validate_serial.jl.
# ============================================================

using Random, Statistics, Printf
using JACC
JACC.@init_backend

# ---- Species constants (mirrors serial section 1/2) ----------

const SPECIES_NAMES = ["C. neoformans", "D. radiodurans", "C. sphaerospermum",
    "B. subtilis", "A. niger", "S. oneidensis", "O. intermedium AM7"]
const N_SPECIES = 7
const CN = 1; const DR = 2; const CS = 3; const BS = 4
const AN = 5; const SO = 6; const OI = 7

const BETA_ION = Float32[-5e-5, 2.5e-5, -5e-5, 3e-3, 2.5e-4, 7.5e-2, 1e-2]
const ALPHA_M  = Float32[0.10, 0.0, 0.14, 0.0, 0.065, 0.0, 0.0]
const UPTAKE   = Float32[0.01, 0.02, 0.01, 0.03, 0.01, 0.03, 0.02]
# 0.5 melanin coupling for the RADIOTROPIC species (CN, CS), else 0 —
# a spatial preference for melanin-rich sites, not energy transduction
const MEL_COEF = Float32[0.5, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0]

function build_J_matrix()
    J = fill(12.0f0, N_SPECIES + 1, N_SPECIES + 1)
    for i in 1:(N_SPECIES + 1)
        J[i, i] = 0.0f0
    end
    for s in 2:(N_SPECIES + 1)
        J[1, s] = 16.0f0; J[s, 1] = 16.0f0
    end
    for (a, b) in [(CN, DR), (CN, CS), (CS, AN), (SO, OI), (DR, BS)]
        J[a + 1, b + 1] = 4.0f0; J[b + 1, a + 1] = 4.0f0
    end
    for (a, b) in [(CN, SO), (CS, SO), (BS, AN)]
        J[a + 1, b + 1] = 20.0f0; J[b + 1, a + 1] = 20.0f0
    end
    return J
end

# ---- Host initialization (mirrors serial init_state; same RNG
#      stream ⇒ identical initial lattice for a given seed) ------

@inline radial_dist(x, y, N) = sqrt((x - N / 2.0)^2 + (y - N / 2.0)^2)

function place_cell!(lattice, cell_id, cx, cy, cz, radius, N)
    count = 0
    for dz in -radius:radius, dy in -radius:radius, dx in -radius:radius
        if dx^2 + dy^2 + dz^2 <= radius^2
            x = cx + dx; y = cy + dy; z = cz + dz
            if 1 <= x <= N && 1 <= y <= N && 1 <= z <= N &&
               lattice[x, y, z] == 0   # wall sites are -1, so this is the interior test
                lattice[x, y, z] = cell_id
                count += 1
            end
        end
    end
    return count
end

function init_host(N::Int, n_cells_per_species::Int, seed::Int, I0, κ, C_wall)
    rng = MersenneTwister(seed)
    lattice = zeros(Int32, N, N, N)
    R = N / 2.0
    for z in 1:N, y in 1:N, x in 1:N
        radial_dist(x, y, N) > R && (lattice[x, y, z] = Int32(-1))
    end

    spec = Int32[]
    vols = Int32[]
    next_id = 1
    for s in 1:N_SPECIES, _ in 1:n_cells_per_species
        placed = false; attempts = 0
        while !placed && attempts < 100
            cx = rand(rng, 5:(N - 4)); cy = rand(rng, 5:(N - 4)); cz = rand(rng, 5:(N - 4))
            if radial_dist(cx, cy, N) < R - 4
                radius = rand(rng, 2:3)
                vol = place_cell!(lattice, Int32(next_id), cx, cy, cz, radius, N)
                if vol > 0
                    push!(spec, Int32(s)); push!(vols, Int32(vol))
                    next_id += 1; placed = true
                end
            end
            attempts += 1
        end
    end

    rad = Array{Float32, 3}(undef, N, N, N)
    nut = Array{Float32, 3}(undef, N, N, N)
    for z in 1:N, y in 1:N, x in 1:N
        r = radial_dist(x, y, N)
        rad[x, y, z] = Float32(I0 * exp(-κ * r / N))
        nut[x, y, z] = Float32(C_wall * clamp(r / R, 0.0, 1.0))
    end
    mel = zeros(Float32, N, N, N)
    return lattice, spec, vols, rad, mel, nut
end

# ---- Device-side pure helpers --------------------------------

@inline function splitmix64(x::UInt64)
    x += 0x9e3779b97f4a7c15
    x = (x ⊻ (x >> 30)) * 0xbf58476d1ce4e5b9
    x = (x ⊻ (x >> 27)) * 0x94d049bb133111eb
    return x ⊻ (x >> 31)
end

@inline u01(r::UInt64) = Float32(r >> 40) * Float32(2.0^-24)

# n ∈ 0:25 → Moore-26 offset (base-3 decode, skipping the center 13)
@inline function nb26(n)
    m = n < 13 ? n : n + 1
    return (m % 3 - 1, (m ÷ 3) % 3 - 1, m ÷ 9 - 1)
end

# Species of a lattice value: 0 for medium/wall, else spec[σ]
@inline spidx(σ::Int32, spec) = σ > Int32(0) ? Int(spec[σ]) : 0

# Adhesion energy of site (x,y,z) assuming its spin is σc (pure —
# replaces the serial write-trial-spin-then-restore, which races)
@inline function site_adh(lat, spec, J, x, y, z, σc::Int32, N)
    spc = spidx(σc, spec)
    E = 0.0f0
    for dz in -1:1, dy in -1:1, dx in -1:1
        (dx == 0 && dy == 0 && dz == 0) && continue
        nx = x + dx; ny = y + dy; nz = z + dz
        if 1 <= nx <= N && 1 <= ny <= N && 1 <= nz <= N
            σn = lat[nx, ny, nz]
            if σn != σc
                E += J[spc + 1, spidx(σn, spec) + 1]
            end
        else
            E += J[spc + 1, 1]   # out of bounds = medium
        end
    end
    return E
end

# ΔH for copying σs into target (tx,ty,tz) currently holding σt.
# Mirrors serial compute_delta_H (adhesion + volume + radiation +
# melanin); volume term algebraically simplified: (V±1-Vt)²-(V-Vt)².
@inline function delta_H(lat, vols, spec, J, βv, melc, rad, mel,
        tx, ty, tz, σs::Int32, σt::Int32, N, λV::Float32, Vt::Int32)
    ΔH = site_adh(lat, spec, J, tx, ty, tz, σs, N) -
         site_adh(lat, spec, J, tx, ty, tz, σt, N)
    I = rad[tx, ty, tz]
    M = mel[tx, ty, tz]
    if σs > Int32(0)
        # ponytail: vols read is stale within a color pass (bounded,
        # unbiased); refresh per color if ensemble volumes ever drift
        d = Int(vols[σs]) - Int(Vt)
        sp = Int(spec[σs])
        ΔH += λV * Float32(2d + 1) + βv[sp] * I - melc[sp] * M
    end
    if σt > Int32(0)
        d = Int(vols[σt]) - Int(Vt)
        sp = Int(spec[σt])
        ΔH += λV * Float32(-2d + 1) - βv[sp] * I + melc[sp] * M
    end
    return ΔH
end

# ---- Kernels -------------------------------------------------

# One checkerboard color pass: thread (i,j,k) owns target site
# (2i-1+ox, 2j-1+oy, 2k-1+oz); all same-color targets are mutually
# non-adjacent under Moore-26, so lattice reads are race-free.
function cpm_color!(i, j, k, lat, vols, spec, J, βv, melc, rad, mel,
        ox, oy, oz, seed::UInt64, step::UInt64, N, λV::Float32,
        Vt::Int32, T::Float32)
    tx = 2 * (Int(i) - 1) + Int(ox) + 1
    ty = 2 * (Int(j) - 1) + Int(oy) + 1
    tz = 2 * (Int(k) - 1) + Int(oz) + 1
    σt = lat[tx, ty, tz]
    σt == Int32(-1) && return nothing

    lin = UInt64((tx - 1) + N * ((ty - 1) + N * (tz - 1)))
    r1 = splitmix64(splitmix64(seed ⊻ step) ⊻ lin)
    r2 = splitmix64(r1)

    dx, dy, dz = nb26(Int(r1 % UInt64(26)))
    sx = tx + dx; sy = ty + dy; sz = tz + dz
    (1 <= sx <= N && 1 <= sy <= N && 1 <= sz <= N) || return nothing
    σs = lat[sx, sy, sz]
    (σs == Int32(-1) || σs == σt || (σs <= 0 && σt <= 0)) && return nothing

    ΔH = delta_H(lat, vols, spec, J, βv, melc, rad, mel,
        tx, ty, tz, σs, σt, N, λV, Vt)
    if ΔH <= 0.0f0 || u01(r2) < exp(-ΔH / T)
        σt > Int32(0) && (JACC.@atomic vols[σt] -= Int32(1))
        σs > Int32(0) && (JACC.@atomic vols[σs] += Int32(1))
        lat[tx, ty, tz] = σs
    end
    return nothing
end

# 6-neighbor Laplacian, Neumann BC (out-of-bounds contributes 0)
@inline function lap6(F, x, y, z, N)
    v = F[x, y, z]
    L = 0.0f0
    x > 1 && (L += F[x - 1, y, z] - v); x < N && (L += F[x + 1, y, z] - v)
    y > 1 && (L += F[x, y - 1, z] - v); y < N && (L += F[x, y + 1, z] - v)
    z > 1 && (L += F[x, y, z - 1] - v); z < N && (L += F[x, y, z + 1] - v)
    return L
end

# ∂M/∂t = D_M ∇²M + α_M(species) · I_γ   (wall sites keep old value)
function melanin_k!(i, j, k, M, Mn, lat, spec, αv, rad, N, dt::Float32, D_M::Float32)
    σ = lat[i, j, k]
    if σ == Int32(-1)
        Mn[i, j, k] = M[i, j, k]
        return nothing
    end
    α = σ > Int32(0) ? αv[Int(spec[σ])] : 0.0f0
    Mn[i, j, k] = max(0.0f0, M[i, j, k] +
        dt * (D_M * lap6(M, i, j, k, N) + α * rad[i, j, k]))
    return nothing
end

# ∂C/∂t = D_C ∇²C - uptake(species)   (wall pinned to m·C_wall)
function nutrient_k!(i, j, k, C, Cn, lat, spec, upt, N, dt::Float32,
        D_C::Float32, C_wall_eff::Float32)
    σ = lat[i, j, k]
    if σ == Int32(-1)
        Cn[i, j, k] = C_wall_eff
        return nothing
    end
    u = σ > Int32(0) ? upt[Int(spec[σ])] : 0.0f0
    Cn[i, j, k] = max(0.0f0, C[i, j, k] +
        dt * (D_C * lap6(C, i, j, k, N) - u))
    return nothing
end

# ---- 1D radiolysis solver (host, copied from serial §12) -----

Base.@kwdef struct RadiolysisParams
    Nr::Int = 40
    D_eff::Float64 = 1e-3
    k_ads::Float64 = 0.05
    k_red::Float64 = 0.02
    k_des::Float64 = 0.005
    k_loss::Float64 = 0.001
    X_total::Float64 = 1.0
    X_red::Float64 = 0.3

    # Basis gate acknowledgement.  false refuses any X_total != 1.0.
    # Bool rather than a Symbol because this field is HDF5-serialised by
    # export_checkpoint.jl and must survive a restart; the exemption is binary
    # anyway, and WHICH sites hold it is pinned by the census test rather than
    # by a symbol name.  See _assert_basis_gate.
    basis_gate_ack::Bool = false
    P0::Float64 = 0.01
    alpha_P::Float64 = 0.02
    k_dam::Float64 = 0.005
    Ddot_R::Float64 = 1.0
    c_ext::Float64 = 1.0
    dt_rd::Float64 = 0.5
end

mutable struct RadiolysisState
    r_grid::Vector{Float64}
    c::Vector{Float64}
    s::Vector{Float64}
    m::Float64
    t::Float64
    params::RadiolysisParams
end

function init_radiolysis(rp::RadiolysisParams; R::Float64 = 1.0)
    r_grid = collect(range(0.0, R, length = rp.Nr))
    RadiolysisState(r_grid, zeros(rp.Nr), zeros(rp.Nr), 1.0, 0.0, rp)
end

"""
Advance the radiodialysis PDE system by one time step dt.

MIRRORS biofilms_potts.jl:1243 EXACTLY, and must continue to. Without this
wrapper the JACC port took a single raw forward-Euler step at dt_rd = 0.5 while
the serial port substepped, so the two would silently disagree the moment the
explicit-Euler diffusion bound was exceeded — and validate_serial.jl exists to
compare them.

Today both ports are saved by the same thing: every call site passes
R = N/2 in LATTICE units, giving dr = 0.513 and dt_stable = 52.6 at N = 40, so
n_sub = 1 and no trajectory moves. That safety is an artifact of the documented
r-in-lattice-units error (RADIODIALYSIS: BLOCKED), not of the time step being
small. Correct the units to R = 1.0 cm and dt_stable becomes 0.132, n_sub = 4 —
at which point a port without this guard diverges.
"""
function step_radiolysis!(rd::RadiolysisState, dt::Float64)
    _assert_basis_gate(rd.params.X_total, rd.params.basis_gate_ack)
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
function _assert_basis_gate(X_total::Float64, ack::Bool = false)
    X_total == 1.0 && return nothing
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

"""One forward-Euler step; the guard lives in step_radiolysis! above."""
function _step_radiolysis_euler!(rd::RadiolysisState, dt::Float64)
    rp = rd.params
    Nr = rp.Nr
    dr = rd.r_grid[2] - rd.r_grid[1]
    c = rd.c; s = rd.s; m = rd.m; t = rd.t
    P_eff = rp.P0 * exp(rp.alpha_P * rp.Ddot_R * t)
    dm_dt = -rp.k_dam * rp.Ddot_R * m
    uptake = rp.k_ads * rp.X_total + rp.k_red * rp.X_red
    dc = zeros(Nr); ds = zeros(Nr)
    dc[1] = rp.D_eff * 2.0 * (c[2] - c[1]) / dr^2 + (-uptake * c[1] + rp.k_des * s[1])
    for i in 2:(Nr - 1)
        r_i = rd.r_grid[i]
        diff_cyl = rp.D_eff * ((r_i + 0.5dr) * (c[i + 1] - c[i]) -
                               (r_i - 0.5dr) * (c[i] - c[i - 1])) / (r_i * dr^2)
        dc[i] = diff_cyl + (-uptake * c[i] + rp.k_des * s[i])
    end
    let i = Nr
        r_i = rd.r_grid[i]
        c_ghost = c[Nr - 1] - 2.0 * dr * P_eff * (c[Nr] - rp.c_ext) / rp.D_eff
        diff_cyl = rp.D_eff * ((r_i + 0.5dr) * (c_ghost - c[Nr]) -
                               (r_i - 0.5dr) * (c[Nr] - c[Nr - 1])) / (r_i * dr^2)
        dc[i] = diff_cyl + (-uptake * c[Nr] + rp.k_des * s[Nr])
    end
    @. ds = uptake * c - (rp.k_des + rp.k_loss) * s
    @. rd.c = max(0.0, c + dt * dc)
    @. rd.s = max(0.0, s + dt * ds)
    rd.m = max(0.0, m + dt * dm_dt)
    rd.t = t + dt
end

# Radial biomass from a host lattice copy (serial compute_radial_biomass)
function radial_biomass(lat_h::Array{Int32, 3}, spec_h::Vector{Int32}, N::Int, Nr::Int)
    dr = (N / 2.0) / (Nr - 1)
    total = zeros(Nr); red = zeros(Nr); counts = zeros(Int, Nr)
    for z in 1:N, y in 1:N, x in 1:N
        σ = lat_h[x, y, z]
        σ == Int32(-1) && continue
        idx = clamp(round(Int, radial_dist(x, y, N) / dr) + 1, 1, Nr)
        counts[idx] += 1
        if σ > 0
            total[idx] += 1.0
            spec_h[σ] == Int32(SO) && (red[idx] += 1.0)
        end
    end
    X_total = [counts[i] > 0 ? total[i] / counts[i] : 0.0 for i in 1:Nr]
    X_red = [counts[i] > 0 ? red[i] / counts[i] : 0.0 for i in 1:Nr]
    return X_total, X_red
end

# ---- Snapshot / reporting ------------------------------------

function species_stats(lat_h, mel_h, vols_h, spec_h)
    vol = zeros(Int, N_SPECIES); ncells = zeros(Int, N_SPECIES)
    for (id, v) in enumerate(vols_h)
        v > 0 || continue
        s = Int(spec_h[id])
        ncells[s] += 1; vol[s] += Int(v)
    end
    mel_sum = zeros(N_SPECIES); sites = zeros(Int, N_SPECIES)
    for idx in eachindex(lat_h)
        σ = lat_h[idx]
        if σ > 0 && vols_h[σ] > 0
            s = Int(spec_h[σ])
            mel_sum[s] += mel_h[idx]; sites[s] += 1
        end
    end
    mean_mel = [sites[s] > 0 ? mel_sum[s] / sites[s] : 0.0 for s in 1:N_SPECIES]
    return vol, ncells, mean_mel
end

function print_stats(mcs, vol, ncells, mean_mel, rd)
    @printf("  MCS %-5d  m=%.4f\n", mcs, rd.m)
    for s in 1:N_SPECIES
        @printf("    %-22s  vol=%6d  cells=%3d  mel=%8.4f\n",
            SPECIES_NAMES[s], vol[s], ncells[s], mean_mel[s])
    end
end

# ---- Driver --------------------------------------------------

function run_coupled(; N::Int = 40, n_cells_per_species::Int = 6,
        n_mcs::Int = 100, seed::Int = 42, snapshot_interval::Int = 20,
        T_cpm = 5.0f0, λ_V = 10.0f0, V_target = Int32(120),
        I0 = 1.0, κ = 2.0, D_M = 0.1f0, dt_field = 0.5f0,
        D_C = 0.2f0, C_wall = 1.0f0, rp = RadiolysisParams(), verbose = true)
    @assert iseven(N) "checkerboard requires even N"

    lat_h, spec_h, vols_h, rad_h, mel_h, nut_h =
        init_host(N, n_cells_per_species, seed, I0, κ, Float64(C_wall))

    lat = JACC.array(lat_h)
    vols = JACC.array(vols_h)
    spec = JACC.array(spec_h)
    J = JACC.array(build_J_matrix())
    βv = JACC.array(BETA_ION)
    αv = JACC.array(ALPHA_M)
    upt = JACC.array(UPTAKE)
    melc = JACC.array(MEL_COEF)
    rad = JACC.array(rad_h)
    mel = JACC.array(mel_h); mel2 = JACC.array(zeros(Float32, N, N, N))
    nut = JACC.array(nut_h); nut2 = JACC.array(zeros(Float32, N, N, N))

    rd = init_radiolysis(rp; R = N / 2.0)
    Nh = N ÷ 2
    gseed = splitmix64(UInt64(seed))

    for mcs in 1:n_mcs
        for c in 0:7
            JACC.parallel_for((Nh, Nh, Nh), cpm_color!,
                lat, vols, spec, J, βv, melc, rad, mel,
                Int32(c & 1), Int32((c >> 1) & 1), Int32((c >> 2) & 1),
                gseed, UInt64(mcs * 8 + c), Int32(N), λ_V, V_target, T_cpm)
        end

        if mcs % 10 == 1
            X_tot, X_red = radial_biomass(JACC.to_host(lat), spec_h, N, rd.params.Nr)
            rp0 = rd.params
            # Keyword form deliberately: this was positional, so adding a
            # field silently shifted every argument past it by one.
            rd.params = RadiolysisParams(
                Nr = rp0.Nr, D_eff = rp0.D_eff, k_ads = rp0.k_ads,
                k_red = rp0.k_red, k_des = rp0.k_des, k_loss = rp0.k_loss,
                X_total = mean(X_tot), X_red = mean(X_red),
                basis_gate_ack = rp0.basis_gate_ack,
                P0 = rp0.P0, alpha_P = rp0.alpha_P, k_dam = rp0.k_dam,
                Ddot_R = rp0.Ddot_R, c_ext = rp0.c_ext, dt_rd = rp0.dt_rd)
        end
        step_radiolysis!(rd, rd.params.dt_rd)

        JACC.parallel_for((N, N, N), melanin_k!,
            mel, mel2, lat, spec, αv, rad, Int32(N), dt_field, D_M)
        mel, mel2 = mel2, mel
        JACC.parallel_for((N, N, N), nutrient_k!,
            nut, nut2, lat, spec, upt, Int32(N), dt_field, D_C,
            Float32(rd.m) * C_wall)
        nut, nut2 = nut2, nut

        if verbose && (mcs % snapshot_interval == 0 || mcs == n_mcs)
            vol, ncells, mean_mel = species_stats(
                JACC.to_host(lat), JACC.to_host(mel), JACC.to_host(vols), spec_h)
            print_stats(mcs, vol, ncells, mean_mel, rd)
        end
    end

    lat_final = JACC.to_host(lat)
    mel_final = JACC.to_host(mel)
    vols_final = JACC.to_host(vols)
    return lat_final, mel_final, vols_final, spec_h, rd
end

function print_csv(seed, vol, ncells, mean_mel, rd)
    for s in 1:N_SPECIES
        @printf("CSV,%d,%d,%d,%d,%.5f,%d\n",
            seed, s, vol[s], ncells[s], mean_mel[s], ncells[s] > 0 ? 1 : 0)
    end
    @printf("CSVTOT,%d,%d,%.5f\n", seed, sum(ncells), rd.m)
end

function main(seeds::Vector{Int})
    println("JACC backend: $(JACC.backend)  array_type: $(JACC.array_type())")
    for seed in seeds
        t0 = time()
        lat_h, mel_h, vols_h, spec_h, rd = run_coupled(; seed)
        vol, ncells, mean_mel = species_stats(lat_h, mel_h, vols_h, spec_h)
        @printf("# seed %d done in %.1f s\n", seed, time() - t0)
        print_csv(seed, vol, ncells, mean_mel, rd)
    end
end

# ---- Selftest: ΔH and one PDE step vs the serial reference ---

# Load the serial script (minus its CairoMakie figure section) into a
# sandbox module so nothing collides. Must run at top level (world age).
function load_serial()
    src = read(joinpath(@__DIR__, "biofilms_potts.jl"), String)
    src = split(src, "#  13. Figure export")[1]
    SerialRef = Module(:SerialRef)
    Base.eval(SerialRef, :(using LinearAlgebra, Statistics, Random, Printf))
    Base.include_string(SerialRef, src, "biofilms_potts.jl")
    return SerialRef
end

function selftest(SerialRef::Module)
    p = SerialRef.CPMParams(N = 40, n_cells_per_species = 6)
    st = SerialRef.init_state(p; seed = 42)
    N = p.N

    # Matching port-side host state (same seed ⇒ same lattice)
    lat_h, spec_h, vols_h, rad_h, mel_h, nut_h =
        init_host(N, 6, 42, p.I0, p.κ, p.C_wall)
    @assert lat_h == st.lattice "initial lattices differ"
    @assert all(vols_h[id] == st.cells[id].volume for id in keys(st.cells))

    # Give melanin a nonzero profile so the ΔH melanin term is exercised
    for idx in eachindex(mel_h)
        mel_h[idx] = Float32(0.01 * (idx % 7))
        st.melanin[idx] = 0.01 * (idx % 7)
    end
    J32 = build_J_matrix()

    # 1. delta_H vs serial compute_delta_H over sampled site pairs
    rng = MersenneTwister(1)
    tested = 0
    while tested < 200
        sx, sy, sz = rand(rng, 2:(N - 1)), rand(rng, 2:(N - 1)), rand(rng, 2:(N - 1))
        dx, dy, dz = rand(rng, -1:1), rand(rng, -1:1), rand(rng, -1:1)
        (dx == 0 && dy == 0 && dz == 0) && continue
        tx, ty, tz = sx + dx, sy + dy, sz + dz
        σs = lat_h[sx, sy, sz]; σt = lat_h[tx, ty, tz]
        (σs == -1 || σt == -1 || σs == σt || (σs <= 0 && σt <= 0)) && continue
        want = SerialRef.compute_delta_H(st, sx, sy, sz, tx, ty, tz)
        got = delta_H(lat_h, vols_h, spec_h, J32, BETA_ION, MEL_COEF,
            rad_h, mel_h, tx, ty, tz, σs, σt, N, 10.0f0, Int32(120))
        @assert isapprox(got, want; rtol = 1e-4, atol = 1e-4) "ΔH mismatch: $got vs $want at ($sx,$sy,$sz)→($tx,$ty,$tz)"
        tested += 1
    end
    println("selftest: delta_H matches serial on $tested site pairs")

    # 2. One melanin + nutrient step vs serial (through JACC kernels)
    mel = JACC.array(mel_h); mel2 = JACC.array(zeros(Float32, N, N, N))
    nut = JACC.array(nut_h); nut2 = JACC.array(zeros(Float32, N, N, N))
    lat = JACC.array(lat_h)
    spec = JACC.array(spec_h)
    JACC.parallel_for((N, N, N), melanin_k!, mel, mel2, lat, spec,
        JACC.array(ALPHA_M), JACC.array(rad_h), Int32(N), 0.5f0, 0.1f0)
    JACC.parallel_for((N, N, N), nutrient_k!, nut, nut2, lat, spec,
        JACC.array(UPTAKE), Int32(N), 0.5f0, 0.2f0, Float32(0.9))

    copyto!(st.nutrient, nut_h)  # ensure identical starting field
    SerialRef.update_melanin!(st)
    SerialRef.update_nutrient_coupled!(st, 0.9)
    mel2_h = JACC.to_host(mel2); nut2_h = JACC.to_host(nut2)
    # serial melanin leaves wall values untouched; kernel copies them — same values
    @assert maximum(abs.(mel2_h .- Float32.(st.melanin))) < 1e-4 "melanin step mismatch"
    @assert maximum(abs.(nut2_h .- Float32.(st.nutrient))) < 1e-4 "nutrient step mismatch"
    println("selftest: melanin/nutrient kernels match serial one-step update")
    println("selftest: PASS")
end

if abspath(PROGRAM_FILE) == @__FILE__
    if "--selftest" in ARGS
        Base.invokelatest(selftest, load_serial())
    else
        seeds = isempty(ARGS) ? [42, 43, 44] : parse.(Int, ARGS)
        main(seeds)
    end
end
