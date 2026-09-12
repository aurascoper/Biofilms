"""
    BindingBenchmark

A closed diffusion-and-binding benchmark on a frozen lattice geometry.

This module makes **no chemistry and no isotope claim**. It integrates

    B(x)  = O(x) · [B0 + ΔB · h(A(x))],   h(A) = Aⁿ / (Kⁿ + Aⁿ)
    r     = k_on · c · (B − b) − k_off · b
    ∂c/∂t = D_c ∇²c − r − λ c + q_ext
    ∂b/∂t = r − λ b

where `c` is a dissolved scalar, `b` a bound scalar, `B` a finite volumetric capacity
carried by occupied sites, and `A` a frozen scalar read off the geometry. A *closed*
benchmark has `q_ext = 0` and no flux through the domain boundary, so the only thing
that leaves the ledger is decay — and decay is accumulated, not inferred.

Units are **declared diagnostic units**: lengths in lattice sites, time in `dtu`, rates
per `dtu`, `c` and `b` in an arbitrary concentration unit. Nothing here converts to
seconds or metres; that conversion is blocked on D-PITCH and D-TIMESERIES, and the
half-life illustration's days are a different diagnostic's units, not these.

Naming follows the specification deliberately: `lambda` is the loss rate applied to both
pools, `B`/`B0`/`dB` are capacity (never a loss rate), `q_ext` is the external source
(never the binding relaxation rate `q_bind = k_on·c + k_off + λ`).
"""
module BindingBenchmark

using HDF5

export Geometry, Params, State, Ledger
export load_geometry, neighbour_fraction, hill, capacity_field, matched_uniform_capacity
export make_state, diffusion_rate, stability_rate, laplacian!, step!, run_to
export inventory, closure_residual, bound_fractions

# One adjacency rule is used for everything in this diagnostic: 6-connected faces.
# The diffusion stencil, the capacity modulator `A`, and the no-flux boundary all use
# it, so there is no second rule to confuse with the first.
const FACE_OFFSETS = ((1, 0, 0), (0, 1, 0), (0, 0, 1))

# ---------------------------------------------------------------------------
# Frozen geometry
# ---------------------------------------------------------------------------

"""
Labels and mask read from one frozen snapshot. Nothing here is advanced in time: the
benchmark holds the geometry fixed and moves only `c` and `b`.
"""
struct Geometry
    occupied::BitArray{3}
    interior::BitArray{3}
    spacing::NTuple{3, Float64}
    mcs::Int
    run_id::String
    label_state_hash::String
    mask_sha256::String
end

"""
    load_geometry(path; spacing = (1.0, 1.0, 1.0))

Read `lattice/cell_id` and `lattice/interior_mask` from a producer snapshot.

Occupancy is `cell_id > 0`: the background sentinel (0) and the wall sentinel (-1) are
both unoccupied, and conflating the wall with occupied medium would put capacity inside
the wall. The snapshot's own sentinel declarations are read and checked rather than
assumed, and occupancy is required to lie inside the interior — an occupied site outside
the mask would carry capacity that the transport operator never visits.
"""
function load_geometry(path::AbstractString; spacing::NTuple{3, Float64} = (1.0, 1.0, 1.0))
    all(>(0), spacing) || throw(ArgumentError("spacing must be positive"))
    h5open(path, "r") do f
        a = attributes(f)
        background = Int(read(a["cell_id_background"]))
        wall = Int(read(a["cell_id_wall"]))
        (background == 0 && wall == -1) ||
            throw(ArgumentError("unexpected cell_id sentinels: background=$background wall=$wall"))
        read(a["logical_axis_order"]) == "xyz" ||
            throw(ArgumentError("snapshot is not in xyz logical axis order"))
        cell_id = read(f["lattice/cell_id"])
        mask = read(f["lattice/interior_mask"])
        size(cell_id) == size(mask) ||
            throw(ArgumentError("cell_id and interior_mask disagree on shape"))
        occupied = BitArray(cell_id .> 0)
        interior = BitArray(mask .== 1)
        any(occupied .& .!interior) &&
            throw(ArgumentError("occupied sites outside the interior mask"))
        Geometry(occupied, interior, spacing, Int(read(a["mcs"])),
                 String(read(a["run_id"])), String(read(a["label_state_hash"])),
                 String(read(a["mask_sha256"])))
    end
end

"""
    neighbour_fraction(occupied)

`A(x)` — the fraction of `x`'s six face neighbours that are occupied, in `[0, 1]`.

This is the capacity modulator, and it is derived from the frozen labels rather than
imported from any field, so the benchmark reads exactly one thing off the geometry.
Neighbours outside the grid count as unoccupied, which is the same convention the
no-flux boundary uses.
"""
function neighbour_fraction(occupied::BitArray{3})
    nx, ny, nz = size(occupied)
    A = zeros(Float64, nx, ny, nz)
    @inbounds for (di, dj, dk) in FACE_OFFSETS, k in 1:nz, j in 1:ny, i in 1:nx
        i2, j2, k2 = i + di, j + dj, k + dk
        (i2 <= nx && j2 <= ny && k2 <= nz) || continue
        occupied[i2, j2, k2] && (A[i, j, k] += 1 / 6)
        occupied[i, j, k] && (A[i2, j2, k2] += 1 / 6)
    end
    A
end

"Hill response `Aⁿ / (Kⁿ + Aⁿ)`, zero at `A = 0` for every positive `n`."
hill(A::Real, K::Real, n::Real) = A <= 0 ? 0.0 : A^n / (K^n + A^n)

# ---------------------------------------------------------------------------
# Parameters and capacity
# ---------------------------------------------------------------------------

"""
Declared benchmark constants. These are *chosen* to exercise the scheme, not measured
from anything, and the receipt says so. `q_ext = 0` is what makes the benchmark closed.
"""
struct Params
    D_c::Float64      # lattice² / dtu
    lambda::Float64   # 1 / dtu, applied to both pools
    k_on::Float64     # 1 / (concentration · dtu)
    k_off::Float64    # 1 / dtu
    B0::Float64       # capacity floor on an occupied site
    dB::Float64       # capacity added at full modulation
    K::Float64        # Hill half-saturation, in units of A
    n::Float64        # Hill exponent
    q_ext::Float64    # external source; zero for a closed benchmark
end

"""
    capacity_field(geo, p; A = neighbour_fraction(geo.occupied))

`B = O · [B0 + ΔB · h(A)]`. Capacity lives only on occupied sites and is *volumetric* —
one capacity per site, independent of how much of the site any label claims.

Refuses a modulator that is constant over the occupied set: with a constant `A` the
matched-total-capacity control compares a field against itself and cannot fail, so a
degenerate input would silently turn the benchmark's central control into a tautology.
"""
function capacity_field(geo::Geometry, p::Params; A = neighbour_fraction(geo.occupied))
    size(A) == size(geo.occupied) || throw(ArgumentError("modulator shape mismatch"))
    p.dB == 0 || begin
        vals = A[geo.occupied]
        isempty(vals) && throw(ArgumentError("no occupied sites"))
        maximum(vals) > minimum(vals) ||
            throw(ArgumentError("capacity modulator is constant over the occupied set; " *
                                "the matched-capacity control would be vacuous"))
    end
    B = zeros(Float64, size(A))
    @inbounds for ix in eachindex(B)
        geo.occupied[ix] && (B[ix] = p.B0 + p.dB * hill(A[ix], p.K, p.n))
    end
    B
end

"""
    matched_uniform_capacity(B, occupied)

The uniform capacity `B̄` carrying the *same integrated capacity* as `B`, so that a
difference between the two runs is attributable to the spatial arrangement of capacity
and not to there being more of it. Returns the field and `B̄`.
"""
function matched_uniform_capacity(B::Array{Float64, 3}, occupied::BitArray{3})
    n = count(occupied)
    n > 0 || throw(ArgumentError("no occupied sites"))
    total = 0.0
    @inbounds for ix in eachindex(B)
        occupied[ix] && (total += B[ix])
    end
    Bbar = total / n
    U = zeros(Float64, size(B))
    @inbounds for ix in eachindex(U)
        occupied[ix] && (U[ix] = Bbar)
    end
    U, Bbar
end

# ---------------------------------------------------------------------------
# State, ledger
# ---------------------------------------------------------------------------

mutable struct State
    c::Array{Float64, 3}
    b::Array{Float64, 3}
    B::Array{Float64, 3}
    lap::Array{Float64, 3}
end

"""
    make_state(geo, B; c0, b0)

Uniform initial pools over the interior. `b0` is clipped to the local capacity, and a
`b0` that would exceed capacity anywhere is refused rather than quietly clipped.
"""
function make_state(geo::Geometry, B::Array{Float64, 3}; c0::Float64, b0::Float64 = 0.0)
    c = zeros(Float64, size(B))
    b = zeros(Float64, size(B))
    @inbounds for ix in eachindex(c)
        geo.interior[ix] || continue
        c[ix] = c0
        # The capacity guard is applied where the pool is placed. Checking `b0 <= B` on
        # every interior site instead would refuse any positive `b0`, because capacity is
        # zero off the occupied set by construction.
        geo.occupied[ix] || continue
        b0 <= B[ix] || throw(ArgumentError(
            "initial bound pool $b0 exceeds capacity $(B[ix]) at " *
            "$(Tuple(CartesianIndices(B)[ix]))"))
        b[ix] = b0
    end
    State(c, b, copy(B), zeros(Float64, size(B)))
end

"""
Every term is accumulated from the term itself as the step computes it. Nothing here is
obtained by differencing the fields, so the closure identity has two independent sides
and a mistake in the update shows up as a residual instead of cancelling.

`dropped` is non-zero only under deliberate defect injection and is deliberately absent
from [`closure_residual`](@ref): a ledger that books the material a bug destroys is a
ledger that cannot report the bug.
"""
mutable struct Ledger
    t::Float64
    steps::Int
    initial_dissolved::Float64
    initial_bound::Float64
    decayed_from_dissolved::Float64
    decayed_from_bound::Float64
    bound_transfer::Float64
    released::Float64
    external_input::Float64
    chemostat_input::Float64
    transport_residual::Float64
    dropped::Float64
    min_dissolved::Float64
    min_bound::Float64
    max_overshoot::Float64
end

"""
The running extrema are seeded from the initial state rather than from sentinels. An
unseeded ledger reports `Inf` for a run of zero steps, and — worse — excludes the initial
condition from the positivity record, so an initial state that already violated a bound
would be reported clean.
"""
function Ledger(st::State, geo::Geometry)
    Id, Ib, _ = inventory(st, geo)
    cmin = Inf; bmin = Inf; over = -Inf
    @inbounds for ix in eachindex(st.c)
        geo.interior[ix] || continue
        cmin = min(cmin, st.c[ix])
        bmin = min(bmin, st.b[ix])
        over = max(over, st.b[ix] - st.B[ix])
    end
    Ledger(0.0, 0, Id, Ib, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, cmin, bmin, over)
end

"Dissolved inventory, bound inventory, and integrated capacity, over the interior."
function inventory(st::State, geo::Geometry)
    V = prod(geo.spacing)
    sc = 0.0; sb = 0.0; sB = 0.0
    @inbounds for ix in eachindex(st.c)
        geo.interior[ix] || continue
        sc += st.c[ix]; sb += st.b[ix]; sB += st.B[ix]
    end
    (sc * V, sb * V, sB * V)
end

"""
    closure_residual(st, geo, led)

`(I_c + I_b) − [I_c(0) + I_b(0) + inputs − decayed]`.

`transport_residual` is **not** a term here. Diffusion is an internal transfer and must
cancel exactly; folding its floating-point residue into the identity would let a real
boundary leak close the books.
"""
function closure_residual(st::State, geo::Geometry, led::Ledger)
    Id, Ib, _ = inventory(st, geo)
    expected = led.initial_dissolved + led.initial_bound +
               led.external_input + led.chemostat_input -
               led.decayed_from_dissolved - led.decayed_from_bound
    (Id + Ib) - expected
end

"""
    bound_fractions(st, geo)

Both denominators, separately, because they answer different questions and a single
"bound fraction" hides which one was meant: `of_inventory = I_b / (I_c + I_b)` is how
much of what is left is bound; `of_capacity = I_b / I_B` is how full the capacity is.
"""
function bound_fractions(st::State, geo::Geometry)
    Id, Ib, IB = inventory(st, geo)
    ((Id + Ib) > 0 ? Ib / (Id + Ib) : 0.0, IB > 0 ? Ib / IB : 0.0)
end

# ---------------------------------------------------------------------------
# Transport, stability, the step
# ---------------------------------------------------------------------------

"""
    laplacian!(out, c, interior, D, spacing)

`D ∇²c` in **flux form**: each interior–interior face is visited once and contributes
`+f` to one cell and `−f` to the other. A face touching a non-interior site contributes
nothing, which is the discrete no-flux condition.

Written this way the operator is antisymmetric by construction, so it conserves the
dissolved inventory on a ragged boundary as exactly as floating-point addition allows.
A site-centred stencil with boundary special-cases does not: it leaks at precisely the
sites where the mask is ragged, and on this geometry that is 6,766 of 50,200 sites.
"""
function laplacian!(out::Array{Float64, 3}, c::Array{Float64, 3},
                    interior::BitArray{3}, D::Float64, spacing::NTuple{3, Float64})
    fill!(out, 0.0)
    nx, ny, nz = size(c)
    @inbounds for (d, (di, dj, dk)) in enumerate(FACE_OFFSETS)
        coef = D / spacing[d]^2
        coef == 0 && continue
        for k in 1:nz, j in 1:ny, i in 1:nx
            i2, j2, k2 = i + di, j + dj, k + dk
            (i2 <= nx && j2 <= ny && k2 <= nz) || continue
            (interior[i, j, k] && interior[i2, j2, k2]) || continue
            f = coef * (c[i2, j2, k2] - c[i, j, k])
            out[i, j, k] += f
            out[i2, j2, k2] -= f
        end
    end
    out
end

"""
    diffusion_rate(D, spacing)

`2 D Σ_j h_j⁻²` — the magnitude of the diffusion operator's diagonal at a site with all
six faces interior, and therefore the explicit-Euler restriction `Δt ≤ 1 / (2 D Σ_j h_j⁻²)`.

Derived from *this* stencil, not imported. In 3-D with uniform `h` that is `h²/(6D)`,
not the one-dimensional `h²/(2D)`: at `h = D = 1, Δt = 0.25` a unit centre with six zero
neighbours goes to `1 − 0.25·6 = −0.5`, which the 1-D bound admits and this one refuses.
"""
diffusion_rate(D::Float64, spacing::NTuple{3, Float64}) = 2 * D * sum(1 ./ (spacing .^ 2))

"""
    stability_rate(st, geo, p)

The largest negative diagonal of the unsplit right-hand side, over both equations:

    dissolved: 2 D Σ_j h_j⁻² + k_on · max(B − b) + λ
    bound:     k_on · max(c) + k_off + λ

`Δt · stability_rate ≤ 1` is the restriction the implemented scheme actually needs when
reaction and decay are taken together with transport rather than split off. It is
evaluated against the current state on every step, because `max(c)` and `max(B − b)`
both move, so a run can leave the admissible region after starting inside it.
"""
function stability_rate(st::State, geo::Geometry, p::Params)
    free = 0.0; cmax = 0.0
    @inbounds for ix in eachindex(st.c)
        geo.interior[ix] || continue
        free = max(free, st.B[ix] - st.b[ix])
        cmax = max(cmax, st.c[ix])
    end
    max(diffusion_rate(p.D_c, geo.spacing) + p.k_on * free + p.lambda,
        p.k_on * cmax + p.k_off + p.lambda)
end

"""
    step!(st, geo, p, dt, led; release_to_dissolved = true, clamp_dissolved = false)

One explicit forward-Euler step of the unsplit system, then the capacity-overflow
release, then the optional dissolved clamp.

`release_to_dissolved = false` is a **deliberate defect**: the overflow is removed from
`b` and credited nowhere, so [`closure_residual`](@ref) goes negative by exactly the
amount lost. It exists so the conservative-release control has something to be red
against, and the production run never sets it.

`clamp_dissolved = true` is **a different scheme**, not a diagnostic view of this one: it
holds `c` at its pre-step value, which makes the system open. The material the clamp puts
back is booked as `chemostat_input` so the ledger still closes, and only the analytical
limit control uses it.
"""
function step!(st::State, geo::Geometry, p::Params, dt::Float64, led::Ledger;
               release_to_dissolved::Bool = true, clamp_dissolved::Bool = false)
    dt > 0 || throw(ArgumentError("dt must be positive"))
    rate = stability_rate(st, geo, p)
    dt * rate <= 1 || throw(ErrorException(
        "forward-Euler restriction violated at t = $(led.t): dt·rate = $(dt * rate) > 1 " *
        "(dt = $dt, rate = $rate); reduce dt below $(1 / rate)"))

    csave = clamp_dissolved ? copy(st.c) : st.c
    laplacian!(st.lap, st.c, geo.interior, p.D_c, geo.spacing)

    V = prod(geo.spacing)
    dec_c = 0.0; dec_b = 0.0; trans = 0.0; inp = 0.0; lapsum = 0.0; over = -Inf
    @inbounds for ix in eachindex(st.c)
        geo.interior[ix] || continue
        c = st.c[ix]; b = st.b[ix]; B = st.B[ix]
        r = p.k_on * c * (B - b) - p.k_off * b
        dec_c += p.lambda * c
        dec_b += p.lambda * b
        trans += r
        inp += p.q_ext
        lapsum += st.lap[ix]
        st.c[ix] = c + dt * (st.lap[ix] - r - p.lambda * c + p.q_ext)
        st.b[ix] = b + dt * (r - p.lambda * b)
        over = max(over, st.b[ix] - B)
    end
    led.decayed_from_dissolved += dec_c * dt * V
    led.decayed_from_bound += dec_b * dt * V
    led.bound_transfer += trans * dt * V
    led.external_input += inp * dt * V
    led.transport_residual += lapsum * dt * V
    led.max_overshoot = max(led.max_overshoot, over)

    rel = 0.0; drop = 0.0
    @inbounds for ix in eachindex(st.b)
        geo.interior[ix] || continue
        ex = st.b[ix] - st.B[ix]
        ex > 0 || continue
        st.b[ix] = st.B[ix]
        if release_to_dissolved
            st.c[ix] += ex
            rel += ex
        else
            drop += ex
        end
    end
    led.released += rel * V
    led.dropped += drop * V

    if clamp_dissolved
        flux = 0.0
        @inbounds for ix in eachindex(st.c)
            geo.interior[ix] || continue
            flux += csave[ix] - st.c[ix]
            st.c[ix] = csave[ix]
        end
        led.chemostat_input += flux * V
    end

    cmin = Inf; bmin = Inf
    @inbounds for ix in eachindex(st.c)
        geo.interior[ix] || continue
        cmin = min(cmin, st.c[ix]); bmin = min(bmin, st.b[ix])
    end
    led.min_dissolved = min(led.min_dissolved, cmin)
    led.min_bound = min(led.min_bound, bmin)
    led.t += dt
    led.steps += 1
    st
end

"""
    run_to(st, geo, p, dt, nsteps, led; observe = nothing, kwargs...)

`nsteps` steps of [`step!`](@ref). `observe(step, led, st)` is called before the first
step and after every step, so a caller records from the integrator's own state rather
than re-deriving it.
"""
function run_to(st::State, geo::Geometry, p::Params, dt::Float64, nsteps::Integer,
                led::Ledger; observe = nothing, kwargs...)
    isnothing(observe) || observe(0, led, st)
    for s in 1:nsteps
        step!(st, geo, p, dt, led; kwargs...)
        isnothing(observe) || observe(s, led, st)
    end
    st
end

end # module
