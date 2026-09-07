module InertSignal

export SignalParams, update_signal!, endpoint, slab_halfwidth

"Declared lattice-unit parameters; production is per occupied voxel, not per parcel."
struct SignalParams
    diffusion::Float64
    decay::Float64
    production::NTuple{7,Float64}
    spacing::Float64
    max_dt::Float64
    threshold::Float64
    function SignalParams(; diffusion=0.2, decay=0.1, production=ones(7),
                          spacing=1.0, max_dt=0.5, threshold=5.0)
        length(production) == 7 || throw(ArgumentError("seven production rates required"))
        all(isfinite, (diffusion, decay, production..., spacing, max_dt, threshold)) ||
            throw(ArgumentError("parameters must be finite"))
        minimum((diffusion, decay, production..., threshold)) >= 0 ||
            throw(ArgumentError("rates and threshold must be nonnegative"))
        spacing > 0 && max_dt > 0 || throw(ArgumentError("spacing and max_dt must be positive"))
        new(diffusion, decay, Tuple(Float64.(production)), spacing, max_dt, threshold)
    end
end

"""
Advance an independent signal by `duration` MCS using fixed, read-only species/mask
arrays. Forward Euler, face-neighbour diffusion, zero Dirichlet at every masked
wall and no flux across outer cube faces. Positive coefficients are enforced by substepping;
there is no clipping. Nothing here holds a CPM state or uses an RNG.

The returned mass ledger integrates production, decay and washout over all
substeps. Source occupancy is held constant over this call, not reconstructed
as physical sub-MCS timing from the accepted-copy record.
"""
function update_signal!(A::Array{Float64,D}, species::AbstractArray{<:Integer,D},
                        mask::AbstractArray{Bool,D}, p::SignalParams; duration=1.0) where D
    size(A) == size(species) == size(mask) || throw(ArgumentError("field/label/mask shapes differ"))
    isfinite(duration) && duration >= 0 || throw(ArgumentError("invalid duration"))
    all(isfinite, A) && minimum(A) >= 0 || throw(ArgumentError("invalid initial signal"))
    all(iszero, A[.!mask]) || throw(ArgumentError("wall signal must be zero"))
    all(s -> 0 <= s <= 7, species[mask]) || throw(ArgumentError("unknown interior species"))
    all(s -> s <= 0, species[.!mask]) || throw(ArgumentError("occupied wall"))
    rate = 2D * p.diffusion / p.spacing^2 + p.decay
    dt_limit = min(p.max_dt, rate == 0 ? Inf : 0.9 / rate)
    nsteps = max(1, ceil(Int, duration / dt_limit))
    dt = duration / nsteps
    B = similar(A)
    source_mass = decay_mass = washout_mass = 0.0
    initial_mass = sum(A)
    offsets = ntuple(d -> CartesianIndex(ntuple(k -> k == d ? 1 : 0, D)), D)
    for _ in 1:nsteps
        fill!(B, 0)
        for I in CartesianIndices(A)
            mask[I] || continue
            a = A[I]
            neighbour_sum = 0.0
            wall_faces = 0
            for off in offsets, sign in (-1, 1)
                J = I + sign * off
                if !checkbounds(Bool, A, J)
                    neighbour_sum += a # omitted flux, like the nutrient stencil
                elseif mask[J]
                    neighbour_sum += A[J]
                else
                    wall_faces += 1
                end
            end
            sp = species[I]
            prod = sp == 0 ? 0.0 : p.production[sp]
            diff = p.diffusion / p.spacing^2 * (neighbour_sum - 2D * a)
            B[I] = a + dt * (diff + prod - p.decay * a)
            B[I] >= 0 || error("negative signal despite positivity substeps")
            source_mass += dt * prod
            decay_mass += dt * p.decay * a
            washout_mass += dt * p.diffusion / p.spacing^2 * wall_faces * a
        end
        copyto!(A, B)
    end
    (; initial_mass, final_mass=sum(A), source_mass, decay_mass, washout_mass, nsteps, dt)
end

"Threshold endpoint: >= threshold among currently occupied interior voxels. Empty denominator -> missing."
function endpoint(A, species, mask, p::SignalParams)
    occupied = mask .& (species .> 0)
    denominator = count(occupied)
    numerator = count(occupied .& (A .>= p.threshold))
    (; above_threshold=numerator, occupied_sites=denominator,
       fraction=denominator == 0 ? nothing : numerator / denominator)
end

"Steady uniformly producing slab, half-width L, A(+-L)=0; not a 3D radius."
function slab_halfwidth(diffusion, decay, production, threshold)
    all(isfinite, (diffusion, decay, production, threshold)) &&
        diffusion > 0 && decay >= 0 && production > 0 && threshold >= 0 ||
        throw(ArgumentError("invalid slab parameters"))
    threshold == 0 && return 0.0
    decay == 0 && return sqrt(2diffusion * threshold / production)
    decay * threshold >= production && return Inf
    sqrt(diffusion / decay) * acosh(inv(1 - decay * threshold / production))
end

end
