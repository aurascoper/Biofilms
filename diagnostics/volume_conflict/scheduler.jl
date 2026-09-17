# A conflict-serializing scheduler, and the proof obligation it has to meet.
#
# THE ENABLING PROPERTY. Within one colour class the proposal SET is fully determined by
# the pass-start lattice, independently of any acceptance decision in that pass:
#
#   - target sites are the class members, fixed by (ox,oy,oz);
#   - the RNG is counter-based on (seed, step, linear index), so each site's draw and its
#     chosen neighbour do not depend on execution order;
#   - two distinct class members differ by >= 2 in some axis, and a site's source is within
#     +-1 of itself, so one member's source can never BE another member. The only lattice
#     write a proposal makes is to its own target site, which no other member reads.
#
# So (donor, recipient) for every proposal is knowable before evaluating any of them. The
# schedule can therefore be built up front. Only `vols` is contended.
#
# THE SCHEDULE. Layer proposals by the conflict DAG, in reference order:
#
#   batch(i) = 1 + max{ batch(j) : j < i, ids(j) ∩ ids(i) ≠ ∅ }      (max over ∅ = 0)
#
# ids(p) is BOTH positive endpoints, donor and recipient. A growth of parcel A and a
# shrinkage of parcel A conflict exactly as two growths do; the contention survey shows
# mixed-role pairs are roughly a third of all conflicts, so a one-sided rule would miss them.
#
# This gives two properties at once:
#   - within a batch no two proposals share a parcel id, so evaluating them against one
#     common `vols` is identical to evaluating them one at a time; and
#   - every conflicting pair keeps its reference relative order, because a conflicting
#     successor lands in a strictly later batch.
#
# The second is not optional. A merely conflict-free partition can reorder conflicting
# pairs and change the trajectory; being concurrency-safe is not the same as being
# semantics-preserving.
#
# The batch count is the longest chain of pairwise-conflicting proposals, i.e. the critical
# path of that DAG -- which is what bounds the available parallelism.

include(joinpath(@__DIR__, "oracle.jl"))

"""Proposal geometry and endpoints only -- everything knowable before evaluation."""
struct Site
    tx::Int; ty::Int; tz::Int
    sx::Int; sy::Int; sz::Int
    donor::Int32; recipient::Int32
    r2::UInt64
end

"""Enumerate a colour class's proposals from the pass-start lattice. No evaluation."""
function propose(lat, ox, oy, oz, seed::UInt64, step::UInt64, N)
    out = Site[]
    half = div(N, 2)
    for k in 1:half, j in 1:half, i in 1:half
        tx = 2 * (i - 1) + ox + 1
        ty = 2 * (j - 1) + oy + 1
        tz = 2 * (k - 1) + oz + 1
        (tx <= N && ty <= N && tz <= N) || continue
        σt = lat[tx, ty, tz]
        σt == Int32(-1) && continue
        lin = UInt64((tx - 1) + N * ((ty - 1) + N * (tz - 1)))
        r1 = splitmix64(splitmix64(seed ⊻ step) ⊻ lin)
        r2 = splitmix64(r1)
        dx, dy, dz = nb26(Int(r1 % UInt64(26)))
        sx = tx + dx; sy = ty + dy; sz = tz + dz
        (1 <= sx <= N && 1 <= sy <= N && 1 <= sz <= N) || continue
        σs = lat[sx, sy, sz]
        (σs == Int32(-1) || σs == σt || (σs <= 0 && σt <= 0)) && continue
        push!(out, Site(tx, ty, tz, sx, sy, sz, σs, σt, r2))
    end
    return out
end

ids(s::Site) = (s.donor > 0, s.recipient > 0) == (false, false) ? Int32[] :
               s.donor > 0 && s.recipient > 0 ? Int32[s.donor, s.recipient] :
               s.donor > 0 ? Int32[s.donor] : Int32[s.recipient]

"""Layer proposals into conflict-free batches preserving the reference order."""
function schedule(sites::Vector{Site})
    batch = zeros(Int, length(sites))
    last_of = Dict{Int32,Int}()
    for (i, s) in enumerate(sites)
        b = 0
        for id in ids(s)
            b = max(b, get(last_of, id, 0))
        end
        batch[i] = b + 1
        for id in ids(s)
            last_of[id] = batch[i]
        end
    end
    return batch
end

"""
Execute one colour pass batch by batch. Every proposal in a batch reads the SAME `vols`,
which is what a concurrent launch would do -- and is safe here precisely because the batch
is id-disjoint. Updates are applied after the batch, in reference order.
"""
function scheduled_pass!(lat, vols, spec, J, βv, melc, rad, mel,
                         ox, oy, oz, seed::UInt64, step::UInt64, N,
                         λV::Float32, Vt::Int32, T::Float32)
    sites = propose(lat, ox, oy, oz, seed, step, N)
    batch = schedule(sites)
    nb = isempty(batch) ? 0 : maximum(batch)
    out = Proposal[]
    order = sortperm(batch)          # stable: reference order preserved within a batch
    frozen = copy(vols)
    cur = 0
    for idx in order
        s = sites[idx]
        if batch[idx] != cur         # new batch: publish the previous batch's updates
            frozen = copy(vols)
            cur = batch[idx]
        end
        vd = s.donor > 0 ? frozen[s.donor] : Int32(0)
        vr = s.recipient > 0 ? frozen[s.recipient] : Int32(0)
        ΔH = delta_H(lat, frozen, spec, J, βv, melc, rad, mel,
                     s.tx, s.ty, s.tz, s.donor, s.recipient, N, λV, Vt)
        u = u01(s.r2)
        acc = ΔH <= 0.0f0 || u < exp(-ΔH / T)
        if acc
            s.recipient > Int32(0) && (vols[s.recipient] -= Int32(1))
            s.donor > Int32(0) && (vols[s.donor] += Int32(1))
            lat[s.tx, s.ty, s.tz] = s.donor
        end
        push!(out, Proposal(s.tx, s.ty, s.tz, s.sx, s.sy, s.sz,
                            s.donor, s.recipient, vd, vr, ΔH, u, acc))
    end
    return out, nb, length(sites)
end
