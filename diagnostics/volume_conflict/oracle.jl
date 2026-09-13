# Deterministic host oracle for one checkerboard colour pass.
#
# WHY THIS EXISTS. `cpm_color!` reads `vols[σs]` and `vols[σt]` to evaluate the volume
# penalty, and atomically mutates the same entries on acceptance. The eight-colour
# decomposition separates target sites under the Moore-26 stencil, which makes the `lat`
# reads race-free -- but it does NOT make the donor and recipient PARCEL IDs disjoint.
# Two sites in one colour class can be spatially far apart and still belong to one cell.
# An atomic increment protects the increment; it does not make (read, read, evaluate,
# decide, update) one transaction.
#
# WHAT THIS IS NOT. It is not the random-sequential serial CPM in biofilms_potts.jl --
# that is a different update schedule and comparing against it would conflate two
# differences. This replays the SAME proposal stream the kernel generates, in a declared
# sequential order, so the only variable is when `vols` is read.
#
# NOT A PROPOSED FIX. `snapshot=true` evaluates every proposal against the sweep-start
# volumes. That is included as a DIAGNOSTIC ARM ONLY. It is not semantics-preserving:
# for H(V)=λ(V-Vt)^2 the joint change after two same-parcel growths exceeds the sum of
# two changes evaluated at the starting volume by exactly 2λ, with the opposite sign for
# a growth paired with a shrinkage. A worked case at λ=10, V=Vt=120, T=5, draw=0.01:
# sequential gives ΔH 10 (accept) then 30 (reject); snapshot gives 10 and 10 and accepts
# both. Removing the concurrent read by freezing changes which moves are accepted.

include(joinpath(@__DIR__, "..", "..", "biofilms_potts_jacc.jl"))

"""One proposal, recorded in full so every term can be compared arm to arm."""
struct Proposal
    tx::Int; ty::Int; tz::Int
    sx::Int; sy::Int; sz::Int
    donor::Int32          # σs, the label copied FROM (grows)
    recipient::Int32      # σt, the label copied INTO (shrinks)
    vol_donor::Int32      # vols[σs] as READ at evaluation time
    vol_recip::Int32      # vols[σt] as READ at evaluation time
    dH::Float32
    draw::Float32
    accepted::Bool
end

"""
    colour_pass!(lat, vols, ...; snapshot=false)

Replay one colour class in a declared sequential order. `snapshot=false` reads `vols`
live, which is the reference semantics. `snapshot=true` reads a sweep-start copy --
the diagnostic arm, not a fix.

Returns the proposals in evaluation order.
"""
function colour_pass!(lat, vols, spec, J, βv, melc, rad, mel,
                      ox, oy, oz, seed::UInt64, step::UInt64, N,
                      λV::Float32, Vt::Int32, T::Float32; snapshot::Bool=false)
    frozen = snapshot ? copy(vols) : vols
    out = Proposal[]
    half = div(N, 2)
    for k in 1:half, j in 1:half, i in 1:half      # declared order: x fastest
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

        vd = σs > 0 ? frozen[σs] : Int32(0)
        vr = σt > 0 ? frozen[σt] : Int32(0)
        ΔH = delta_H(lat, frozen, spec, J, βv, melc, rad, mel,
                     tx, ty, tz, σs, σt, N, λV, Vt)
        u = u01(r2)
        acc = ΔH <= 0.0f0 || u < exp(-ΔH / T)
        if acc
            σt > Int32(0) && (vols[σt] -= Int32(1))
            σs > Int32(0) && (vols[σs] += Int32(1))
            lat[tx, ty, tz] = σs
        end
        push!(out, Proposal(tx, ty, tz, sx, sy, sz, σs, σt, vd, vr, ΔH, u, acc))
    end
    return out
end

"""Independent bookkeeping check: recompute volumes from the lattice itself."""
function id_histogram(lat, nlabels)
    h = zeros(Int32, nlabels)
    for v in lat
        v > 0 && v <= nlabels && (h[v] += Int32(1))
    end
    return h
end

"""
    contention(props)

Every pair of proposals in one pass that touches a shared parcel id, on EITHER side.
Conflict detection must cover donor and recipient: a growth and a shrinkage of the same
label conflict just as two growths do.
"""
function contention(props::Vector{Proposal})
    pairs = Tuple{Int,Int,Int32,Symbol}[]
    for a in 1:length(props), b in (a+1):length(props)
        pa, pb = props[a], props[b]
        for (ida, rolea) in ((pa.donor, :donor), (pa.recipient, :recipient)),
            (idb, roleb) in ((pb.donor, :donor), (pb.recipient, :recipient))
            if ida > 0 && ida == idb
                push!(pairs, (a, b, ida, rolea === roleb ? rolea : :mixed))
            end
        end
    end
    return pairs
end
