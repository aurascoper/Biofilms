"""Per-parcel face contacts; this reader neither changes the lattice nor runs a model.

The six faces define surface contact, not the model's 26-neighbor adhesion energy.
Same-ID faces never count. Walls, medium and the array exterior enter the boundary
denominator. Species averages weight parcels equally, not by volume or surface.
"""
function parcel_contacts(lattice::AbstractArray{<:Integer,3}, species_of::AbstractDict{Int,Int}, n_species::Int)
    n_species > 0 || error("n_species must be positive")
    faces = ((1,0,0), (-1,0,0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1))
    # Each parcel retains its own [volume, boundary faces, same-species faces].
    counts = Dict{Int,Vector{Int}}()
    for position in CartesianIndices(lattice)
        id = Int(lattice[position])
        id > 0 || continue
        haskey(species_of, id) || error("orphan parcel $id at $position")
        species = species_of[id]
        1 <= species <= n_species || error("invalid species $species for parcel $id")
        c = get!(counts, id, [0, 0, 0])
        c[1] += 1
        for delta in faces
            neighbor = position + CartesianIndex(delta)
            other = checkbounds(Bool, lattice, neighbor) ? Int(lattice[neighbor]) : 0
            other == id && continue
            c[2] += 1
            if other > 0
                haskey(species_of, other) || error("orphan neighbor parcel $other at $neighbor")
                1 <= species_of[other] <= n_species || error("invalid species for parcel $other")
                species_of[other] == species && (c[3] += 1)
            end
        end
    end
    return [begin
        ids = sort([id for id in keys(counts) if species_of[id] == species])
        volumes = [counts[id][1] for id in ids]
        boundary = [counts[id][2] for id in ids]
        same = [counts[id][3] for id in ids]
        all(>(0), boundary) || error("parcel with no boundary faces")
        contact = isempty(ids) ? NaN : sum(same ./ boundary) / length(ids)
        (; ids, volumes, boundary, same, contact)
    end for species in 1:n_species]
end
