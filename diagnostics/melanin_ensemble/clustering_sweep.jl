#!/usr/bin/env julia
# Atomic, restartable shards of the registered full/no_mel contact experiment.
# julia clustering_sweep.jl /external/contact_42-57.csv --seeds 42:57
using Printf, SHA

module ClusterRadial
include(joinpath(@__DIR__, "radial_sweep.jl"))
end
include(joinpath(@__DIR__, "clustering_observable.jl"))

const CLUSTER_REPO = normpath(joinpath(@__DIR__, "..", ".."))
const CLUSTER_REFERENCE = joinpath(@__DIR__, "radial_n40_p6_seeds42-297.csv")
const CLUSTER_PREDICTION = joinpath(@__DIR__, "PREDICTION_clustering_n40_p6.md")
const CLUSTER_ARMS = ("full", "no_mel")
const CLUSTER_COLUMNS = ("seed", "arm", "mcs", "species", "species_name", "ncells", "volume",
    "mean_r", "mean_r_over_R", "mean_melanin", "mean_contact_fraction", "parcel_ids",
    "parcel_volumes", "parcel_boundary_faces", "parcel_same_faces", "placement_sha256")
const BASELINE_FIELDS = ("ncells", "volume", "mean_r", "mean_r_over_R", "mean_melanin")

"Read the producer's plain CSV dialect, with explicit metadata and unique row keys."
function cluster_table(path)
    meta = Dict{String,String}()
    rows = Dict{Tuple{Int,String,Int,Int},Dict{String,String}}()
    columns = String[]
    for line in eachline(path)
        isempty(strip(line)) && continue
        if startswith(line, "#")
            for m in eachmatch(r"(\w+)=(\S+)", line)
                key, value = String(m[1]), String(m[2])
                haskey(meta, key) && meta[key] != value && error("$path: conflicting metadata $key")
                meta[key] = value
            end
            continue
        end
        fields = String.(split(line, ","))
        if isempty(columns)
            columns = fields
            allunique(columns) || error("$path: duplicate column")
            all(in(columns), ("seed", "arm", "mcs", "species")) || error("$path: incomplete header")
            continue
        end
        length(fields) == length(columns) || error("$path: row has wrong field count")
        row = Dict(zip(columns, fields))
        key = (parse(Int, row["seed"]), row["arm"], parse(Int, row["mcs"]), parse(Int, row["species"]))
        haskey(rows, key) && error("$path: duplicate row $key")
        rows[key] = row
    end
    isempty(rows) && error("$path: no data")
    (; meta, rows, columns)
end

function check_cluster_baseline(row, reference)
    key = (parse(Int, row["seed"]), row["arm"], parse(Int, row["mcs"]), parse(Int, row["species"]))
    haskey(reference.rows, key) || error("baseline lacks $key")
    for field in BASELINE_FIELDS
        row[field] == reference.rows[key][field] || error(
            "baseline mismatch $key $field: $(row[field]) versus $(reference.rows[key][field])")
    end
end

function cluster_source_record(smoke, reference)
    git(args...) = strip(read(`git -C $CLUSTER_REPO $args`, String))
    dirty = !isempty(git("status", "--porcelain", "--untracked-files=no"))
    !smoke && dirty && error("production requires a clean committed preregistration and source tree")
    record = Dict("source_commit" => git("rev-parse", "HEAD"),
                  "source_state" => dirty ? "dirty" : "clean",
                  "mode" => smoke ? "smoke" : "production",
                  "julia_version" => string(VERSION), "machine" => string(Sys.MACHINE))
    for (label, path) in (
        "preregistration" => CLUSTER_PREDICTION, "producer" => @__FILE__,
        "observable" => joinpath(@__DIR__, "clustering_observable.jl"),
        "analysis" => joinpath(@__DIR__, "clustering_analyse.jl"),
        "tests" => joinpath(@__DIR__, "test_clustering.jl"),
        "radial_sweep" => joinpath(@__DIR__, "radial_sweep.jl"),
        "model" => joinpath(CLUSTER_REPO, "biofilms_potts.jl"), "reference" => reference)
        record[label * "_sha256"] = bytes2hex(sha256(read(path)))
    end
    # An untracked preregistration must not masquerade as a committed protocol.
    if !smoke
        for path in (CLUSTER_PREDICTION, @__FILE__, joinpath(@__DIR__, "clustering_observable.jl"),
                     joinpath(@__DIR__, "clustering_analyse.jl"), joinpath(@__DIR__, "test_clustering.jl"))
            committed = read(`git -C $CLUSTER_REPO show HEAD:$(relpath(path, CLUSTER_REPO))`)
            committed == read(path) || error("uncommitted study source $path")
        end
    end
    record
end

function placement_digest(state)
    io = IOBuffer()
    print(io, size(state.lattice), ':', join(vec(state.lattice), ','), ';')
    for id in sort(collect(keys(state.cells)))
        print(io, id, ':', state.cells[id].species, ';')
    end
    bytes2hex(sha256(take!(io)))
end

function clustering_sweep(args)
    isempty(args) && error("usage: clustering_sweep.jl /external/out.csv --seeds 42:57 [--reference radial.csv]")
    out = abspath(args[1])
    startswith(out, CLUSTER_REPO * "/") && error("write shards outside the source checkout")
    ispath(out) && error("destination already exists: $out")
    partial = out * ".partial"
    ispath(partial) && error("incomplete prior shard exists: $partial; inspect/remove it before retrying")
    allowed = Set(("--seeds", "--n", "--parcels", "--at", "--reference"))
    i = 2
    while i <= length(args)
        if args[i] == "--smoke"
            i += 1
        else
            args[i] in allowed && i < length(args) || error("unknown or incomplete option $(args[i])")
            i += 2
        end
    end
    getarg(flag, default) = ClusterRadial.getopt(args, flag, default)
    seeds = ClusterRadial.parse_seeds(getarg("--seeds", "42:297"))
    N = parse(Int, getarg("--n", "40"))
    parcels = parse(Int, getarg("--parcels", "6"))
    at = parse(Int, getarg("--at", "100"))
    smoke = "--smoke" in args
    reference_path = abspath(getarg("--reference", CLUSTER_REFERENCE))
    !isempty(seeds) && allunique(seeds) || error("seeds must be nonempty and unique")
    if smoke
        10 <= N <= 16 && parcels == 1 && 1 <= at <= 4 && all(s -> !(42 <= s <= 297), seeds) ||
            error("smoke requires N 10:16, one parcel, 1:4 MCS and seeds outside 42:297")
    else
        (N, parcels, at) == (40, 6, 100) && all(s -> 42 <= s <= 297, seeds) ||
            error("production configuration is fixed at N40,p6,MCS100,seeds42:297")
        reference_path == CLUSTER_REFERENCE || error("production must use the committed radial reference")
    end
    reference = cluster_table(reference_path)
    (reference.meta["N"], reference.meta["parcels_per_species"], reference.meta["read_at_mcs"]) ==
        string.((N, parcels, at)) || error("reference configuration differs")
    source_record = cluster_source_record(smoke, reference_path)
    SR = ClusterRadial.load_serial()
    params = Dict(arm => ClusterRadial.arm_params(SR, arm, N, parcels, at) for arm in CLUSTER_ARMS)
    checked = 0
    mkpath(dirname(out))
    open(partial, "w") do io
        println(io, "# clustering contact sweep")
        println(io, "# N=$N parcels_per_species=$parcels read_at_mcs=$at arms=full,no_mel")
        println(io, "# declared_seeds=", join(seeds, ';'))
        for key in sort(collect(keys(source_record)))
            println(io, "# $key=", source_record[key])
        end
        println(io, join(CLUSTER_COLUMNS, ','))
        for seed in seeds
            placement = nothing
            started = time()
            for arm in CLUSTER_ARMS
                initial, _ = redirect_stdout(devnull) do
                    Base.invokelatest(SR.run_simulation, params[arm], 0; seed=seed)
                end
                digest = placement_digest(initial)
                if isnothing(placement)
                    placement = digest
                else
                    digest == placement || error("seed $seed: arms do not share the same initial lattice")
                end
                state, trajectory = redirect_stdout(devnull) do
                    Base.invokelatest(SR.run_simulation, params[arm], at; seed=seed)
                end
                for (mcs, lattice_state) in ((0, initial), (at, state))
                    hits = filter(s -> s.mcs == mcs, trajectory)
                    length(hits) == 1 || error("missing/duplicate snapshot at MCS $mcs")
                    snapshot = only(hits)
                    registry_r = ClusterRadial.registry_mean_r(SR, lattice_state, N)
                    species = Dict(id => cell.species for (id, cell) in lattice_state.cells)
                    metrics = parcel_contacts(lattice_state.lattice, species, SR.N_SPECIES)
                    for sd in snapshot.species_data
                        sp = sd.species
                        registry_r[sp] == sd.mean_r || error("registry/snapshot radius mismatch")
                        metric = metrics[sp]
                        length(metric.ids) == sd.n_cells && sum(metric.volumes) == sd.total_volume ||
                            error("parcel metric/snapshot registry mismatch")
                        isfinite(metric.contact) || error("undefined contact endpoint: seed $seed arm $arm species $sp")
                        arm == "no_mel" && sd.mean_melanin != 0 && error("no_mel has nonzero melanin")
                        row = Dict("seed" => string(seed), "arm" => arm, "mcs" => string(mcs),
                            "species" => string(sp), "species_name" => replace(SR.SPECIES_NAMES[sp], ',' => ' '),
                            "ncells" => string(sd.n_cells), "volume" => string(sd.total_volume),
                            "mean_r" => @sprintf("%.6f", sd.mean_r),
                            "mean_r_over_R" => @sprintf("%.6f", sd.mean_r / (N / 2)),
                            "mean_melanin" => @sprintf("%.6f", sd.mean_melanin),
                            "mean_contact_fraction" => @sprintf("%.17g", metric.contact),
                            "parcel_ids" => join(metric.ids, ';'), "parcel_volumes" => join(metric.volumes, ';'),
                            "parcel_boundary_faces" => join(metric.boundary, ';'),
                            "parcel_same_faces" => join(metric.same, ';'), "placement_sha256" => digest)
                        check_cluster_baseline(row, reference)
                        checked += 1
                        println(io, join((row[key] for key in CLUSTER_COLUMNS), ','))
                    end
                end
            end
            flush(io)
            @printf(stderr, "clustering seed %d complete, %.1f s\n", seed, time() - started)
        end
        checked == length(seeds) * 28 || error("baseline check did not cover every row")
        cluster_source_record(smoke, reference_path) == source_record || error("study source changed during the shard")
        println(io, "# baseline_checked_rows=$checked completed=true")
    end
    mv(partial, out; force=false)
    println("wrote $out; $checked rows checked against the radial reference")
end

if abspath(PROGRAM_FILE) == @__FILE__
    clustering_sweep(ARGS)
end
