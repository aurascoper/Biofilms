#!/usr/bin/env julia
# julia clustering_analyse.jl /external/clustering_*.csv
# --validate-only checks complete declared shards but emits no scientific verdict.
using Printf, Statistics, SHA
module ClusterInput
include(joinpath(@__DIR__, "clustering_sweep.jl"))
end

const CONTACT_ALPHA = 0.05
const CONTACT_MEOI = 0.05
const CONTACT_RADIUS_BOUND = 0.05
const CONTACT_Z90 = 1.6448536269514722
const CONTACT_PRIMARY_SPECIES = (1, 3)

function validate_contact_row(row)
    integer_list(field) = parse.(Int, split(row[field], ';'))
    ids = integer_list("parcel_ids")
    volumes = integer_list("parcel_volumes")
    boundary = integer_list("parcel_boundary_faces")
    same = integer_list("parcel_same_faces")
    n = parse(Int, row["ncells"])
    n > 0 && length(ids) == length(volumes) == length(boundary) == length(same) == n ||
        error("missing parcels or inconsistent parcel list lengths")
    issorted(ids) && allunique(ids) && all(>(0), ids) || error("invalid or duplicate parcel ids")
    all(>(0), volumes) && sum(volumes) == parse(Int, row["volume"]) || error("parcel volume mismatch")
    all(i -> 0 < boundary[i] <= 6 * volumes[i] && 0 <= same[i] <= boundary[i], eachindex(ids)) ||
        error("invalid contact counts")
    iseven(sum(same)) || error("same-species face contacts must count both parcels")
    n == 1 && any(>(0), same) && error("a single parcel cannot contact a distinct same-species parcel")
    contact = sum(same ./ boundary) / n
    stored = parse(Float64, row["mean_contact_fraction"])
    isfinite(stored) && stored == contact || error("contact mean disagrees with per-parcel face counts")
    for field in ("mean_r", "mean_r_over_R", "mean_melanin")
        value = parse(Float64, row[field])
        isfinite(value) && value >= 0 || error("undefined or invalid $field")
    end
    occursin(r"^[0-9a-f]{64}$", row["placement_sha256"]) || error("invalid placement digest")
    (; ids, contact, radius=parse(Float64, row["mean_r_over_R"]))
end

"Validate the full Cartesian product before any pairing or hypothesis calculation."
function read_clustering(paths; complete=true, smoke=false, reference_path=ClusterInput.CLUSTER_REFERENCE)
    isempty(paths) && error("no input shards")
    reference = ClusterInput.cluster_table(reference_path)
    rows = Dict{Tuple{Int,String,Int,Int},Dict{String,String}}()
    canonical = nothing
    allseeds = Int[]
    for path in paths
        table = ClusterInput.cluster_table(path)
        table.columns == collect(ClusterInput.CLUSTER_COLUMNS) || error("$path: unexpected schema")
        meta = table.meta
        get(meta, "completed", "false") == "true" || error("$path: incomplete shard")
        get(meta, "arms", "") == "full,no_mel" || error("$path: wrong arms")
        get(meta, "mode", "") == (smoke ? "smoke" : "production") || error("$path: incompatible run mode")
        !smoke && get(meta, "source_state", "") != "clean" && error("$path: production source was not clean")
        occursin(r"^[0-9a-f]{40}$", get(meta, "source_commit", "")) || error("$path: missing source commit")
        for label in ("julia_version", "machine")
            !isempty(get(meta, label, "")) || error("$path: missing $label")
        end
        for (label, source) in (
            "preregistration" => ClusterInput.CLUSTER_PREDICTION,
            "producer" => joinpath(@__DIR__, "clustering_sweep.jl"),
            "observable" => joinpath(@__DIR__, "clustering_observable.jl"),
            "analysis" => @__FILE__, "tests" => joinpath(@__DIR__, "test_clustering.jl"),
            "radial_sweep" => joinpath(@__DIR__, "radial_sweep.jl"),
            "model" => joinpath(ClusterInput.CLUSTER_REPO, "biofilms_potts.jl"),
            "reference" => reference_path)
            get(meta, label * "_sha256", "") == bytes2hex(sha256(read(source))) ||
                error("$path: $label hash mismatch")
        end
        config = Tuple(parse(Int, meta[key]) for key in ("N", "parcels_per_species", "read_at_mcs"))
        config == Tuple(parse(Int, reference.meta[key]) for key in ("N", "parcels_per_species", "read_at_mcs")) ||
            error("$path: baseline configuration differs")
        !smoke && config != (40, 6, 100) && error("$path: configuration differs from registration")
        seeds = parse.(Int, split(meta["declared_seeds"], ';'))
        !isempty(seeds) && allunique(seeds) || error("$path: duplicate/empty seed declaration")
        !smoke && !all(s -> 42 <= s <= 297, seeds) && error("$path: seed outside registration")
        at = config[3]
        expected = Set((seed, arm, mcs, sp) for seed in seeds for arm in ClusterInput.CLUSTER_ARMS
                       for mcs in (0, at) for sp in 1:7)
        Set(keys(table.rows)) == expected || error("$path: missing or unexpected seed/arm/time/species row")
        parse(Int, meta["baseline_checked_rows"]) == length(expected) || error("$path: incomplete baseline check")
        common = Dict(k => v for (k, v) in meta if !(k in ("declared_seeds", "baseline_checked_rows")))
        if isnothing(canonical)
            canonical = common
        else
            common == canonical || error("$path: incompatible shard metadata")
        end
        for (key, row) in table.rows
            haskey(rows, key) && error("duplicate pairing across shards: $key")
            validate_contact_row(row)
            ClusterInput.check_cluster_baseline(row, reference)
            rows[key] = row
        end
        append!(allseeds, seeds)
        for seed in seeds
            digests = Set(table.rows[(seed, arm, mcs, sp)]["placement_sha256"]
                          for arm in ClusterInput.CLUSTER_ARMS for mcs in (0, at) for sp in 1:7)
            length(digests) == 1 || error("seed $seed: arms do not share initial placement")
            for sp in 1:7
                a, b = table.rows[(seed, "full", 0, sp)], table.rows[(seed, "no_mel", 0, sp)]
                all(key -> key == "arm" || a[key] == b[key], table.columns) ||
                    error("seed $seed species $sp: unequal initial observables")
            end
            for arm in ClusterInput.CLUSTER_ARMS, mcs in (0, at)
                ids = vcat([validate_contact_row(table.rows[(seed, arm, mcs, sp)]).ids for sp in 1:7]...)
                allunique(ids) || error("seed $seed arm $arm: parcel id assigned to multiple species")
            end
        end
    end
    allunique(allseeds) || error("overlapping seed shards")
    complete && sort(allseeds) != collect(42:297) && error("primary analysis requires every seed 42:297 exactly once")
    (; rows, seeds=sort(allseeds), meta=canonical)
end

"Exact one-sided sign probability P[Binomial(n, 1/2) >= k]; ties are excluded upstream."
function contact_sign_probability(k, n)
    0 <= k <= n || error("invalid sign counts")
    sum(binomial(big(n), big(i)) for i in k:n) // big(2)^n
end

function contact_summary(contact_differences, radius_differences)
    length(contact_differences) == length(radius_differences) >= 2 || error("paired vectors differ or are too short")
    all(isfinite, contact_differences) && all(isfinite, radius_differences) || error("nonfinite paired value")
    nonzero = filter(!=(0.0), contact_differences)
    k, n = count(>(0), nonzero), length(nonzero)
    p = contact_sign_probability(k, n)
    effect = mean(contact_differences)
    radial = mean(radius_differences)
    half = CONTACT_Z90 * std(radius_differences) / sqrt(length(radius_differences))
    lo, hi = radial - half, radial + half
    size_pass = effect >= CONTACT_MEOI
    sign_pass = p < CONTACT_ALPHA
    radius_pass = lo > -CONTACT_RADIUS_BOUND && hi < CONTACT_RADIUS_BOUND
    (; effect, sd=std(contact_differences), k, n, ties=length(contact_differences)-n, p,
       radial, radial_sd=std(radius_differences), lo, hi, size_pass, sign_pass, radius_pass,
       passes=size_pass && sign_pass && radius_pass)
end

function clustering_analysis(args)
    validation = "--validate-only" in args
    paths = filter(!=("--validate-only"), args)
    any(startswith("--"), paths) && error("unknown analysis option")
    data = read_clustering(paths; complete=!validation)
    if validation
        println("VALIDATED complete shards for $(length(data.seeds)) seeds; no hypothesis verdict issued")
        return
    end
    println("Registered clustering comparison: full minus no_mel; 256 paired seeds; N40,p6,MCS100")
    println("source_commit=$(data.meta["source_commit"]) preregistration_sha256=$(data.meta["preregistration_sha256"])")
    println("Primary contact criterion: mean >= 0.05, one-sided exact sign p < 0.05; radius 90% CI strictly inside +-0.05 r/R")
    results = Bool[]
    for sp in 1:7
        contact_differences, radius_differences = Float64[], Float64[]
        for seed in data.seeds
            full = validate_contact_row(data.rows[(seed, "full", 100, sp)])
            off = validate_contact_row(data.rows[(seed, "no_mel", 100, sp)])
            push!(contact_differences, full.contact - off.contact)
            push!(radius_differences, full.radius - off.radius)
        end
        r = contact_summary(contact_differences, radius_differences)
        name = data.rows[(first(data.seeds), "full", 100, sp)]["species_name"]
        primary = sp in CONTACT_PRIMARY_SPECIES
        println(primary ? "PRIMARY: $name" : "DESCRIPTIVE: $name")
        @printf("  contact mean %+.8f, paired SD %.8f; %d/%d positive, %d ties, one-sided p %.8g\n",
                r.effect, r.sd, r.k, r.n, r.ties, Float64(r.p))
        @printf("  radius mean %+.8f, paired SD %.8f, 90%% CI [%+.8f, %+.8f]\n", r.radial, r.radial_sd, r.lo, r.hi)
        if primary
            println("  effect_size=$(r.size_pass) sign_test=$(r.sign_pass) radius_equivalence=$(r.radius_pass)")
            push!(results, r.passes)
        end
    end
    length(results) == 2 || error("both primary species were not evaluated")
    println(all(results) ? "COMPOUND PREDICTION: supported under the registered criteria" :
                          "COMPOUND PREDICTION: did not meet the registered criteria")
    println("The radius values reproduce previously inspected data; this is not an independent radius confirmation.")
    println("The contact endpoint does not measure temporal retention, an adhesion coefficient, or organisms.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    clustering_analysis(ARGS)
end
