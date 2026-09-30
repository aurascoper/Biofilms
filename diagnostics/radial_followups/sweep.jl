#!/usr/bin/env julia
# One seed/configuration per atomic external shard. See PREDICTION.md.
using Printf, SHA
module LegacyRadial
include(joinpath(@__DIR__, "..", "melanin_ensemble", "radial_sweep.jl"))
end
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const BASELINE = joinpath(ROOT, "diagnostics", "melanin_ensemble", "radial_n40_p6_seeds42-297.csv")
const FACTORS = (1, 10, 100, 1000, 10000, 100000, 1000000)
const COLUMNS = "seed,arm,mcs,species,species_name,ncells,volume,mean_r,mean_r_over_R,mean_melanin"
digest(path) = bytes2hex(sha256(read(path)))

function options(args)
    isempty(args) && error("usage: sweep.jl OUT.csv --axis beta|alpha|horizon --factor F --arm full|null --seed S --expected-revision SHA [--test --n N --parcels P --at MCS]")
    out, opts = abspath(args[1]), Dict{String,String}()
    allowed = ("--axis", "--factor", "--arm", "--seed", "--expected-revision", "--n", "--parcels", "--at", "--baseline-root")
    i = 2
    while i <= length(args)
        key = args[i]
        haskey(opts, key) && error("duplicate option $key")
        if key == "--test"
            opts[key] = "true"; i += 1; continue
        end
        key in allowed || error("unknown option $key")
        i < length(args) || error("missing value for $key")
        opts[key] = args[i + 1]; i += 2
    end
    axis = get(opts, "--axis", ""); axis in ("beta", "alpha", "horizon", "baseline") || error("unknown axis $axis")
    factor = parse(Int, get(opts, "--factor", "1"))
    arm = get(opts, "--arm", "full"); arm in ("full", "null") || error("unknown arm $arm")
    axis in ("horizon", "baseline") || arm == "full" || error("coefficient sweep must use full arm")
    !(axis in ("horizon", "baseline")) || factor == 1 || error("horizon/baseline uses factor 1")
    axis != "baseline" || arm == "null" || error("baseline axis is the null reproduction")
    seed = parse(Int, get(opts, "--seed", "0"))
    test = haskey(opts, "--test")
    n = parse(Int, get(opts, "--n", "40")); parcels = parse(Int, get(opts, "--parcels", "6"))
    at = parse(Int, get(opts, "--at", axis == "horizon" ? "400" : "100"))
    factor > 0 && n >= 12 && parcels > 0 && at > 0 || error("invalid numerical option")
    if test
        n <= 16 && parcels == 1 && at <= 4 && seed >= 9000 || error("--test is restricted to tiny, nonregistered configurations")
    else
        (n, parcels, at) == (40, 6, axis == "horizon" ? 400 : 100) || error("production configuration is fixed")
        seed in 42:297 && factor in FACTORS || error("production seed or factor outside registered grid")
        revision = get(opts, "--expected-revision", "")
        occursin(r"^[0-9a-f]{40}$", revision) || error("production requires full --expected-revision")
        strip(read(`git -C $ROOT rev-parse HEAD`, String)) == revision || error("wrong production revision")
        isempty(strip(read(`git -C $ROOT status --porcelain`, String))) || error("production requires a clean worktree")
        startswith(out, ROOT * "/") && error("production output must be outside the repository")
    end
    ispath(out) && error("destination exists: $out; validate before resuming")
    (; out, axis, factor, arm, seed, test, n, parcels, at,
       revision = get(opts, "--expected-revision", "test"), baseline_root = get(opts, "--baseline-root", ""))
end

function require_baseline_pair(cfg)
    cfg.test && return
    (cfg.axis == "alpha" || (cfg.axis == "beta" && cfg.factor > 1)) || return
    isempty(cfg.baseline_root) && error("scaled runs require --baseline-root with full and null reproductions")
    for (arm, rel) in (("full", joinpath("beta", "factor_1")), ("null", joinpath("baseline", "null")))
        path = joinpath(cfg.baseline_root, rel, "seed_$(cfg.seed).csv")
        isfile(path) || error("missing baseline reproduction $path")
        lines = readlines(path)
        metadata = Dict{String,String}()
        for line in lines
            startswith(line, "# ") || continue
            k, v = split(line[3:end], '='; limit = 2); metadata[k] = v
        end
        for (k, v) in (("mode", "production"), ("revision", cfg.revision),
                       ("baseline_sha256", digest(BASELINE)), ("baseline_checked", "14"),
                       ("model_sha256", digest(joinpath(ROOT, "biofilms_potts.jl"))),
                       ("producer_sha256", digest(@__FILE__)),
                       ("prediction_sha256", digest(joinpath(@__DIR__, "PREDICTION.md"))))
            get(metadata, k, "") == v || error("invalid baseline $k in $path")
        end
        rows = filter(l -> startswith(l, "$(cfg.seed),$arm,"), lines)
        length(rows) == 14 || error("baseline reproduction must contain exactly 14 rows")
        baseline_check(rows, merge(cfg, (; arm, axis = "baseline", factor = 1)))
    end
end

function parameters(SR, cfg)
    base = LegacyRadial.arm_params(SR, cfg.arm, cfg.n, cfg.parcels, cfg.test ? cfg.at : 100)
    beta, alpha = copy(base.β_ion), copy(base.α_M_species)
    cfg.axis == "beta" && (beta .*= cfg.factor)
    cfg.axis == "alpha" && (alpha .*= cfg.factor)
    Base.invokelatest(SR.CPMParams; N = cfg.n, n_cells_per_species = cfg.parcels,
        snapshot_interval = cfg.test ? cfg.at : 100, β_ion = beta, α_M_species = alpha)
end

function snapshot_lines(SR, cfg, snapshots)
    lines = String[]
    for snap in snapshots, sd in snap.species_data
        vals = (sd.mean_r, sd.mean_r / (cfg.n / 2), sd.mean_melanin)
        all(isfinite, vals) || error("nonfinite snapshot at seed $(cfg.seed), MCS $(snap.mcs)")
        push!(lines, @sprintf("%d,%s,%d,%d,%s,%d,%d,%.6f,%.6f,%.6f", cfg.seed,
            cfg.arm, snap.mcs, sd.species, replace(SR.SPECIES_NAMES[sd.species], "," => ""),
            sd.n_cells, sd.total_volume, vals...))
    end
    lines
end

function baseline_check(lines, cfg, baseline = BASELINE)
    reference = Dict{Tuple{Int,Int},String}()
    for line in eachline(baseline)
        startswith(line, "$(cfg.seed),$(cfg.arm),") || continue
        f = split(line, ','); reference[(parse(Int, f[3]), parse(Int, f[4]))] = line
    end
    checked = 0; seen = Set{Tuple{Int,Int}}()
    for line in lines
        f = split(line, ','); mcs, species = parse(Int, f[3]), parse(Int, f[4])
        must_match = mcs == 0 || ((cfg.axis == "horizon" || cfg.factor == 1) && mcs == 100)
        must_match || continue
        (mcs, species) in seen && error("duplicate baseline row $mcs, $species")
        push!(seen, (mcs, species))
        get(reference, (mcs, species), "") == line || error("baseline mismatch: seed $(cfg.seed), arm $(cfg.arm), MCS $mcs, species $species")
        checked += 1
    end
    expected = (cfg.axis == "horizon" || cfg.factor == 1) ? 14 : 7
    checked == expected || error("baseline check covered $checked rows, expected $expected")
    checked
end

function run_one(cfg)
    require_baseline_pair(cfg)
    SR = LegacyRadial.load_serial()
    # The dynamically loaded module's bindings must be read in their definition world.
    Base.invokelatest(run_loaded, SR, cfg)
end

function run_loaded(SR, cfg)
    p = parameters(SR, cfg); t0 = time()
    state, trajectory = redirect_stdout(devnull) do
        Base.invokelatest(SR.run_simulation, p, cfg.at; seed = cfg.seed)
    end
    wanted = cfg.test ? (0, cfg.at) : cfg.axis == "horizon" ? (0, 100, 400) : (0, 100)
    snapshots = [only(filter(s -> s.mcs == at, trajectory)) for at in wanted]
    recomputed = LegacyRadial.registry_mean_r(SR, state, cfg.n)
    all(sd.mean_r == recomputed[sd.species] for sd in snapshots[end].species_data) || error("registry/snapshot mismatch")
    cfg.arm != "null" || all(sd.mean_melanin == 0 for s in snapshots for sd in s.species_data) || error("null produces melanin")
    lines = snapshot_lines(SR, cfg, snapshots)
    checked = cfg.test ? 0 : baseline_check(lines, cfg)
    elapsed = time() - t0
    meta = ["format=radial_followups_v1", "mode=$(cfg.test ? "test" : "production")",
        "revision=$(cfg.revision)", "axis=$(cfg.axis)", "factor=$(cfg.factor)", "arm=$(cfg.arm)",
        "seed=$(cfg.seed)", "N=$(cfg.n)", "parcels=$(cfg.parcels)", "at=$(cfg.at)",
        "julia_version=$(VERSION)", "machine=$(Sys.MACHINE)", "baseline_checked=$checked",
        "baseline_sha256=$(digest(BASELINE))", "model_sha256=$(digest(joinpath(ROOT, "biofilms_potts.jl")))",
        "producer_sha256=$(digest(@__FILE__))", "prediction_sha256=$(digest(joinpath(@__DIR__, "PREDICTION.md")))",
        "beta=" * join(p.β_ion, ';'), "alpha=" * join(p.α_M_species, ';'), "elapsed_seconds=$elapsed"]
    mkpath(dirname(cfg.out)); temporary = cfg.out * ".partial.$(getpid())"
    try
        open(temporary, "w") do io
            foreach(m -> println(io, "# ", m), meta)
            println(io, COLUMNS); foreach(line -> println(io, line), lines)
        end
        mv(temporary, cfg.out; force = false)
    finally
        isfile(temporary) && rm(temporary)
    end
    println("completed $(cfg.axis) factor=$(cfg.factor) arm=$(cfg.arm) seed=$(cfg.seed) seconds=$elapsed")
    cfg.out
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_one(options(ARGS))
end
