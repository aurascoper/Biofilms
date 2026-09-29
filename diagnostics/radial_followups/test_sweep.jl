#!/usr/bin/env julia
using Test
include(joinpath(@__DIR__, "sweep.jl"))

@testset "coefficient/horizon producer contracts; tiny tests only" begin
    directory = mktempdir()
    cli = [joinpath(directory, "tiny.csv"), "--axis", "beta", "--factor", "3",
           "--seed", "9001", "--test", "--n", "12", "--parcels", "1", "--at", "2"]
    cfg = options(cli)
    SR = LegacyRadial.load_serial()
    base = Base.invokelatest(SR.CPMParams)
    beta = parameters(SR, cfg)
    @test beta.β_ion == 3 .* base.β_ion
    @test beta.α_M_species == base.α_M_species
    alpha = parameters(SR, merge(cfg, (; axis = "alpha")))
    @test alpha.α_M_species == 3 .* base.α_M_species
    @test alpha.β_ion == base.β_ion
    @test base.β_ion == Base.invokelatest(SR.CPMParams).β_ion
    nul = parameters(SR, merge(cfg, (; axis = "horizon", factor = 1, arm = "null")))
    @test all(iszero, nul.β_ion) && all(iszero, nul.α_M_species)
    @test_throws ErrorException options([cli; "--axis"; "alpha"])
    @test_throws ErrorException options(replace(cli, "beta" => "unknown"))
    @test_throws ErrorException options(replace(cli, "9001" => "42"))
    @test_throws ErrorException options(replace(cli, "12" => "40"))
    @test_throws ErrorException options([joinpath(directory, "bad.csv"), "--axis", "beta", "--seed", "42"])
    @test !isfile(joinpath(directory, "bad.csv"))

    # Real producer versus the original producer at a nonregistered test configuration.
    one = merge(cfg, (; factor = 1))
    run_one(one)
    original = joinpath(directory, "original.csv")
    redirect_stdout(devnull) do
        LegacyRadial.sweep_main([original, "--n", "12", "--parcels", "1", "--at", "2",
                                "--arms", "full", "--seeds", "9001"])
    end
    data(path) = filter(l -> startswith(l, "9001,"), readlines(path))
    @test length(data(one.out)) == 14
    @test data(one.out) == data(original)
    @test occursin("# mode=test", read(one.out, String))
    @test isempty(filter(p -> occursin(".partial.", p), readdir(directory)))
    @test_throws ErrorException options(cli)  # refuses to overwrite a complete shard

    # These are existing archive rows, not a new run at any production coefficient.
    historical = filter(l -> startswith(l, "42,full,"), readlines(BASELINE))
    historical_cfg = merge(cfg, (; seed = 42, factor = 1, test = false))
    @test baseline_check(historical, historical_cfg) == 14
    @test_throws ErrorException baseline_check(historical[1:end-1], historical_cfg)
    @test_throws ErrorException baseline_check([historical[1:end-1]; historical[1]], historical_cfg)
    damaged = copy(historical)
    cells = split(damaged[end], ','); cells[8] = "99.000000"; damaged[end] = join(cells, ',')
    @test_throws ErrorException baseline_check(damaged, historical_cfg)
    @test_throws ErrorException require_baseline_pair(merge(historical_cfg, (; factor = 10)))

    # A complete historical pair permits execution; missing/stale companions refuse it.
    pair_cfg = merge(historical_cfg, (; factor = 10, baseline_root = directory))
    for (arm, rel) in (("full", joinpath("beta", "factor_1")), ("null", joinpath("baseline", "null")))
        path = joinpath(directory, rel, "seed_42.csv"); mkpath(dirname(path))
        rows = filter(l -> startswith(l, "42,$arm,"), readlines(BASELINE))
        metadata = ["mode=production", "revision=$(pair_cfg.revision)", "baseline_checked=14",
                    "baseline_sha256=$(digest(BASELINE))", "producer_sha256=$(digest(joinpath(@__DIR__, "sweep.jl")))",
                    "model_sha256=$(digest(joinpath(ROOT, "biofilms_potts.jl")))",
                    "prediction_sha256=$(digest(joinpath(@__DIR__, "PREDICTION.md")))"]
        write(path, join(["# " .* metadata; COLUMNS; rows], "\n") * "\n")
    end
    @test require_baseline_pair(pair_cfg) === nothing
    null_path = joinpath(directory, "baseline", "null", "seed_42.csv")
    original_null = read(null_path, String)
    write(null_path, replace(original_null, "baseline_checked=14" => "baseline_checked=0"))
    @test_throws ErrorException require_baseline_pair(pair_cfg)
    rm(null_path)
    @test_throws ErrorException require_baseline_pair(pair_cfg)
end
