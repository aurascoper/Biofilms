#!/usr/bin/env julia
# Data-free checks on the configuration reader, the destination guards, the shared step
# count, the receipt hashes and the control verdicts: everything setup.jl and controls.jl
# decide before or after the integrator runs.
#
#   julia --project=diagnostics/binding_benchmark diagnostics/binding_benchmark/test_setup.jl
#
# Nothing here reads the frozen geometry. Each case writes its own configuration into a
# scratch directory from the shipped benchmark.toml, so the cases stay in step with the
# schema the reader enforces.
using Test, TOML
include(joinpath(@__DIR__, "setup.jl"))

# controls.jl includes setup.jl itself, so it is loaded in its own module rather than a
# second time into Main, where two copies of BindingBenchmark would make `Params` ambiguous.
# Its main() is guarded by PROGRAM_FILE and does not run here.
module CtlTest
include(joinpath(@__DIR__, "controls.jl"))
end

const SHIPPED = TOML.parsefile(joinpath(@__DIR__, "benchmark.toml"))

"Write the shipped configuration with `edit!` applied, and return its path."
function config_with(edit!)
    cfg = deepcopy(SHIPPED)
    edit!(cfg)
    path = joinpath(mktempdir(), "benchmark.toml")
    open(path, "w") do io; TOML.print(io, cfg); end
    path
end

@testset "binding benchmark setup" begin

@testset "read_config refuses what --config could otherwise smuggle in" begin
    cfg, p = read_config(joinpath(@__DIR__, "benchmark.toml"))
    @test p isa Params && p.lambda == 0.01
    # A TOML integer is a configuration, not a MethodError.
    _, p0 = read_config(config_with(c -> (c["params"]["D_c"] = 0; c["params"]["B0"] = 1)))
    @test p0.D_c === 0.0 && p0.B0 === 1.0
    # Signs: a negative loss rate is growth; the closed-decay claim does not survive it.
    @test_throws ArgumentError read_config(config_with(c -> c["params"]["lambda"] = -0.01))
    @test_throws ArgumentError read_config(config_with(c -> c["params"]["k_on"] = -1.0))
    @test_throws ArgumentError read_config(config_with(c -> c["params"]["k_on"] = 0.0))   # K3 would be order(0, 0)
    @test_throws ArgumentError read_config(config_with(c -> c["params"]["D_c"] = -0.1))
    @test_throws ArgumentError read_config(config_with(c -> c["params"]["B0"] = -1.0))
    # Hill parameters: K = 0 divides by zero at A = 0, n = 0 is not a response.
    @test_throws ArgumentError read_config(config_with(c -> c["params"]["K"] = 0.0))
    @test_throws ArgumentError read_config(config_with(c -> c["params"]["n"] = 0.0))
    # No capacity anywhere is a benchmark of nothing.
    @test_throws ArgumentError read_config(config_with(c -> (c["params"]["B0"] = 0.0; c["params"]["dB"] = 0.0)))
    # NaN passes every comparison written as `x > 0` and fails every one written as this.
    @test_throws ArgumentError read_config(config_with(c -> c["params"]["lambda"] = NaN))
    @test_throws ArgumentError read_config(config_with(c -> c["dt"] = NaN))
    # Initial pools: negative is refused, and an empty inventory makes 0/0 residuals.
    @test_throws ArgumentError read_config(config_with(c -> c["c0"] = -1.0))
    @test_throws ArgumentError read_config(config_with(c -> (c["c0"] = 0.0; c["b0"] = 0.0)))
    _, pb = read_config(config_with(c -> (c["c0"] = 0.0; c["b0"] = 0.5)))
    @test pb isa Params
end

@testset "the destination is created exclusively, and never inside the parent" begin
    root = mktempdir()
    d = joinpath(root, "out")
    @test fresh_destination(d) == d && isdir(d)
    @test_throws ArgumentError fresh_destination(d)            # the first check
    # The gap between the check and the creation: a directory that appears inside it
    # must be refused by the creation, not accepted by it. With `mkpath` this test
    # returned normally and the run would have written into the intruder's directory.
    d2 = joinpath(root, "out2")
    @test_throws Base.IOError fresh_destination(d2; race = () -> mkdir(d2))
    # A failed run leaves nothing behind: the parent is made, the leaf is not.
    d3 = joinpath(root, "deep", "out3")
    @test_throws Base.IOError fresh_destination(d3; race = () -> mkdir(d3))
    @test isdir(joinpath(root, "deep"))

    parent = mkdir(joinpath(root, "parent"))
    @test_throws ArgumentError refuse_inside(parent, joinpath(parent, "bench"))
    @test_throws ArgumentError refuse_inside(parent, parent)
    @test_throws ArgumentError refuse_inside(parent, joinpath(parent, "a", "..", "b"))
    @test isnothing(refuse_inside(parent, joinpath(root, "parent2")))  # a sibling, not a child
    @test isnothing(refuse_inside(parent, joinpath(root, "elsewhere", "bench")))
    # Through a symlink from outside: the lexical prefix differs, the resolved one does not.
    link = joinpath(root, "link"); symlink(parent, link)
    @test_throws ArgumentError refuse_inside(parent, joinpath(link, "bench"))
    # A root parent: the old prefix test built "//" and accepted everything.
    @test_throws ArgumentError refuse_inside("/", "/tmp/out")
    # A child whose name begins with ".." is a child, not an ancestor.
    @test_throws ArgumentError refuse_inside(parent, joinpath(parent, "..bench"))

    # The scope check is repeated on the created leaf: an intermediate component swapped
    # for a symlink into the parent after the first check is caught, and nothing is left
    # inside the parent.
    mid = mkdir(joinpath(root, "mid"))
    target = joinpath(mid, "out4")
    @test isnothing(refuse_inside(parent, target))
    swap = () -> (rm(mid); symlink(parent, mid))
    @test_throws ArgumentError fresh_destination(target; parent = parent, race = swap)
    @test isempty(readdir(parent))
end

@testset "the step count is shared, and a non-integral horizon is refused" begin
    @test nsteps_for(50.0, 0.5) == 100
    @test nsteps_for(50.0, 0.1) == 500                         # 0.1 is inexact; the tolerance absorbs it
    @test_throws ArgumentError nsteps_for(50.1, 0.5)          # what controls.jl used to round to 100
    @test_throws ArgumentError nsteps_for(1.0, 0.3)
    @test_throws ArgumentError nsteps_for(1e-10, 1.0)         # rounded to zero steps, within tolerance
    # Both entry points route through it: neither may round on its own.
    for f in ("run.jl", "controls.jl")
        src = read(joinpath(@__DIR__, f), String)
        @test occursin("nsteps_for(total, dt)", src)
        @test !occursin("round(Int, total / dt)", src)
    end
end

@testset "the receipt hashes bind the configuration that was read, and the environment" begin
    custom = config_with(c -> c["params"]["lambda"] = 0.02)
    h = code_hashes(custom)
    @test h["config"] == sha256_file(custom)
    @test h["config"] != sha256_file(joinpath(@__DIR__, "benchmark.toml"))
    @test h["config_file"] == "benchmark.toml"      # the name alone would not have told them apart
    for f in ("Project.toml", "Manifest.toml", "test_setup.jl", "setup.jl")
        @test h[f] == sha256_file(joinpath(@__DIR__, f))
    end
end

@testset "control verdicts and the expected-red list" begin
    k6_verdict = CtlTest.k6_verdict; Ledger = CtlTest.Ledger; EXPECTED_RED = CtlTest.EXPECTED_RED
    # A ledger whose pre-release record shows a transient overshoot, released before T:
    # the final state is clean and the old predicate said PASS.
    clean = Ledger(0.0, 0, 1.0, 0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1, 0.0, 0.0)
    @test k6_verdict(clean, 0.0) == "PASS"
    transient = Ledger(0.0, 0, 1.0, 0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1, 0.0, 1e-9)
    @test k6_verdict(transient, 0.0) == "FAIL"
    negative = Ledger(0.0, 0, 1.0, 0.0, 0, 0, 0, 0, 0, 0, 0, 0, -1e-12, 0.0, 0.0)
    @test k6_verdict(negative, 0.0) == "FAIL"
    @test k6_verdict(clean, 1e-9) == "FAIL"

    # Both red-by-design controls are declared, by the names controls.jl gives them.
    src = read(joinpath(@__DIR__, "controls.jl"), String)
    @test length(EXPECTED_RED) == 2
    @test any(startswith("K8 "), EXPECTED_RED) && any(startswith("K12 "), EXPECTED_RED)
    for n in EXPECTED_RED
        @test occursin("add!(\"$n\",", src)
    end
    # Every control has one expected verdict, and a red control that passes is unexpected.
    ev = CtlTest.expected_verdict
    @test ev(EXPECTED_RED[1]) == "FIRES" && ev(EXPECTED_RED[2]) == "FIRES"
    @test ev("K1 matched total capacity") == "MEASURED"
    @test ev("K10 both bound-fraction denominators") == "MEASURED"
    @test ev("K11 ledger closure from independent accumulators") == "PASS"
    @test ev("K7 conservative release under a forced capacity decrease") == "FIRES"   # green is FIRES here
    # Every add! call in controls.jl names a control the table classifies without falling
    # through to PASS by accident: the FIRES-vocabulary controls are exactly K7, K8, K12.
    names = [m.captures[1] for m in eachmatch(r"add!\(\"(K\d+ [^\"]+)\"", src)]
    @test length(names) == 12
    @test [n[1:findfirst(' ', n)-1] for n in names if ev(n) == "FIRES"] == ["K7", "K8", "K12"]
    @test "PASS" != ev(EXPECTED_RED[1])
end

@testset "the committed receipts were produced by this code, in this environment" begin
    # A receipt binds a run to the code and environment that made it. Both committed
    # receipts drifted while every function they call kept passing: they held seven hash
    # keys after code_hashes grew to eleven, and one expected-red control after controls.jl
    # declared two. Nothing read the committed bytes; this does, so a producer edit without
    # a regeneration against the parent bundle is red here, not discovered in review.
    want = code_hashes(joinpath(@__DIR__, "benchmark.toml"))
    for f in ("benchmark_receipt.json", "control_verification.json")
        r = JSON3.read(read(joinpath(@__DIR__, f), String))
        got = Dict(String(k) => String(v) for (k, v) in pairs(r[:code_sha256]))
        @test sort(collect(keys(got))) == sort(collect(keys(want)))
        for (k, v) in want
            @test get(got, k, "<absent from $f>") == v
        end
        @test String(r[:parent][:manifest_sha256]) == SHIPPED["parent_manifest_sha256"]
        @test String(r[:parent][:snapshot]) == SHIPPED["parent_snapshot"]
    end
    ctl = JSON3.read(read(joinpath(@__DIR__, "control_verification.json"), String))
    @test collect(String.(ctl[:expected_red])) == CtlTest.EXPECTED_RED
    # The mutation receipt is data-free and binds only the module and the suite it mutated.
    mut = JSON3.read(read(joinpath(@__DIR__, "mutation_verification.json"), String))
    @test String(mut[:module_sha256]) == sha256_file(joinpath(@__DIR__, "BindingBenchmark.jl"))
    @test String(mut[:suite_sha256]) == sha256_file(joinpath(@__DIR__, "test_numerics.jl"))
end

end # setup
