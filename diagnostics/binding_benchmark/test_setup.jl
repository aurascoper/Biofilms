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
end

@testset "the step count is shared, and a non-integral horizon is refused" begin
    @test nsteps_for(50.0, 0.5) == 100
    @test nsteps_for(50.0, 0.1) == 500                         # 0.1 is inexact; the tolerance absorbs it
    @test_throws ArgumentError nsteps_for(50.1, 0.5)          # what controls.jl used to round to 100
    @test_throws ArgumentError nsteps_for(1.0, 0.3)
    # Both entry points route through it: neither may round on its own.
    for f in ("run.jl", "controls.jl")
        src = read(joinpath(@__DIR__, f), String)
        @test occursin("nsteps_for(total, dt)", src)
        @test !occursin("round(Int, total / dt)", src)
    end
end

end # setup
