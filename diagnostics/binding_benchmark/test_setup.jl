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

end # setup
