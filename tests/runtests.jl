#!/usr/bin/env julia
# Test entry point:  julia --project=. tests/runtests.jl
# (The repo is not a Julia package, so there is no Pkg.test path.)

using Test

const REPO = dirname(@__DIR__)

# Same loading trick as validate_serial.jl: the serial monolith minus its
# CairoMakie figure section, evaluated in a sandbox module. The split marker
# "#  13. Figure export" is load-bearing — see docs/audit_biofilms_potts.md §9.
function load_serial()
    src = read(joinpath(REPO, "biofilms_potts.jl"), String)
    src = split(src, "#  13. Figure export")[1]
    M = Module(:SerialRef)
    Base.eval(M, :(using LinearAlgebra, Statistics, Random, Printf))
    Base.include_string(M, src, "biofilms_potts.jl")
    return M
end

const SR = load_serial()

@testset "Biofilms serial contract" begin
    include("contract_csv.jl")
end

@testset "Deterministic radiation response" begin
    include("deterministic_radiation.jl")
end

@testset "Lifecycle, dose contract, windowed API" begin
    include("genealogy_tests.jl")
end

@testset "Exchange schemas (snapshot + restart)" begin
    include("checkpoint_io_tests.jl")
end
