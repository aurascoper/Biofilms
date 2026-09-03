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

include("radiodialysis_basis_gate.jl")

@testset "Lifecycle, dose contract, windowed API" begin
    include("genealogy_tests.jl")
end

@testset "Exchange schemas (snapshot + restart)" begin
    include("checkpoint_io_tests.jl")
end

# The JACC port had no automated execution at all until this tier. It runs on
# whichever backend JACC selects — threads where there is no GPU, which is what
# a CI runner sees — so the portability layer's portability is itself tested.
@testset "JACC port kernels versus the serial reference" begin
    include("jacc_port_tests.jl")
end

@testset "Console report honesty" begin
    include("console_report_tests.jl")
end

@testset "Manuscript claims" begin
    include("manuscript_claims_tests.jl")
end
