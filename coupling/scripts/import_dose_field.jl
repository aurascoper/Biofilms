#!/usr/bin/env julia
# The return leg: a PHYSICAL dose field, produced by OpenMC, installed into CPM
# state. This is the last step of the synthetic one-way chain, and it is where
# the chain STOPS.
#
#   julia --project=. coupling/scripts/import_dose_field.jl <dose_field.h5>
#
# What it deliberately does NOT do: pass `hamiltonian_transform` or
# `melanin_transform`, so the Metropolis radiation term and the melanin drive
# keep their legacy dimensionless values and the dynamics are untouched. And it
# checks that `accrue_dose!` still REFUSES, because `seconds_per_mcs` is NaN
# until a temporal calibration exists. Per-parcel dose in Gy needs a clock that
# the spatial gate has not produced; the mass-weighted attribution on the Python
# side needs none, which is why that is where attribution happens.
#
# So the strongest statement this run supports is: a physical dose-rate field
# crossed the boundary in the correct orientation. Not that anything responded
# to it.

using HDF5
using Test

# The same loading trick as tests/runtests.jl and validate_serial.jl: the serial
# monolith MINUS its CairoMakie figure section, in a sandbox module. The split
# marker "#  13. Figure export" is load-bearing (docs/audit_biofilms_potts.md
# §9) and nothing above it may use a non-stdlib package — a plain `include`
# pulls in CairoMakie and fails.
const REPO = dirname(dirname(@__DIR__))

function load_serial()
    src = read(joinpath(REPO, "biofilms_potts.jl"), String)
    src = split(src, "#  13. Figure export")[1]
    M = Module(:SerialRef)
    Base.eval(M, :(using LinearAlgebra, Statistics, Random, Printf))
    Base.include_string(M, src, "biofilms_potts.jl")
    return M
end

const SR = load_serial()

function main(path::AbstractString)
    mean_Gy_s, sd_Gy_s, attrs = h5open(path, "r") do f
        # Written logical (x,y,z); h5py stores dims reversed, so the Julia
        # reader gets (x,y,z) back directly. The orientation contract is in
        # docs/exchange_schema.md.
        read(f["mesh/dose_rate_mean_Gy_s"]),
        read(f["mesh/dose_rate_sd_Gy_s"]),
        Dict(k => read_attribute(f, k) for k in keys(attributes(f)))
    end

    println("dose field      : ", size(mean_Gy_s), " ", attrs["logical_axis_order"])
    println("reference system: ", attrs["reference_system_id"])
    println("execution class : ", attrs["execution_class"],
            "  target_calibration=", attrs["target_calibration"])
    println("range           : ", minimum(mean_Gy_s), " .. ", maximum(mean_Gy_s), " Gy/s")

    N = size(mean_Gy_s, 1)
    @testset "one-way dose import (synthetic)" begin
        @test size(mean_Gy_s) == (N, N, N)
        @test attrs["logical_axis_order"] == "xyz"
        # The label has to survive the field being read somewhere else.
        @test attrs["target_calibration"] == false

        st = SR.init_state(SR.CPMParams(N = N, n_cells_per_species = 2); seed = 42)
        rad0 = copy(st.radiation)
        mel0 = copy(st.melanin_drive)

        SR.import_dose_field!(st, mean_Gy_s, sd_Gy_s)

        @test st.dose_active
        @test st.dose_rate_mean_Gy_s == mean_Gy_s
        @test st.dose_rate_sd_Gy_s == sd_Gy_s
        # No transform supplied => the biological signals are untouched. This is
        # the line between "a dose field is present" and "the model responded".
        @test st.radiation == rad0
        @test st.melanin_drive == mel0
        # And with no clock, no dose may be accrued to any parcel.
        @test isnan(st.params.seconds_per_mcs)
        @test_throws ErrorException SR.accrue_dose!(st)
    end
    println("\nthe field is installed; the dynamics were not advanced and no ",
            "biological\nresponse was applied. THE CHAIN STOPS HERE.")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("usage: import_dose_field.jl <dose_field.h5>")
    exit(main(ARGS[1]))
end
