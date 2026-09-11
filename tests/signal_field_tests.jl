# Independent diagnostic: no serial CPM edits or basis-gate acknowledgement.
#
# Four test files sit in diagnostics/inert_signal/; until this patch only
# test_numerics.jl was reachable from here, so three files that exist and are
# maintained had never run in a suite. They were not wired in by oversight:
# test_viewer.jl and test_native_layout.jl take a parent run and its derived
# companion through ARGS and walk all 101 snapshots, and no CI runner has a
# trajectory. Including them unconditionally turns a green suite red everywhere.
#
# So the data-dependent pair is recorded as UNCOVERED when the data is absent,
# never as a skip. A skipped test reads as a passing suite; this repository's
# rule (docs/visualization/lattice_viewer.md:213-218) is that unrun surface is
# stated uncovered surface. Set BIOFILMS_PARENT_RUN and BIOFILMS_DERIVED_RUN to
# a verified bundle to cover them.
module SignalDiagnosticTests
using Test

const D=normpath(joinpath(@__DIR__,"..","diagnostics","inert_signal"))
const PROJ=normpath(joinpath(@__DIR__,".."))
const PARENT=get(ENV,"BIOFILMS_PARENT_RUN","")
const DERIVED=get(ENV,"BIOFILMS_DERIVED_RUN","")
const UNCOVERED=String[]

const have_parent  = !isempty(PARENT)  && isdir(joinpath(PARENT,"snapshots"))
const have_derived = !isempty(DERIVED) && isdir(joinpath(DERIVED,"fields"))

# These three take their inputs through ARGS and install module-level names that
# collide across files -- test_viewer.jl ships a PlotNamespaceControl shim precisely
# because Makie exports `attributes` too, and test_native_layout.jl cd()s into the
# viewer directory. A separate process per file is the isolation they already assume,
# and matches the subprocess idiom test_pipeline.jl uses for the guard CLI.
# `$args`, not `$(args...)`: inside a command literal an interpolated splat is a
# cartesian product, not a shell splat -- `$(["A","B"]...)` becomes the single
# argument "AB". An interpolated Vector{String} is what expands to separate
# arguments. The concatenated form handed test_viewer.jl one fused path and it
# failed with a MethodError that looked like a defect in the test.
run_file(file,args)=success(pipeline(ignorestatus(
    `$(Base.julia_cmd()) --project=$PROJ $(joinpath(D,file)) $args`),
    stdout=stdout,stderr=stderr))

include(joinpath(D,"test_numerics.jl"))
include(joinpath(D,"test_guard_identifiers.jl"))

@testset "Inert-signal pipeline boundaries and known-bad configs" begin
    # Runs standalone; the parent-manifest defect controls need a verified run.
    @test run_file("test_pipeline.jl", have_parent ? [PARENT] : String[])
    have_parent || push!(UNCOVERED,
        "test_pipeline.jl parent-manifest defect controls (missing/duplicate/parcels/hash)")
end

@testset "Native viewer reader, CLI, and Makie layout" begin
    if have_parent && have_derived
        @test run_file("test_viewer.jl",[PARENT,DERIVED])
        @test run_file("test_native_layout.jl",[PARENT,DERIVED])
    else
        push!(UNCOVERED,"test_viewer.jl (101-snapshot reader contract)")
        push!(UNCOVERED,"test_native_layout.jl (Makie objects and callbacks)")
    end
end

if !isempty(UNCOVERED)
    println(stderr,"\nUNCOVERED -- these files exist, are maintained, and did not run here:")
    for u in UNCOVERED; println(stderr,"  - ",u); end
    println(stderr,"Set BIOFILMS_PARENT_RUN (a run with snapshots/ and run_manifest.json)")
    println(stderr,"and BIOFILMS_DERIVED_RUN (a run with fields/) to cover them.\n")
end
end
