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
# The viewer's CLI refusals need no trajectory. They used to sit inside test_viewer.jl,
# which destructures `parent,derived=ARGS` on its third line, so the whole file was gated
# behind the data and these never ran by default -- 2 of 106 assertions, both negative
# controls, reported as UNCOVERED when nothing prevented covering them. They cannot be
# reached by invoking test_viewer.jl argument-less (line 13 would raise BoundsError), so
# they are extracted rather than called differently.
include(joinpath(D,"test_viewer_cli.jl"))

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
        push!(UNCOVERED,"test_viewer.jl 101-snapshot reader contract (its CLI refusals "*
                        "are covered above by test_viewer_cli.jl)")
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

# The inertness census behind section 6.7, which ran once by hand and never again.
module SignalCensusTests
using Test
const ROOT = normpath(joinpath(@__DIR__, ".."))
const DIAG = joinpath(ROOT, "diagnostics", "inert_signal")
include(joinpath(DIAG, "guard.jl"))

@testset "Inert signal census over the CPM sources" begin
    sources = census_sources(ROOT)

    # Which sources are empty is asserted, not assumed. signal_references refuses an
    # empty file on purpose -- an empty file and a clean file grep identically -- so an
    # unnamed zero-byte source would abort the census rather than pass it. If
    # scoby_3d.jl ever gains content it enters the census and this line says so.
    @test basename.(filter(p -> filesize(p) == 0, sources)) == ["scoby_3d.jl"]
    live = filter(p -> filesize(p) > 0, sources)
    @test length(live) == 10

    # Regression, not a control. This passes today because the branch makes no CPM
    # edits, so its entire value is in the future edit it catches; nothing here shows
    # the census can fire. The planted line below is what shows that.
    for path in live
        @test signal_references(path) == Tuple{Int,String}[]
    end

    # The known-bad is recovered from the diagnostic's own driver, never typed here.
    # A needle I wrote would test my idea of the failure: I would have reached for
    # `InertSignal.foo`, which matches the alternation that already worked, and the
    # `update_signal!(` gap this control exists for would have survived the control.
    driver = joinpath(DIAG, "run.jl")
    calls = filter(l -> occursin("update_signal!(", l), readlines(driver))
    # Zero or two-plus means run.jl changed. That is a finding about the driver, not a
    # bug in this test -- fix the control, and do not reach for a remembered string.
    @test length(calls) == 1
    needle = get(calls, 1, "")

    mktempdir() do scratch
        # The plant enters by the production door: a .jl file at a root census_sources
        # walks, in a directory that did not exist when the census was written.
        write(joinpath(scratch, "biofilms_potts.jl"),
              read(joinpath(ROOT, "biofilms_potts.jl"), String) * "\n" * needle * "\n")
        found = census_sources(scratch)
        @test length(found) == 1
        refs = signal_references(only(found))
        @test length(refs) == 1
        @test occursin("update_signal!", last(only(refs)))
    end

    # Comment stripping is normalisation, and normalisation can delete the evidence.
    # All five measured against the real entry point. Three are characterisations of
    # known blind spots, marked as such: they assert the census is currently wrong, so
    # that closing a gap is a visible edit here rather than a silent behaviour change.
    mktempdir() do scratch
        probe(name, text) = signal_references(
            (p = joinpath(scratch, name * ".jl"); write(p, text * "\n"); p))
        @test !isempty(probe("trailing",  needle * " # drive the field"))
        @test  isempty(probe("commented", "# " * needle))
        # BLIND SPOT (false negative, load-bearing): a `#` in a string literal
        # truncates the line and hides the call that follows it.
        @test  isempty(probe("in_string", "msg = \"no # sign here\"; " * needle))
        # BLIND SPOT (false negative): a one-line block comment reads as a comment.
        @test  isempty(probe("block_one", "#= " * needle * " =#"))
        # BLIND SPOT (false positive, opposite direction): an interior line of a
        # multi-line block carries no `#`, so commented-out code is reported.
        @test !isempty(probe("block_many", "#=\n" * needle * "\n=#"))
    end
end
end

# The CairoMakie layout check: menu, colorbar, legend, slider and title callbacks,
# with the OpenGL window excluded. Needs Main.SR, so runtests.jl includes this file
# after load_serial().
module SignalLayoutTests
using Test, HDF5, SHA
const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "export_checkpoint.jl"))   # CLI-guarded; functions only
const SR = Main.SR

"""
    build_fixture(dir) -> (parent, derived)

101 frame pairs for `show_signal`. The parent snapshots come from the real exporter,
so a change to the snapshot schema reaches this fixture automatically. The companions
do not: `run.jl`'s writer needs a verified parent run and cannot be called here, so the
six attributes `signal_grid` reads are restated below. That restatement is the fixture's
one weakness and it is named in docs/visualization/inert_signal_validation.md -- if
`run.jl`'s companion format drifts, this keeps passing against the old shape and only
the manual `test_viewer.jl` against the real run notices.
"""
function build_fixture(dir)
    parent, derived = joinpath(dir, "parent"), joinpath(dir, "derived")
    mkpath(joinpath(parent, "snapshots")); mkpath(joinpath(derived, "fields"))
    p = SR.CPMParams(N = 12, n_cells_per_species = 2, snapshot_interval = 100)
    # No basis_gate_ack and no advance_window!: the initial configuration already has
    # every parcel on the lattice, which is all a layout check needs, and stepping is
    # the only thing the gate blocks. Acknowledging it here would have made this a
    # fourth ack site and failed the census in tests/radiodialysis_basis_gate.jl --
    # correctly, since a figure fixture has no business declaring a blocked basis.
    rp = SR.RadiolysisParams(Nr = 20, Ddot_R = 1.0, c_ext = 1.0)
    sim = SR.init_coupled_simulation(p, rp; seed = 9)
    seed = joinpath(dir, "seed.h5")
    export_transport_snapshot(SR, sim, seed)
    for t in 0:100
        stem = lpad(t, 6, '0')
        snap = joinpath(parent, "snapshots", "snap_mcs$(stem).h5")
        cp(seed, snap)
        h5open(snap, "r+") do f
            delete_attribute(f, "mcs"); HDF5.attributes(f)["mcs"] = t
        end
        mask = h5open(f -> read(f["lattice/interior_mask"]), snap, "r")
        A = zeros(size(mask))   # zero field: this checks figure wiring, not physics
        h5open(joinpath(derived, "fields", "signal_mcs$(stem).h5"), "w") do f
            f["fields/signal"] = A
            f["lattice/interior_mask"] = mask
            a = HDF5.attributes(f)
            a["mcs"] = t
            a["logical_axis_order"] = "xyz"
            a["dataset_axis_order_h5py"] = "zyx"
            a["parent_snapshot_sha256"] = bytes2hex(open(sha256, snap))
            a["signal_sha256"] = bytes2hex(sha256(reinterpret(UInt8, vec(A))))
            a["acceptance_coupling"] = 0
        end
    end
    parent, derived
end

mktempdir() do dir
    parent, derived = build_fixture(dir)
    # test_native_layout.jl reads the global ARGS so that it stays runnable by hand
    # against the real run. ARGS is global and append! mutates it, so restore by
    # snapshot in a finally: a failure mid-test must not leak paths into later files.
    saved = copy(ARGS)
    try
        empty!(ARGS); append!(ARGS, [parent, derived])
        include(joinpath(ROOT, "diagnostics", "inert_signal", "test_native_layout.jl"))
    finally
        empty!(ARGS); append!(ARGS, saved)
    end
end
end
