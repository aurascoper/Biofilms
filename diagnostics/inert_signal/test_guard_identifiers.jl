using Test
include(joinpath(@__DIR__,"guard.jl"))

# guard.jl names five identifiers. Only `.signal` was ever shown to fire -- test_pipeline.jl
# plants `diagnostic_probe = state.signal` into a copy of the serial source. The other four
# were asserted-by-regex and never watched refuse, and one of them did not: a single trailing
# \b after the alternation cannot hold for an alternative ending in `!`, so `update_signal!(state)`
# passed the census with a clean exit. These controls plant every guarded identifier, and
# require both the helper and the production CLI to refuse each one.
#
# Planting appends to the end of a copy of the real serial source. The census is lexical, so
# position carries no meaning; test_pipeline.jl plants inside mcs_step! because it is also
# asserting something about that function, which this file is not.

const SERIAL=normpath(joinpath(@__DIR__,"..","..","biofilms_potts.jl"))

# One realistic call site per guarded identifier.
const GUARDED=["InertSignal.solve!(state)",
               "update_signal!(state)",
               "autoinducer = 0.0",
               "signal_field[1] = 0.0",
               "diagnostic_probe = state.signal"]

# Names the census must NOT claim: near-misses, plurals, and the commented-out form.
const BENIGN=["signalling = true",
              "my_signal_fieldx = 1",
              "x.signals = 2",
              "signal = 1",
              "autoinducers = 2",
              "my_update_signal!(x)",
              "# InertSignal.solve!(state)"]

guard_exit(path)=run(pipeline(ignorestatus(
    `$(Base.julia_cmd()) $(joinpath(@__DIR__,"guard.jl")) $path`),
    stdout=devnull,stderr=devnull)).exitcode

@testset "Every guarded identifier fires the census and the CLI" begin
    source=read(SERIAL,String)
    # The control is only meaningful against a source the guard currently clears.
    @test isempty(signal_references(SERIAL))
    @test guard_exit(SERIAL) == 0
    mktempdir() do tmp
        for (i,stmt) in enumerate(GUARDED)
            path=joinpath(tmp,"planted_$i.jl"); write(path,source*"\n"*stmt*"\n")
            refs=signal_references(path)
            @test length(refs) == 1
            @test only(refs)[2] == stmt
            @test guard_exit(path) != 0
        end
    end
end

@testset "The census does not claim near-misses or commented code" begin
    source=read(SERIAL,String)
    mktempdir() do tmp
        for (i,stmt) in enumerate(BENIGN)
            path=joinpath(tmp,"benign_$i.jl"); write(path,source*"\n"*stmt*"\n")
            @test isempty(signal_references(path))
            @test guard_exit(path) == 0
        end
    end
end

@testset "The census refuses a missing or empty source" begin
    @test_throws ErrorException signal_references(joinpath(@__DIR__,"nonexistent.jl"))
    mktempdir() do tmp
        empty_path=joinpath(tmp,"empty.jl"); write(empty_path,"")
        @test_throws ErrorException signal_references(empty_path)
    end
end
