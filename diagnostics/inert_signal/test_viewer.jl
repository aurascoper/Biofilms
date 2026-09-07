using Test
include(joinpath(@__DIR__,"..","..","viewer","signal_grid.jl"))

@testset "Native viewer reader and CLI without OpenGL" begin
    @test_throws ArgumentError signal_viewer_options(String[])
    @test_throws ArgumentError signal_viewer_options(["/nonexistent-parent","/nonexistent-signal"])
    parent,derived=ARGS
    for t in 0:100
        snapshot=joinpath(parent,"snapshots","snap_mcs$(lpad(t,6,'0')).h5")
        companion=joinpath(derived,"fields","signal_mcs$(lpad(t,6,'0')).h5")
        grid,A,mcs=signal_grid(snapshot,companion)
        @test mcs==t && size(A)==size(grid)==(40,40,40)
    end
    snapshot=joinpath(parent,"snapshots","snap_mcs000000.h5")
    @test_throws ArgumentError signal_grid(snapshot,joinpath(derived,"fields","signal_mcs000001.h5"))
    mktempdir() do tmp
        @test signal_viewer_options([parent,derived]).derived==derived
        source=joinpath(derived,"fields","signal_mcs000000.h5")
        bad=joinpath(tmp,"altered.h5");cp(source,bad)
        h5open(bad,"r+") do f
            @test read(f["fields/signal"])[20,20,20]==0
            f["fields/signal"][20,20,20]=1.
        end
        @test_throws ArgumentError signal_grid(snapshot,bad)
    end
end
