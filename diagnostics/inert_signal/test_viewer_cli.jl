using Test
include(joinpath(@__DIR__,"..","..","viewer","signal_grid.jl"))

# The two refusals in signal_viewer_options need no trajectory: they check argument
# count and directory existence, and read no file. They lived inside test_viewer.jl,
# which destructures `parent,derived=ARGS` on its third line, so the whole file was
# gated behind BIOFILMS_PARENT_RUN and BIOFILMS_DERIVED_RUN and these never ran in a
# default environment -- 2 of 106 assertions, and both of them negative controls,
# reported as UNCOVERED when nothing prevented covering them.
#
# They cannot be reached by calling test_viewer.jl argument-less: lines 11-12 would
# pass and line 13 would then raise BoundsError. Extraction is the fix, not invocation.
#
# The include chain needs HDF5 and SHA only -- no Makie, no OpenGL.

@testset "Viewer CLI refusals (no trajectory required)" begin
    @test_throws ArgumentError signal_viewer_options(String[])
    @test_throws ArgumentError signal_viewer_options(["/nonexistent-parent"])
    @test_throws ArgumentError signal_viewer_options(["/nonexistent-parent","/nonexistent-signal"])
    @test_throws ArgumentError signal_viewer_options([(@__DIR__), "/nonexistent-signal"])
    @test_throws ArgumentError signal_viewer_options(["/nonexistent-parent", (@__DIR__)])
    # Two existing directories is the accepted shape; the refusals above must not be
    # firing for some unrelated reason.
    here = @__DIR__
    opts = signal_viewer_options([here, here])
    @test opts.parent == here
    @test opts.derived == here
end
