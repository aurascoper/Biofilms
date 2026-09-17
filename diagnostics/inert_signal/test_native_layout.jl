# Exercise the native figure's Makie objects/callbacks without an OpenGL window.
# This substitutes only the backend import and suppresses display, not the reader,
# geometry, transfer function, menu or timeline implementation.
using Test, CairoMakie
module NativeLayoutUnderTest end
path=normpath(joinpath(@__DIR__,"..","..","viewer","visualize_signal.jl"))
source=read(path,String)
@test count("using GLMakie",source)==1 && count("display(fig)",source)==1
adapted=replace(source,"using GLMakie"=>"using CairoMakie","display(fig)"=>"nothing")
cd(dirname(path)) do
    Base.include_string(NativeLayoutUnderTest,adapted,path)
end
Base.invokelatest() do
    fig=NativeLayoutUnderTest.show_signal(ARGS...)
    menu=only(filter(x->x isa Menu,fig.content))
    slider=only(filter(x->x isa Slider,fig.content))
    bar=only(filter(x->x isa Colorbar,fig.content))
    legend=only(filter(x->x isa Legend,fig.content))
    axis=only(filter(x->x isa Axis3,fig.content))
    @testset "Native Makie layout and callbacks; OpenGL window excluded" begin
        menu.selection[]="Species labels"
        @test !bar.blockscene.visible[] && legend.blockscene.visible[] && legend.scene.visible[]
        menu.selection[]="Signal volume"
        @test bar.blockscene.visible[] && !legend.blockscene.visible[] && !legend.scene.visible[]
        set_close_to!(slider,100)
        @test axis.title[]=="Inert generic signal - MCS 100"
        set_close_to!(slider,0)
        @test axis.title[]=="Inert generic signal - MCS 0"
    end
end
