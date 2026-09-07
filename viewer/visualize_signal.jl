#!/usr/bin/env julia
# Native GLMakie view; OpenGL window execution must be checked on the target Mac/Linux host.
using GLMakie
include("signal_grid.jl")

function show_signal(parent,derived)
    frames=[signal_grid(joinpath(parent,"snapshots","snap_mcs$(lpad(t,6,'0')).h5"),
                        joinpath(derived,"fields","signal_mcs$(lpad(t,6,'0')).h5")) for t in 0:100]
    fig=Figure(size=(1100,800))
    ax=Axis3(fig[1,1];aspect=:data,xlabel="x (sites)",ylabel="y (sites)",zlabel="z (sites)",
             limits=(0,40,0,40,0,40),title="Inert generic signal - MCS 0")
    species=Observable(frames[1][1]);signal=Observable(Float32.(frames[1][2]))
    labelplot=voxels!(ax,0..40,0..40,0..40,species;
        color=parse.(Makie.Colorant,COLORS),is_air=(==(0x00)))
    transfer=[RGBAf(c.r,c.g,c.b,Float32(i==1 ? 0 : .015+.12*(i-1)/255))
              for (i,c) in enumerate(Makie.resample_cmap(:viridis,256))]
    fieldplot=volume!(ax,0..40,0..40,0..40,signal;algorithm=:absorption,
                     colormap=transfer,colorrange=(0,10))
    Colorbar(fig[1,2],fieldplot;label="A (declared arbitrary units)")
    slider=Slider(fig[2,1:2];range=0:100,startvalue=0)
    menu=Menu(fig[3,1];options=["Signal volume","Species labels"],default="Signal volume")
    labelplot.visible[]=false
    on(menu.selection) do value
        labelplot.visible[]=value=="Species labels"
        fieldplot.visible[]=value=="Signal volume"
    end
    on(slider.value) do t
        species[]=frames[t+1][1];signal[]=Float32.(frames[t+1][2])
        ax.title[]="Inert generic signal - MCS $t"
    end
    Label(fig[4,1:2],"42 initial parcels; seed 42. No physical time conversion or biological QS response.")
    display(fig)
    fig
end

if abspath(PROGRAM_FILE)==@__FILE__
    options=signal_viewer_options(ARGS)
    show_signal(options.parent,options.derived)
    wait(GLMakie.Screen())
end
