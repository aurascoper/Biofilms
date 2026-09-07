#!/usr/bin/env julia
using Test, ReadVTK
include("run.jl")

function verify_derived(parent, output)
    path=joinpath(output,"derived_manifest.json")
    require(first(split(read(joinpath(output,"derived_manifest.sha256"),String))) == sha(path),"derived receipt mismatch")
    d=JSON3.read(read(path,String),Dict{String,Any})
    m=verify_parent(parent,d["parent_manifest_sha256"])
    for (rel,hash) in d["artifacts"]
        require(!isabspath(rel) && !(".." in splitpath(rel)),"unsafe artifact path")
        require(sha(joinpath(output,rel)) == hash,"derived artifact hash mismatch: $rel")
    end
    _,p=read_config(joinpath(output,"provenance","config.toml"))
    A=zeros(40,40,40)
    previous_species=nothing
    metrics=JSON3.read(read(joinpath(output,"metrics.json"),String))
    @test length(metrics)==101
    @testset "101 HDF5 and VTI frames: hashes, solver replay, axes and endpoints" begin
        for row in m["snapshots"]
            t=row["mcs"]; name="signal_mcs$(lpad(t,6,'0'))"
            ids,sp,mask=h5open(joinpath(parent,row["path"]),"r") do f
                read(f["lattice/cell_id"]),read(f["lattice/species_id"]),Bool.(read(f["lattice/interior_mask"]))
            end
            t>0 && update_signal!(A,previous_species,mask,p)
            previous_species=sp
            h5open(joinpath(output,"fields",name*".h5"),"r") do f
                a=attributes(f); B=read(f["fields/signal"])
                @test read(a["mcs"])==t
                @test read(a["parent_snapshot_sha256"])==row["sha256"]
                @test read(a["logical_axis_order"])=="xyz" && read(a["dataset_axis_order_h5py"])=="zyx"
                @test Bool.(read(f["lattice/interior_mask"]))==mask
                @test read(a["signal_sha256"])==bytes2hex(sha256(reinterpret(UInt8,vec(B))))
                @test B==A
                @test all(iszero,B[.!mask]) && minimum(B)>=0
            end
            vtk=VTKFile(joinpath(output,"paraview",name*".vti"))
            data=get_cell_data(vtk)
            @test reshape(get_data(data["signal"]),size(A))==A
            @test reshape(get_data(data["cell_id"]),size(A))==ids
            @test reshape(get_data(data["occupied_above_threshold"]),size(A))==UInt8.(mask .& (sp .> 0) .& (A .>= p.threshold))
            @test get_data(get_field_data(vtk)["mcs"])[1]==t
            e=endpoint(A,sp,mask,p)
            @test metrics[t+1].above_threshold==e.above_threshold
            @test metrics[t+1].occupied_sites==e.occupied_sites
            @test metrics[t+1].fraction==e.fraction
        end
    end
    println("Verified immutable parent and all 101 derived HDF5/VTI frames.")
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==2 || error("usage: verify.jl PARENT_RUN DERIVED_RUN")
    verify_derived(ARGS...)
end
