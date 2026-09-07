using Test
include("run.jl")
include("guard.jl")

@testset "Production boundaries and known-bad inputs" begin
    config=joinpath(@__DIR__,"demo.toml")
    @test read_config(config)[2].threshold == 5
    mktempdir() do tmp
        c=TOML.parsefile(config)
        # Unknown direct acknowledgements enter through the real config loader.
        c["basis_gate_ack"]=true
        path=joinpath(tmp,"bad.toml")
        open(io -> TOML.print(io,c),path,"w")
        @test_throws ArgumentError read_config(path)
        delete!(c,"basis_gate_ack"); c["acceptance_coupling"]=true
        open(io -> TOML.print(io,c),path,"w")
        @test_throws ArgumentError read_config(path)
        c["acceptance_coupling"]=false; c["physical_conversion"]="seconds"
        open(io -> TOML.print(io,c),path,"w")
        @test_throws ArgumentError read_config(path)

        serial=joinpath(@__DIR__,"..","..","biofilms_potts.jl")
        @test isempty(signal_references(serial))
        source=read(serial,String)
        start=findfirst("function mcs_step!",source)
        @test start !== nothing
        at=findnext("    N = ",source,last(start))
        @test at !== nothing
        planted=source[1:first(at)-1]*"    diagnostic_probe = state.signal\n"*source[first(at):end]
        @test planted != source
        @test occursin("diagnostic_probe = state.signal",split(planted,"function mcs_step!")[2])
        badserial=joinpath(tmp,"biofilms_potts.jl"); write(badserial,planted)
        refs=signal_references(badserial)
        @test length(refs) == 1
        @test only(refs)[2] == "diagnostic_probe = state.signal"
        # The guard's CLI is the production refusal, not just a helper assertion.
        proc=run(pipeline(ignorestatus(`$(Base.julia_cmd()) $(joinpath(@__DIR__,"guard.jl")) $badserial`),stdout=devnull,stderr=devnull))
        @test proc.exitcode != 0
    end

    if !isempty(ARGS)
        parent=only(ARGS); parent_sha=sha(joinpath(parent,"run_manifest.json"))
        @test length(verify_parent(parent,parent_sha)["snapshots"]) == 101
        mktempdir() do tmp
            # Symlink immutable artifacts; mutate only private manifest bytes.
            for name in readdir(parent)
                name == "run_manifest.json" && continue
                symlink(joinpath(parent,name),joinpath(tmp,name))
            end
            original=JSON3.read(read(joinpath(parent,"run_manifest.json"),String),Dict{String,Any})
            for defect in (:missing,:duplicate,:parcels,:hash)
                m=deepcopy(original)
                if defect == :missing
                    deleteat!(m["snapshots"],31)
                elseif defect == :duplicate
                    m["snapshots"][31]["mcs"]=29
                elseif defect == :parcels
                    m["configuration"]["initial_parcels"]=14
                else
                    victim=m["snapshots"][31]["path"]
                    @test haskey(m["artifacts"],victim)
                    m["artifacts"][victim]=repeat("0",64)
                end
                path=joinpath(tmp,"run_manifest.json"); write(path,JSON3.write(m))
                @test_throws ArgumentError verify_parent(tmp,sha(path))
            end
            @test_throws ArgumentError verify_parent(parent,repeat("0",64))
            @test_throws ArgumentError run_diagnostic(parent,tmp,parent_sha,config)
        end
    end
end
