#!/usr/bin/env julia
# Synthetic metric controls plus two off-plan N16/p1/2MCS seeds through the real producer.
using Test, Printf, Statistics
module Producer
include(joinpath(@__DIR__, "clustering_sweep.jl"))
end
module Analysis
include(joinpath(@__DIR__, "clustering_analyse.jl"))
end

function protocol_mismatches(text)
    declared = Dict(m[1] => parse(Float64,m[2]) for m in eachmatch(r"\| (CONTACT_\w+) \| ([0-9.]+) \|",text))
    expected = Dict("CONTACT_ALPHA"=>Analysis.CONTACT_ALPHA, "CONTACT_MEOI"=>Analysis.CONTACT_MEOI,
                    "CONTACT_RADIUS_BOUND"=>Analysis.CONTACT_RADIUS_BOUND, "CONTACT_Z90"=>Analysis.CONTACT_Z90)
    Set(keys(declared)) == Set(keys(expected)) || return Set(["missing/extra constants"])
    Set(key for (key,value) in expected if declared[key] != value)
end

@testset "decision constants match the preregistration" begin
    protocol = read(Producer.CLUSTER_PREDICTION,String)
    @test isempty(protocol_mismatches(protocol))
    victim = only(collect(eachmatch(r"\| CONTACT_MEOI \| ([0-9.]+) \|",protocol)))
    mutation = replace(protocol, victim.match=>"| CONTACT_MEOI | 0.01 |"; count=1)
    @test mutation != protocol
    @test protocol_mismatches(mutation) == Set(["CONTACT_MEOI"])
end

@testset "contact observable separates parcel identity from species and shape" begin
    lattice = zeros(Int, 7, 7, 7)
    lattice[3:4, 3, 3] .= 1
    one = Producer.parcel_contacts(lattice, Dict(1=>1), 2)
    @test one[1].boundary == [10]
    @test one[1].same == [0]
    @test one[1].contact == 0
    @test isnan(one[2].contact) # an absent species is not zero contact

    lattice[4, 3, 3] = 2
    contact = Producer.parcel_contacts(lattice, Dict(1=>1, 2=>1), 1)[1]
    @test contact.ids == [1, 2]
    @test contact.boundary == [6, 6]
    @test contact.same == [1, 1]
    @test contact.contact == 1/6
    different = Producer.parcel_contacts(lattice, Dict(1=>1, 2=>2), 2)
    @test [x.contact for x in different] == [0, 0]
    @test [x.boundary for x in different] == [[6], [6]]

    lattice[4, 3, 3] = 0
    lattice[4, 4, 3] = 2 # touching only along an edge: a 26-neighbor metric would count it
    diagonal = Producer.parcel_contacts(lattice, Dict(1=>1, 2=>1), 1)[1]
    @test diagonal.same == [0, 0]
    @test diagonal.contact == 0

    lattice .= 0
    lattice[2, 3, 3] = 1
    lattice[3:4, 3, 3] .= 2
    unequal = Producer.parcel_contacts(lattice, Dict(1=>1, 2=>1), 1)[1]
    @test unequal.volumes == [1, 2]
    @test unequal.boundary == [6, 10]
    @test unequal.same == [1, 1]
    @test unequal.contact ≈ 2/15
    @test unequal.contact != sum(unequal.same)/sum(unequal.boundary)
    refined = Producer.parcel_contacts(repeat(lattice; inner=(2,2,2)), Dict(1=>1, 2=>1), 1)[1]
    @test refined.boundary == 4 .* unequal.boundary
    @test refined.same == 4 .* unequal.same
    @test refined.volumes == 8 .* unequal.volumes
    @test refined.contact == unequal.contact

    lattice .= 0
    lattice[1, 1, 1] = 1
    lattice[2, 1, 1] = -1
    edge = Producer.parcel_contacts(lattice, Dict(1=>1), 1)[1]
    @test edge.boundary == [6] # exterior and wall both contribute to exposure
    @test edge.same == [0]
    lattice[1, 1, 2] = 9
    @test_throws r"orphan" Producer.parcel_contacts(lattice, Dict(1=>1), 1)
    @test_throws r"invalid species" Producer.parcel_contacts(lattice, Dict(1=>1, 9=>3), 2)
end

@testset "each registered decision condition can independently refuse" begin
    @test Analysis.contact_sign_probability(16,16) == 1//big(2)^16
    @test Analysis.contact_sign_probability(15,16) == 17//big(2)^16
    @test Analysis.contact_sign_probability(0,0) == 1
    @test Analysis.contact_sign_probability(0,16) == 1
    supported = Analysis.contact_summary(fill(.06,16), repeat([-.001,.001],8))
    @test (supported.size_pass, supported.sign_pass, supported.radius_pass, supported.passes) == (true,true,true,true)
    small = Analysis.contact_summary(fill(.01,16), zeros(16))
    @test (small.size_pass, small.sign_pass, small.radius_pass, small.passes) == (false,true,true,false)
    inconsistent = Analysis.contact_summary(vcat(fill(.2,7),fill(-.001,9)), zeros(16))
    @test (inconsistent.size_pass, inconsistent.sign_pass, inconsistent.radius_pass, inconsistent.passes) == (true,false,true,false)
    shifted = Analysis.contact_summary(fill(.06,16), fill(.08,16))
    @test (shifted.size_pass, shifted.sign_pass, shifted.radius_pass, shifted.passes) == (true,true,false,false)
    boundary = Analysis.contact_summary(fill(.06,16), fill(.05,16))
    @test !boundary.radius_pass
    ties = Analysis.contact_summary(zeros(16), zeros(16))
    @test (ties.k,ties.n,ties.ties,ties.p,ties.passes) == (0,0,16,1,false)
    @test_throws r"nonfinite" Analysis.contact_summary([NaN,0.1], zeros(2))
end

"Mutate the actual producer row identified by its key, and prove only that row changed."
function changed_csv(source, dest, key, field, value)
    lines = readlines(source)
    header = findfirst(line -> startswith(line,"seed,arm,mcs,species,"), lines)
    columns = split(lines[header], ',')
    col = findfirst(==(field), columns)
    isnothing(col) && error("mutation field absent")
    victims = findall(i -> begin
        startswith(lines[i],"#") && return false
        fields = split(lines[i], ',')
        length(fields) == length(columns) || return false
        join(fields[1:4], ',') == key
    end, eachindex(lines))
    length(victims) == 1 || error("mutation key missing or ambiguous")
    victim = only(victims)
    fields = String.(split(lines[victim], ','))
    fields[col] != value || error("mutation is a no-op")
    fields[col] = value
    mutated = copy(lines)
    mutated[victim] = join(fields, ',')
    @test findall(mutated .!= lines) == [victim]
    write(dest, join(mutated, '\n') * "\n")
end

@testset "real producer and reader, off-plan smoke only" begin
    mktempdir() do dir
        reference = joinpath(dir,"radial.csv")
        # These are NOT planned study seeds or the study configuration.
        Producer.ClusterRadial.sweep_main([reference,"--seeds","501,502","--n","16",
            "--parcels","1","--at","2","--arms","full,no_mel"])
        out = joinpath(dir,"contact.csv")
        args = [out,"--smoke","--seeds","501,502","--n","16","--parcels","1",
                "--at","2","--reference",reference]
        Producer.clustering_sweep(args)
        @test isfile(out) && !isfile(out * ".partial")
        reader(path; complete=false) = Analysis.read_clustering([path]; complete, smoke=true, reference_path=reference)
        data = reader(out)
        @test data.seeds == [501,502]
        @test length(data.rows) == 56
        @test data.meta["mode"] == "smoke"
        @test_throws r"every seed 42:297" reader(out; complete=true)
        @test_throws r"destination already exists" Producer.clustering_sweep(args)
        @test_throws r"incompatible run mode" Analysis.read_clustering([out]; complete=false, reference_path=reference)
        @test_throws r"duplicate pairing across shards" Analysis.read_clustering([out,out]; complete=false, smoke=true, reference_path=reference)

        original = readlines(out)
        victim = findfirst(line -> startswith(line,"501,no_mel,2,1,"), original)
        @test !isnothing(victim)
        missing = joinpath(dir,"missing.csv")
        shortened = original[setdiff(eachindex(original),[victim])]
        @test length(shortened) == length(original)-1
        write(missing,join(shortened,'\n')*"\n")
        @test_throws r"missing or unexpected" reader(missing)

        duplicate = joinpath(dir,"duplicate.csv")
        write(duplicate,join(vcat(original,original[victim]),'\n')*"\n")
        @test_throws r"duplicate row" reader(duplicate)

        badseed = joinpath(dir,"badseed.csv")
        changed_csv(out,badseed,"501,no_mel,2,1","seed","503")
        @test_throws r"missing or unexpected" reader(badseed)
        placement = joinpath(dir,"placement.csv")
        changed_csv(out,placement,"501,no_mel,0,1","placement_sha256",repeat("f",64))
        @test_throws r"do not share initial placement" reader(placement)
        radius = joinpath(dir,"radius.csv")
        original_radius = parse(Float64,data.rows[(501,"full",2,1)]["mean_r_over_R"])
        changed_csv(out,radius,"501,full,2,1","mean_r_over_R",@sprintf("%.6f",original_radius+.01))
        @test_throws r"baseline mismatch" reader(radius)
        wrongmean = joinpath(dir,"wrongmean.csv")
        changed_csv(out,wrongmean,"501,full,2,1","mean_contact_fraction","0.125")
        @test_throws r"contact mean disagrees" reader(wrongmean)
        nonfinite = joinpath(dir,"nonfinite.csv")
        changed_csv(out,nonfinite,"501,full,2,1","mean_contact_fraction","NaN")
        @test_throws r"contact mean disagrees" reader(nonfinite)
        incomplete = joinpath(dir,"incomplete.csv")
        @test startswith(last(original),"# baseline_checked_rows=")
        write(incomplete,join(original[1:end-1],'\n')*"\n")
        @test_throws r"incomplete shard" reader(incomplete)
    end
end
