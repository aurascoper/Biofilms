using Test, LinearAlgebra, Random
include("InertSignal.jl")
using .InertSignal

@testset "Inert signal numerical contract" begin
    @test_throws ArgumentError SignalParams(diffusion=-1)
    @test_throws ArgumentError SignalParams(production=[1,2])
    @test_throws ArgumentError SignalParams(decay=NaN)
    @test_throws ArgumentError SignalParams(max_dt=0)
    mask = trues(7,7,7); mask[1,:,:] .= false
    species = zeros(Int,7,7,7); species[3,3,3] = 2
    mask_before, sp_before = copy(mask), copy(species)
    rng_before = copy(Random.default_rng())
    A = zeros(7,7,7)
    p = SignalParams(production=[0,2,0,0,0,0,0], diffusion=0, decay=0, max_dt=10)
    b = update_signal!(A, species, mask, p)
    @test A[3,3,3] == 2 && sum(A) == 2
    @test species == sp_before && mask == mask_before
    @test rand(copy(rng_before), 12) == rand(copy(Random.default_rng()), 12)
    @test b.source_mass == 2 && b.washout_mass == b.decay_mass == 0
    @test_throws ArgumentError update_signal!(A, species, mask, p; duration=-1)
    bad = copy(species); bad[3,3,3] = 8
    @test_throws ArgumentError update_signal!(A, bad, mask, p)
    badA = copy(A); badA[1,1,1] = 1
    @test_throws ArgumentError update_signal!(badA, species, mask, p)
    A .= 0
    p = SignalParams(production=zeros(7))
    update_signal!(A, species, mask, p; duration=100)
    @test all(iszero,A)

    # First-order time convergence to an independent analytic source-decay ODE.
    errors = Float64[]
    for dt in (0.2,0.1,0.05)
        a = zeros(1); s = [1]; m = trues(1)
        p = SignalParams(diffusion=0, decay=0.3, max_dt=dt)
        update_signal!(a,s,m,p; duration=2)
        push!(errors, abs(a[1] - (1-exp(-0.6))/0.3))
    end
    @test all(1.9 .< errors[1:2] ./ errors[2:3] .< 2.2)
    # Uniform concentration with Neumann outer faces: diffusion must vanish.
    # Seven distinct source rates exercise every lookup, unlike the uniform demo.
    for s in 1:7
        a=zeros(3,3,3); ids=fill(s,3,3,3)
        rates=collect(1.0:7.0)
        p=SignalParams(diffusion=0.2,decay=0.3,production=rates,max_dt=0.1)
        update_signal!(a,ids,trues(3,3,3),p; duration=1)
        expected=rates[s]/0.3*(1-(1-0.3*0.1)^10)
        @test all(isapprox.(a,expected;atol=2e-14))
    end

    # A deliberately large requested step must subcycle. The mass identity forbids
    # unaccounted source, sink or reflecting-wall changes, independent of clipping.
    a = zeros(7,7,7); a[2,3,3] = 3
    p = SignalParams(diffusion=2, decay=1, max_dt=100, production=fill(0.3,7))
    b = update_signal!(a,species,mask,p; duration=2)
    @test b.nsteps > 20
    @test minimum(a) >= 0 && all(iszero,a[.!mask])
    @test b.washout_mass > 0 && b.decay_mass > 0
    @test isapprox(b.final_mass-b.initial_mass,
                   b.source_mass-b.decay_mass-b.washout_mass; atol=2e-12)
    @test maximum(a) <= 3

    # Spatial convergence on the SAME uniformly producing slab [-L,L]. Boundaries
    # are explicit zero nodes; transverse dimensions are absent, not resized.
    L=3.0; D=0.7; gamma=0.4; prod=1.2
    errors = Float64[]
    for intervals in (12,24,48)
        h=2L/intervals
        x=collect(range(-L,L,length=intervals+1))
        a=zeros(length(x)); s=ones(Int,length(x)); m=trues(length(x))
        m[[1,end]].=false; s[[1,end]].=0
        p=SignalParams(diffusion=D,decay=gamma,production=fill(prod,7),spacing=h,max_dt=0.1)
        update_signal!(a,s,m,p; duration=70)
        truth=prod/gamma .* (1 .- cosh.(sqrt(gamma/D).*x)./cosh(sqrt(gamma/D)*L))
        push!(errors,maximum(abs.(a.-truth)))
    end
    @test all(3.8 .< errors[1:2]./errors[2:3] .< 4.2)
    threshold=prod/gamma*(1-1/cosh(sqrt(gamma/D)*L))
    @test slab_halfwidth(D,gamma,prod,threshold) ≈ L
    @test slab_halfwidth(D,gamma,prod,prod/gamma) == Inf
    @test slab_halfwidth(D,0,prod,threshold) ≈ sqrt(2D*threshold/prod)
    @test_throws ArgumentError slab_halfwidth(-1,gamma,prod,threshold)
    @test endpoint([5.,8.,10.],[1,0,-1],Bool[1,1,0],SignalParams()).fraction == 1
    @test endpoint([0.],[0],trues(1),SignalParams()).fraction === nothing
end

@testset "Production-source mutations are load-bearing" begin
    source=read(joinpath(@__DIR__,"InertSignal.jl"),String)
    function mutant(needle,replacement,name)
        @test count(needle,source)==1
        text=replace(source,needle=>replacement)
        @test text!=source && occursin(replacement,text)
        m=Module(name); Base.include_string(m,text,"mutated_InertSignal.jl")
        Base.invokelatest(getfield,m,:InertSignal)
    end
    # Bypass adaptive substeps and its positivity assertion in a private source
    # copy: an initial impulse genuinely goes negative under an oversized step.
    unstable=replace(source,"nsteps = max(1, ceil(Int, duration / dt_limit))"=>"nsteps = 1",
        "B[I] >= 0 || error(\"negative signal despite positivity substeps\")"=>"nothing")
    @test unstable!=source && occursin("nsteps = 1",unstable)
    m=Module(:NoStability);Base.include_string(m,unstable,"unstable.jl");bad=m.InertSignal
    a=zeros(3,3,3);a[2,2,2]=1
    bad.update_signal!(a,zeros(Int,3,3,3),trues(3,3,3),bad.SignalParams(diffusion=2,decay=1,production=zeros(7),max_dt=100))
    @test minimum(a)<0
    # Omitted decay must fail the same independent slab accuracy bound.
    wrong=mutant("diff + prod - p.decay * a","diff + prod",:NoDecay)
    function slab_error(mod)
        n=25;L=3.;h=2L/(n-1);x=range(-L,L,length=n)
        a=zeros(n);s=ones(Int,n);mask=trues(n);mask[[1,end]].=false;s[[1,end]].=0
        p=mod.SignalParams(diffusion=.7,decay=.4,production=fill(1.2,7),spacing=h,max_dt=.1)
        mod.update_signal!(a,s,mask,p;duration=70)
        maximum(abs.(a .- 3 .* (1 .- cosh.(sqrt(.4/.7).*x)./cosh(sqrt(.4/.7)*L))))
    end
    @test slab_error(InertSignal)<.004
    @test slab_error(wrong)>.5
end
