using Test, Random, Statistics, Serialization, TOML
include("model.jl")
const F=FitnessPhase
const CONFIG=joinpath(@__DIR__,"config.toml")
c=F.read_config(CONFIG)
function altered(changes)
    d=TOML.parsefile(CONFIG);merge!(d,changes)
    p=tempname();open(io->TOML.print(io,d),p,"w");F.read_config(p)
end

@testset "restart after division preserves removal event order" begin
    # Reproduces the failed production seed, not a hand-made surrogate. At the
    # next MCS several parcels vanish mechanically. Serialization changes Dict
    # bucket order; enabled cleanup must not let that reorder the event log.
    a=F.new_run(c;seed=44,arm="baseline",growth=true)
    for _ in 1:100;F.step!(a);end
    checkpoint=tempname();F.save_checkpoint(checkpoint,a;source="regression")
    b=F.restore_checkpoint(checkpoint;source="regression",config_hash=c.hash)
    previous=length(a.state.death_log)
    for _ in 1:3
        F.step!(a);F.step!(b);@test F.state_hash(a;phase=true)==F.state_hash(b;phase=true)
    end
    @test length(a.state.death_log)>previous
end
@testset "configuration refusal before mutation" begin
    for (k,v) in [("phase_step_mcs",0.),("seconds_per_mcs",NaN),("axis_rate_Gy_s",-1.),("alpha_Gy_inv",[.1]),("backend","jacc"),("graph_rule","periodic"),("phase_noise",.1),("radiation_coupling","dose"),("density",zeros(7))]
        @test_throws Exception altered(Dict(k=>v))
    end
    r=F.new_run(c;seed=42,arm="combined");r.state.dose_rate_sd_Gy_s[1]=NaN;h=F.state_hash(r;phase=true)
    @test_throws Exception F.step!(r);@test h==F.state_hash(r;phase=true)
    r=F.new_run(c;seed=42,arm="combined");delete!(r.phase.theta,1);h=F.state_hash(r)
    @test_throws Exception F.step!(r);@test h==F.state_hash(r)
end
@testset "analytic hazard; incremental exposure; nonvacuous recovery" begin
    rng=MersenneTwister(9001);n=100000
    survivors=count(_->!F.hazard_draw(rng,log(2.),1.),1:n);phat=survivors/n;z=3.290526731
    center=(phat+z*z/(2n))/(1+z*z/n);half=z*sqrt(phat*(1-phat)/n+z*z/(4n*n))/(1+z*z/n)
    @test center-half<=0.5<=center+half
    @test !(center-half<=0.25<=center+half)
    r=F.new_run(c;seed=9001,arm="survival_only",fixed=true)
    for rate in (0.,0.001,0.002)
        r.state.dose_rate_mean_Gy_s.=rate;F.step!(r)
    end
    @test all(q->isapprox(q.dose,.18;atol=1e-12),values(r.parcels))
    @test all(q->q.status!=:alive||isapprox(q.hazard,c.alpha_Gy_inv[q.species]*.18;atol=1e-12),values(r.parcels))
    @test all(q->q.status!=:alive||isapprox(q.increment,.12;atol=1e-12),values(r.parcels))
    r0=F.new_run(altered(Dict("alpha_Gy_inv"=>zeros(7)));seed=42,arm="survival_only")
    F.step!(r0);@test isempty(r0.events)
end
@testset "phase analytic controls and bad graph rejection" begin
    A=[0. 1.;1. 0.];K=ones(2,2);theta=[.2,-.3];omega=[1.,-.1];dt=.01
    next=F.oscillator_step(theta,omega,A,K,dt)
    @test next[1]-next[2]≈theta[1]-theta[2]+dt*(1.1-2sin(.5))
    @test F.oscillator_step(theta,omega,zeros(2,2),K,dt)≈theta.+dt.*omega
    @test_throws Exception F.oscillator_step(theta,omega,[1. 1.;1. 0.],K,dt)
    @test_throws Exception F.oscillator_step(theta,omega,[0. 1.;0. 0.],K,dt)
    star=[0. 1. 1.;1. 0. 0.;1. 0. 0.];th=[0.,pi/2,pi/2]
    @test F.oscillator_step(th,zeros(3),star,ones(3,3),dt)[1]≈dt
    function evolve(h,om)
        x=[0.,0.]
        for _ in 1:round(Int,40/h);x=F.oscillator_step(x,om,A,K,h);end
        x[1]-x[2]
    end
    @test abs(evolve(.01,[.5,-.5])-asin(.5))<1e-8
    @test abs(evolve(.02,[.5,-.5])-evolve(.01,[.5,-.5]))<1e-8
    @test evolve(.01,[2.,-2.])>2pi
    function transient(h)
        x=[.5,-.5]
        for _ in 1:round(Int,1/h);x=F.oscillator_step(x,zeros(2),A,K,h);end
        x[1]-x[2]
    end
    exact=2atan(exp(-2)*tan(.5))
    @test 0<abs(transient(.025)-exact)<abs(transient(.05)-exact)<.02
    @test F.coherence([1,2,3],Dict(1=>0.,2=>0.,3=>pi))≈1/3
    @test F.coherence([1,2],Dict(1=>0.,2=>0.,3=>pi))≈1
    @test F.components(zeros(3,3))==([1,2,3],[1,1,1])
end
@testset "canonical paired paths; export and restart" begin
    for growth in (false,true), pair in (("baseline","phase_only"),("survival_only","combined"))
        a=F.new_run(c;seed=42,arm=pair[1],growth=growth)
        b=F.new_run(c;seed=42,arm=pair[2],growth=growth)
        for _ in 1:12
            F.step!(a);F.step!(b);@test F.state_hash(a)==F.state_hash(b)
        end
        if pair[1]=="survival_only";@test any(x->x.kind=="hazard_commitment",a.events);end
        broken=deepcopy(b);rand(broken.hazard_rng);@test F.state_hash(a)!=F.state_hash(broken)
        broken=deepcopy(b);broken.state.accumulated_dose_Gy[1]+=1;@test F.state_hash(a)!=F.state_hash(broken)
        path=tempname();F.save_checkpoint(path,b;source="test")
        resumed=F.restore_checkpoint(path;source="test",config_hash=c.hash)
        @test_throws Exception F.restore_checkpoint(path;source="wrong",config_hash=c.hash)
        for _ in 1:3
            F.step!(b);F.step!(resumed);@test F.state_hash(b;phase=true)==F.state_hash(resumed;phase=true)
        end
        before=F.state_hash(b;phase=true);F.export_frame(mktempdir(),b;source="test",run_id="test")
        @test before==F.state_hash(b;phase=true)
        if growth
            budget=F.budgets(b);@test abs(budget.nutrient_residual)<1e-7;@test abs(budget.biomass_residual)<1e-7
        end
    end
end
@testset "legacy mechanics and per-parcel energy" begin
    r=F.new_run(c;seed=42,arm="baseline");s=deepcopy(r.state);rng=deepcopy(r.movement_rng)
    for _ in 1:10
        F.CPM.mcs_step!(s,rng);F.CPM.update_melanin!(s);F.CPM.update_nutrient!(s);F.CPM.update_centers_of_mass!(s)
        F.step!(r)
        @test s.lattice==r.state.lattice && s.melanin==r.state.melanin && s.nutrient==r.state.nutrient && rng==r.movement_rng
    end
    r=F.new_run(c;seed=43,arm="survival_only");s=r.state;id=first(sort(collect(keys(s.cells))));r.parcels[id].target=0.
    source=Tuple(findfirst(==(id),s.lattice));target=Tuple(findfirst(==(0),s.lattice))
    before=sum(s.params.λ_V*(x.volume-r.parcels[i].target)^2 for (i,x) in s.cells)
    terms=F.energy_terms(r,s,source...,target...)
    after=sum(s.params.λ_V*(x.volume+(i==id ? 1 : 0)-r.parcels[i].target)^2 for (i,x) in s.cells)
    @test terms.vol==after-before
end
@testset "growth transfers and lifecycle" begin
    C=fill(2.,3,4,5);C[1,1,1]=8.;mask=trues(size(C));next=F.diffuse_closed(C,mask,.1,.1,(1.,2.,3.))
    @test sum(next)≈sum(C);@test minimum(next)>=0
    @test_throws Exception F.diffuse_closed(C,mask,100.,1.,(1.,1.,1.))
    after=F.nutrient_reaction(2.,.3,.5);@test after+.3*after/(.5+after)≈2.
    @test F.nutrient_reaction(0.,.3,.5)==0.
    r=F.new_run(c;seed=42,arm="combined",growth=true);r.state.current_mcs=10
    id=first(sort(collect(F.alive(r))));p=r.parcels[id];p.biomass=250.;p.target=250.;before=p.biomass
    event=F.divide!(r,id)
    @test event isa F.CPM.DivisionEvent
    if event isa F.CPM.DivisionEvent
        a,b=event.daughter_a_id,event.daughter_b_id
        @test r.parcels[a].biomass+r.parcels[b].biomass==before
        @test r.phase.theta[a]==r.phase.theta[b]==r.phase.theta[id]
        @test r.parcels[a].ancestral_dose==p.dose
        @test F.validate(r)===r
        r.parcels[a].status=:dying;r.parcels[a].target=0.;h=F.state_hash(r;phase=true)
        @test F.divide!(r,a)==:nonviable_parent;@test h==F.state_hash(r;phase=true)
    end
end

@testset "selection, fixed cohorts, slips and extinction" begin
    zero=altered(Dict("phase_K"=>0.))
    r=F.new_run(zero;seed=9002,arm="combined",fixed=true)
    ids=sort(collect(F.alive(r)));victim=ids[end]
    for id in ids
        r.phase.theta[id]=id==victim ? pi : 0.;r.phase.replay[id]=r.phase.theta[id];r.phase.omega[id]=0.
    end
    r.state.dose_rate_mean_Gy_s.=0.
    r.state.dose_rate_mean_Gy_s[r.state.lattice.==victim].=1e6
    records=NamedTuple[F.observation(r)]
    ph=F.step!(r);push!(records,F.observation(r;phase=ph))
    @test r.parcels[victim].status==:dying
    @test ph.selection>0 && ph.evolution==0 && ph.coherence≈1
    outcome=F.window_metrics(records;first_mcs=0,last_mcs=1,drift_bound=.01,range_bound=.5)
    @test outcome.count==length(ids)-1 && !outcome.coupling_supported
    r.state.dose_rate_mean_Gy_s.=1e6;ph=F.step!(r)
    @test isempty(ph.ids) && ph.coherence===nothing && ph.selection===nothing
    x=F.projection(r);@test all(iszero,x["phase_valid"]) && all(isnan,x["phase_rad"])
    fake=[(mcs=t,phase=(ids=[1,2],theta=[0.,2pi*t],replay=[0.,2pi*t])) for t in 0:4]
    outcome=F.window_metrics(fake;first_mcs=0,last_mcs=4,drift_bound=.01,range_bound=.5)
    @test !outcome.locking && outcome.coupled.slips[1]==4
    @test_throws Exception F.window_metrics(fake[1:end-1];first_mcs=0,last_mcs=4,drift_bound=.01,range_bound=.5)
end
