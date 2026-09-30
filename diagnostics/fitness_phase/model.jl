module FitnessPhase
using Random, Statistics, LinearAlgebra, Serialization, SHA, TOML, JSON3
const ROOT = normpath(joinpath(@__DIR__, "../.."))
module CPM
    const source_path = normpath(joinpath(@__DIR__, "../../biofilms_potts.jl"))
    include_string(@__MODULE__, split(read(source_path, String), "#  13. Figure export")[1], source_path)
end
include("../synthetic_response/config.jl")

struct Config
    values::NamedTuple
    hash::String
end
Base.getproperty(c::Config, k::Symbol) = k in (:values, :hash) ? getfield(c,k) : getproperty(getfield(c,:values),k)
function read_config(path)
    d = TOML.parsefile(path)
    required = ["schema","N","parcels_per_species","species_order","backend","boundary","length_unit","spacing","origin","seconds_per_mcs","axis_rate_Gy_s","alpha_Gy_inv","physical_field","graph_rule","radiation_coupling","phase_noise","phase_step_mcs","phase_omega_mean","phase_omega_sd","phase_K","coupling_mode","nutrient_initial","nutrient_diffusion","nutrient_step_mcs","uptake_max","monod_K","yield","maintenance","density","division_mass","minimum_age_mcs","minimum_daughter_voxels","parameter_status"]
    Set(keys(d)) == Set(required) || error("configuration keys differ from schema")
    d["schema"] == "cpm-sp-config-1" || error("unsupported schema")
    d["species_order"] == CPM.SPECIES_NAMES || error("species order mismatch")
    d["backend"] == "serial" || error("extension has no verified JACC backend")
    d["boundary"] == "legacy_nonperiodic_cylinder" || error("unsupported topology")
    d["radiation_coupling"] == "none" && d["phase_noise"] == 0 || error("forced/noisy phase is unsupported")
    d["graph_rule"] in ("contact26", "face6") || error("unsupported graph")
    d["coupling_mode"] == "nonnegative" || error("signed coupling unsupported")
    d["physical_field"] in ("gradient", "uniform") || error("unsupported physical field")
    d["parameter_status"] == "synthetic" || error("no biological calibration is established")
    d["length_unit"] == "synthetic_length" || error("no physical pitch is calibrated")
    for k in ("N", "parcels_per_species", "minimum_age_mcs", "minimum_daughter_voxels")
        d[k] isa Integer && d[k] > 0 || error("$k must be a positive integer")
    end
    d["N"] >= 16 || error("legacy placement requires N >= 16")
    for k in ("spacing", "origin")
        length(d[k]) == 3 && all(isfinite, d[k]) || error("invalid $k")
    end
    all(>(0), d["spacing"]) || error("nonpositive spacing")
    for k in ("alpha_Gy_inv","uptake_max","monod_K","yield","maintenance","density","division_mass")
        length(d[k]) == CPM.N_SPECIES && all(x -> isfinite(x) && x >= 0, d[k]) || error("invalid species vector $k")
    end
    for k in ("monod_K","yield","density","division_mass")
        all(>(0), d[k]) || error("$k must be positive")
    end
    for k in ("seconds_per_mcs","axis_rate_Gy_s","phase_noise","phase_step_mcs","phase_omega_mean","phase_omega_sd","phase_K","nutrient_initial","nutrient_diffusion","nutrient_step_mcs")
        d[k] isa Real && isfinite(d[k]) && d[k] >= 0 || error("invalid $k")
    end
    for k in ("seconds_per_mcs","phase_step_mcs","nutrient_step_mcs")
        d[k] > 0 || error("nonpositive $k")
    end
    for k in ("phase_step_mcs","nutrient_step_mcs")
        n = round(Int, 1/d[k]); n >= 1 && isapprox(n*d[k],1; atol=1e-12) || error("$k must divide one MCS")
    end
    2*d["nutrient_diffusion"]*d["nutrient_step_mcs"]*sum(1 ./ d["spacing"].^2) <= 1 || error("diffusion positivity bound exceeded")
    maximum(d["maintenance"])*d["nutrient_step_mcs"] < 1 || error("maintenance positivity bound exceeded")
    SyntheticResponse.ResponseConfig(alpha_Gy_inv=d["alpha_Gy_inv"],seconds_per_mcs=d["seconds_per_mcs"],axis_rate_Gy_s=d["axis_rate_Gy_s"])
    ks = sort(collect(keys(d)))
    immutable(x) = x isa AbstractVector ? Tuple(immutable.(x)) : x
    Config(NamedTuple{Tuple(Symbol.(ks))}(Tuple(immutable(d[k]) for k in ks)),bytes2hex(sha256(read(path))))
end

mutable struct Parcel
    species::Int
    parent::Int
    lineage::Int
    birth::Int
    status::Symbol
    death::Int
    removal::Int
    initial_r::Float64
    biomass::Float64
    target::Float64
    dose::Float64
    increment::Float64
    hazard::Float64
    commitment_dose::Float64
    ancestral_dose::Float64
    uptake::Float64
    growth_rate::Float64
    hazard_rate::Float64
end
mutable struct Phase
    theta::Dict{Int,Float64} # unwrapped; wrap only when projecting for display
    omega::Dict{Int,Float64}
    replay::Dict{Int,Float64}
end
mutable struct Run
    state::CPM.CPMState
    config::Config
    seed::Int
    arm::String
    growth::Bool
    fixed::Bool
    movement_rng::MersenneTwister
    hazard_rng::MersenneTwister
    parcels::Dict{Int,Parcel}
    phase::Union{Nothing,Phase}
    events::Vector{NamedTuple}
    budget::Dict{String,Float64}
end
alive(r::Run) = Set(id for (id,p) in r.parcels if p.status == :alive)
survival(r::Run) = r.arm in ("survival_only","combined")
model_id(r::Run) = r.growth ? "cpm-sp-growth-v1" : "cpm-sp-v1"
function stable_seed(seed, tag)
    b=sha256("cpm-sp-seed-v1:$seed:$tag")
    sum(UInt64(b[i]) << (8*(i-1)) for i in 1:8)
end
function new_run(c::Config; seed::Int, arm::String, growth::Bool=false, fixed::Bool=false, tropism::Bool=true)
    arm in ("baseline","survival_only","phase_only","combined") || error("unknown arm")
    fixed && growth && error("fixed-position control excludes growth")
    beta = tropism ? CPM.CPMParams().β_ion : zeros(CPM.N_SPECIES)
    p=CPM.CPMParams(N=c.N,n_cells_per_species=c.parcels_per_species,seconds_per_mcs=c.seconds_per_mcs,β_ion=beta)
    s=CPM.init_state(p; seed=seed); CPM.update_centers_of_mass!(s)
    length(s.cells) == CPM.N_SPECIES*c.parcels_per_species || error("initialisation failed to place all parcels")
    rate=c.physical_field == "gradient" ? c.axis_rate_Gy_s .* s.radiation ./ p.I0 : fill(c.axis_rate_Gy_s,size(s.radiation))
    CPM.import_dose_field!(s,rate,zeros(size(rate)))
    s.physical_time_s=0.0
    growth && (s.nutrient .= ifelse.(s.interior,c.nutrient_initial,0.0))
    volume=prod(c.spacing)
    parcels=Dict(id=>Parcel(x.species,x.parent_id,x.lineage_id,x.birth_mcs,:alive,-1,-1,
        CPM.radial_dist(x.com[1],x.com[2],c.N)/(c.N/2),
        growth ? p.V_target*c.density[x.species]*volume : NaN,Float64(p.V_target),0.,0.,0.,NaN,0.,0.,0.,0.) for (id,x) in s.cells)
    ph=nothing
    if arm in ("phase_only","combined")
        prng=MersenneTwister(stable_seed(seed,"phase-initialisation"))
        ids=sort(collect(keys(parcels)))
        theta=Dict(id=>2pi*rand(prng) for id in ids)
        omega=Dict(id=>c.phase_omega_mean+c.phase_omega_sd*randn(prng) for id in ids)
        ph=Phase(theta,omega,copy(theta))
    end
    b=Dict("initial_nutrient"=>sum(s.nutrient[s.interior])*volume,"initial_biomass"=>growth ? sum(x.biomass for x in values(parcels)) : 0.,"uptake"=>0.,"production"=>0.,"maintenance"=>0.,"removed_biomass"=>0.)
    r=Run(s,c,seed,arm,growth,fixed,MersenneTwister(seed),MersenneTwister(stable_seed(seed,"mortality")),parcels,ph,NamedTuple[],b)
    refresh_rates!(r); validate(r); r
end

function validate(r::Run)
    s=r.state; c=r.config; shape=(c.N,c.N,c.N)
    size(s.lattice)==shape && size(s.interior)==shape || error("shape mismatch")
    for f in (:radiation,:melanin,:nutrient,:dose_rate_mean_Gy_s,:dose_rate_sd_Gy_s,:dose_increment_Gy,:accumulated_dose_Gy,:melanin_drive)
        a=getfield(s,f); size(a)==shape && all(x->isfinite(x)&&x>=0,a) || error("invalid field $f")
    end
    CPM._check_dose_accrual(s)
    isfinite(maximum(s.dose_rate_mean_Gy_s)*c.seconds_per_mcs*maximum(c.alpha_Gy_inv)) || error("hazard increment overflow")
    counts=Dict(id=>0 for id in keys(s.cells))
    for k in eachindex(s.lattice)
        id=Int(s.lattice[k])
        ((s.interior[k] && id>=0) || (!s.interior[k] && id == -1)) || error("illegal mask or ID")
        if id>0
            haskey(counts,id) || error("orphan ID $id")
            counts[id]+=1
        end
    end
    for (id,x) in s.cells
        counts[id]==x.volume && x.volume>0 || error("registry volume mismatch")
        haskey(r.parcels,id) || error("missing extension ID")
        q=r.parcels[id]
        q.status in (:alive,:dying) && q.species==x.species || error("inconsistent parcel")
        q.status==:dying && q.target!=0 && error("dying parcel has nonzero target")
        all(isfinite,(q.target,q.dose,q.increment,q.hazard)) && min(q.target,q.dose,q.increment,q.hazard)>=0 || error("invalid parcel values")
        r.growth && !(isfinite(q.biomass)&&q.biomass>0) && error("invalid biomass")
    end
    for (id,q) in r.parcels
        q.status in (:alive,:dying) && !haskey(s.cells,id) && error("missing occupied parcel")
    end
    s.next_id>maximum(keys(r.parcels);init=0) || error("ID reuse")
    if r.phase!==nothing
        Set(keys(r.phase.theta))==Set(keys(r.parcels))==Set(keys(r.phase.omega))==Set(keys(r.phase.replay)) || error("phase ID mismatch")
        all(isfinite,values(r.phase.theta)) && all(isfinite,values(r.phase.omega)) && all(isfinite,values(r.phase.replay)) || error("invalid phase")
    end
    r
end

function energy_terms(r::Run, s, sx,sy,sz,tx,ty,tz)
    t=CPM.compute_delta_H_terms(s,sx,sy,sz,tx,ty,tz)
    source=Int(s.lattice[sx,sy,sz]); target=Int(s.lattice[tx,ty,tz])
    v=0.0
    for (id,change) in ((source,1),(target,-1))
        id>0 || continue
        x=s.cells[id]; goal=r.parcels[id].target
        v+=s.params.λ_V*((x.volume+change-goal)^2-(x.volume-goal)^2)
    end
    (adh=t.adh,vol=v,rad=t.rad,mel=t.mel)
end
function event!(r,kind,id; detail="")
    push!(r.events,(mcs=r.state.current_mcs,kind=String(kind),id=id,detail=detail))
end
function hazard_draw(rng,alpha,increment)
    all(isfinite,(alpha,increment)) && min(alpha,increment)>=0 || error("invalid hazard inputs")
    dh=alpha*increment;isfinite(dh)||error("hazard overflow")
    dh>0 && rand(rng)<-expm1(-dh)
end
function record_loss!(r,id,cause)
    q=r.parcels[id]
    q.status=:removed; q.removal=r.state.current_mcs
    r.growth && (r.budget["removed_biomass"]+=q.biomass)
    event!(r,cause,id)
end

include("phase.jl")
include("growth.jl")

function refresh_rates!(r::Run)
    s=r.state;c=r.config
    for q in values(r.parcels);q.hazard_rate=0.;q.growth_rate=0.;end
    for k in eachindex(s.lattice)
        id=Int(s.lattice[k]);id>0 || continue
        q=r.parcels[id];q.status==:alive || continue
        sp=q.species;v=s.cells[id].volume
        survival(r) && (q.hazard_rate+=c.alpha_Gy_inv[sp]*s.dose_rate_mean_Gy_s[k]*c.seconds_per_mcs/v)
        r.growth && (q.growth_rate+=c.yield[sp]*c.uptake_max[sp]*s.nutrient[k]/(c.monod_K[sp]+s.nutrient[k])/v)
    end
    if r.growth
        for q in values(r.parcels);q.status==:alive && (q.growth_rate-=c.maintenance[q.species]);end
    end
end

function step!(r::Run)
    validate(r) # invalid input must fail before any mutation
    s=r.state; c=r.config; oldids=alive(r)
    oldtheta=r.phase===nothing ? Dict{Int,Float64}() : copy(r.phase.theta)
    if r.fixed
        s.current_mcs+=1
    elseif survival(r) || r.growth
        CPM.mcs_step!(s,r.movement_rng;energy_terms=(args...)->energy_terms(r,args...))
    else
        CPM.mcs_step!(s,r.movement_rng)
    end
    for id in sort(collect(keys(r.parcels)))
        q=r.parcels[id]
        q.status in (:alive,:dying) && !haskey(s.cells,id) && record_loss!(r,id,q.status==:dying ? "resorbed" : "mechanical_loss")
    end
    live=alive(r)
    CPM.update_melanin!(s;active=survival(r) ? live : nothing)
    if r.growth
        growth_step!(r)
    else
        CPM.update_nutrient!(s;active=survival(r) ? live : nothing)
    end
    # Store the current occupied-site increment independently of cumulative dose.
    increments=Dict(id=>0.0 for id in keys(s.cells))
    for k in eachindex(s.lattice)
        id=Int(s.lattice[k]); id>0 || continue
        increments[id]+=s.dose_rate_mean_Gy_s[k]*c.seconds_per_mcs/s.cells[id].volume
    end
    CPM.accrue_dose!(s); s.physical_time_s=s.current_mcs*c.seconds_per_mcs
    for id in sort(collect(keys(s.cells)))
        q=r.parcels[id]; q.dose=s.cells[id].lifetime_dose_Gy; q.increment=increments[id]
        q.status==:alive || continue
        dh=survival(r) ? c.alpha_Gy_inv[q.species]*q.increment : 0.0
        q.hazard_rate=dh # per MCS
        q.hazard+=dh
        if hazard_draw(r.hazard_rng,survival(r) ? c.alpha_Gy_inv[q.species] : 0.,q.increment)
            q.status=:dying; q.death=s.current_mcs; q.target=0.; q.commitment_dose=q.dose
            event!(r,"hazard_commitment",id;detail="dose=$(q.dose);increment=$(q.increment)")
        end
    end
    if !r.fixed
        for id in sort(collect(keys(s.cells)))
            q=r.parcels[id]
            if q.status==:dying && s.cells[id].volume<2
                x=s.cells[id]; s.lattice[s.lattice.==id].=0
                push!(s.archive,CPM.ArchivedCell(id,x.species,x.lineage_id,x.parent_id,x.generation,x.birth_mcs,s.current_mcs,:died,x.lifetime_dose_Gy))
                push!(s.death_log,CPM.DeathEvent(s.current_mcs,s.physical_time_s,id,x.species,x.lineage_id,x.generation,x.lifetime_dose_Gy))
                delete!(s.cells,id); record_loss!(r,id,"resorbed")
            end
        end
    end
    CPM.update_centers_of_mass!(s)
    afterloss=alive(r)
    selection=r.phase===nothing ? nothing : coherence(afterloss,oldtheta)-coherence(oldids,oldtheta)
    # Extinction is missing, never an arithmetic zero.
    (isempty(oldids)||isempty(afterloss)) && (selection=nothing)
    r.growth && schedule_divisions!(r)
    # End-boundary rates use the current ownership, including newborn daughters.
    refresh_rates!(r)
    phasestats=phase_step!(r,oldtheta,afterloss,selection)
    validate(r)
    phasestats
end

function canonical(r::Run; phase=false)
    s=r.state
    state=Tuple(f==:cells ? [(id,Tuple(getfield(s.cells[id],k) for k in fieldnames(CPM.CellInfo))) for id in sort(collect(keys(s.cells)))] : getfield(s,f) for f in fieldnames(CPM.CPMState))
    extra=[(id,Tuple(getfield(r.parcels[id],k) for k in fieldnames(Parcel))) for id in sort(collect(keys(r.parcels)))]
    p=phase && r.phase!==nothing ? [(id,r.phase.theta[id],r.phase.omega[id],r.phase.replay[id]) for id in sort(collect(keys(r.parcels)))] : nothing
    (r.config.values,r.seed,r.growth,r.fixed,state,extra,r.movement_rng,r.hazard_rng,r.events,sort(collect(r.budget)),p)
end
function digest(x)
    io=IOBuffer(); serialize(io,x); bytes2hex(sha256(take!(io)))
end
state_hash(r;phase=false)=digest(canonical(r;phase=phase))
function save_checkpoint(path,r::Run; source)
    validate(r); mkpath(dirname(path))
    ispath(path) && error("checkpoint already exists")
    tmp=path*".partial"
    open(tmp,"w") do io
        serialize(io,(schema="cpm-sp-restart-1",julia=string(VERSION),source=source,config=r.config.hash,hash=state_hash(r;phase=true),run=r))
    end
    mv(tmp,path); path
end
function restore_checkpoint(path; config_hash,source)
    x=open(deserialize,path)
    x.schema=="cpm-sp-restart-1" && x.julia==string(VERSION) && x.source==source && x.config==config_hash || error("checkpoint identity mismatch")
    x.run.config.hash==config_hash || error("checkpoint configuration mismatch")
    validate(x.run); state_hash(x.run;phase=true)==x.hash || error("checkpoint contents mismatch")
    x.run
end
include("export.jl")
end
