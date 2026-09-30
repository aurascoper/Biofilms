using WriteVTK
const EXPORT_SCHEMA="cpm-sp-vti-1"
const UNITS=Dict("dose_rate"=>"Gy/s (synthetic exposure)","ambient_cum_dose"=>"Gy", "parcel_cum_dose"=>"Gy",
 "cumulative_hazard"=>"dimensionless", "conditional_survival"=>"probability", "biomass"=>"synthetic mass; parcel total repeated per voxel",
 "biomass_density"=>"synthetic mass/synthetic length^3", "target_volume"=>"synthetic length^3", "growth_rate"=>"MCS^-1",
 "mortality_hazard"=>"MCS^-1", "net_biomass_rate_proxy"=>"MCS^-1", "phase_rad"=>"radian", "phase_cos"=>"dimensionless",
 "phase_sin"=>"dimensionless", "natural_frequency"=>"radian/MCS", "local_coherence"=>"dimensionless", "contact_strength"=>"Boolean degree for contact26; face area for face6")
function atomic_json(path,object)
    tmp=path*".partial"; open(io->JSON3.write(io,object),tmp,"w");mv(tmp,path;force=true)
end
function projection(r::Run)
    s=r.state;shape=size(s.lattice);d=Dict{String,Any}()
    for key in ("cell_id","species_id","parent_id","lineage_id","graph_degree")
        d[key]=zeros(Int64,shape)
    end
    for key in ("interior_mask","occupied","alive","phase_valid","biomass_valid","local_coherence_valid")
        d[key]=zeros(UInt8,shape)
    end
    d["growth_enabled"]=fill(UInt8(r.growth),shape)
    d["survival_enabled"]=fill(UInt8(survival(r)),shape)
    for key in keys(UNITS); d[key]=fill(NaN,shape);end
    d["interior_mask"].=UInt8.(s.interior)
    d["dose_rate"]=copy(s.dose_rate_mean_Gy_s);d["ambient_cum_dose"]=copy(s.accumulated_dose_Gy)
    d["nutrient"]=copy(s.nutrient);d["melanin"]=copy(s.melanin);d["radiation_cpm"]=copy(s.radiation)
    ids,A=contact_graph(r);indices=Dict(id=>i for (i,id) in enumerate(ids))
    for k in eachindex(s.lattice)
        id=Int(s.lattice[k]);id>0||continue
        p=r.parcels[id];d["cell_id"][k]=id;d["species_id"][k]=p.species;d["occupied"][k]=1
        d["parent_id"][k]=p.parent;d["lineage_id"][k]=p.lineage;d["alive"][k]=UInt8(p.status==:alive)
        d["parcel_cum_dose"][k]=p.dose;d["cumulative_hazard"][k]=p.hazard;d["conditional_survival"][k]=exp(-p.hazard)
        d["target_volume"][k]=p.target*prod(r.config.spacing)
        d["growth_rate"][k]=p.status==:alive ? p.growth_rate : 0.
        d["mortality_hazard"][k]=p.status==:alive ? p.hazard_rate : 0.
        d["net_biomass_rate_proxy"][k]=d["growth_rate"][k]-d["mortality_hazard"][k]
        if r.growth
            d["biomass_valid"][k]=1;d["biomass"][k]=p.biomass
            d["biomass_density"][k]=p.biomass/(s.cells[id].volume*prod(r.config.spacing))
        end
        if r.phase!==nothing && p.status==:alive
            i=indices[id];theta=r.phase.theta[id];degree=count(>(0),A[i,:]);strength=sum(A[i,:])
            d["phase_valid"][k]=1;d["phase_rad"][k]=mod(theta,2pi);d["phase_cos"][k]=cos(theta);d["phase_sin"][k]=sin(theta)
            d["natural_frequency"][k]=r.phase.omega[id];d["graph_degree"][k]=degree;d["contact_strength"][k]=strength
            if strength>0
                d["local_coherence_valid"][k]=1
                d["local_coherence"][k]=abs(sum(A[i,j]*cis(r.phase.theta[ids[j]]) for j in eachindex(ids)))/strength
            end
        end
    end
    d
end

function write_image(stem,arrays; origin,spacing)
    shape=size(arrays["cell_id"])
    length(shape)==3 && all(size(x)==shape for x in values(arrays)) || error("voxel array shape mismatch")
    axes=ntuple(i->range(Float64(origin[i]),step=Float64(spacing[i]),length=shape[i]+1),3)
    vtk=vtk_grid(stem,axes...;compress=true)
    for key in sort(collect(keys(arrays)));vtk[key,VTKCellData()]=arrays[key];end
    vtk_save(vtk)
end
function observation(r::Run; phase=phase_observation(r))
    s=r.state;species=NamedTuple[];parcelrows=NamedTuple[]
    for id in sort(collect(keys(r.parcels)))
        p=r.parcels[id];x=get(s.cells,id,nothing)
        radius=x===nothing ? nothing : CPM.radial_dist(x.com[1],x.com[2],r.config.N)/(r.config.N/2)
        push!(parcelrows,(id=id,species=p.species,status=String(p.status),parent=p.parent,lineage=p.lineage,birth=p.birth,death=p.death,removal=p.removal,
            volume=x===nothing ? 0 : x.volume,com=x===nothing ? nothing : copy(x.com),initial_r=p.initial_r,radius=radius,
            biomass=nullable(p.biomass),target=p.target,dose=p.dose,increment=p.increment,hazard=p.hazard,commitment_dose=nullable(p.commitment_dose),
            ancestral_dose=p.ancestral_dose,uptake=p.uptake,growth_rate=p.growth_rate,hazard_rate=p.hazard_rate))
    end
    for sp in 1:CPM.N_SPECIES
        rows=filter(x->x.species==sp && x.status=="alive",parcelrows);n=length(rows)
        push!(species,(species=sp,alive=n,founder_count=r.config.parcels_per_species,
          founder_survivors=count(x->x.species==sp&&x.parent==0&&x.status=="alive",parcelrows),
          mean_radius=n==0 ? nothing : mean(x.radius for x in rows),initial_radius=n==0 ? nothing : mean(x.initial_r for x in rows),
          displacement=n==0 ? nothing : mean(x.radius-x.initial_r for x in rows),
          biomass=r.growth ? sum((x.biomass for x in rows);init=0.) : nothing,
          reason=n==0 ? "zero_viable" : r.growth ? "radius_relative_to_each_parcels_birth" : "survivor_conditioned"))
    end
    (mcs=s.current_mcs,species=species,parcels=parcelrows,phase=phase,budget=budgets(r),nonphase_hash=state_hash(r),full_hash=state_hash(r;phase=true))
end
function export_frame(dir,r::Run; source,run_id)
    validate(r);before=state_hash(r;phase=true);mkpath(dir)
    file="frame_"*lpad(r.state.current_mcs,6,'0')*".vti"
    ispath(joinpath(dir,file))&&error("duplicate export time")
    arrays=projection(r);stem=joinpath(dir,"pending_"*lpad(r.state.current_mcs,6,'0'))
    write_image(stem,arrays;origin=r.config.origin,spacing=r.config.spacing)
    mv(stem*".vti",joinpath(dir,file))
    hashes=Dict(k=>bytes2hex(sha256(reinterpret(UInt8,vec(v)))) for (k,v) in arrays)
    types=Dict(k=>string(eltype(v)) for (k,v) in arrays)
    receipt=(file=file,mcs=r.state.current_mcs,source=source,config=r.config.hash,run_id=run_id,
        schema=EXPORT_SCHEMA,shape=size(r.state.lattice),origin=r.config.origin,spacing=r.config.spacing,
        sha256=bytes2hex(sha256(read(joinpath(dir,file)))),arrays=hashes,dtypes=types,nonphase_hash=state_hash(r))
    before==state_hash(r;phase=true)||error("export mutated simulation")
    atomic_json(joinpath(dir,file*".json"),receipt);receipt
end
function write_collection(dir,frames)
    ts=[f.mcs for f in frames];issorted(ts)&&length(unique(ts))==length(ts)||error("invalid time collection")
    open(joinpath(dir,"series.pvd.partial"),"w") do io
        println(io,"<?xml version=\"1.0\"?><VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\"><Collection>")
        for f in frames
            isfile(joinpath(dir,f.file))||error("missing frame")
            println(io,"<DataSet timestep=\"$(f.mcs)\" group=\"\" part=\"0\" file=\"$(f.file)\"/>")
        end
        println(io,"</Collection></VTKFile>")
    end
    mv(joinpath(dir,"series.pvd.partial"),joinpath(dir,"series.pvd");force=true)
end
