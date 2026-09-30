"""No-flux finite-volume diffusion plus positive implicit Monod uptake; one transfer amount."""
function nutrient_reaction(C, a, K)
    all(isfinite,(C,a,K)) && C>=0 && a>=0 && K>0 || error("invalid reaction input")
    C==0 && return 0.0
    b=K+a-C; root=sqrt(b*b+4*C*K)
    b>=0 ? 2*C*K/(b+root) : (root-b)/2
end
function diffuse_closed(C,mask,D,dt,spacing)
    size(C)==size(mask) && all(isfinite,C) && all(>=(0),C) || error("invalid diffusion field")
    isfinite(D)&&D>=0&&isfinite(dt)&&dt>0&&all(x->isfinite(x)&&x>0,spacing) || error("invalid diffusion parameters")
    2*D*dt*sum(1 ./ collect(spacing).^2)<=1 || error("diffusion CFL")
    result=copy(C); nx,ny,nz=size(C)
    for z in 1:nz,y in 1:ny,x in 1:nx
        mask[x,y,z] || continue
        for (axis,(xx,yy,zz)) in enumerate(((x+1,y,z),(x,y+1,z),(x,y,z+1)))
            xx<=nx&&yy<=ny&&zz<=nz&&mask[xx,yy,zz] || continue
            flux=D*dt*(C[xx,yy,zz]-C[x,y,z])/spacing[axis]^2
            result[x,y,z]+=flux; result[xx,yy,zz]-=flux
        end
    end
    minimum(result)>=0 && all(isfinite,result) || error("diffusion positivity failure; no clipping")
    result
end
function growth_step!(r::Run)
    s=r.state;c=r.config;dv=prod(c.spacing);dt=c.nutrient_step_mcs
    for q in values(r.parcels);q.uptake=0.;q.growth_rate=0.;end
    for _ in 1:round(Int,1/dt)
        next=diffuse_closed(s.nutrient,s.interior,c.nutrient_diffusion,dt,c.spacing)
        used=Dict(id=>0. for id in alive(r))
        for k in eachindex(s.lattice)
            id=Int(s.lattice[k]);haskey(used,id)||continue
            q=r.parcels[id];sp=q.species
            a=dt*c.uptake_max[sp]*q.biomass/(s.cells[id].volume*dv)
            after=nutrient_reaction(next[k],a,c.monod_K[sp])
            used[id]+=(next[k]-after)*dv;next[k]=after
        end
        for id in sort(collect(keys(used)))
            q=r.parcels[id];sp=q.species
            made=c.yield[sp]*used[id];loss=dt*c.maintenance[sp]*q.biomass
            q.biomass+=made-loss;q.uptake+=used[id]
            q.biomass>0 || error("biomass positivity failure")
            q.target=q.biomass/(c.density[sp]*dv)
            r.budget["uptake"]+=used[id];r.budget["production"]+=made;r.budget["maintenance"]+=loss
        end
        s.nutrient.=next
    end
    # End-state instantaneous conditional growth rate; integrated uptake remains separate.
    for k in eachindex(s.lattice)
        id=Int(s.lattice[k]);id>0 || continue
        q=r.parcels[id];q.status==:alive||continue
        sp=q.species;C=s.nutrient[k]
        q.growth_rate+=c.yield[sp]*c.uptake_max[sp]*C/(c.monod_K[sp]+C)/s.cells[id].volume
    end
    for q in values(r.parcels);q.status==:alive && (q.growth_rate-=c.maintenance[q.species]);end
end
function divide!(r::Run,id::Int)
    validate(r)
    r.growth || return :growth_disabled
    haskey(r.parcels,id)&&r.parcels[id].status==:alive || return :nonviable_parent
    parent=deepcopy(r.parcels[id]);s=r.state;c=r.config
    s.current_mcs-parent.birth>=c.minimum_age_mcs || return :below_minimum_age
    parent.biomass>=c.division_mass[parent.species] || return :below_division_mass
    result=CPM.divide_cell!(s,id;min_daughter_volume=c.minimum_daughter_voxels)
    if result isa Symbol
        event!(r,"division_deferred",id;detail=String(result));return result
    end
    r.parcels[id].status=:divided
    da=result.daughter_a_id;db=result.daughter_b_id
    mass_a=parent.biomass*result.daughter_a_volume/result.parent_volume
    for (daughter,mass) in ((da,mass_a),(db,parent.biomass-mass_a))
        x=s.cells[daughter]
        r.parcels[daughter]=Parcel(x.species,id,x.lineage_id,s.current_mcs,:alive,-1,-1,
          CPM.radial_dist(x.com[1],x.com[2],c.N)/(c.N/2),mass,mass/(c.density[x.species]*prod(c.spacing)),0.,0.,0.,NaN,parent.ancestral_dose+parent.dose,0.,parent.growth_rate,parent.hazard_rate)
        if r.phase!==nothing
            ph=r.phase;ph.theta[daughter]=ph.theta[id];ph.omega[daughter]=ph.omega[id];ph.replay[daughter]=ph.replay[id]
        end
    end
    event!(r,"division",id;detail="$da;$db")
    result
end
function schedule_divisions!(r)
    for id in sort(collect(alive(r)))
        q=r.parcels[id]
        r.state.current_mcs-q.birth>=r.config.minimum_age_mcs && q.biomass>=r.config.division_mass[q.species] || continue
        divide!(r,id)
    end
end
function budgets(r)
    r.growth||return nothing
    nutrient=sum(r.state.nutrient[r.state.interior])*prod(r.config.spacing)
    biomass=sum((q.biomass for q in values(r.parcels) if q.status in (:alive,:dying));init=0.)
    b=r.budget
    (nutrient=nutrient,biomass=biomass,nutrient_residual=nutrient-b["initial_nutrient"]+b["uptake"],
     biomass_residual=biomass-b["initial_biomass"]-b["production"]+b["maintenance"]+b["removed_biomass"],ledger=copy(b))
end
