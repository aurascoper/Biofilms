function contact_graph(r::Run)
    s=r.state; ids=sort(collect(alive(r))); index=Dict(id=>i for (i,id) in enumerate(ids))
    A=zeros(length(ids),length(ids)); n=s.params.N
    offsets=r.config.graph_rule=="contact26" ? CPM.NEIGHBORS_26 : CPM.NEIGHBORS_6
    for z in 1:n,y in 1:n,x in 1:n
        i=get(index,Int(s.lattice[x,y,z]),0); i==0 && continue
        for (dx,dy,dz) in offsets
            xx=x+dx;yy=y+dy;zz=z+dz
            1<=xx<=n && 1<=yy<=n && 1<=zz<=n || continue
            j=get(index,Int(s.lattice[xx,yy,zz]),0); j>i || continue
            if r.config.graph_rule=="contact26"
                A[i,j]=A[j,i]=1.
            else
                axis=dx!=0 ? 1 : dy!=0 ? 2 : 3
                area=prod(r.config.spacing)/r.config.spacing[axis]
                A[i,j]+=area; A[j,i]+=area
            end
        end
    end
    ids,A
end

function oscillator_step(theta,omega,A,K,dt)
    n=length(theta)
    size(A)==(n,n) && size(K)==(n,n) && length(omega)==n || error("phase shapes mismatch")
    all(isfinite,theta) && all(isfinite,omega) && all(isfinite,A) && all(isfinite,K) || error("nonfinite phase inputs")
    isfinite(dt)&&dt>0 && all(>=(0),A) && all(>=(0),K) || error("invalid phase parameters")
    A==transpose(A) && all(iszero,diag(A)) || error("invalid undirected graph")
    next=similar(theta)
    for i in 1:n
        d=sum(@view A[i,:]); v=omega[i]
        if d>0
            v+=sum(A[i,j]*K[i,j]*sin(theta[j]-theta[i]) for j in 1:n)/d
        end
        next[i]=theta[i]+dt*v
    end
    all(isfinite,next) || error("nonfinite phase result")
    next
end
coherence(ids,theta) = isempty(ids) ? NaN : abs(sum(cis(theta[id]) for id in sort(collect(ids))))/length(ids)
nullable(x)=isfinite(x) ? x : nothing
function components(A)
    n=size(A,1); membership=zeros(Int,n); sizes=Int[]
    for i in 1:n
        membership[i]==0 || continue
        push!(sizes,0); k=length(sizes); stack=[i]; membership[i]=k
        while !isempty(stack)
            j=pop!(stack); sizes[k]+=1
            for v in 1:n
                if A[j,v]>0 && membership[v]==0
                    membership[v]=k;push!(stack,v)
                end
            end
        end
    end
    membership,sizes
end
function phase_observation(r::Run)
    r.phase===nothing && return nothing
    ids,A=contact_graph(r); n=length(ids); membership,sizes=components(A)
    component_r=[coherence(ids[findall(==(j),membership)],r.phase.theta) for j in eachindex(sizes)]
    edges=[(ids[i],ids[j],A[i,j]) for i in 1:n for j in i+1:n if A[i,j]>0]
    weights=r.growth ? [r.parcels[id].biomass for id in ids] : Float64[]
    rb=isempty(weights) ? nothing : abs(sum(weights[i]*cis(r.phase.theta[id]) for (i,id) in enumerate(ids)))/sum(weights)
    (ids=ids,theta=[r.phase.theta[id] for id in ids],replay=[r.phase.replay[id] for id in ids],
     omega=[r.phase.omega[id] for id in ids],coherence=nullable(coherence(ids,r.phase.theta)),biomass_coherence=rb,
     edges=edges,components=membership,component_sizes=sizes,component_coherence=component_r,
     degree=[count(>(0),A[i,:]) for i in 1:n],strength=vec(sum(A;dims=2)),
     isolates=count(i->all(iszero,A[i,:]),1:n),largest_component=n==0 ? nothing : maximum(sizes)/n,
     reason=n==0 ? "zero_active_phase" : n==1 ? "singleton_no_locking_inference" : "")
end
function phase_step!(r::Run,oldtheta,afterloss,selection)
    r.phase===nothing && return nothing
    ids,A=contact_graph(r); ph=r.phase
    before=coherence(ids,ph.theta)
    birth=isempty(afterloss)||isempty(ids) ? nothing : before-coherence(afterloss,oldtheta)
    theta=[ph.theta[id] for id in ids]; omega=[ph.omega[id] for id in ids]
    K=fill(r.config.phase_K,length(ids),length(ids)); dt=r.config.phase_step_mcs
    for _ in 1:round(Int,1/dt)
        theta=oscillator_step(theta,omega,A,K,dt)
    end
    for (i,id) in enumerate(ids)
        ph.theta[id]=theta[i]; ph.replay[id]+=omega[i] # exact zero-coupling integral for one MCS
    end
    merge(phase_observation(r),(selection=selection,birth=nullable(something(birth,NaN)),evolution=nullable(coherence(ids,ph.theta)-before)))
end

"""Fixed cohort, all ID pairs including changing/nonedges; endpoint drift and path range."""
function window_metrics(records; first_mcs,last_mcs,drift_bound,range_bound)
    rows=filter(x->first_mcs<=x.mcs<=last_mcs,records)
    [x.mcs for x in rows]==collect(first_mcs:last_mcs) || error("missing/duplicate window records")
    all(x->x.phase!==nothing,rows) || error("phase OFF has no phase outcome")
    ids=sort(collect(reduce(intersect,(Set(x.phase.ids) for x in rows))))
    n=length(ids)
    n>=2 || return (count=n,ids=ids,reason=n==0 ? "zero_fixed_cohort" : "singleton_no_locking_inference",locking=nothing)
    function paths(field)
        [Float64[getproperty(row.phase,field)[findfirst(==(id),row.phase.ids)] for id in ids] for row in rows]
    end
    t=paths(:theta); free=paths(:replay); elapsed=last_mcs-first_mcs; elapsed>0 || error("empty analysis window")
    function summary(v)
        drift=Float64[]; ranges=Float64[]; slips=Int[]
        for i in 1:n,j in i+1:n
            delta=[x[i]-x[j] for x in v]
            push!(drift,(delta[end]-delta[1])/elapsed);push!(ranges,maximum(delta)-minimum(delta))
            # Full 2pi windings from the first relative phase, not branch crossings at +/-pi.
            push!(slips,sum(abs(trunc(Int,(delta[k]-delta[1])/(2pi))-trunc(Int,(delta[k-1]-delta[1])/(2pi))) for k in 2:length(delta)))
        end
        freq=(v[end].-v[1])./elapsed
        (drift=drift,ranges=ranges,slips=slips,frequency_spread=std(freq;corrected=false),coherence=[abs(sum(cis,x))/n for x in v])
    end
    coupled=summary(t); replay=summary(free)
    locking=maximum(abs,coupled.drift)<=drift_bound && maximum(coupled.ranges)<=range_bound && all(iszero,coupled.slips)
    reduction=replay.frequency_spread>0 ? 1-coupled.frequency_spread/replay.frequency_spread : nothing
    (count=n,ids=ids,reason="survivor_conditioned",coupled=coupled,replay=replay,locking=locking,
     coupling_supported=locking && reduction!==nothing && reduction>=0.5,frequency_spread_reduction=reduction)
end
