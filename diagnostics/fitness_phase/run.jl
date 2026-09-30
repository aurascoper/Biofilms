include("model.jl")
using .FitnessPhase, JSON3, SHA, Dates, TOML
const F=FitnessPhase
const ARMS=("baseline","survival_only","phase_only","combined")
function source_revision(;development=false)
    source=strip(read(`git -C $(F.ROOT) rev-parse HEAD`,String))
    dirty=!isempty(strip(read(`git -C $(F.ROOT) status --porcelain`,String)))
    dirty&&!development&&error("production requires a clean committed source")
    source,dirty
end
function run_seed(config,out,seed,horizon;growth=false,dense=false,development=false)
    source,dirty=source_revision(;development=development)
    c=F.read_config(config);horizon>=2||error("horizon too short")
    prediction=joinpath(@__DIR__,"PREDICTION.md")
    !development&&!isfile(prediction)&&error("missing prediction")
    protocol=TOML.parsefile(joinpath(@__DIR__,"protocol.toml"))
    plan=protocol[growth ? "growth" : "replication"]
    if !development
        c.hash==protocol["config_sha256"] || error("configuration differs from frozen plan")
        seed in plan["first_seed"]:plan["last_seed"] && horizon==plan["horizon"] || error("seed or horizon differs from frozen plan")
        dense==(seed in plan["dense_seeds"]) || error("export schedule differs from frozen plan")
    end
    root=joinpath(out,"growth_$(Int(growth))_seed_$seed")
    ispath(root)&&error("refusing to overwrite run $root");mkpath(root)
    runs=[F.new_run(c;seed=seed,arm=arm,growth=growth) for arm in ARMS]
    ids=["$(F.model_id(r))_$(r.arm)_seed$(seed)" for r in runs]
    dirs=[joinpath(root,arm) for arm in ARMS];foreach(mkpath,dirs)
    histories=[NamedTuple[] for _ in ARMS];frames=[NamedTuple[] for _ in ARMS]
    t0=time();mid=fld(horizon,2);restarted=Any[]
    for t in 0:horizon
        for (i,r) in enumerate(runs)
            phase=t==0 ? F.phase_observation(r) : F.step!(r)
            push!(histories[i],F.observation(r;phase=phase))
            if t>mid
                F.step!(restarted[i]);F.state_hash(restarted[i];phase=true)==F.state_hash(r;phase=true)||error("restart divergence at $t $(r.arm)")
            end
            if dense || t in (0,mid,horizon)
                push!(frames[i],F.export_frame(dirs[i],r;source=source,run_id=ids[i]));F.write_collection(dirs[i],frames[i])
            end
            if t==mid
                checkpoint=joinpath(dirs[i],"midpoint.jls")
                F.save_checkpoint(checkpoint,r;source=source)
                push!(restarted,F.restore_checkpoint(checkpoint;config_hash=c.hash,source=source))
            end
        end
        for (a,b) in ((1,3),(2,4))
            histories[a][end].nonphase_hash==histories[b][end].nonphase_hash||error("paired invariant failed at $t $(ARMS[a])/$(ARMS[b])")
        end
    end
    for (i,r) in enumerate(runs)
        F.save_checkpoint(joinpath(dirs[i],"terminal.jls"),r;source=source)
        F.atomic_json(joinpath(dirs[i],"metrics.json"),histories[i])
        F.atomic_json(joinpath(dirs[i],"events.json"),r.events)
        window=r.phase===nothing ? nothing : F.window_metrics(histories[i];first_mcs=mid,last_mcs=horizon,drift_bound=0.01,range_bound=0.5)
        F.atomic_json(joinpath(dirs[i],"phase_window.json"),window)
        F.atomic_json(joinpath(dirs[i],"replay_receipt.json"),(run_id=ids[i]*"_zeroK_replay",parent_run_id=ids[i],source=source,
            config_hash=c.hash,implemented=r.phase!==nothing,method="analytic zero-coupling evolution with identical membership and inherited daughter phase; no forcing or phase noise",
            window=window,phase_feedback=false))
        F.atomic_json(joinpath(dirs[i],"run_manifest.json"),(schema=F.EXPORT_SCHEMA,model=F.model_id(r),run_id=ids[i],source=source,dirty=dirty,
          config_hash=c.hash,config=c.values,seed=seed,arm=r.arm,growth=growth,horizon=horizon,
          time_unit="MCS",exposure_clock_status="synthetic 60 seconds/MCS; not biological calibration",length_unit=c.length_unit,
          association="CellData",tuple_order="x-fastest",units=merge(F.UNITS,Dict("nutrient"=>growth ? "synthetic nutrient amount/synthetic volume" : "legacy nutrient score")),
          frames=frames[i],expected_frame_times=dense ? collect(0:horizon) : [0,mid,horizon],species=F.CPM.SPECIES_NAMES,julia=string(VERSION),threads=Threads.nthreads(),
          prediction_sha256=isfile(prediction) ? bytes2hex(sha256(read(prediction))) : nothing,development=development,
          simulation_complete=true,export_readback="pending",paraview="pending"))
    end
    F.atomic_json(joinpath(root,"invariants.json"),(source=source,config_hash=c.hash,seed=seed,growth=growth,steps=horizon+1,
      pairs=["baseline/phase_only","survival_only/combined"],exact=true,restart_exact=true,restart_at=mid,elapsed_seconds=time()-t0,
      note="Each paired MCS includes canonical nonphase RNG, parcel, dose, field and event state; restart compares complete same-arm state."))
    println("completed seed=$seed growth=$growth horizon=$horizon seconds=$(round(time()-t0;digits=3))")
    root
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)>=6||error("usage: run.jl CONFIG OUTPUT FIRST_SEED LAST_SEED HORIZON GROWTH_0_OR_1 [--development] [--dense]")
    config,out=abspath(ARGS[1]),abspath(ARGS[2]);first,last,horizon=parse.(Int,ARGS[3:5]);growth=parse(Int,ARGS[6])
    growth in (0,1)||error("growth must be 0 or 1")
    all(x->x in ("--development","--dense"),ARGS[7:end])||error("unknown option")
    for seed in first:last
        run_seed(config,out,seed,horizon;growth=growth==1,dense=("--dense" in ARGS)||seed==42,development="--development" in ARGS)
    end
end
