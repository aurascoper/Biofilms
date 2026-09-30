include("run.jl")
length(ARGS)==2||error("usage: controls.jl CONFIG NEW_OUTPUT")
source,dirty=source_revision();c=F.read_config(ARGS[1]);out=abspath(ARGS[2])
ispath(out)&&error("output exists");mkpath(out)
protocol=TOML.parsefile(joinpath(@__DIR__,"protocol.toml"));c.hash==protocol["config_sha256"]||error("config differs from freeze")
p=protocol["controls"];results=NamedTuple[]
for seed in p["first_seed"]:p["last_seed"]
    cases=[("fixed_baseline",F.new_run(c;seed=seed,arm="baseline",fixed=true)),
           ("fixed_survival",F.new_run(c;seed=seed,arm="survival_only",fixed=true)),
           ("moving_tropism_on",F.new_run(c;seed=seed,arm="survival_only")),
           ("moving_tropism_off",F.new_run(c;seed=seed,arm="survival_only",tropism=false))]
    for (label,r) in cases
        for _ in 1:p["horizon"];F.step!(r);end
        obs=F.observation(r)
        if r.fixed
            all(x->x.displacement===nothing||x.displacement==0.,obs.species)||error("fixed control moved")
        end
        push!(results,(seed=seed,case=label,observation=obs,events=r.events))
    end
end
F.atomic_json(joinpath(out,"controls.json"),(source=source,config_hash=c.hash,plan=p,results=results,
    interpretation="Fixed-position selection and direct-tropism diagnostics; all eight seeds retained."))
println("completed $(length(results)) control runs")
