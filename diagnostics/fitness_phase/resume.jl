include("model.jl")
using .FitnessPhase, JSON3
const F=FitnessPhase
length(ARGS)==4||error("usage: resume.jl RUN_DIRECTORY CHECKPOINT_MCS TERMINAL_MCS RECEIPT.json")
dir=abspath(ARGS[1]);at=parse(Int,ARGS[2]);last=parse(Int,ARGS[3]);output=abspath(ARGS[4])
manifest=JSON3.read(read(joinpath(dir,"run_manifest.json"),String))
source=strip(read(`git -C $(F.ROOT) rev-parse HEAD`,String))
source==manifest.source||error("resume requires the original source revision")
isempty(strip(read(`git -C $(F.ROOT) status --porcelain`,String)))||error("resume requires clean source")
r=F.restore_checkpoint(joinpath(dir,"midpoint.jls");source=source,config_hash=String(manifest.config_hash))
r.state.current_mcs==at && last==manifest.horizon || error("resume interval differs from declared run")
history=JSON3.read(read(joinpath(dir,"metrics.json"),String))
for t in at+1:last
    F.step!(r);F.state_hash(r;phase=true)==history[t+1].full_hash||error("reference continuation mismatch at $t")
end
F.atomic_json(output,(source=source,run_id=manifest.run_id,first=at,last=last,status="passed",exact_complete_state=true))
println(output)
