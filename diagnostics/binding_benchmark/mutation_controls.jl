#!/usr/bin/env julia
# Does the data-free suite bite? Seven single-line defects, each planted in a scratch copy
# of the module, each expected to turn `test_numerics.jl` red.
#
#   julia --project=diagnostics/binding_benchmark \
#         diagnostics/binding_benchmark/mutation_controls.jl <out_dir>
#
# An 81-assertion suite that went green the first time it ran is a suite nobody has seen
# fail. This is the receipt that it can.
#
# Every mutation is checked for having applied EXACTLY once before the suite is run. A
# textual patch that quietly matched nothing is indistinguishable, in the report, from a
# defect the suite failed to catch -- and that is not a hypothetical failure mode: a
# colour-table control elsewhere in this repository spent a whole run rewriting a digit
# inside an XML attribute name instead of the value it was aimed at, and reported a pass.
using Dates, JSON3, SHA

const MODULE = "BindingBenchmark.jl"

struct Mutation
    name::String
    defect::String
    from::String
    to::String
end

const MUTATIONS = [
    Mutation("laplacian loses flux form",
             "the neighbour's half of each face flux is not applied, so transport " *
             "creates material at every interior face",
             "            out[i, j, k] += f\n            out[i2, j2, k2] -= f",
             "            out[i, j, k] += f"),
    Mutation("one-dimensional diffusion bound",
             "the stencil's diagonal is halved, admitting the timestep that sends a " *
             "unit centre to -0.5",
             "diffusion_rate(D::Float64, spacing::NTuple{3, Float64}) = 2 * D * sum(1 ./ (spacing .^ 2))",
             "diffusion_rate(D::Float64, spacing::NTuple{3, Float64}) = D * sum(1 ./ (spacing .^ 2))"),
    Mutation("stability ignores the binding term",
             "reaction is treated as if it were split off, so the unsplit scheme is run " *
             "outside its own restriction",
             "    max(diffusion_rate(p.D_c, geo.spacing) + p.k_on * free + p.lambda,",
             "    max(diffusion_rate(p.D_c, geo.spacing) + p.lambda,"),
    Mutation("stability ignores decay",
             "the loss term is dropped from the diagonal",
             "        p.k_on * cmax + p.k_off + p.lambda)",
             "        p.k_on * cmax + p.k_off)"),
    Mutation("capacity overflow is destroyed, not transferred",
             "material above capacity leaves the bound pool and is credited nowhere, " *
             "which is the defect K8 exists to expose",
             "            st.c[ix] += ex\n            rel += ex",
             "            rel += ex"),
    Mutation("the closure identity books the dropped material",
             "the ledger absorbs the loss it exists to report, so the red control turns " *
             "green",
             "               led.external_input + led.chemostat_input -",
             "               led.external_input + led.chemostat_input + led.dropped -"),
    Mutation("the degenerate-modulator guard is removed",
             "a constant modulator is accepted, making the matched-capacity control a " *
             "comparison of a field with itself",
             "        maximum(vals) > minimum(vals) ||\n            throw(ArgumentError(",
             "        false &&\n            throw(ArgumentError("),
]

function apply_mutation(src::String, m::Mutation)
    n = length(collect(eachmatch(Regex(escape_string_for_regex(m.from)), src)))
    n == 1 || return nothing, n
    replace(src, m.from => m.to), n
end

escape_string_for_regex(s) = replace(s, r"([\\^$.|?*+()\[\]{}])" => s"\\\1")

function main(out::String)
    ispath(out) && error("destination already exists: $out")
    here = @__DIR__
    src = read(joinpath(here, MODULE), String)
    proj = here
    results = Vector{Dict{String, Any}}()

    # The unmutated suite first: if this is not green, nothing below means anything.
    base = mktempdir()
    cp(joinpath(here, MODULE), joinpath(base, MODULE))
    cp(joinpath(here, "test_numerics.jl"), joinpath(base, "test_numerics.jl"))
    baseline = success(pipeline(`julia --project=$proj $(joinpath(base, "test_numerics.jl"))`;
                                stdout = devnull, stderr = devnull))

    for m in MUTATIONS
        mutated, hits = apply_mutation(src, m)
        if isnothing(mutated)
            push!(results, Dict{String, Any}(
                "mutation" => m.name, "defect" => m.defect, "verdict" => "DID-NOT-APPLY",
                "match_count" => hits,
                "note" => "the patch matched $hits times, not once; no conclusion about " *
                          "the suite can be drawn from this row"))
            continue
        end
        dir = mktempdir()
        write(joinpath(dir, MODULE), mutated)
        cp(joinpath(here, "test_numerics.jl"), joinpath(dir, "test_numerics.jl"))
        log = joinpath(dir, "out.txt")
        ok = success(pipeline(`julia --project=$proj $(joinpath(dir, "test_numerics.jl"))`;
                              stdout = log, stderr = log))
        failing = [strip(l) for l in eachline(log)
                   if occursin("Test Failed", l) || occursin("Error During Test", l)]
        push!(results, Dict{String, Any}(
            "mutation" => m.name, "defect" => m.defect,
            "verdict" => ok ? "SUITE-STAYED-GREEN" : "CAUGHT",
            "match_count" => hits,
            "failing_assertions" => length(failing),
            "first_failures" => first(failing, 3)))
    end

    mkpath(out)
    doc = Dict{String, Any}(
        "diagnostic" => "mutation controls for the binding-benchmark numerics suite",
        "question" => "can test_numerics.jl fail?",
        "baseline_suite_green" => baseline,
        "module_sha256" => bytes2hex(open(SHA.sha256, joinpath(here, MODULE))),
        "suite_sha256" => bytes2hex(open(SHA.sha256, joinpath(here, "test_numerics.jl"))),
        "mutations" => results,
        "julia_version" => string(VERSION),
        "created_utc" => string(now(UTC)) * "Z")
    open(joinpath(out, "mutation_verification.json"), "w") do io
        JSON3.pretty(io, JSON3.write(doc)); println(io)
    end

    println("baseline suite green: ", baseline)
    for r in results
        println(rpad(r["verdict"], 20), r["mutation"],
                haskey(r, "failing_assertions") ? "  ($(r["failing_assertions"]) assertions)" : "")
    end
    bad = [r["mutation"] for r in results if r["verdict"] != "CAUGHT"]
    println(isempty(bad) && baseline ? "every planted defect was caught" :
            "not caught: " * join(bad, "; "))
    println("wrote ", joinpath(out, "mutation_verification.json"))
    (baseline && isempty(bad)) || exit(1)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("usage: mutation_controls.jl <out_dir>")
    main(ARGS[1])
end
