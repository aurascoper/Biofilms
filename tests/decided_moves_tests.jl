# ---------- Table 4 is re-executable, and this is what executes it ----------
#
# The table was published from a configuration no entry point reproduced:
# `run_simulation` and `run_simulation_coupled` both seed MersenneTwister and
# give 14281 accepted moves at 100 MCS against the table's 16037. The run came
# from the idiom in delta_h_decomposition.jl -- `mcs_step!` by hand with
# Random.Xoshiro(seed) and `update_melanin!` each sweep -- and nothing recorded
# it, so for two weeks the table could not be checked against anything.
# PP-62-11 records that; decided_moves.jl discharges it; this tier is what makes
# the discharge automatic rather than something a person has to remember to run.
#
# ONE SOURCE FOR THE NUMBERS. The expected counts are read out of the script's
# own PUBLISHED table rather than restated here. A second copy would let the
# script and the suite drift and each look green.

include(joinpath(REPO, "decided_moves.jl"))   # defines functions only; its
                                              # PROGRAM_FILE guard does not fire

@testset "Table 4 reproduces from the shipped entry point" begin
    @test isfile(joinpath(REPO, "decided_moves.jl"))
    @test !isempty(PUBLISHED)

    for ((seed, n_mcs), expected) in sort(collect(PUBLISHED))
        got = decided_moves(SR, seed, n_mcs)
        @test got == expected
    end

    # THE TWO NUMBERS THE VERSION 1.2 CORRECTION MOVED, pinned directly, because
    # they are the ones a regression would quietly restore. `rad` is the fourth
    # label and 206042 is the denominator PP-62-04 published.
    rad_idx = findfirst(==(:rad), SR.DRIVER_LABELS)
    @test rad_idx == 4
    at400 = [PUBLISHED[(s, 400)] for s in (42, 43, 44)]
    @test sum(sum(c) for c in at400) == 206_042
    @test sum(c[rad_idx] for c in at400) == 1      # NOT zero -- seed 44 carries it

    # ...and the control: the comparison must be able to fail. MersenneTwister
    # is what every shipped run_simulation* path uses, and it is precisely what
    # does NOT reproduce the table.
    p  = SR.CPMParams(N = 40, n_cells_per_species = 6)
    st = SR.init_state(p; seed = 42)
    d  = SR.DriverCounts(p.N)
    rng = Random.MersenneTwister(42)
    for _ in 1:100
        SR.mcs_step!(st, rng; driver = d)
        SR.update_melanin!(st)
    end
    @test sum(vec(sum(d.counts, dims = (1, 2, 3)))) != sum(PUBLISHED[(42, 100)])
end
