using Serialization
include(joinpath(REPO, "diagnostics", "synthetic_response", "config.jl"))

function dose_test_bytes(value)
    io = IOBuffer()
    serialize(io, value)
    take!(io)
end

function dose_test_state(; dt = 2.0)
    SR.init_state(SR.CPMParams(N = 12, n_cells_per_species = 1,
                              seconds_per_mcs = dt); seed = 9001)
end

function dose_refusal_unchanged(f, state)
    before = dose_test_bytes(state)
    @test_throws Exception f()
    @test dose_test_bytes(state) == before
end

@testset "Dose imports reject before writing" begin
    shape = (12, 12, 12)
    for active in (false, true), field in (:mean, :sd), bad in (
        fill(0.1, 12, 6, 24), fill(0.1, 12, 12, 11), fill(0.1, 12, 12, 13),
        fill(-0.1, shape), fill(NaN, shape), fill(Inf, shape), fill(-Inf, shape))
        st = dose_test_state()
        active && SR.import_dose_field!(st, fill(0.3, shape), fill(0.02, shape))
        mean = field == :mean ? bad : fill(0.5, shape)
        sd = field == :sd ? bad : fill(0.01, shape)
        dose_refusal_unchanged(() -> SR.import_dose_field!(st, mean, sd), st)
    end

    st = dose_test_state()
    rad, mel = copy(st.radiation), copy(st.melanin_drive)
    SR.import_dose_field!(st, zeros(shape), zeros(shape))
    @test st.dose_active && all(iszero, st.dose_rate_mean_Gy_s)
    @test st.radiation == rad && st.melanin_drive == mel
    for transform in (x -> error("transform failed"), x -> Inf, x -> NaN)
        dose_refusal_unchanged(() -> SR.import_dose_field!(st, fill(0.5, shape),
            zeros(shape); hamiltonian_transform = x -> 2x,
            melanin_transform = transform), st)
    end
    SR.import_dose_field!(st, fill(0.5, shape), fill(0.01, shape);
        hamiltonian_transform = x -> 2x, melanin_transform = x -> x / 2)
    @test all(==(1.0), st.radiation)
    @test all(==(0.25), st.melanin_drive)
    # Inputs may alias existing state buffers. Install the original values of both.
    SR.import_dose_field!(st, st.dose_rate_sd_Gy_s, st.dose_rate_mean_Gy_s)
    @test all(==(0.01), st.dose_rate_mean_Gy_s)
    @test all(==(0.5), st.dose_rate_sd_Gy_s)
end

@testset "Dose accrual is guarded and atomic" begin
    shape = (12, 12, 12)
    for dt in (NaN, 0.0, -1.0, Inf, -Inf)
        st = dose_test_state(; dt)
        SR.import_dose_field!(st, fill(0.5, shape), zeros(shape))
        dose_refusal_unchanged(() -> SR.accrue_dose!(st), st)
    end
    st = dose_test_state(dt = NaN)
    before = dose_test_bytes(st)
    @test SR.accrue_dose!(st) === st
    @test dose_test_bytes(st) == before

    for field in (:dose_rate_mean_Gy_s, :dose_rate_sd_Gy_s, :accumulated_dose_Gy),
        bad in (-1.0, NaN, Inf, -Inf)
        st = dose_test_state()
        SR.import_dose_field!(st, fill(0.5, shape), zeros(shape))
        getfield(st, field)[end] = bad
        dose_refusal_unchanged(() -> SR.accrue_dose!(st), st)
    end
    for field in (:dose_rate_mean_Gy_s, :dose_rate_sd_Gy_s,
                  :dose_increment_Gy, :accumulated_dose_Gy)
        st = dose_test_state()
        SR.import_dose_field!(st, fill(0.5, shape), zeros(shape))
        setfield!(st, field, zeros(12, 6, 24))
        dose_refusal_unchanged(() -> SR.accrue_dose!(st), st)
    end
    for bad in (-1.0, NaN, Inf)
        st = dose_test_state()
        SR.import_dose_field!(st, fill(0.5, shape), zeros(shape))
        first(values(st.cells)).lifetime_dose_Gy = bad
        dose_refusal_unchanged(() -> SR.accrue_dose!(st), st)
    end
    # Finite inputs can overflow; no voxel or parcel may receive a partial update.
    for kind in (:increment, :accumulated, :parcel)
        st = dose_test_state()
        rate = kind == :increment ? floatmax(Float64) : floatmax(Float64) / 4
        SR.import_dose_field!(st, fill(rate, shape), zeros(shape))
        kind == :accumulated && fill!(st.accumulated_dose_Gy, floatmax(Float64))
        kind == :parcel && (first(values(st.cells)).lifetime_dose_Gy = floatmax(Float64))
        dose_refusal_unchanged(() -> SR.accrue_dose!(st), st)
    end

    st = dose_test_state()
    rate = [Float64(x) / 100 for x in 1:12, y in 1:12, z in 1:12]
    SR.import_dose_field!(st, rate, zeros(shape))
    expected = Dict(id => sum(rate[st.lattice .== id]) * 2 / c.volume
                    for (id, c) in st.cells)
    SR.accrue_dose!(st)
    @test st.dose_increment_Gy == rate .* 2
    for (id, c) in st.cells
        @test c.lifetime_dose_Gy ≈ expected[id]
    end
    SR.accrue_dose!(st)
    @test st.accumulated_dose_Gy == rate .* 4
    for (id, c) in st.cells
        @test c.lifetime_dose_Gy ≈ 2expected[id]
    end
end

@testset "A window rejects physical inputs before advancing" begin
    for dt in (NaN, 0.0, -1.0, Inf), active in (false, true)
        isnan(dt) && !active && continue # the legacy configuration is valid
        sim = SR.init_coupled_simulation(
            SR.CPMParams(N = 12, n_cells_per_species = 1, seconds_per_mcs = dt),
            SR.RadiolysisParams(Nr = 12); seed = 9001)
        active && SR.import_dose_field!(sim.state, fill(0.1, 12, 12, 12), zeros(12, 12, 12))
        dose_refusal_unchanged(() -> SR.advance_window!(sim, 1), sim)
    end
    for field in (:dose_rate_mean_Gy_s, :dose_rate_sd_Gy_s)
        sim = SR.init_coupled_simulation(
            SR.CPMParams(N = 12, n_cells_per_species = 1, seconds_per_mcs = 2.0),
            SR.RadiolysisParams(Nr = 12); seed = 9001)
        SR.import_dose_field!(sim.state, fill(0.1, 12, 12, 12), zeros(12, 12, 12))
        getfield(sim.state, field)[end] = NaN
        dose_refusal_unchanged(() -> SR.advance_window!(sim, 1), sim)
    end
end

@testset "Explicit, immutable response inputs" begin
    make(alpha; dt = 60.0, axis = 0.01) = SyntheticResponse.ResponseConfig(
        alpha_Gy_inv = alpha, seconds_per_mcs = dt, axis_rate_Gy_s = axis)
    alpha = zeros(7)
    cfg = make(alpha)
    @test cfg.alpha_Gy_inv == Tuple(alpha)
    alpha[1] = -1
    @test cfg.alpha_Gy_inv[1] == 0
    @test make(fill(0.2, 7); axis = 0).axis_rate_Gy_s == 0
    for bad in (zeros(6), zeros(8), zeros(1, 7), [-1.0; zeros(6)],
                [NaN; zeros(6)], [Inf; zeros(6)], ["bad"; zeros(6)])
        @test_throws ArgumentError make(bad)
    end
    for dt in (NaN, 0, -1, Inf, -Inf)
        @test_throws ArgumentError make(zeros(7); dt)
    end
    for axis in (NaN, -1, Inf, -Inf)
        @test_throws ArgumentError make(zeros(7); axis)
    end
    @test_throws MethodError SyntheticResponse.ResponseConfig(Tuple(zeros(7)), 0.0, 0.01)
end
