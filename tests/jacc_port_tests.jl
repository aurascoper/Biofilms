# ---------- JACC port: the tier that was missing --------------------------
#
# `biofilms_potts_jacc.jl` carried a `--selftest` entry point that nothing ever
# called. Neither `runtests.jl` nor `validate_serial.jl` touched it — both load
# the SERIAL monolith through the split-marker trick — so the port could drift
# arbitrarily, or be broken outright, while every suite stayed green. A
# radiodialysis substep guard was added to it before this tier existed and no
# automated path would have caught a mistake in it.
#
# WHY THIS COMPARES KERNELS AND NOT TRAJECTORIES. The port is an 8-colour
# checkerboard Metropolis sweep: sites are visited in a different order from the
# serial reference, so the RNG stream differs and the two CANNOT agree
# trajectory-by-trajectory. `tests/fixtures/serial_seed42.csv` pins the serial
# stream byte-for-byte; demanding the same of a parallel port would be demanding
# that it not be parallel. The declared equivalence criterion is therefore
# per-kernel agreement on identical inputs — ΔH over sampled site pairs, and one
# melanin/nutrient step — at Float32 tolerance. That is a weaker claim than
# bitwise identity and it is the strongest claim a parallel port can support.
#
# WHY IT RUNS ANYWHERE. JACC selects its backend from `LocalPreferences.toml`,
# which is UNTRACKED. With no preference file, JACC defaults to threads, so a
# runner with no GPU executes these kernels for real rather than skipping them.
# On a machine whose preferences name a GPU backend the identical assertions run
# there instead. The backend actually used is reported, because "PASS" means
# something different on each and the log should say which.

const JACC_PORT_PATH = joinpath(REPO, "biofilms_potts_jacc.jl")

@testset "JACC port" begin
    @test isfile(JACC_PORT_PATH)

    # Load the port as a module. Its top-level driver is guarded by
    # `abspath(PROGRAM_FILE) == @__FILE__`, so including it defines the kernels
    # without running `main`.
    JaccPort = Module(:JaccPort)
    Base.eval(JaccPort, :(using Random, Statistics, Printf))
    included = try
        Base.include_string(JaccPort, read(JACC_PORT_PATH, String),
                            "biofilms_potts_jacc.jl")
        true
    catch err
        @error "JACC port failed to load" exception = err
        false
    end
    @test included

    if included
        backend = try
            string(nameof(typeof(Base.eval(JaccPort, :(JACC.default_backend())))))
        catch
            "unknown"
        end
        @info "JACC port selftest running" backend
        @test backend != "unknown"

        # The port's own selftest, which asserts internally with informative
        # messages. Reaching the end without throwing IS the contract; the
        # assertions inside name the site pair or the field that disagreed.
        ran = try
            Base.invokelatest(Base.eval(JaccPort, :selftest), SR)
            true
        catch err
            @error "JACC selftest failed" exception = (err, catch_backtrace())
            false
        end
        @test ran

        # The substep guard must exist in BOTH ports and use the same rule.
        # Without it the JACC port took a single raw forward-Euler step while
        # the serial port substepped, so they would diverge silently the moment
        # the explicit-Euler diffusion bound was exceeded — which is exactly the
        # comparison this repository keeps a serial validator for.
        port_src = read(JACC_PORT_PATH, String)
        serial_src = read(joinpath(REPO, "biofilms_potts.jl"), String)
        rule = "0.4 * dr^2 / (2.0 * rd.params.D_eff)"
        @test occursin(rule, port_src)
        @test occursin(rule, serial_src)
        @test occursin("_step_radiolysis_euler!", port_src)

        # And the guard must be a no-op at the geometry every call site uses,
        # or adding it would have moved a published trajectory.
        for (N, Nr, D_eff, dt_rd) in ((40, 40, 1e-3, 0.5), (20, 40, 1e-3, 0.5))
            R = N / 2.0
            dr = R / (Nr - 1)
            n_sub = max(1, ceil(Int, dt_rd / (0.4 * dr^2 / (2.0 * D_eff))))
            @test n_sub == 1
        end
    end
end
