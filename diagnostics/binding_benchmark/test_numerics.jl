#!/usr/bin/env julia
# Data-free numerics for the closed diffusion-and-binding benchmark.
#
#   julia --project=diagnostics/binding_benchmark diagnostics/binding_benchmark/test_numerics.jl
#
# Nothing here reads the frozen geometry, takes an argument, or consults an environment
# variable: every case builds its own lattice. That is deliberate. A suite whose
# assertions are reachable only when a data bundle happens to be mounted is a suite that
# reports "passed" on a machine where it never ran.
#
# Three cases below are *exact*, not approximate, and they are the load-bearing ones:
# a discrete Neumann eigenmode of the flux-form Laplacian, the geometric decay sequence,
# and the capacity-overflow transfer. Where only a rate of convergence is available, the
# observed order is measured and asserted, never assumed.
using Test, Random
include(joinpath(@__DIR__, "BindingBenchmark.jl"))
using .BindingBenchmark

synthetic(occ, int; h = (1.0, 1.0, 1.0)) =
    Geometry(BitArray(occ), BitArray(int), h, 0, "synthetic", "", "")

zero_params(; kw...) = Params(get(kw, :D_c, 0.0), get(kw, :lambda, 0.0), get(kw, :k_on, 0.0),
                              get(kw, :k_off, 0.0), get(kw, :B0, 0.0), get(kw, :dB, 0.0),
                              get(kw, :K, 0.5), get(kw, :n, 2.0), get(kw, :q_ext, 0.0))

order(e_coarse, e_fine) = log2(e_coarse / e_fine)

@testset "binding benchmark numerics" begin

@testset "capacity modulator and Hill response" begin
    @test hill(0.0, 0.5, 2.0) == 0.0
    @test hill(0.5, 0.5, 2.0) ≈ 0.5
    @test hill(1e6, 0.5, 2.0) ≈ 1.0 atol = 1e-9
    @test hill(0.3, 0.5, 2.0) < hill(0.4, 0.5, 2.0) < hill(0.5, 0.5, 2.0)

    # A lone occupied site has no occupied neighbours; each of its six neighbours has one.
    occ = falses(3, 3, 3); occ[2, 2, 2] = true
    A = neighbour_fraction(BitArray(occ))
    @test A[2, 2, 2] == 0.0
    @test A[1, 2, 2] ≈ 1 / 6
    @test A[3, 2, 2] ≈ 1 / 6
    @test A[2, 1, 2] ≈ 1 / 6
    @test A[2, 2, 3] ≈ 1 / 6
    @test A[1, 1, 1] == 0.0          # a face-diagonal neighbour is not a face neighbour
    @test count(!iszero, A) == 6

    # Each member of an occupied pair sees the other, and only the other.
    occ2 = falses(3, 3, 3); occ2[2, 2, 2] = true; occ2[2, 2, 3] = true
    A2 = neighbour_fraction(BitArray(occ2))
    @test A2[2, 2, 2] ≈ 1 / 6
    @test A2[2, 2, 3] ≈ 1 / 6

    # A modulator that does not vary makes the matched-capacity control vacuous, so a
    # configuration that would produce one is refused rather than run.
    solid = trues(3, 3, 3)
    geo = synthetic(solid, solid)
    flat = fill(0.5, 3, 3, 3)
    p = zero_params(B0 = 1.0, dB = 2.0)
    @test_throws ArgumentError capacity_field(geo, p; A = flat)
    @test capacity_field(geo, zero_params(B0 = 1.0, dB = 0.0); A = flat) == fill(1.0, 3, 3, 3)
end

@testset "matched uniform capacity carries the same integral" begin
    # A solid 2x2x2 block is degenerate: every site in it has exactly three occupied face
    # neighbours, so A = 1/2 everywhere and the guard refuses it. One site hung off a face
    # breaks the symmetry -- 1/6, 3/6 and 4/6 all appear.
    cube = falses(4, 4, 4); cube[2:3, 2:3, 2:3] .= true
    @test_throws ArgumentError capacity_field(synthetic(cube, trues(4, 4, 4)),
                                              zero_params(B0 = 1.0, dB = 2.0))
    occ = copy(cube); occ[2, 2, 4] = true
    geo = synthetic(occ, trues(4, 4, 4))
    A = neighbour_fraction(geo.occupied)
    @test sort(unique(round.(A[geo.occupied] .* 6))) == [1.0, 3.0, 4.0]
    B = capacity_field(geo, zero_params(B0 = 1.0, dB = 2.0))
    U, Bbar = matched_uniform_capacity(B, geo.occupied)
    @test sum(U) ≈ sum(B) rtol = 1e-13
    @test count(!iszero, U) == count(geo.occupied)
    @test all(iszero, U[.!geo.occupied])
    @test minimum(B[geo.occupied]) < Bbar < maximum(B[geo.occupied])
end

@testset "flux-form Laplacian conserves on a ragged mask" begin
    int = trues(6, 6, 6)
    int[1, :, :] .= false; int[:, 6, :] .= false; int[3, 3, 3] = false   # ragged interior
    c = randn(Xoshiro(20260911), 6, 6, 6)
    out = zeros(6, 6, 6)
    laplacian!(out, c, BitArray(int), 0.3, (1.0, 1.0, 1.0))
    @test abs(sum(out[int])) < 1e-12 * sum(abs, out[int])
    @test all(iszero, out[.!int])   # nothing is written outside the interior
end

@testset "discrete Neumann eigenmode is reproduced exactly" begin
    # cos(θ_d (x_d − ½)) is an exact eigenvector of this stencil under no flux, with
    # eigenvalue −4D Σ_d sin²(θ_d/2)/h_d². Forward Euler must therefore reproduce the
    # geometric sequence (1 + Δt μ)ⁿ to floating-point precision — this checks the
    # operator, the boundary condition and the update together, not just conservation.
    N = (8, 6, 5); m = (1, 2, 1); h = (1.0, 1.5, 0.75); D = 0.05
    θ = ntuple(d -> π * m[d] / N[d], 3)
    μ = -4D * sum(sin(θ[d] / 2)^2 / h[d]^2 for d in 1:3)
    c0 = [cos(θ[1] * (i - 0.5)) * cos(θ[2] * (j - 0.5)) * cos(θ[3] * (k - 0.5))
          for i in 1:N[1], j in 1:N[2], k in 1:N[3]]
    geo = synthetic(trues(N...), trues(N...); h = h)
    p = zero_params(D_c = D)
    st = State(copy(c0), zeros(N...), zeros(N...), zeros(N...))
    led = Ledger(st, geo)
    dt = 0.9 / diffusion_rate(D, h)
    nsteps = 40
    run_to(st, geo, p, dt, nsteps, led)
    @test maximum(abs, st.c .- (1 + dt * μ)^nsteps .* c0) < 1e-13
    @test abs(1 + dt * μ) < 1                       # the admitted step is a contraction
    @test abs(led.transport_residual) < 1e-12
end

@testset "the timestep bound is the 3-D one, and it is enforced" begin
    # h = D = 1: the one-dimensional bound h²/(2D) admits Δt = 0.25, and at Δt = 0.25 a
    # unit centre with six zero neighbours lands at −0.5. The implemented bound refuses it.
    geo = synthetic(trues(3, 3, 3), trues(3, 3, 3))
    p = zero_params(D_c = 1.0)
    @test diffusion_rate(1.0, (1.0, 1.0, 1.0)) == 6.0

    c = zeros(3, 3, 3); c[2, 2, 2] = 1.0
    lap = zeros(3, 3, 3)
    laplacian!(lap, c, geo.interior, 1.0, (1.0, 1.0, 1.0))
    @test lap[2, 2, 2] ≈ -6.0
    @test c[2, 2, 2] + 0.25 * lap[2, 2, 2] ≈ -0.5     # what the 1-D bound would allow

    st = State(copy(c), zeros(3, 3, 3), zeros(3, 3, 3), zeros(3, 3, 3))
    led = Ledger(st, geo)
    @test_throws ErrorException step!(st, geo, p, 0.25, led)
    @test led.steps == 0                              # the refusal happens before any write
    @test st.c[2, 2, 2] == 1.0

    step!(st, geo, p, 1 / 6, led)                     # equality is admitted
    @test st.c[2, 2, 2] ≈ 0.0 atol = 1e-15
    @test led.steps == 1

    # Reaction and decay are unsplit, so they enter the same restriction. Each of the
    # three terms is pinned separately: a bound that dropped any one of them would still
    # look plausible against a single combined number.
    p2 = Params(1.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.5, 2.0, 0.0)
    st2 = State(copy(c), zeros(3, 3, 3), zeros(3, 3, 3), zeros(3, 3, 3))
    @test stability_rate(st2, geo, p2) ≈ 6.5               # 2DΣh⁻² + λ
    @test_throws ErrorException step!(st2, geo, p2, 1 / 6, Ledger(st2, geo))

    p3 = Params(1.0, 0.5, 0.25, 0.0, 0.0, 0.0, 0.5, 2.0, 0.0)
    st3 = State(zeros(3, 3, 3), zeros(3, 3, 3), fill(4.0, 3, 3, 3), zeros(3, 3, 3))
    @test stability_rate(st3, geo, p3) ≈ 6.0 + 0.25 * 4.0 + 0.5   # + k_on·max(B−b)
    st3.b .= 1.0
    @test stability_rate(st3, geo, p3) ≈ 6.0 + 0.25 * 3.0 + 0.5   # tracks B−b, not B
    st4 = State(fill(10.0, 3, 3, 3), zeros(3, 3, 3), zeros(3, 3, 3), zeros(3, 3, 3))
    p4 = Params(0.0, 0.5, 0.25, 0.75, 0.0, 0.0, 0.5, 2.0, 0.0)
    @test stability_rate(st4, geo, p4) ≈ 0.25 * 10.0 + 0.75 + 0.5  # the bound-pool row
end

@testset "pure decay: the exact sequence, and why Δt ≤ 1/λ is not accuracy" begin
    geo = synthetic(trues(2, 2, 2), trues(2, 2, 2))
    λ = 0.4
    p = zero_params(lambda = λ)
    for dt in (0.1, 0.5, 1 / λ)
        st = State(fill(1.0, 2, 2, 2), zeros(2, 2, 2), zeros(2, 2, 2), zeros(2, 2, 2))
        led = Ledger(st, geo)
        n = 5
        run_to(st, geo, p, dt, n, led)
        @test all(x -> x ≈ (1 - λ * dt)^n, st.c)          # exact discrete solution
        @test led.min_dissolved >= 0                       # positivity holds at the bound
        @test abs(closure_residual(st, geo, led)) < 1e-12
    end
    # At Δt = 1/λ the scheme is positive and wrong: it returns 0 where the continuum
    # returns e⁻¹. Positivity is a constraint on the sign, not on the value.
    st = State(fill(1.0, 2, 2, 2), zeros(2, 2, 2), zeros(2, 2, 2), zeros(2, 2, 2))
    run_to(st, geo, p, 1 / λ, 1, Ledger(st, geo))
    @test st.c[1, 1, 1] == 0.0
    @test abs(st.c[1, 1, 1] - exp(-1)) > 0.36

    errs = Float64[]
    for dt in (0.05, 0.025, 0.0125)
        st = State(fill(1.0, 2, 2, 2), zeros(2, 2, 2), zeros(2, 2, 2), zeros(2, 2, 2))
        n = round(Int, 2.0 / dt)
        run_to(st, geo, p, dt, n, Ledger(st, geo))
        push!(errs, abs(st.c[1, 1, 1] - exp(-λ * 2.0)))
    end
    @test order(errs[1], errs[2]) ≈ 1.0 atol = 0.1
    @test order(errs[2], errs[3]) ≈ 1.0 atol = 0.1
end

@testset "analytical binding limit under a clamped dissolved pool" begin
    # With c held at c0 the bound equation is linear: b(t) = b_eq + (b(0) − b_eq)·e^{−q_bind t},
    # q_bind = k_on·c0 + k_off + λ. This is a different scheme from the closed benchmark —
    # the clamp makes the system open — and the material the clamp returns is booked.
    c0 = 0.8; k_on = 0.3; k_off = 0.12; λ = 0.05; B0 = 2.0; b0 = 0.25
    q_bind = k_on * c0 + k_off + λ
    b_eq = k_on * c0 * B0 / q_bind
    geo = synthetic(trues(3, 3, 3), trues(3, 3, 3))
    p = Params(0.0, λ, k_on, k_off, B0, 0.0, 0.5, 2.0, 0.0)
    B = capacity_field(geo, p)
    @test all(==(B0), B)
    @test b_eq < B0                                  # so the overflow release never fires

    T = 6.0; errs = Float64[]; released = Float64[]
    for dt in (0.2, 0.1, 0.05)
        st = make_state(geo, B; c0 = c0, b0 = b0)
        led = Ledger(st, geo)
        run_to(st, geo, p, dt, round(Int, T / dt), led; clamp_dissolved = true)
        exact = b_eq + (b0 - b_eq) * exp(-q_bind * T)
        push!(errs, abs(st.b[2, 2, 2] - exact))
        push!(released, led.released)
        @test all(x -> x ≈ c0, st.c)                  # the clamp held
        @test abs(closure_residual(st, geo, led)) < 1e-11
    end
    @test all(iszero, released)
    @test order(errs[1], errs[2]) ≈ 1.0 atol = 0.15
    @test order(errs[2], errs[3]) ≈ 1.0 atol = 0.15
end

@testset "capacity overflow is a transfer, and omitting it is visible" begin
    # Every rate is zero, so the Euler update moves nothing and the release is the only
    # thing that happens. The released amount is then known in closed form.
    occ = falses(4, 4, 4); occ[2:3, 2:3, 2:3] .= true
    geo = synthetic(occ, trues(4, 4, 4))
    p = zero_params(B0 = 2.0)
    B = capacity_field(geo, p)
    st = make_state(geo, B; c0 = 0.5, b0 = 2.0)
    @test sum(st.b) ≈ 2.0 * count(occ)

    st.B .*= 0.25                                     # a forced capacity decrease
    led = Ledger(st, geo)
    before = sum(st.c) + sum(st.b)
    step!(st, geo, p, 1.0, led)
    @test led.released ≈ 1.5 * count(occ)
    @test led.dropped == 0.0
    @test sum(st.c) + sum(st.b) ≈ before              # a transfer, not a source
    @test abs(closure_residual(st, geo, led)) < 1e-12
    @test maximum(st.b .- st.B) <= 0
    @test led.max_overshoot ≈ 1.5                     # measured before the release

    # The same run with the transfer omitted. The ledger must not absorb it.
    st2 = make_state(geo, B; c0 = 0.5, b0 = 2.0)
    st2.B .*= 0.25
    led2 = Ledger(st2, geo)
    step!(st2, geo, p, 1.0, led2; release_to_dissolved = false)
    @test led2.released == 0.0
    @test led2.dropped ≈ 1.5 * count(occ)
    @test closure_residual(st2, geo, led2) ≈ -led2.dropped rtol = 1e-12
    @test abs(closure_residual(st2, geo, led2)) > 1.0  # unmistakably red, not a rounding tail
end

@testset "make_state refuses an initial pool above capacity" begin
    occ = falses(3, 3, 3); occ[2, 2, 2] = true
    geo = synthetic(occ, trues(3, 3, 3))
    B = capacity_field(geo, zero_params(B0 = 1.0))
    @test_throws ArgumentError make_state(geo, B; c0 = 0.1, b0 = 1.5)
    st = make_state(geo, B; c0 = 0.1, b0 = 1.0)
    @test st.b[2, 2, 2] == 1.0
    @test st.b[1, 2, 2] == 0.0                        # unoccupied sites carry no bound pool
    @test st.c[1, 2, 2] == 0.1                        # but they do carry dissolved material
end

@testset "both bound-fraction denominators, reported apart" begin
    occ = falses(4, 4, 4); occ[2:3, 2:3, 2:3] .= true
    geo = synthetic(occ, trues(4, 4, 4))
    B = capacity_field(geo, zero_params(B0 = 2.0))
    st = make_state(geo, B; c0 = 1.0, b0 = 0.5)
    of_inv, of_cap = bound_fractions(st, geo)
    @test of_inv ≈ (0.5 * 8) / (1.0 * 64 + 0.5 * 8)
    @test of_cap ≈ (0.5 * 8) / (2.0 * 8)
    @test !(of_inv ≈ of_cap)     # one number cannot stand for both
end

end
