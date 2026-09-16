# Lifecycle, dose-contract, and windowed-API tests for the commit-3 refactor.
# The byte-level serial contract (contract_csv.jl) separately proves the
# legacy path is untouched; these tests exercise the NEW machinery.

using Statistics: mean

# ---------- divide_cell!: successful binary fission ----------

let
    p = SR.CPMParams(N = 24, n_cells_per_species = 2)
    st = SR.init_state(p; seed = 11)
    # largest cell divides most robustly
    id = argmax(Dict(k => c.volume for (k, c) in st.cells))
    parent = st.cells[id]
    pvol, pspecies = parent.volume, parent.species
    plineage, pgen = parent.lineage_id, parent.generation
    next0 = st.next_id

    ev = SR.divide_cell!(st, id; min_daughter_volume = 4)
    @test ev isa SR.DivisionEvent

    # parent identity ended
    @test !haskey(st.cells, id)
    @test st.archive[end].id == id
    @test st.archive[end].fate == :divided
    # two fresh daughters
    da, db = ev.daughter_a_id, ev.daughter_b_id
    @test (da, db) == (next0, next0 + 1)
    @test st.next_id == next0 + 2
    @test haskey(st.cells, da) && haskey(st.cells, db)
    # volume conservation and bookkeeping
    @test ev.parent_volume == pvol
    @test ev.daughter_a_volume + ev.daughter_b_volume == pvol
    @test st.cells[da].volume + st.cells[db].volume == pvol
    for d in (da, db)
        c = st.cells[d]
        @test c.species == pspecies
        @test c.lineage_id == plineage
        @test c.parent_id == id
        @test c.generation == pgen + 1
        @test c.lifetime_dose_Gy == 0.0     # daughters start fresh
    end
    # lattice consistency: no orphan sites, counts match volumes
    counts = Dict{Int, Int}()
    for v in st.lattice
        v > 0 && (counts[Int(v)] = get(counts, Int(v), 0) + 1)
    end
    @test !haskey(counts, id)
    @test counts[da] == st.cells[da].volume
    @test counts[db] == st.cells[db].volume
    @test st.division_log[end] == ev
end

# ---------- divide_cell!: failure paths mutate nothing ----------

let
    p = SR.CPMParams(N = 24, n_cells_per_species = 1)
    st = SR.init_state(p; seed = 11)
    id = first(keys(st.cells))
    lattice0 = copy(st.lattice)
    ncells0 = length(st.cells)

    @test SR.divide_cell!(st, 99999) == :no_such_cell
    @test SR.divide_cell!(st, id; min_daughter_volume = 10^6) == :below_min_volume
    @test st.lattice == lattice0
    @test length(st.cells) == ncells0
    @test isempty(st.division_log)
    @test isempty(st.archive)
end

# ---------- divide_cell!: disconnected daughter is refused ----------

let
    p = SR.CPMParams(N = 32, n_cells_per_species = 1)
    st = SR.init_state(p; seed = 3)
    for i in eachindex(st.lattice)
        st.lattice[i] > 0 && (st.lattice[i] = Int32(0))
    end
    empty!(st.cells)
    # one cell id 1 painted as THREE separated 2x2x2 blobs along x:
    # any principal-axis plane split leaves a daughter in two pieces.
    for xs in ((6, 7), (15, 16), (24, 25)), x in xs, y in 15:16, z in 15:16
        st.lattice[x, y, z] = Int32(1)
    end
    st.cells[1] = SR.CellInfo(1, 24, Float64[15.5, 15.5, 15.5])
    lattice0 = copy(st.lattice)

    @test SR.divide_cell!(st, 1; min_daughter_volume = 4) == :disconnected_daughter
    @test st.lattice == lattice0
    @test haskey(st.cells, 1)
end

# ---------- divide_cell!: deterministic (RNG-free) ----------

let
    split_lattice(seed) = begin
        p = SR.CPMParams(N = 24, n_cells_per_species = 2)
        st = SR.init_state(p; seed = seed)
        id = argmax(Dict(k => c.volume for (k, c) in st.cells))
        ev = SR.divide_cell!(st, id; min_daughter_volume = 4)
        (ev, copy(st.lattice))
    end
    ev1, lat1 = split_lattice(11)
    ev2, lat2 = split_lattice(11)
    @test ev1 == ev2
    @test lat1 == lat2
end

# ---------- dose contract: guarded, attributed per-MCS, signals separate ----------

let
    # unset seconds_per_mcs must refuse to accrue
    st = SR.init_state(SR.CPMParams(N = 16, n_cells_per_species = 1); seed = 5)
    @test isnan(st.params.seconds_per_mcs)
    @test st.melanin_drive == st.radiation      # legacy: identical signals
    dose = fill(0.5, 16, 16, 16)
    SR.import_dose_field!(st, dose, zeros(16, 16, 16))
    @test st.dose_active
    @test_throws ErrorException SR.accrue_dose!(st)

    # configured clock: ΔD = 0.5 Gy/s × 2 s = 1 Gy per MCS
    p = SR.CPMParams(N = 16, n_cells_per_species = 1, seconds_per_mcs = 2.0)
    st = SR.init_state(p; seed = 5)
    rad0 = copy(st.radiation)
    SR.import_dose_field!(st, dose, zeros(16, 16, 16))
    @test st.radiation == rad0                  # no transform ⇒ dynamics untouched
    SR.accrue_dose!(st)
    @test all(st.accumulated_dose_Gy .≈ 1.0)
    @test all(st.dose_increment_Gy .≈ 1.0)
    for (_, c) in st.cells
        @test c.lifetime_dose_Gy ≈ 1.0          # uniform field ⇒ site-mean = 1 Gy
    end

    # explicit transforms are the ONLY path into the biological signals
    SR.import_dose_field!(st, dose, zeros(16, 16, 16);
                          hamiltonian_transform = d -> 2.0 * d,
                          melanin_transform = d -> 0.5 * d)
    @test all(st.radiation .≈ 1.0)
    @test all(st.melanin_drive .≈ 0.25)
end

# ---------- windowed API ≡ legacy coupled driver ----------

let
    p = SR.CPMParams(N = 20, n_cells_per_species = 2, snapshot_interval = 100)
    # basis_gate_ack: this compares blocked output against blocked output (two code
    # paths, or a round trip) and asserts no magnitude, so it needs the gated basis
    # to step but makes no claim about it. Enumerated in the ack census in
    # tests/radiodialysis_basis_gate.jl -- adding a fourth site fails that test.
    rp = SR.RadiolysisParams(Nr = 20, Ddot_R = 1.0, c_ext = 1.0,
                             basis_gate_ack = true)
    n = 21   # crosses the mcs % 10 == 1 biomass-rebuild boundary twice

    st_legacy, rd_legacy, _, _ = redirect_stdout(devnull) do
        SR.run_simulation_coupled(p, rp, n; seed = 7)
    end

    sim = SR.init_coupled_simulation(p, rp; seed = 7)
    SR.advance_window!(sim, 12)     # split across windows on purpose
    SR.advance_window!(sim, n - 12)

    @test sim.state.lattice == st_legacy.lattice
    @test sim.state.melanin == st_legacy.melanin
    @test sim.rd.m == rd_legacy.m
    @test sim.rd.c == rd_legacy.c
    @test sim.rd.s == rd_legacy.s
    @test sim.mcs == n == sim.state.current_mcs
end

# ---------- online feedback: the fail-closed side of the authorization ------
#
# The authorization DECISION lives in Python (biofilm_openmc.feedback_gate).
# Julia's obligation is narrower and more important: it must be structurally
# unable to act on a decision that was never made. These assert the properties
# that make "ONLINE_FEEDBACK = DISABLED" true of the code rather than of a
# configuration flag somebody could flip.
let
    N = 12
    st = SR.init_state(SR.CPMParams(N = N, n_cells_per_species = 1); seed = 7)
    dose = fill(0.25, N, N, N)
    sd = fill(0.01, N, N, N)

    rad0 = copy(st.radiation)
    mel0 = copy(st.melanin_drive)
    lat0 = copy(st.lattice)
    acc0 = copy(st.accumulated_dose_Gy)

    # A physical dose field may ALWAYS be imported — it is a diagnostic until
    # something transforms it. Importing is not authorizing.
    SR.import_dose_field!(st, dose, sd)
    @test st.dose_active
    @test st.dose_rate_mean_Gy_s == dose

    # No transform supplied => no biological signal moved, and no cell moved.
    @test st.radiation == rad0
    @test st.melanin_drive == mel0
    @test st.lattice == lat0
    @test st.accumulated_dose_Gy == acc0

    # And with no clock, no dose may be accrued to any parcel, so a
    # time-integrated response cannot be assembled by accident.
    @test isnan(st.params.seconds_per_mcs)
    @test_throws ErrorException SR.accrue_dose!(st)

    # Re-importing a DIFFERENT field still changes no dynamics. The only path
    # into the biological signals is an explicit transform, which is what the
    # online gate authorizes and what nothing here supplies.
    SR.import_dose_field!(st, fill(9.0, N, N, N), sd)
    @test st.radiation == rad0
    @test st.melanin_drive == mel0
    @test st.lattice == lat0
end
