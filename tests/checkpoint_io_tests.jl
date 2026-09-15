# Exchange-schema tests: transport snapshot integrity (probes, conventions)
# and restart checkpoint bit-exact resume (RNG included).

using HDF5

include(joinpath(REPO, "export_checkpoint.jl"))   # CLI-guarded; functions only

tmp = mktempdir()
p = SR.CPMParams(N = 20, n_cells_per_species = 2, snapshot_interval = 100)
# basis_gate_ack: this compares blocked output against blocked output (two code
# paths, or a round trip) and asserts no magnitude, so it needs the gated basis
# to step but makes no claim about it. Enumerated in the ack census in
# tests/radiodialysis_basis_gate.jl -- adding a fourth site fails that test.
rp = SR.RadiolysisParams(Nr = 20, Ddot_R = 1.0, c_ext = 1.0,
                         basis_gate_ack = true)

# ---------- a checkpoint with no provenance is UNKNOWN, not ungated ----------

let
    # Codex on f49fe94: defaulting missing provenance to false is rule 3 inside
    # a gate. Simulate a pre-provenance file by deleting the dataset.
    rpx = SR.RadiolysisParams(Nr = 20, Ddot_R = 1.0, c_ext = 1.0,
                              basis_gate_ack = true)
    sim = SR.init_coupled_simulation(p, rpx; seed = 11)
    SR.advance_window!(sim, 2)
    legacy = joinpath(tmp, "legacy_no_provenance.h5")
    export_restart_checkpoint(SR, sim, legacy)
    h5open(legacy, "r+") do f
        delete_object(f, "rd/basis_from_occupancy")
    end
    @test !h5open(g -> haskey(g, "rd/basis_from_occupancy"), legacy, "r")

    # refuses rather than assuming
    err = try
        restore_restart_checkpoint(SR, legacy); nothing
    catch e
        sprint(showerror, e)
    end
    @test err !== nothing
    @test occursin("UNKNOWN", err)

    # and both explicit declarations are honoured
    simT = restore_restart_checkpoint(SR, legacy; declare_basis_from_occupancy = true)
    @test simT.rd.basis_from_occupancy === true
    simF = restore_restart_checkpoint(SR, legacy; declare_basis_from_occupancy = false)
    @test simF.rd.basis_from_occupancy === false

    # a CURRENT checkpoint still restores with no declaration needed, or the
    # refusal above would just be blocking everything.
    current = joinpath(tmp, "current_provenance.h5")
    export_restart_checkpoint(SR, sim, current)
    @test restore_restart_checkpoint(SR, current).rd.basis_from_occupancy === true
end

# ---------- the exported file declares a gated basis (rule 4) ----------

let
    # A consumer reading rd/c cannot tell the basis was gated by looking at the
    # array, so the producer says it. Both directions, because a marker that is
    # always true and one that is always false are equally uninformative.
    # The `false` case must NOT advance: advance_window! reconstructs the params
    # with X_total = mean(compute_radial_biomass(...)), so every coupled run
    # ends up on the gated basis by construction. That is the finding, not a
    # quirk -- an advanced coupled sim can never export an ungated basis.
    for (mcs, expected) in ((0, false), (2, true))
        rpx = SR.RadiolysisParams(Nr = 20, Ddot_R = 1.0, c_ext = 1.0,
                                  basis_gate_ack = true)
        sim = SR.init_coupled_simulation(p, rpx; seed = 5)
        mcs > 0 && SR.advance_window!(sim, mcs)
        @test (sim.rd.params.X_total != 1.0) == expected
        # export_restart_checkpoint, not the transport snapshot: only the
        # restart file carries rd/c and rd/s. The transport snapshot exports
        # lattice, fields and dose, none of which touch the gated basis.
        path = joinpath(tmp, "gate_$(mcs).h5")
        export_restart_checkpoint(SR, sim, path)
        h5open(path, "r") do f
            @test read(f["rd/basis_gate_blocked"]) == expected
            note = read(f["rd/basis_gate_note"])
            @test occursin(expected ? "RADIODIALYSIS: BLOCKED" :
                                      "the basis gate does not apply", note)
        end
    end
end

# ---------- transport snapshot: conventions + probe integrity ----------

let
    sim = SR.init_coupled_simulation(p, rp; seed = 9)
    SR.advance_window!(sim, 5)
    snap = joinpath(tmp, "snap.h5")
    export_transport_snapshot(SR, sim, snap)

    h5open(snap, "r") do f
        a = attributes(f)
        @test read(a["schema_version"]) == 1
        @test read(a["logical_axis_order"]) == "xyz"
        @test read(a["dataset_axis_order_h5py"]) == "zyx"
        @test read(a["coordinate_index_base"]) == 0
        @test read(a["cell_id_wall"]) == -1
        @test read(a["material_class_source"]) == "absent"
        @test length(read(a["label_state_hash"])) == 64

        cell_id = read(f["lattice/cell_id"])
        @test eltype(cell_id) == Int32
        @test size(cell_id) == (20, 20, 20)
        @test read(f["config_toml"]) == ""

        # every probe row must match the array at its 0-based coordinates
        probes = read(f["orientation_probes"])
        @test size(probes, 1) >= 3
        for i in axes(probes, 1)
            x, y, z, expect = probes[i, :]
            @test Int64(cell_id[x + 1, y + 1, z + 1]) == expect
        end

        # label arrays are consistent with the live-cell registry
        species = read(f["lattice/species_id"])
        ids = read(f["cells/id"])
        @test all(>(0), ids)
        @test sum(species .> 0) == sum(read(f["cells/volume"]))
    end
end

# ---------- restart checkpoint: resumed ≡ unbroken (bit-exact) ----------

let
    simA = SR.init_coupled_simulation(p, rp; seed = 13)
    SR.advance_window!(simA, 15)
    ckpt = joinpath(tmp, "restart.h5")
    export_restart_checkpoint(SR, simA, ckpt)

    SR.advance_window!(simA, 10)                    # unbroken run continues

    simB = restore_restart_checkpoint(SR, ckpt)     # resumed run
    @test simB.mcs == 15
    SR.advance_window!(simB, 10)

    @test simB.state.lattice == simA.state.lattice  # RNG-exact continuation
    @test simB.state.melanin == simA.state.melanin
    @test simB.state.nutrient == simA.state.nutrient
    @test simB.state.current_mcs == simA.state.current_mcs
    @test simB.rd.c == simA.rd.c
    @test simB.rd.s == simA.rd.s
    @test simB.rd.m == simA.rd.m
    @test simB.mcs == simA.mcs
    @test length(simB.state.archive) == length(simA.state.archive)

    # live-cell registries agree exactly
    @test sort(collect(keys(simB.state.cells))) == sort(collect(keys(simA.state.cells)))
    for (id, cA) in simA.state.cells
        cB = simB.state.cells[id]
        @test (cB.species, cB.volume, cB.lineage_id, cB.generation) ==
              (cA.species, cA.volume, cA.lineage_id, cA.generation)
    end
end
