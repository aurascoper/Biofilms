# Exchange-schema tests: transport snapshot integrity (probes, conventions)
# and restart checkpoint bit-exact resume (RNG included).

using HDF5

include(joinpath(REPO, "export_checkpoint.jl"))   # CLI-guarded; functions only

tmp = mktempdir()
p = SR.CPMParams(N = 20, n_cells_per_species = 2, snapshot_interval = 100)
rp = SR.RadiolysisParams(Nr = 20, Ddot_R = 1.0, c_ext = 1.0)

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
