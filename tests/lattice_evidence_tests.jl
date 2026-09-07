using Test, Random, HDF5, Printf
isdefined(@__MODULE__, :LatticeEvidence) || include(joinpath(REPO, "lattice_evidence.jl"))
const LE = LatticeEvidence
const LD = LE.LabelDynamics

@testset "label histories use identity, order and a fixed mask" begin
    sequence(xs) = [reshape(Int32[x], 1, 1, 1) for x in xs]
    mask = trues(1, 1, 1)
    sp = ones(Int32, 1, 1, 1)
    a, b = sequence([1, 1, 1, 2, 2, 2]), sequence([1, 2, 1, 2, 1, 2])
    @test isequal(LD.species_occupancy(a, mask), LD.species_occupancy(b, mask))
    ma, mb = LD.sampled_history(a, mask, sp), LD.sampled_history(b, mask, sp)
    @test only(ma.transitions) == 1
    @test only(mb.transitions) == 5
    @test only(ma.returns) == 0
    @test only(mb.returns) == 2
    @test !only(ma.persistent) && !only(mb.persistent)
    constant = LD.sampled_history(sequence([1, 1, 1]), mask, sp)
    @test only(constant.transitions) == only(constant.returns) == 0
    @test only(constant.persistent)
    hidden = sequence([1, 2, 1])
    fine = LD.sampled_history(hidden, mask, sp)
    coarse = LD.sampled_history(hidden, mask, sp; stride = 2)
    @test only(fine.transitions) == 2 && only(fine.returns) == 1
    @test only(coarse.transitions) == only(coarse.returns) == 0
    @test only(coarse.hidden) == 1
    @test only(coarse.persistent) && !only(fine.persistent)
    @test only(fine.endpoint_identity)
    # Same-species parcel handoff disappears at species resolution.
    parcel = LD.sampled_history(sequence([11, 12, 11]), mask, sp)
    species = LD.sampled_history(sequence([1, 1, 1]), mask, sp)
    @test only(parcel.transitions) == 2 && only(parcel.returns) == 1
    @test only(species.transitions) == only(species.returns) == 0
    # Relabelling is bijective and changes no identity statistic.
    permuted = LD.sampled_history(sequence([7, 3, 7, 3, 7, 3]), mask, fill(Int32(7), 1, 1, 1))
    for field in (:transitions, :returns, :hidden, :persistent, :endpoint_identity)
        @test getproperty(permuted, field) == getproperty(mb, field)
    end
    p_occ = LD.species_occupancy(sequence([7, 7, 7, 3, 3, 3]), mask)
    @test p_occ[1, 1, 1, 8] == p_occ[1, 1, 1, 4] == 0.5
    # Empty interior space participates; a changing wall does not.
    mask2 = reshape(Bool[true, false], 2, 1, 1)
    labels = [reshape(Int32[0, -1], 2, 1, 1), reshape(Int32[1, 99], 2, 1, 1), reshape(Int32[0, -1], 2, 1, 1)]
    result = LD.sampled_history(labels, mask2, zeros(Int32, 2, 1, 1))
    @test vec(result.transitions) == [2, 0]
    @test vec(result.returns) == [1, 0]
    @test result.rows[1]["denominator_sites"] == 1
    @test result.rows[2]["stratum"] == "initial_species"
    @test result.rows[2]["denominator_sites"] == 1
    @test result.rows[3]["denominator_sites"] == 0
    @test all(isnothing, result.rows[3]["persistence_fraction"])
    occ = LD.species_occupancy(labels, mask2)
    @test occ[1, 1, 1, 1] == 2 / 3
    @test isnan(occ[2, 1, 1, 1])
    @test_throws ErrorException LD.sampled_history(labels, falses(2, 1, 1), zeros(Int32, 2, 1, 1))
end

@testset "canonical export and inert accepted copies" begin
    sim = LE.manuscript_trajectory(SR)
    @test sim.state.params.N == 40
    @test sim.state.params.n_cells_per_species == 6
    @test sim.state.params.T_cpm == 5.0
    @test sim.rd.params.basis_gate_ack
    @test sim.rng isa MersenneTwister
    @test LE.require_initial_inventory(LE.inventory(SR, sim)) === nothing
    bad = deepcopy(LE.inventory(SR, sim))
    bad["live_registry_count"] = 41
    @test_throws ErrorException LE.require_initial_inventory(bad)
    wrong_species = deepcopy(LE.inventory(SR, sim))
    wrong_species["per_species_parcel_counts"] = [5, 7, 6, 6, 6, 6, 6]
    @test_throws ErrorException LE.require_initial_inventory(wrong_species)

    log = LE.AcceptedCopies()
    checked_before = Ref(0)
    observe = e -> begin
        # This closure inspects but never mutates state; the production logger
        # receives only the immutable event. This proves callback timing.
        @test sim.state.lattice[e.donor_site] == e.donor_id
        @test sim.state.lattice[e.recipient_site] == e.recipient_id
        @test e.adh + e.vol + e.rad + e.mel == e.delta_h
        @test e.delta_h <= 0 ? isnan(e.draw) : 0 <= e.draw < exp(-e.delta_h / 5)
        log(e)
        checked_before[] += 1
    end
    bare = LE.manuscript_trajectory(SR)
    perturbed = LE.perturb_uptake!(SR, LE.manuscript_trajectory(SR))
    @test perturbed.rd.params.k_ads == 10bare.rd.params.k_ads
    @test perturbed.rd.params.k_red == 10bare.rd.params.k_red
    mktempdir() do dir
        mkdir(joinpath(dir, "snapshots"))
        mkdir(joinpath(dir, "perturbed"))
        for mcs in 0:2
            if mcs > 0
                SR.advance_window!(sim, 1; on_accepted = observe)
                SR.advance_window!(bare, 1)
                SR.advance_window!(perturbed, 1)
            end
            name = @sprintf("snap_mcs%06d.h5", mcs)
            path = joinpath(dir, "snapshots", name)
            pertpath = joinpath(dir, "perturbed", name)
            LE.write_snapshot(SR, sim, path, "fixture")
            LE.write_snapshot(SR, perturbed, pertpath, "fixture")
            # Through the production writer/reader, not just in-memory arrays.
            actual, pert = LE.read_snapshot(path), LE.read_snapshot(pertpath)
            @test actual.attrs["label_state_hash"] == pert.attrs["label_state_hash"]
            @test actual.info == pert.info
        end
        @test checked_before[] > 0
        @test LE.complete_equal(sim, bare)
        @test rand(copy(sim.rng), UInt64, 20) == rand(copy(bare.rng), UInt64, 20)
        LE.write_events(joinpath(dir, "accepted_copies.h5"), log, "fixture")
        frames = LE.verify_snapshots(dir; expected_mcs = collect(0:2))
        replay = LE.replay_events(dir, frames)
        @test replay["accepted_copy_count"] == checked_before[]
        @test replay["replayed_mcs"] == [0, 1, 2]
        path = joinpath(dir, "snapshots", "snap_mcs000001.h5")
        saved = read(path)
        # Missing MCS: rename outside production's snapshot discovery root.
        mv(path, joinpath(dir, "removed.h5"))
        @test_throws ErrorException LE.verify_snapshots(dir; expected_mcs = collect(0:2))
        mv(joinpath(dir, "removed.h5"), path)
        h5open(path, "r+") do f
            HDF5.delete_attribute(f, "mcs")
            attributes(f)["mcs"] = 0
        end
        @test LE.read_snapshot(path).attrs["mcs"] == 0
        @test_throws ErrorException LE.verify_snapshots(dir; expected_mcs = collect(0:2))
        write(path, saved)
        h5open(path, "r+") do f
            HDF5.delete_attribute(f, "label_state_hash")
            attributes(f)["label_state_hash"] = repeat("0", 64)
        end
        @test read(path) != saved
        @test_throws ErrorException LE.verify_snapshots(dir; expected_mcs = collect(0:2))
        write(path, saved)
        h5open(path, "r+") do f
            HDF5.delete_attribute(f, "live_registry_count")
            attributes(f)["live_registry_count"] = 41
        end
        @test_throws ErrorException LE.verify_snapshots(dir; expected_mcs = collect(0:2))
        write(path, saved)
        @test length(LE.verify_snapshots(dir; expected_mcs = collect(0:2))) == 3
        # Ordered event corruption enters the real replay reader.
        eventpath = joinpath(dir, "accepted_copies.h5")
        eventbytes = read(eventpath)
        h5open(eventpath, "r+") do f
            ds = f["accepted/recipient_id"]
            ds[1] = Int64(999)
        end
        @test read(eventpath) != eventbytes
        @test_throws ErrorException LE.replay_events(dir, frames)
        write(eventpath, eventbytes)
        @test LE.replay_events(dir, frames)["accepted_copy_count"] == checked_before[]
        mkdir(joinpath(dir, "already_exists"))
        @test_throws ErrorException LE.produce_run(SR, "already_exists"; root = dir)
    end
end
