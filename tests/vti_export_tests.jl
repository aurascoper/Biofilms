# .vti exporter: the file says what the snapshot said (round trip), VTK x is Julia
# axis 1 (orientation, checked against the snapshot's own probes), sentinels survive,
# and a physical spacing is refused without a declaration. Everything here is in
# lattice units; the one physical array is carried under the schema's own label.

using HDF5, WriteVTK, ReadVTK

isdefined(Main, :export_transport_snapshot) || include(joinpath(REPO, "export_checkpoint.jl"))
include(joinpath(REPO, "export_vti.jl"))

# String field data is written appended and raw (compress = false): a UInt64 byte
# count then the bytes, at the array's offset after the "_" that opens the block.
# ReadVTK reads numeric arrays only, so the strings are read back by hand here.
function vti_field_string(path, name)
    bytes = read(path)
    txt = String(copy(bytes))
    m = match(Regex("<Array type=\"String\" Name=\"$name\"[^>]*offset=\"(\\d+)\""), txt)
    m === nothing && return nothing
    off = parse(Int, m.captures[1])
    i = findfirst("<AppendedData encoding=\"raw\">", txt)
    start = findnext(==(UInt8('_')), bytes, i[end]) + 1 + off
    n = Int(reinterpret(UInt64, bytes[start:start + 7])[1])
    return rstrip(String(bytes[start + 8:start + 7 + n]), '\0')
end

tmp = mktempdir()
p = SR.CPMParams(N = 12, n_cells_per_species = 2, snapshot_interval = 100)
# basis_gate_ack: the snapshot fixture only; nothing here reads rd/c or asserts a
# magnitude. Enumerated in the ack census in tests/radiodialysis_basis_gate.jl.
rp = SR.RadiolysisParams(Nr = 20, Ddot_R = 1.0, c_ext = 1.0, basis_gate_ack = true)
sim = SR.init_coupled_simulation(p, rp; seed = 9)
SR.advance_window!(sim, 2)
snap = joinpath(tmp, "snap.h5")
export_transport_snapshot(SR, sim, snap)
N = (12, 12, 12)

@testset "round trip: every site array reads back equal, in lattice units" begin
    stem = joinpath(tmp, "rt")
    export_vti(snap, stem)
    path = stem * ".vti"
    f = VTKFile(path)
    @test collect(get_spacing(f)) == [1.0, 1.0, 1.0]
    @test collect(get_origin(f)) == [0.0, 0.0, 0.0]
    # The exporter's default is the uncompressed appended form, and this file was written
    # with that default: the path the tests read is the path a user gets.
    @test !occursin("compressor=", read(path, String)) && occursin("encoding=\"raw\"", read(path, String))
    @test vti_field_string(path, "units") == "lattice"
    @test occursin("zeros until a dose was imported", vti_field_string(path, "accumulated_dose_Gy_units"))
    @test vti_field_string(path, "logical_axis_order") == "xyz"
    fd = get_field_data(f)
    @test get_data(fd["mcs"])[1] == 2.0
    @test get_data(fd["cell_id_wall"])[1] == -1.0
    # The snapshot above was written with no seed, so the .vti states the absent
    # declaration. NaN pairs with it, as physical_time_s already does here.
    @test vti_field_string(path, "cpm_seed_source") == "absent"
    @test isnan(get_data(fd["cpm_seed"])[1])

    @testset "a snapshot predating the seed attributes still exports, and says so" begin
        # Every file written before 2026-09-17 has no cpm_seed_source at all.
        # Requiring it would make this exporter throw on all of them, and
        # skipping it would leave the .vti silent about which case it is.
        old = joinpath(tmp, "old.h5")
        cp(snap, old)
        h5open(old, "r+") do g
            delete_attribute(g, "cpm_seed_source")
        end
        stem_old = joinpath(tmp, "old_rt")
        export_vti(old, stem_old)
        pold = stem_old * ".vti"
        @test vti_field_string(pold, "cpm_seed_source") == "unrecorded"
        @test isnan(get_data(get_field_data(VTKFile(pold))["cpm_seed"])[1])
    end

    @testset "a declared seed reaches the .vti as a number" begin
        # The other half of the control: absent, unrecorded and declared must
        # produce three different field-data values, not one.
        seeded = joinpath(tmp, "seeded.h5")
        export_transport_snapshot(SR, sim, seeded; cpm_seed = 7)
        stem_seed = joinpath(tmp, "seeded_rt")
        export_vti(seeded, stem_seed)
        pseed = stem_seed * ".vti"
        @test vti_field_string(pseed, "cpm_seed_source") == "declared"
        @test get_data(get_field_data(VTKFile(pseed))["cpm_seed"])[1] == 7.0
    end
    h5open(snap, "r") do h
        cd = get_cell_data(f)
        @test all(iszero, read(h["lattice/generation"]))
        @test all(iszero, read(h["dose/accumulated_Gy"]))
        @test !("generation" in keys(cd))
        @test !("accumulated_dose_Gy" in keys(cd))
        for (ds, name) in SITE_ARRAYS
            name in ("generation", "accumulated_dose_Gy") && continue
            src = read(h[ds])
            back = reshape(get_data(cd[name]), N)
            @test eltype(back) == eltype(src)
            @test back == src
        end
        species = reshape(get_data(cd["species"]), N)
        @test eltype(species) == UInt8
        cid = read(h["lattice/cell_id"])
        @test all(species[cid .<= 0] .== 0)          # medium and wall are air
        @test species[cid .> 0] == UInt8.(read(h["lattice/species_id"])[cid .> 0])
        @test Set(unique(species)) ⊇ Set(0x01:0x07)  # the fixture places all seven
    end

    @testset "non-zero generation and accumulated_dose_Gy are written and round-trip (issue #34)" begin
        # Negative control for OMIT_WHEN_ZERO: a blanket removal of generation and
        # accumulated_dose_Gy from SITE_ARRAYS would pass the all-zero check above and
        # silently drop real division or imported-dose data when present.
        nz_snap = joinpath(tmp, "snap_nonzero.h5")
        cp(snap, nz_snap; force = true)
        gen_expected = zeros(Int32, N)
        gen_expected[2, 3, 4] = Int32(2)
        dose_expected = zeros(Float64, N)
        dose_expected[5, 6, 7] = 1.25
        h5open(nz_snap, "r+") do g
            g["lattice/generation"][:, :, :] = gen_expected
            g["dose/accumulated_Gy"][:, :, :] = dose_expected
        end
        nz_stem = joinpath(tmp, "rt_nonzero")
        export_vti(nz_snap, nz_stem)
        cd_nz = get_cell_data(VTKFile(nz_stem * ".vti"))
        @test "generation" in keys(cd_nz)
        @test "accumulated_dose_Gy" in keys(cd_nz)
        gen_back = reshape(get_data(cd_nz["generation"]), N)
        dose_back = reshape(get_data(cd_nz["accumulated_dose_Gy"]), N)
        @test eltype(gen_back) == Int32 && gen_back == gen_expected
        @test eltype(dose_back) == Float64 && dose_back == dose_expected

        # Per-array independence: non-zero generation with all-zero accumulated_dose_Gy
        # writes generation and still omits accumulated_dose_Gy.
        half_snap = joinpath(tmp, "snap_half_nonzero.h5")
        cp(snap, half_snap; force = true)
        h5open(half_snap, "r+") do g
            g["lattice/generation"][:, :, :] = gen_expected
        end
        half_stem = joinpath(tmp, "rt_half_nonzero")
        export_vti(half_snap, half_stem)
        cd_half = get_cell_data(VTKFile(half_stem * ".vti"))
        @test "generation" in keys(cd_half)
        @test !("accumulated_dose_Gy" in keys(cd_half))

        # Symmetric independence: all-zero generation with non-zero accumulated_dose_Gy
        # writes accumulated_dose_Gy and still omits generation.
        dose_only_snap = joinpath(tmp, "snap_dose_only_nonzero.h5")
        cp(snap, dose_only_snap; force = true)
        h5open(dose_only_snap, "r+") do g
            g["dose/accumulated_Gy"][:, :, :] = dose_expected
        end
        dose_only_stem = joinpath(tmp, "rt_dose_only_nonzero")
        export_vti(dose_only_snap, dose_only_stem)
        cd_dose_only = get_cell_data(VTKFile(dose_only_stem * ".vti"))
        @test !("generation" in keys(cd_dose_only))
        @test "accumulated_dose_Gy" in keys(cd_dose_only)
    end
end

@testset "sentinels: background 0 and wall -1 survive as stored" begin
    path = joinpath(tmp, "rt.vti")
    cid = reshape(get_data(get_cell_data(VTKFile(path))["cell_id"]), N)
    @test eltype(cid) == Int32
    @test count(==(Int32(-1)), cid) > 0
    @test count(==(Int32(0)), cid) > 0
    @test minimum(cid) == -1
end

@testset "orientation: VTK x is Julia axis 1, by a valued array and by the snapshot's probes" begin
    # A 3x4x5 array valued 100i+10j+k: any axis swap or flip moves a value.
    A = Int32[100i + 10j + k for i in 1:3, j in 1:4, k in 1:5]
    stem = joinpath(tmp, "orient")
    vtk = vtk_grid(stem, 0:3, 0:4, 0:5; compress = false)
    vtk["v", VTKCellData()] = A
    vtk_save(vtk)
    back = reshape(get_data(get_cell_data(VTKFile(stem * ".vti"))["v"]), 3, 4, 5)
    @test back == A
    @test back[2, 3, 4] == 234
    # The exporter's own output against the snapshot's orientation probes
    # (0-based x, y, z, value rows written by export_checkpoint.jl).
    cid = reshape(get_data(get_cell_data(VTKFile(joinpath(tmp, "rt.vti")))["cell_id"]), N)
    probes = h5open(h -> read(h["orientation_probes"]), snap, "r")
    @test size(probes, 1) >= 3
    for r in 1:size(probes, 1)
        x, y, z, v = probes[r, :]
        @test cid[x + 1, y + 1, z + 1] == v
    end
end

@testset "a physical spacing is refused without a declaration" begin
    @test_throws ArgumentError export_vti(snap, joinpath(tmp, "phys"); spacing = 0.012)
    @test !isfile(joinpath(tmp, "phys.vti"))
    stem = joinpath(tmp, "declared")
    export_vti(snap, stem; spacing = 0.012,
               declared_pitch = (value = 0.012, unit = "cm", source = "test value, not a measurement"))
    f = VTKFile(stem * ".vti")
    @test collect(get_spacing(f)) == [0.012, 0.012, 0.012]
    u = vti_field_string(stem * ".vti", "units")
    @test startswith(u, "declared:") && occursin("test value", u)
    # and a snapshot whose axis order is not the one this exporter maps is refused
    odd = joinpath(tmp, "odd.h5"); cp(snap, odd)
    h5open(odd, "r+") do h; delete_attribute(h, "logical_axis_order"); attributes(h)["logical_axis_order"] = "zyx"; end
    @test_throws ArgumentError export_vti(odd, joinpath(tmp, "odd"))
end

@testset "a Gy/s field is labelled from the transport result's own qualifiers, or refused" begin
    # A transport result the shape results.py writes: mesh/dose_rate_mean_Gy_s plus attrs.
    function fake_result(path; attrs...)
        h5open(path, "w") do g
            g["mesh/dose_rate_mean_Gy_s"] = fill(0.25, N)
            for (k, v) in attrs
                attributes(g)[string(k)] = v
            end
        end
        path
    end
    synthetic = fake_result(joinpath(tmp, "tr_synth.h5"); target_calibration = 0,
                            source_rate_photons_per_s = 3.7e9, logical_axis_order = "xyz")
    stem = joinpath(tmp, "with_dose")
    export_vti(snap, stem; dose = synthetic)
    f = VTKFile(stem * ".vti")
    @test reshape(get_data(get_cell_data(f)["dose_rate_mean_Gy_s"]), N) == fill(0.25, N)
    u = vti_field_string(stem * ".vti", "dose_rate_mean_Gy_s_units")
    @test occursin("synthetic source rate, not a physical target", u)
    @test occursin("3.7e9", u) && occursin("mesh/dose_rate_mean_Gy_s", u)
    target = fake_result(joinpath(tmp, "tr_target.h5"); target_calibration = 1,
                         source_rate_photons_per_s = 1.0e8, logical_axis_order = "xyz")
    export_vti(snap, joinpath(tmp, "with_target"); dose = target)
    @test occursin("target_calibration = true", vti_field_string(joinpath(tmp, "with_target.vti"), "dose_rate_mean_Gy_s_units"))
    bare = fake_result(joinpath(tmp, "tr_bare.h5"); source_rate_photons_per_s = 1.0e8, logical_axis_order = "xyz")
    @test_throws ArgumentError export_vti(snap, joinpath(tmp, "bare"); dose = bare)
    @test !isfile(joinpath(tmp, "bare.vti"))
    wrong = joinpath(tmp, "tr_wrong.h5")
    h5open(wrong, "w") do g
        g["mesh/dose_rate_mean_Gy_s"] = zeros(4, 4, 4)
        attributes(g)["target_calibration"] = 0; attributes(g)["source_rate_photons_per_s"] = 1.0
        attributes(g)["logical_axis_order"] = "xyz"
    end
    @test_throws ArgumentError export_vti(snap, joinpath(tmp, "wrong"); dose = wrong)
end

@testset "series: one .vti per snapshot, .pvd keyed by mcs, duplicates refused" begin
    d = joinpath(tmp, "series"); mkpath(d)
    cp(snap, joinpath(d, "a.h5"))
    sim2 = SR.init_coupled_simulation(p, rp; seed = 9)
    SR.advance_window!(sim2, 5)
    export_transport_snapshot(SR, sim2, joinpath(d, "b.h5"))
    written = export_series(d, joinpath(tmp, "run"))
    @test length(written) == 2 && all(isfile, written)
    pvd = read(joinpath(tmp, "run.pvd"), String)
    @test occursin("timestep=\"2.0\"", pvd) && occursin("timestep=\"5.0\"", pvd)
    cp(snap, joinpath(d, "c.h5"))   # a second snapshot at mcs 2
    @test_throws ArgumentError export_series(d, joinpath(tmp, "dup"))
end

@testset "a restart from another time is refused; the same time is accepted" begin
    same = joinpath(tmp, "restart_same.h5")
    export_restart_checkpoint(SR, sim, same)          # sim is at mcs 2, like snap
    export_vti(snap, joinpath(tmp, "with_nutrient"); restart = same)
    f = VTKFile(joinpath(tmp, "with_nutrient.vti"))
    @test "nutrient" in keys(get_cell_data(f))
    later = SR.init_coupled_simulation(p, rp; seed = 9)
    SR.advance_window!(later, 5)
    other = joinpath(tmp, "restart_other.h5")
    export_restart_checkpoint(SR, later, other)
    @test_throws ArgumentError export_vti(snap, joinpath(tmp, "wrong_time"); restart = other)
    @test !isfile(joinpath(tmp, "wrong_time.vti"))
end

@testset "a dose file declaring another axis order is refused" begin
    zyx = joinpath(tmp, "tr_zyx.h5")
    h5open(zyx, "w") do g
        g["mesh/dose_rate_mean_Gy_s"] = fill(0.25, N)
        attributes(g)["target_calibration"] = 0
        attributes(g)["source_rate_photons_per_s"] = 1.0
        attributes(g)["logical_axis_order"] = "zyx"
    end
    @test_throws ArgumentError export_vti(snap, joinpath(tmp, "zyx"); dose = zyx)
    @test !isfile(joinpath(tmp, "zyx.vti"))
end

# ---------- the viewer's renderer-free half, in CI ----------
include(joinpath(REPO, "viewer", "lattice_grid.jl"))

@testset "viewer grid: air exactly where the file's sentinels say, ids elsewhere" begin
    grid, mcs = species_grid(snap)
    @test mcs == 2 && size(grid) == N && eltype(grid) == UInt8
    h5open(snap, "r") do h
        cid = read(h["lattice/cell_id"]); sp = read(h["lattice/species_id"])
        @test all(grid[cid .<= 0] .== 0x00)
        @test grid[cid .> 0] == UInt8.(sp[cid .> 0])
    end
    @test length(COLORS) == 7 && length(LABELS) == 7
    @test_throws ArgumentError species_grid(joinpath(tmp, "odd.h5"))   # logical_axis_order zyx
end

@testset "viewer CLI contract" begin
    o = viewer_options(["s.h5"])
    @test o.snapshot == "s.h5" && o.still === nothing && o.record_to === nothing && o.frames == 120
    o = viewer_options(["s.h5", "--still", "a.png", "--frames", "12"])
    @test o.still == "a.png" && o.frames == 12
    @test_throws ArgumentError viewer_options(String[])
    @test_throws ArgumentError viewer_options(["s.h5", "--bogus", "x"])
    @test_throws ArgumentError viewer_options(["s.h5", "--still"])
    @test_throws ArgumentError viewer_options(["s.h5", "--still", "a.png", "--record", "b.mp4"])
    @test_throws ArgumentError viewer_options(["s.h5", "--frames", "0"])
end

# The GL half cannot run here. What can be decided statically is decided statically: the
# file parses, and every name it uses resolves in Base, HDF5, lattice_grid.jl, a local
# binding, or the declared list of GLMakie names it relies on. A first-render UndefVarError
# is the class this catches (AGENTS.md rule 2).
GLMAKIE_NAMES = Set([:GLMakie, :Makie, :Figure, :Axis3, :voxels!, :Legend, :PolyElement,
                     :save, :record, :display, :Colorant, :Screen, :wait, :..])
function free_names(ex)
    defined = Set{Symbol}([:ARGS, :PROGRAM_FILE, :__DIR__, :__FILE__])
    used = Set{Symbol}()
    function define!(lhs)
        lhs isa Symbol && push!(defined, lhs)
        lhs isa Expr && lhs.head in (:tuple, :parameters, :kw, :(::), :(=), :call) && foreach(define!, lhs.args)
    end
    function walk(e, ctx)
        if e isa Symbol
            push!(used, e)
        elseif e isa QuoteNode
            return
        elseif e isa Expr
            h = e.head
            if h === :function
                define!(e.args[1]); foreach(a -> walk(a, ctx), e.args[2:end])
            elseif h === :(=) && !(e.args[1] isa Expr && e.args[1].head === :.)
                define!(e.args[1]); walk(e.args[2], ctx)
            elseif h === :kw
                walk(e.args[2], ctx)               # keyword name is not a reference
            elseif h === :for
                define!(e.args[1].args[1]); walk(e.args[1].args[2], ctx); walk(e.args[2], ctx)
            elseif h === :-> || h === :do
                define!(e.args[1]); foreach(a -> walk(a, ctx), e.args[2:end])
            elseif h === :.
                walk(e.args[1], ctx)               # module.name: the module is the reference
            elseif h === :macrocall
                foreach(a -> walk(a, ctx), e.args[2:end])
            elseif h in (:using, :import, :line)
                return
            else
                foreach(a -> walk(a, ctx), e.args)
            end
        end
    end
    walk(ex, nothing)
    return setdiff(used, defined)
end
function unresolved_viewer_names(src::AbstractString)
    ex = Meta.parseall(src)
    names = free_names(ex)
    filter(n -> !(n in GLMAKIE_NAMES) && !isdefined(Base, n) && !isdefined(Core, n) &&
                !isdefined(HDF5, n) && !isdefined(Main, n) && !startswith(string(n), "@"), collect(names))
end

@testset "viewer GL half: parses, and every name it uses resolves statically" begin
    src = read(joinpath(REPO, "viewer", "visualize_lattice.jl"), String)
    @test Meta.parseall(src) isa Expr
    bad = unresolved_viewer_names(src)
    isempty(bad) || println("unresolved names in visualize_lattice.jl: ", bad)
    @test isempty(bad)
    # control: an undefined name planted in the source is reported by the same walker
    planted = replace(src, "save(still, fig)" => "save(still, fig); polish_frame(fig)")
    @test planted != src
    @test :polish_frame in unresolved_viewer_names(planted)
    # and the declared GLMakie list is not a blanket: a bare unknown symbol is still reported
    @test :not_a_makie_thing in unresolved_viewer_names(src * "\nnot_a_makie_thing(1)\n")
end

@testset "one run, a snapshot every k MCS, keyed for the series" begin
    sim3 = SR.init_coupled_simulation(p, rp; seed = 9)
    d = joinpath(tmp, "every")
    paths = export_transport_series(SR, sim3, d; every = 2, n_mcs = 7)
    @test basename.(paths) == ["snap_mcs000002.h5", "snap_mcs000004.h5", "snap_mcs000006.h5", "snap_mcs000007.h5"]
    @test [h5open(g -> Int(read(attributes(g)["mcs"])), q, "r") for q in paths] == [2, 4, 6, 7]
    written = export_series(d, joinpath(tmp, "every_run"))
    @test length(written) == 4
    pvd = read(joinpath(tmp, "every_run.pvd"), String)
    @test all(occursin("timestep=\"$t\"", pvd) for t in ("2.0", "4.0", "6.0", "7.0"))
    @test_throws ArgumentError export_transport_series(SR, sim3, d; every = 0, n_mcs = 4)
    @test_throws ArgumentError export_transport_series(SR, sim3, d; every = 5, n_mcs = 4)
end

@testset "restart companions: same state, same basename, nutrient carried exactly" begin
    simr = SR.init_coupled_simulation(p, rp; seed = 9)
    d = joinpath(tmp, "companion_snaps"); r = joinpath(tmp, "companion_restarts")
    paths = export_transport_series(SR, simr, d; every = 2, n_mcs = 4, restart_dir = r)
    @test basename.(paths) == ["snap_mcs000002.h5", "snap_mcs000004.h5"]
    @test sort(readdir(r)) == basename.(paths)
    for q in paths
        rq = joinpath(r, basename(q))
        @test h5open(g -> Int(read(attributes(g)["mcs"])), rq, "r") ==
              h5open(g -> Int(read(attributes(g)["mcs"])), q, "r")
    end

    written = export_series(d, joinpath(tmp, "companion_run"); restart_dir = r)
    @test length(written) == 2
    for (q, w) in zip(paths, written)
        cd = get_cell_data(VTKFile(w))
        @test "nutrient" in keys(cd)
        nut_vti = reshape(get_data(cd["nutrient"]), N)
        @test nut_vti == h5open(g -> read(g["fields/nutrient"]), joinpath(r, basename(q)), "r")
        @test vti_field_string(w, "nutrient_uptake_model") == "Monod"
        fd = get_field_data(VTKFile(w))
        @test only(get_data(fd["K_C_monod"])) == p.K_C_monod
        @test only(get_data(fd["C_wall"])) == p.C_wall
        @test occursin("not a physical concentration", vti_field_string(w, "nutrient_basis"))
    end

    # The two directories cannot be one: companions share basenames.
    @test_throws ArgumentError export_transport_series(SR, SR.init_coupled_simulation(p, rp; seed = 9),
                                                       joinpath(tmp, "one_dir"); every = 2, n_mcs = 2,
                                                       restart_dir = joinpath(tmp, "one_dir", "."))
end

@testset "a restart from another run at the same mcs and size is refused" begin
    # The mcs and size checks alone pass this pair; only the state can tell.
    other_run = SR.init_coupled_simulation(p, rp; seed = 10)
    SR.advance_window!(other_run, 2)
    wrong = joinpath(tmp, "restart_wrong_run.h5")
    export_restart_checkpoint(SR, other_run, wrong)
    @test h5open(g -> Int(read(attributes(g)["mcs"])), wrong, "r") ==
          h5open(g -> Int(read(attributes(g)["mcs"])), snap, "r")
    err = try
        export_vti(snap, joinpath(tmp, "wrong_run"); restart = wrong); nothing
    catch e
        sprint(showerror, e)
    end
    @test err !== nothing && occursin("another run", err)
    @test !isfile(joinpath(tmp, "wrong_run.vti"))
end

@testset "a pre-Monod restart is labelled unrecorded, never inferred as Monod" begin
    same = joinpath(tmp, "restart_same_premonod.h5")
    export_restart_checkpoint(SR, sim, same)
    h5open(g -> delete_object(g, "cpm_params/K_C_monod"), same, "r+")
    export_vti(snap, joinpath(tmp, "premonod"); restart = same)
    w = joinpath(tmp, "premonod.vti")
    @test startswith(vti_field_string(w, "nutrient_uptake_model"), "unrecorded")
    @test isnan(only(get_data(get_field_data(VTKFile(w))["K_C_monod"])))
end

@testset "without restart_dir the series writes snapshots only, as before" begin
    d = joinpath(tmp, "no_companions")
    paths = export_transport_series(SR, SR.init_coupled_simulation(p, rp; seed = 9), d;
                                    every = 2, n_mcs = 4)
    @test sort(readdir(d)) == basename.(paths)
    written = export_series(d, joinpath(tmp, "no_companions_run"))
    @test !("nutrient" in keys(get_cell_data(VTKFile(first(written)))))
    @test vti_field_string(first(written), "nutrient_uptake_model") === nothing
end
