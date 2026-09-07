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
    h5open(snap, "r") do h
        for (ds, name) in SITE_ARRAYS
            src = read(h[ds])
            back = reshape(get_data(get_cell_data(f)[name]), N)
            @test eltype(back) == eltype(src)
            @test back == src
        end
        species = reshape(get_data(get_cell_data(f)["species"]), N)
        @test eltype(species) == UInt8
        cid = read(h["lattice/cell_id"])
        @test all(species[cid .<= 0] .== 0)          # medium and wall are air
        @test species[cid .> 0] == UInt8.(read(h["lattice/species_id"])[cid .> 0])
        @test Set(unique(species)) ⊇ Set(0x01:0x07)  # the fixture places all seven
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
                         source_rate_photons_per_s = 1.0e8)
    export_vti(snap, joinpath(tmp, "with_target"); dose = target)
    @test occursin("target_calibration = true", vti_field_string(joinpath(tmp, "with_target.vti"), "dose_rate_mean_Gy_s_units"))
    bare = fake_result(joinpath(tmp, "tr_bare.h5"); source_rate_photons_per_s = 1.0e8)
    @test_throws ArgumentError export_vti(snap, joinpath(tmp, "bare"); dose = bare)
    @test !isfile(joinpath(tmp, "bare.vti"))
    wrong = joinpath(tmp, "tr_wrong.h5")
    h5open(wrong, "w") do g
        g["mesh/dose_rate_mean_Gy_s"] = zeros(4, 4, 4)
        attributes(g)["target_calibration"] = 0; attributes(g)["source_rate_photons_per_s"] = 1.0
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
