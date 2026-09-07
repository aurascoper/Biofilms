using HDF5, SHA
include("lattice_grid.jl")

"Read a signal companion only when its source bytes, MCS, axes and mask agree."
function signal_grid(snapshot, companion)
    grid,mcs=species_grid(snapshot)
    h5open(companion,"r") do f
        a=HDF5.attributes(f)
        read(a["mcs"])==mcs || throw(ArgumentError("signal MCS differs from snapshot"))
        read(a["logical_axis_order"])=="xyz" && read(a["dataset_axis_order_h5py"])=="zyx" ||
            throw(ArgumentError("unknown signal axis order"))
        read(a["parent_snapshot_sha256"])==bytes2hex(open(sha256,snapshot)) ||
            throw(ArgumentError("signal belongs to different snapshot bytes"))
        signal=read(f["fields/signal"])
        size(signal)==size(grid) || throw(ArgumentError("signal grid shape mismatch"))
        mask=h5open(g->read(g["lattice/interior_mask"]),snapshot,"r")
        read(f["lattice/interior_mask"])==mask || throw(ArgumentError("signal mask mismatch"))
        read(a["signal_sha256"])==bytes2hex(sha256(reinterpret(UInt8,vec(signal)))) ||
            throw(ArgumentError("signal array hash mismatch"))
        read(a["acceptance_coupling"])==0 || throw(ArgumentError("signal is not inert"))
        all(isfinite,signal) && minimum(signal)>=0 || throw(ArgumentError("invalid signal"))
        return grid,signal,mcs
    end
end

function signal_viewer_options(args)
    length(args)==2 || throw(ArgumentError("usage: visualize_signal.jl PARENT_RUN DERIVED_RUN"))
    all(isdir,args) || throw(ArgumentError("both run directories must exist"))
    (;parent=args[1],derived=args[2])
end
