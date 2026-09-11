"""Minimal read-only .vti reader: uncompressed appended binary, no VTK needed.

The exporter writes format="appended" with header_type="UInt64" and no compression,
so each array is a UInt64 byte count followed by that many raw little-endian bytes,
at `_marker + offset`. numpy is the only dependency, which is the point -- a portable
audit should not need the renderer's toolchain to check the renderer's output.
"""
import re, numpy as np

_DT = {"UInt8": np.uint8, "Int8": np.int8, "UInt16": np.uint16, "Int16": np.int16,
       "UInt32": np.uint32, "Int32": np.int32, "UInt64": np.uint64, "Int64": np.int64,
       "Float32": np.float32, "Float64": np.float64}


def read_vti(path):
    raw = open(path, "rb").read()
    head = raw[:raw.index(b"<AppendedData")].decode("utf-8", "replace")
    ext = [int(v) for v in re.search(r'WholeExtent="([^"]+)"', head).group(1).split()]
    nx, ny, nz = ext[1] - ext[0], ext[3] - ext[2], ext[5] - ext[4]   # CELL counts
    marker = raw.index(b"_", raw.index(b"<AppendedData")) + 1

    arrays = {}
    for m in re.finditer(r'<DataArray type="(\w+)" Name="(\w+)"[^>]*offset="(\d+)"', head):
        dtype, name, off = _DT[m.group(1)], m.group(2), int(m.group(3))
        n = int(np.frombuffer(raw, np.uint64, count=1, offset=marker + off)[0])
        buf = np.frombuffer(raw, dtype, count=n // dtype().itemsize, offset=marker + off + 8)
        # VTK writes x fastest; order="F" gives [x, y, z] indexing.
        arrays[name] = buf.reshape((nx, ny, nz), order="F") if buf.size == nx*ny*nz else buf

    # FieldData strings are <Array type="String"> -- NOT <DataArray> -- and in this
    # exporter they are appended, not inline. Matching only <DataArray> returned {} and
    # lost every provenance string without saying so.
    fields = {}
    fd = re.search(r"<FieldData>(.*?)</FieldData>", head, re.S)
    if fd:
        for m in re.finditer(r'<Array type="String" Name="(\w+)"[^>]*offset="(\d+)"', fd.group(1)):
            off = int(m.group(2))
            n = int(np.frombuffer(raw, np.uint64, count=1, offset=marker + off)[0])
            # payload is NUL-terminated and n includes the terminator
            fields[m.group(1)] = raw[marker + off + 8: marker + off + 8 + n].rstrip(b"\x00").decode("utf-8", "replace")
        for m in re.finditer(r'<DataArray type="(\w+)" Name="(\w+)"[^>]*offset="(\d+)"', fd.group(1)):
            dtype, off = _DT[m.group(1)], int(m.group(3))
            n = int(np.frombuffer(raw, np.uint64, count=1, offset=marker + off)[0])
            fields[m.group(2)] = np.frombuffer(raw, dtype, count=n // dtype().itemsize,
                                               offset=marker + off + 8)[0]
    return arrays, fields, (nx, ny, nz)
