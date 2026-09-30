"""Independent VTK readback, including source array hashes and parcel projections."""
from pathlib import Path
import argparse
import hashlib
import json
import xml.etree.ElementTree as ET
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def read_image(path, receipt):
    path = Path(path)
    if sha(path) != receipt['sha256']:
        raise ValueError('frame hash mismatch')
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(str(path)); reader.Update()
    image = reader.GetOutput()
    shape = tuple(receipt['shape'])
    if image.GetDimensions() != tuple(x+1 for x in shape):
        raise ValueError('point dimensions do not enclose the declared cells')
    if image.GetNumberOfCells() != np.prod(shape) or image.GetPointData().GetNumberOfArrays():
        raise ValueError('wrong association or tuple count')
    if tuple(image.GetOrigin()) != tuple(receipt['origin']) or tuple(image.GetSpacing()) != tuple(receipt['spacing']):
        raise ValueError('origin/spacing mismatch')
    arrays = {}
    data = image.GetCellData()
    if {data.GetArrayName(i) for i in range(data.GetNumberOfArrays())} != set(receipt['arrays']):
        raise ValueError('array inventory mismatch')
    types = {'Float64': np.dtype('float64'), 'Int64': np.dtype('int64'), 'UInt8': np.dtype('uint8')}
    for name, expected in receipt['arrays'].items():
        a = vtk_to_numpy(data.GetArray(name))
        if len(a) != np.prod(shape) or a.dtype != types[receipt['dtypes'][name]]:
            raise ValueError(f'{name}: dtype or tuple count mismatch')
        if hashlib.sha256(a.tobytes()).hexdigest() != expected:
            raise ValueError(f'{name}: source-value/orientation mismatch')
        arrays[name] = a
    return image, arrays

def validate_run(directory):
    directory = Path(directory)
    manifest = json.loads((directory/'run_manifest.json').read_text())
    if manifest['schema'] != 'cpm-sp-vti-1' or manifest['association'] != 'CellData':
        raise ValueError('unsupported schema')
    frames = manifest['frames']
    entries = ET.parse(directory/'series.pvd').findall('.//DataSet')
    times = [int(e.attrib['timestep']) for e in entries]
    expected = [f['mcs'] for f in frames]
    if expected != manifest.get('expected_frame_times', expected):
        raise ValueError('frame inventory differs from declared schedule')
    if times != expected or times != sorted(set(times)) or times[0] != 0 or times[-1] != manifest['horizon']:
        raise ValueError('missing, duplicate or incorrect time')
    if [e.attrib['file'] for e in entries] != [f['file'] for f in frames]:
        raise ValueError('PVD points to different frames')
    metrics = {x['mcs']: x for x in json.loads((directory/'metrics.json').read_text())}
    for frame in frames:
        path = Path(frame['file'])
        if path.is_absolute() or '..' in path.parts:
            raise ValueError('nonportable frame path')
        receipt = json.loads((directory/(str(path)+'.json')).read_text())
        if receipt != frame:
            raise ValueError('frame sidecar differs from manifest')
        for key in ('source', 'run_id'):
            if frame[key] != manifest[key]:
                raise ValueError('mixed run/source identity')
        if frame['config'] != manifest['config_hash']:
            raise ValueError('mixed configuration')
        _, a = read_image(directory/path, frame)
        occupied = a['cell_id'] > 0
        if not np.array_equal(a['occupied'], occupied) or np.any(occupied & ~a['interior_mask'].astype(bool)):
            raise ValueError('incorrect occupied/domain mask')
        rows = {x['id']: x for x in metrics[frame['mcs']]['parcels'] if x['volume'] > 0}
        if set(np.unique(a['cell_id'][occupied])) != set(rows):
            raise ValueError('parcel inventory mismatch')
        for parcel_id, row in rows.items():
            mask = a['cell_id'] == parcel_id
            if np.count_nonzero(mask) != row['volume'] or not np.all(a['species_id'][mask] == row['species']):
                raise ValueError('parcel volume/species mismatch')
            if not np.all(a['alive'][mask] == (row['status'] == 'alive')):
                raise ValueError('viability mismatch')
            if not np.all(a['parcel_cum_dose'][mask] == row['dose']):
                raise ValueError('dose projection mismatch')
        phase_on = manifest['arm'] in ('phase_only', 'combined')
        if not np.array_equal(a['phase_valid'], a['alive'] if phase_on else np.zeros_like(a['alive'])):
            raise ValueError('invalid phase mask')
        if np.any(np.isfinite(a['phase_rad'][a['phase_valid'] == 0])):
            raise ValueError('undefined phase presented as a value')
        if manifest['growth']:
            from_table = sum(row['biomass'] for row in rows.values())
            integrated = np.nansum(a['biomass_density'])*np.prod(frame['spacing'])
            if not np.isclose(integrated, from_table, rtol=1e-12, atol=1e-9):
                raise ValueError('biomass budget/projection mismatch')
        elif np.any(a['biomass_valid']) or np.any(np.isfinite(a['biomass'])):
            raise ValueError('undefined biomass presented as measured zero')
        if phase_on:
            values = [a['phase_rad'][np.flatnonzero(a['cell_id'] == i)[0]] for i,row in rows.items() if row['status']=='alive']
            expected_r = metrics[frame['mcs']]['phase']['coherence']
            actual = abs(np.exp(1j*np.array(values)).mean()) if values else None
            if (actual is None) != (expected_r is None) or (actual is not None and not np.isclose(actual,expected_r,atol=1e-12)):
                raise ValueError('parcel-count coherence mismatch')
    return {'run_id':manifest['run_id'],'source':manifest['source'],'config_hash':manifest['config_hash'],
            'manifest_sha256':sha(directory/'run_manifest.json'),'vtk':vtk.vtkVersion.GetVTKVersion(),
            'frames_checked':len(frames),'times':times,'status':'passed'}

if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('run',type=Path)
    args=parser.parse_args();result=validate_run(args.run)
    (args.run/'export_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
