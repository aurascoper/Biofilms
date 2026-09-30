import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import xml.etree.ElementTree as ET
import numpy as np
import pytest
from validate_exports import read_image, validate_run

ROOT=Path(__file__).resolve().parents[2]

@pytest.fixture(scope='module')
def smoke():
    path=Path(os.environ['CPM_SMOKE_RUN'])
    assert (path/'run_manifest.json').is_file()
    return path

def test_real_run(smoke):
    assert validate_run(smoke)['frames_checked'] == 5

@pytest.mark.parametrize('mutation',['missing','duplicate','mask','decoy','repeated_mass'])
def test_production_reader_refuses_defects(smoke,tmp_path,mutation):
    run=tmp_path/'run';shutil.copytree(smoke,run)
    before=(run/'series.pvd').read_bytes()
    tree=ET.parse(run/'series.pvd');collection=tree.find('Collection')
    if mutation=='missing':
        collection.remove(collection[-1]);tree.write(run/'series.pvd')
    elif mutation=='duplicate':
        collection.append(ET.fromstring(ET.tostring(collection[0])));tree.write(run/'series.pvd')
    else:
        manifest=json.loads((run/'run_manifest.json').read_text())
        frame=manifest['frames'][0]
        if mutation=='decoy': frame['run_id']='another-run'
        else:
            # Change actual VTI values and update hashes, so semantic checks
            # rather than a stale checksum must reject the bad projection.
            import vtk
            from vtk.util.numpy_support import vtk_to_numpy
            reader=vtk.vtkXMLImageDataReader();reader.SetFileName(str(run/frame['file']));reader.Update()
            data=reader.GetOutput()
            name='phase_valid' if mutation=='mask' else 'biomass_density'
            values=vtk_to_numpy(data.GetCellData().GetArray(name));before_values=values.copy()
            if mutation=='mask':values[:]=0
            else:
                assert manifest['growth'], 'mass control requires growth smoke'
                values[:]=vtk_to_numpy(data.GetCellData().GetArray('biomass'))
            assert not np.array_equal(values,before_values,equal_nan=True)
            writer=vtk.vtkXMLImageDataWriter();writer.SetFileName(str(run/frame['file']));writer.SetInputData(data);assert writer.Write()==1
            frame['sha256']=hashlib.sha256((run/frame['file']).read_bytes()).hexdigest()
            frame['arrays'][name]=hashlib.sha256(values.tobytes()).hexdigest()
        (run/'run_manifest.json').write_text(json.dumps(manifest))
        (run/(frame['file']+'.json')).write_text(json.dumps(frame))
        assert frame != json.loads((smoke/'run_manifest.json').read_text())['frames'][0]
    if mutation in ('missing','duplicate'): assert (run/'series.pvd').read_bytes()!=before
    with pytest.raises(ValueError):validate_run(run)

def test_asymmetric_orientation_and_int64(tmp_path):
    stem=tmp_path/'asymmetric'
    code=r'''
include("diagnostics/fitness_phase/model.jl")
using .FitnessPhase, JSON3, SHA
a=reshape(Int64[2^53+1+i+10*j+100*k for i in 1:2,j in 1:3,k in 1:4],2,3,4)
d=Dict("cell_id"=>a)
FitnessPhase.write_image(ARGS[1],d;origin=(1.,2.,3.),spacing=(.5,1.,2.))
meta=(shape=size(a),origin=(1.,2.,3.),spacing=(.5,1.,2.),sha256=bytes2hex(sha256(read(ARGS[1]*".vti"))),arrays=Dict("cell_id"=>bytes2hex(sha256(reinterpret(UInt8,vec(a))))),dtypes=Dict("cell_id"=>"Int64"))
open(io->JSON3.write(io,meta),ARGS[1]*".json","w")
'''
    subprocess.run(['julia','--project=.','-e',code,str(stem)],cwd=ROOT,check=True)
    receipt=json.loads(stem.with_suffix('.json').read_text())
    _,arrays=read_image(stem.with_suffix('.vti'),receipt)
    grid=arrays['cell_id'].reshape((2,3,4),order='F')
    assert grid[1,2,3]==2**53+1+2+30+400
    bad=dict(receipt);bad['shape']=[3,2,4]
    with pytest.raises(ValueError):read_image(stem.with_suffix('.vti'),bad)
    bad=dict(receipt);bad['arrays']={'cell_id':hashlib.sha256(grid.transpose(1,0,2).tobytes(order='F')).hexdigest()}
    assert bad['arrays']!=receipt['arrays']
    with pytest.raises(ValueError):read_image(stem.with_suffix('.vti'),bad)
