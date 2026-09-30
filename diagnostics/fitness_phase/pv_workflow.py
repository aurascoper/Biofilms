"""Run with pvpython. Build four views or reopen with explicit relocation."""
import argparse
import json
import math
from pathlib import Path
from paraview import simple as pv
from paraview import servermanager
from validate_exports import validate_run, sha

def build(run, output, clip=False):
    receipt=validate_run(run)
    manifest=json.loads((run/'run_manifest.json').read_text())
    output.mkdir(parents=True,exist_ok=True)
    reader=pv.PVDReader(registrationName='CPM series',FileName=str(run/'series.pvd'))
    reader.UpdatePipelineInformation()
    times=list(reader.TimestepValues)
    if times != receipt['times']: raise ValueError('ParaView time mismatch')
    scene=pv.GetAnimationScene();scene.UpdateAnimationUsingDataTimeSteps()
    n=manifest['config']['N'];source=reader
    if clip:
        # Cut away the half facing the camera so the axis, where dose peaks, is visible.
        # Crinkle clip keeps whole voxels: IDs and categories are never interpolated.
        source=pv.Clip(registrationName='Axis clip',Input=reader);source.ClipType='Plane'
        source.ClipType.Origin=[n/2,n/2,n/2];source.ClipType.Normal=[1,-1,0];source.Crinkleclip=1;source.Invert=1
        full=servermanager.Fetch(reader).GetNumberOfCells();kept=servermanager.Fetch(source).GetNumberOfCells()
        if not 0<kept<full:raise ValueError(f'clip kept {kept} of {full} cells')
    layout=pv.CreateLayout(name='CPM scientific views')
    layout.SplitHorizontal(0,0.5);layout.SplitVertical(1,0.5);layout.SplitVertical(2,0.5)
    specs=[('species_id','occupied',1,7,'Species / parcel geometry'),
           ('net_biomass_rate_proxy','alive',-.06,.05,'Growth rate − mortality hazard [1/MCS]'),
           ('phase_rad','phase_valid',0,2*math.pi,'Phase [rad]; invalid parcels excluded'),
           ('dose_rate','interior_mask',0,.01,'Synthetic dose rate [Gy/s]')]
    colors=[.25,.7,.45, .93,.69,.13, .12,.65,.86, .8,.42,.3, .68,.48,.76, .89,.53,.77, .6,.6,.6]
    views=[];clocks=[]
    for slot,(name,mask,lo,hi,title) in zip((3,4,5,6),specs):
        view=pv.CreateView('RenderView');pv.AssignViewToLayout(view=view,layout=layout,hint=slot)
        view.Background=[.08,.09,.12];view.OrientationAxesVisibility=0
        filt=pv.Threshold(registrationName=f'{title}: {mask}',Input=source)
        filt.Scalars=['CELLS',mask];filt.LowerThreshold=.5;filt.UpperThreshold=1.5;filt.ThresholdMethod='Between'
        display=pv.Show(filt,view);display.Representation='Surface';pv.ColorBy(display,('CELLS',name))
        lut=pv.GetColorTransferFunction(name);lut.RescaleTransferFunction(lo,hi)
        lut.AutomaticRescaleRangeMode='Never'
        if name=='species_id':
            lut.InterpretValuesAsCategories=1;lut.Annotations=[v for i,s in enumerate(manifest['species'],1) for v in (str(i),s)]
            lut.IndexedColors=colors
        elif name=='phase_rad':
            lut.RGBPoints=[0,1,0,0, math.pi/3,1,1,0,2*math.pi/3,0,1,0,math.pi,0,1,1,4*math.pi/3,0,0,1,5*math.pi/3,1,0,1,2*math.pi,1,0,0]
            lut.ColorSpace='RGB'
        else:
            lut.RGBPoints=[lo,.267,.005,.329,(lo+hi)/2,.128,.567,.551,hi,.993,.906,.144]
            lut.ColorSpace='RGB'
        display.SetScalarBarVisibility(view,True)
        bar=pv.GetScalarBar(lut,view);bar.Title=title;bar.ComponentTitle=''
        text=pv.Text();text.Text=f"{manifest['arm']} | growth={int(manifest['growth'])} | {title}"
        label=pv.Show(text,view);label.FontSize=11
        clock=pv.AnnotateTimeFilter(Input=reader);clock.Format='MCS: {time:.0f}'
        clocks.append(clock)
        stamp=pv.Show(clock,view);stamp.WindowLocation='Lower Left Corner';stamp.FontSize=12
        view.CameraPosition=[2.2*n,-2.2*n,1.7*n];view.CameraFocalPoint=[n/2,n/2,n/2]
        view.CameraViewUp=[0,0,1];view.CameraParallelProjection=1;view.CameraParallelScale=.95*n
        views.append(view)
    # Saved, named toggles share the net-rate view but start hidden.
    for name,hi in [('conditional_survival',1.0),('mortality_hazard',.06)]:
        filt=pv.Threshold(registrationName=f'Toggle: {name}',Input=source)
        filt.Scalars=['CELLS','alive'];filt.LowerThreshold=.5;filt.UpperThreshold=1.5;filt.ThresholdMethod='Between'
        display=pv.Show(filt,views[1]);pv.ColorBy(display,('CELLS',name))
        lut=pv.GetColorTransferFunction(name);lut.RescaleTransferFunction(0,hi);lut.AutomaticRescaleRangeMode='Never'
        pv.Hide(filt,views[1])
    layout.SetSize(1400,1000)
    selected=sorted(set([times[0],times[len(times)//2],times[-1]]))
    images=[]
    for t in selected:
        scene.AnimationTime=t
        for v in views:v.ViewTime=t
        for clock in clocks:
            clock.UpdatePipeline(time=t)
            if servermanager.Fetch(clock).GetColumn(0).GetValue(0)!=f'MCS: {t:.0f}':
                raise ValueError('displayed time does not match the data')
        pv.RenderAllViews()
        target=output/f'frame_{int(t):06}.png'
        pv.SaveScreenshot(str(target),layout,ImageResolution=[1400,1000])
        images.append({'time':t,'file':target.name,'sha256':sha(target)})
    state=output/'review.pvsm';pv.SaveState(str(state))
    record={'run_id':manifest['run_id'],'source':manifest['source'],'manifest_sha256':sha(run/'run_manifest.json'),
            'times':times,'paraview':str(pv.GetParaViewVersion()),'views':specs,'images':images,'state_sha256':sha(state),
            'reopened':False,'pipeline_sha256':sha(__file__),'time_annotations_verified':True,'axis_clip':clip,
            'note':'Only recorded simulation frames; fixed color limits and common camera.'}
    (output/'render_receipt.json').write_text(json.dumps(record,indent=2)+'\n')
    return record

def reopen(run, output):
    receipt=validate_run(run);manifest=json.loads((run/'run_manifest.json').read_text())
    record=json.loads((output/'render_receipt.json').read_text())
    if record['manifest_sha256']!=sha(run/'run_manifest.json') or record['run_id']!=manifest['run_id']:
        raise ValueError('state belongs to another run')
    if sha(output/'review.pvsm')!=record['state_sha256']:raise ValueError('state changed')
    pv.LoadState(str(output/'review.pvsm'),filenames=[{'name':'CPM series','FileName':str(run/'series.pvd')}])
    reader=pv.FindSource('CPM series');reader.UpdatePipelineInformation()
    if list(reader.TimestepValues)!=record['times']:raise ValueError('relocated times changed')
    scene=pv.GetAnimationScene();scene.AnimationTime=record['times'][-1]
    reader.UpdatePipeline(time=scene.AnimationTime)
    data=servermanager.Fetch(reader)
    if data.GetCellData().GetArray('cell_id') is None or data.GetCellData().GetArray('phase_valid') is None:
        raise ValueError('relocated arrays missing')
    views=pv.GetViews()
    if len([v for v in views if v.IsA('vtkSMRenderViewProxy')])!=4:raise ValueError('view count changed')
    for view in views:
        if hasattr(view,'ViewTime'):view.ViewTime=scene.AnimationTime
    clocks=[proxy for proxy in pv.GetSources().values() if proxy.SMProxy.GetXMLName()=='TimeToTextConvertor']
    if len(clocks)!=4:raise ValueError('time annotations missing')
    for clock in clocks:
        clock.UpdatePipeline(time=scene.AnimationTime)
        if servermanager.Fetch(clock).GetColumn(0).GetValue(0)!=f'MCS: {scene.AnimationTime:.0f}':
            raise ValueError('relocated displayed time mismatch')
    pv.RenderAllViews()
    pv.SaveScreenshot(str(output/'relocated_final.png'),pv.GetLayout(),ImageResolution=[1400,1000])
    receipt.update({'reopened':True,'paraview':str(pv.GetParaViewVersion()),'state_sha256':sha(output/'review.pvsm'),
        'render_receipt_sha256':sha(output/'render_receipt.json'),'image_sha256':sha(output/'relocated_final.png')})
    (output/'reopen_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['build','reopen']);parser.add_argument('run',type=Path);parser.add_argument('output',type=Path)
    parser.add_argument('--clip',action='store_true',help='build only: cut away the camera-facing half to show the axis')
    a=parser.parse_args()
    print(json.dumps(build(a.run.resolve(),a.output.resolve(),a.clip) if a.mode=='build' else reopen(a.run.resolve(),a.output.resolve())))
