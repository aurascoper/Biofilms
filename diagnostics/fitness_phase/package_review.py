"""Assemble a review directory; final acceptance remains a separate verified step."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

ARMS=('baseline','survival_only','phase_only','combined')

def assemble(raw,output,source_repo):
    if output.exists():raise ValueError('refusing to overwrite a package')
    plans=[(0,42,297,100),(1,42,57,200)]
    runs=[];sources=set();configs=set()
    for growth,first,last,horizon in plans:
        for seed in range(first,last+1):
            base=raw/f'growth_{growth}_seed_{seed}'
            invariant=json.loads((base/'invariants.json').read_text())
            if not invariant['exact'] or not invariant['restart_exact']:raise ValueError('failed invariant')
            for arm in ARMS:
                run=base/arm;m=json.loads((run/'run_manifest.json').read_text())
                v=json.loads((run/'export_validation.json').read_text())
                if m['development'] or m['dirty'] or m['horizon']!=horizon:raise ValueError('not a frozen outcome run')
                if v['status']!='passed' or v['manifest_sha256']!=hashlib.sha256((run/'run_manifest.json').read_bytes()).hexdigest():
                    raise ValueError('missing/stale independent readback')
                sources.add(m['source']);configs.add(m['config_hash']);runs.append(run)
    if len(sources)!=1 or len(configs)!=1:raise ValueError('mixed source/config')
    source=sources.pop();output.mkdir(parents=True)
    # A shallow checkout preserves the actual commit and lets strict restart and
    # production commands work after relocation without shipping prior history.
    checkout=output/'source'
    subprocess.run(['git','init','-q',str(checkout)],check=True)
    subprocess.run(['git','-C',str(checkout),'fetch','--depth=1',source_repo.resolve().as_uri(),source],check=True)
    subprocess.run(['git','-C',str(checkout),'checkout','--detach','-q','FETCH_HEAD'],check=True)
    shutil.copy2(source_repo/'Manifest.toml',checkout/'Manifest.toml')
    for run in runs:
        rel=run.relative_to(raw)
        if run.parent.name.endswith('_seed_42'):
            target=output/'diagnostic_runs'/rel;target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copytree(run,target)
            shutil.copy2(run.parent/'invariants.json',target.parent/'invariants.json')
        else:
            target=output/'ensemble_records'/rel;target.mkdir(parents=True)
            # Compact records are labeled as such: no PVD or incomplete time
            # series masquerades as a portable diagnostic run.
            for name in ('metrics.json','events.json','phase_window.json','replay_receipt.json'):
                with (run/name).open('rb') as src,gzip.open(target/(name+'.gz'),'wb') as dst:shutil.copyfileobj(src,dst)
            for name in ('run_manifest.json','export_validation.json'):
                shutil.copy2(run/name,target/name)
            shutil.copy2(run.parent/'invariants.json',target.parent/'invariants.json')
    for growth,_,_,_ in plans:
        for arm in ARMS:
            record=json.loads((raw/'views'/f'growth_{growth}_{arm}'/'render_receipt.json').read_text())
            manifest=raw/f'growth_{growth}_seed_42'/arm/'run_manifest.json'
            if record['manifest_sha256']!=hashlib.sha256(manifest.read_bytes()).hexdigest() or not record.get('time_annotations_verified'):
                raise ValueError('missing/stale verified render')
    for name in ('controls','verification','logs','views'):
        shutil.copytree(raw/name,output/name)
    tools=output/'review_tools';tools.mkdir()
    for name in ('analyse.py','validate_exports.py','pv_workflow.py','package_review.py','schema.json','species_registry.csv'):
        shutil.copy2(Path(__file__).parent/name,tools/name)
    for p in raw.glob('*analysis.*'):shutil.copy2(p,output/p.name)
    shutil.copy2(raw/'execution_jobs.json',output/'execution_jobs.json')
    manifest={'source':source,'config_hash':configs.pop(),'primary_runs':len(runs),'diagnostic_runs':8,
              'scope':'Eight complete dense diagnostic runs, plus compressed per-MCS ensemble records. Sparse VTI/checkpoints for other seeds remain in the raw archive.',
              'raw_archive_local':str(raw.resolve()),'complete':False,'pending':['relocated ParaView acceptance','final report','checksums']}
    (output/'package_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return manifest

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('raw',type=Path);p.add_argument('output',type=Path);p.add_argument('source_repo',type=Path)
    a=p.parse_args();print(json.dumps(assemble(a.raw,a.output,a.source_repo)))
