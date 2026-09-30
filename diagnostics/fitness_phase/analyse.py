"""Seed-level analysis of the declared four-arm synthetic study; refuse incomplete sets."""
from pathlib import Path
import argparse
import json
import numpy as np

ARMS=('baseline','survival_only','phase_only','combined')

def interval(values,rng):
    values=np.array(values,dtype=float)
    if not len(values):return {'n':0,'mean':None,'ci95':None}
    means=values[rng.integers(0,len(values),(5000,len(values)))].mean(axis=1)
    return {'n':len(values),'mean':float(values.mean()),'ci95':np.quantile(means,[.025,.975]).tolist()}

def analyse(root,first,last,horizon,growth):
    root=Path(root);rng=np.random.default_rng(9004);groups=[];sources=set();configs=set()
    for seed in range(first,last+1):
        directory=root/f'growth_{int(growth)}_seed_{seed}'
        invariant=json.loads((directory/'invariants.json').read_text())
        if not invariant['exact'] or not invariant['restart_exact'] or invariant['steps']!=horizon+1:
            raise ValueError('invariant receipt incomplete')
        rows={}
        for arm in ARMS:
            path=directory/arm;manifest=json.loads((path/'run_manifest.json').read_text())
            if manifest['seed']!=seed or manifest['arm']!=arm or manifest['growth']!=growth or manifest['horizon']!=horizon or manifest['dirty'] or manifest['development']:
                raise ValueError('run differs from production plan')
            sources.add(manifest['source']);configs.add(manifest['config_hash'])
            metrics=json.loads((path/'metrics.json').read_text())
            if [x['mcs'] for x in metrics]!=list(range(horizon+1)):raise ValueError('incomplete trajectory')
            rows[arm]={'metrics':metrics,'window':json.loads((path/'phase_window.json').read_text()),'events':json.loads((path/'events.json').read_text())}
        for a,b in [('baseline','phase_only'),('survival_only','combined')]:
            if [x['nonphase_hash'] for x in rows[a]['metrics']]!=[x['nonphase_hash'] for x in rows[b]['metrics']]:
                raise ValueError('canonical history mismatch')
        groups.append(rows)
    if len(sources)!=1 or len(configs)!=1:raise ValueError('mixed revisions/configurations')
    result={'source':sources.pop(),'config_hash':configs.pop(),'seeds':[first,last],'seed_count':len(groups),'horizon':horizon,'growth':growth,
            'primary_runs':4*len(groups),'paired_nonphase_identity':True,'restart_identity':True,'species':[],'phase':{}}
    for species in range(7):
        entry={'species_index':species+1}
        for quantity in ('alive','founder_survivors','mean_radius','initial_radius','displacement','biomass'):
            values=[];excluded=0
            for group in groups:
                a=group['survival_only']['metrics'][-1]['species'][species][quantity]
                b=group['baseline']['metrics'][-1]['species'][species][quantity]
                if a is None or b is None:excluded+=1
                else:values.append(a-b)
            entry[quantity]=interval(values,rng)|{'excluded_seeds':excluded,'contrast':'survival_only minus baseline','conditioning':'evaluable paired seeds; no all-seed radius verdict when missing'}
        result['species'].append(entry)
    for arm in ('phase_only','combined'):
        windows=[g[arm]['window'] for g in groups]
        last_phase=[g[arm]['metrics'][-1]['phase'] for g in groups]
        result['phase'][arm]={
            'coherence':interval([x['coherence'] for x in last_phase if x['coherence'] is not None],rng),
            'active_count':interval([len(x['ids']) for x in last_phase],rng),
            'components':interval([len(x['component_sizes']) for x in last_phase],rng),
            'isolates':interval([x['isolates'] for x in last_phase],rng),
            'fixed_cohort_evaluable':sum(w['count']>=2 for w in windows),
            'finite_window_locking':sum(w['locking'] is True for w in windows),
            'coupling_supported':sum(w.get('coupling_supported') is True for w in windows),
            'selection_sum_complete_trajectories':interval([sum(x['phase']['selection'] for x in g[arm]['metrics'][1:]) for g in groups if all(x['phase'].get('selection') is not None for x in g[arm]['metrics'][1:])],rng),
            'interpretation':'finite-window synthetic diagnostics; selection and division can change coherence without locking'}
    result['deaths']={arm:sum(sum(e['kind']=='hazard_commitment' for e in g[arm]['events']) for g in groups) for arm in ARMS}
    matched=[];phase_deltas=[]
    for group in groups:
        pa=group['phase_only']['metrics'][-1];pb=group['combined']['metrics'][-1]
        founders_a={x['id'] for x in pa['parcels'] if x['parent']==0 and x['status']=='alive'}
        founders_b={x['id'] for x in pb['parcels'] if x['parent']==0 and x['status']=='alive'}
        ids=sorted(founders_a & founders_b)
        if len(ids)>=2:
            def coh(row):
                p=row['phase'];angles=[p['theta'][p['ids'].index(i)] for i in ids]
                return float(abs(np.exp(1j*np.array(angles)).mean()))
            matched.append(coh(pb)-coh(pa))
        a=pa['phase']['coherence'];b=pb['phase']['coherence']
        if a is not None and b is not None:phase_deltas.append(b-a)
    result['phase_combined_minus_phase_only']={'all_current_ids':interval(phase_deltas,rng),'same_founder_ids_at_endpoint':interval(matched,rng),
        'interpretation':'Matched membership diagnostic, not a graph-controlled causal contrast; excludes divided founders in growth mode.'}
    result['nonphase_factorial_interaction']='exactly zero for every defined shared endpoint, established by paired full-state identity'
    result['divisions']={arm:sum(sum(e['kind']=='division' for e in g[arm]['events']) for g in groups) for arm in ARMS}
    if growth:
        budgets=[x['budget'] for g in groups for arm in ARMS for x in g[arm]['metrics']]
        result['max_abs_nutrient_residual']=max(abs(x['nutrient_residual']) for x in budgets)
        result['max_abs_biomass_residual']=max(abs(x['biomass_residual']) for x in budgets)
        if result['max_abs_nutrient_residual']>1e-6 or result['max_abs_biomass_residual']>1e-6:raise ValueError('budget tolerance exceeded')
    return result

def plot_diagnostic(root,growth,seed,output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,3,figsize=(13,7),layout='constrained')
    for arm in ARMS:
        rows=json.loads((Path(root)/f'growth_{int(growth)}_seed_{seed}'/arm/'metrics.json').read_text())
        t=[x['mcs'] for x in rows]
        axes[0,0].plot(t,[sum(s['alive'] for s in x['species']) for x in rows],label=arm)
        axes[0,1].plot(t,[x['budget']['biomass'] if x['budget'] else np.nan for x in rows],label=arm)
        if arm in ('phase_only','combined'):
            axes[0,2].plot(t,[x['phase']['coherence'] for x in rows],label=arm)
            axes[1,0].plot(t,[len(x['phase']['component_sizes']) for x in rows],label=arm)
            axes[1,1].plot(t,np.cumsum([0.]+[x['phase'].get('selection') if x['phase'].get('selection') is not None else np.nan for x in rows[1:]]),label=arm)
            # Same fixed IDs for every value in this window; retain unwrapped drift.
            selected=rows[len(rows)//2:];ids=set(selected[0]['phase']['ids'])
            for row in selected:ids.intersection_update(row['phase']['ids'])
            if len(ids)>=2:
                ids=sorted(ids);theta=np.array([[row['phase']['theta'][row['phase']['ids'].index(i)] for i in ids] for row in selected])
                delta=theta[:,1:]-theta[:,[0]]
                axes[1,2].plot([x['mcs'] for x in selected],np.max(np.abs(delta-delta[0]),axis=1),label=arm)
    titles=['Viable parcels','Occupied biomass (synthetic units)','Count-weighted coherence','Graph components','Cumulative immediate selection effect','Max unwrapped relative-phase change\n(fixed survivor cohort, reference ID)']
    for ax,title in zip(axes.flat,titles):ax.set_title(title);ax.set_xlabel('MCS');ax.grid(alpha=.25)
    axes[0,0].legend(fontsize=8);axes[0,2].set_ylim(0,1)
    fig.suptitle(f'Synthetic model: seed {seed}, growth={int(growth)}; diagnostic, not ensemble evidence')
    fig.savefig(output,dpi=160);plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('first',type=int);p.add_argument('last',type=int);p.add_argument('horizon',type=int);p.add_argument('growth',type=int,choices=[0,1]);p.add_argument('output',type=Path)
    a=p.parse_args();result=analyse(a.root,a.first,a.last,a.horizon,bool(a.growth));a.output.write_text(json.dumps(result,indent=2)+'\n')
    plot_diagnostic(a.root,bool(a.growth),a.first,a.output.with_suffix('.png'));print(a.output)
