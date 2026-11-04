"""Collect completed runs without interpreting numerical failures as bounds."""
from common import *
import shutil
import importlib.metadata

if __name__=='__main__':
    reserve();out=ROOT/'evidence';out.mkdir(exist_ok=True);batches=[]
    for name,expected,key in [('patches842',192,'best'),('insertions842',192,'best'),
                              ('break_contacts',144,'best'),('sparse_insertions842',96,'best'),
                              ('direct841',128,'final'),('projected841',64,'minimax'),
                              ('trade_refine',128,'best'),('repair841',16,'final')]:
        r=json.loads((ROOT/'runs'/name/'summary.json').read_text());assert len(r)==expected,(name,len(r))
        counts=sorted({z[key]['n'] for z in r})
        grouped=[]
        for count in counts:
            group=[z for z in r if z[key]['n']==count];best=min(group,key=lambda z:z[key]['max_cos'])
            grouped.append(dict(n=count,trials=len(group),best_max_cos=best[key]['max_cos'],best_index=best['index'],
                                 exact_successes=sum(bool(z.get('exact_certificate',{}).get('passed')) for z in group)))
        row=dict(name=name,trials=len(r),groups=grouped);batches.append(row)
        save_json(out/(name+'_summary.json'),row)
        config=ROOT/'runs'/name/'config.json'
        if config.exists():shutil.copyfile(config,out/(name+'_config.json'))
    control=json.loads((ROOT/'runs'/'control841_cpu'/'summary.json').read_text());assert len(control)==64
    assert all(z['trace'][-1]['stage']==12 for z in control)
    default=json.loads((ROOT/'runs'/'control841_default512'/'summary.json').read_text());assert default['steps']==35000
    controls=dict(paper_batch64=dict(starts=64,full_schedule_steps=35000,
                                    best_max_cos=min(z['best_max'] for z in control),
                                    first_9000_steps='original PyTorch GPU batch',remaining_26000_steps='original PyTorch functions on separate CPU workers, preserved Adam state'),
                  default_batch512_scaling=dict(starts=8,full_schedule_steps=35000,best_max_cos=min(default['best'])),
                  eligible_for_released_polisher=0,save_threshold=.501)
    assert controls['paper_batch64']['best_max_cos']>.501 and min(default['best'])>.501
    save_json(out/'discovery_controls.json',controls)
    native=[]
    for name in ['trades842','trades_extended842']:
        row=json.loads((ROOT/'runs'/name/'summary.json').read_text());row['name']=name;native.append(row)
    result=dict(lower_bound=841,new_842_configuration=False,
                new_near_miss=metrics(np.load(out/'improved842.npy')),
                exact_upper_cosine='500900475/1000000000',
                published_near_miss_exclusion_radius='6121/250000',
                improved_near_miss_exclusion_radius='24409/1000000',
                batches=batches,discovery_controls=controls,native=native,
                total_native_moves=sum(v['moves'] for v in native),
                trials_842=sum(g['trials'] for b in batches for g in b['groups'] if g['n']==842),
                trials_841=sum(g['trials'] for b in batches for g in b['groups'] if g['n']==841)+72,
                completed=True)
    save_json(ROOT/'campaign.json',result)
    versions={name:importlib.metadata.version(name) for name in ['numpy','scipy','numba','torch','cupy-cuda12x','threadpoolctl']}
    save_json(ROOT/'environment.json',dict(packages=versions,reserved_cpu=127,gpu='NVIDIA RTX A4000',
                                         gpu_index=0,maximum_configured_cpu_worker_budget=121))
    # Record current replay sources and the unchanged shared implementations.
    files=sorted(ROOT.glob('*.py'))+sorted(ROOT.glob('*.cpp'))+sorted((ROOT/'snapshots').glob('*.py'))
    deps=[OLD/'search.py',OLD/'kernels.py',OLD/'continuation.py',OLD/'verify.py',
          AUDIT/'audit.py',AUDIT/'alternative_bridges.py',AUDIT/'riesz_search.py',AUDIT/'inputs'/'search_841_riesz.py']
    save_json(ROOT/'source_manifest.json',dict(sources={str(p.relative_to(ROOT)):sha(p) for p in files},
                                              dependencies={str(p.relative_to(ROOT.parent)):sha(p) for p in deps}))
    print(json.dumps(dict(trials_842=result['trials_842'],trials_841=result['trials_841'],
                         native_moves=result['total_native_moves'],controls=controls),indent=2))
