#!/usr/bin/env python3
"""Recheck all saved candidates and aggregate this campaign without promoting numerical data."""
import hashlib,json
from pathlib import Path
import numpy as np
from baselines import DATA
from prepare19 import le0


def main():
    prep19=json.loads((DATA/'preparation19.json').read_text())
    prep21=json.loads((DATA/'independent21.json').read_text())
    assert prep19['baseline']==12460 and prep19['violations']==0
    assert prep21['count']==31560 and prep21['violations']==0
    # Corrupting a ray list by duplicating any of its three ray types must fail.
    for q in (8,4889129,162):
        assert not bool((le0([q],[0])|le0([3*q*q],[0]))[0])
    labels_checks=[];z=np.load(DATA/'component21.npz');N=z['neighbors']
    for path in sorted(DATA.glob('cp21*_labels.npy')):
        labels=np.load(path);assert labels.shape==(1024,)
        assert np.all((labels>=-1)&(labels<21))
        for u,p in enumerate(labels):
            if p<0:continue
            assert labels[N[u,p]]<0
            assert all(labels[v]!=p for v in N[u] if labels[v]>=0)
        total=27720+4*int(np.sum(labels>=0))
        labels_checks.append(dict(file=path.name,total_if_replicated=total))
        if total>31560:raise RuntimeError('A discrete improvement needs full expansion and exact all-pairs certification.')
    numerical=[]
    for path in sorted(DATA.glob('*_candidate.npy')):
        if not path.name.startswith(('focused','polish')):continue
        X=np.load(path);dimension=X.shape[1];fixed=10668 if dimension==19 else 27720
        assert len(X)==(12461 if dimension==19 else 31561)
        assert np.isfinite(X).all()
        norm_error=float(np.abs((X*X).sum(1)-1).max());assert norm_error<1e-12
        worst=.5;base_violations=0;aux_violations=0
        for a in range(fixed,len(X),256):
            G=X[a:a+256]@X.T
            for t in range(len(G)):G[t,a+t]=-1.
            worst=max(worst,float(G.max()))
            base_violations+=int(np.sum(G[:,:fixed]>.5+1e-10))
            aux_violations+=int(np.sum(G[:,fixed:]>.5+1e-10))
        numerical.append(dict(file=path.name,dimension=dimension,count=len(X),max_dot=worst,
                              norm_error=norm_error,violating_pairs=base_violations+aux_violations//2,
                              valid_improvement=False))
        if worst<=.5+1e-9:raise RuntimeError('Potential numerical improvement requires an exact certificate.')
    results=[]
    for path in sorted(DATA.glob('*_result.json')):
        results.append(dict(file=path.name,**json.loads(path.read_text())))
    pending=[name for name in ('polish21','focused19','focused21') if not (DATA/(name+'_result.json')).exists()]
    summary=dict(date='2026-11-04',outcome='No improved lower bound found.' if not pending else 'Search incomplete.',
                 pending_runs=pending,best_verified_lower_bounds={'19':12460,'21':31560},
                 discrete_replacement_attempts=sum(r.get('iterations',0) for r in results),
                 gpu_multistarts=sum(r.get('restarts',0) for r in results if 'fixed_points' in r),
                 numerical_candidate_checks=numerical,line_label_checks=labels_checks,
                 baseline19=prep19,baseline21=prep21,completed_runs=results,
                 limitations=['No unrestricted upper bound is claimed.',
                              'Symmetry-restricted infeasibility says nothing about asymmetric labelings.',
                              'The 19D finite pool samples rational polar vertices; it is not the full continuous admissible region.',
                              'All numerical candidates remain invalid and are not lower-bound certificates.'])
    (DATA/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    sources=Path(__file__).resolve().parent
    manifest=dict(date='2026-11-04',primary_source='https://hexagonmath.org/pdf/2610.00248v1',
                  source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(sources.iterdir()) if p.suffix in ('.py','.cpp')},
                  input_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (DATA/'baseline19.npz',DATA/'baseline21.npz',DATA/'polar19_rays.npy')},
                  summary_sha256=hashlib.sha256((DATA/'summary.json').read_bytes()).hexdigest())
    (DATA/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ('completed_runs','baseline19','baseline21')}),flush=True)


if __name__=='__main__':main()
