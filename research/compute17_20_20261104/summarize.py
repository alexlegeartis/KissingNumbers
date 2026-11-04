#!/usr/bin/env python3
"""Gather completed runs and independently rescan the saved GPU candidates."""
import hashlib
import json
import numpy as np
from orbits17 import DATA


def main():
    verification=json.loads((DATA/'independent_verification.json').read_text())
    certificate=json.loads((DATA/'profile_ceiling_certificate.json').read_text())
    assert certificate['status']=='PASS'
    assert hashlib.sha256((DATA/'profiles20.json').read_bytes()).hexdigest()==certificate['profile_file_sha256']
    gpu_checks=[]
    for dimension in (17,20):
        result=DATA/f'gpu{dimension}_result.json'
        if not result.exists():continue
        R=json.loads(result.read_text());X=np.load(DATA/f'gpu{dimension}_candidate.npy')
        assert X.shape==(R['count'],dimension) and np.isfinite(X).all()
        norm_error=float(np.abs((X*X).sum(1)-1).max())
        assert norm_error<1e-12
        worst=-1.;violations=0
        for a in range(0,len(X),256):
            G=X[a:a+256]@X.T
            for t in range(len(G)):G[t,a+t]=-1.
            worst=max(worst,float(G.max()))
            violations+=int(np.sum(G>.5+1e-10))
        assert abs(worst-R['max_dot'])<1e-11
        gpu_checks.append(dict(dimension=dimension,count=len(X),max_dot=worst,
                               norm_error=norm_error,violating_pairs=violations//2,
                               certified=False))
        if worst<=.5+1e-9:
            raise RuntimeError('Potential improvement found: exact certification is required before summarizing failure.')
    runs=[dict(file=p.name,**json.loads(p.read_text())) for p in sorted(DATA.glob('*_result.json'))]
    local=[r for r in runs if 'iterations' in r]
    pending=[f'gpu{d}' for d in (17,20) if not (DATA/f'gpu{d}_result.json').exists()]
    summary=dict(date='2026-11-04',outcome=('Search incomplete; no improvement found so far.' if pending else 'No improved lower bound found.'),
                 pending_runs=pending,
                 best_verified_lower_bounds={'17':verification['baseline17']['count'],
                                             '20':verification['baseline20']['count']},
                 local_search_attempts=sum(r['iterations'] for r in local),
                 local_search_attempts_by_dimension={str(d):sum(r['iterations'] for r in local if r['dimension']==d) for d in (17,20)},
                 gpu_verification=gpu_checks,
                 new_exact_finite_family_certificate={k:certificate[k] for k in ('profiles','sign_words','dense_candidate_points','fixed_base','exact_family_ceiling','scope')},
                 completed_runs=runs,
                 limitations=['The searches do not prove optimality of the unrestricted configurations.',
                              'The integer profile certificate applies only to the saved 267 profiles.',
                              'Orbit models exclude configurations breaking their imposed sign subgroup.',
                              'Numerical candidates are not certificates.'])
    (DATA/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='completed_runs'}),flush=True)


if __name__=='__main__':main()
