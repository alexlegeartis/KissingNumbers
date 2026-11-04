#!/usr/bin/env python3
"""LP-ascent insertion directions over the current, independently rebuilt records."""
import argparse,json,time,warnings
import numpy as np
from scipy.optimize import linprog
from baselines import baseline,DATA

p=argparse.ArgumentParser();p.add_argument('--dimension',type=int,choices=(19,21),required=True)
p.add_argument('--restarts',type=int,default=100);p.add_argument('--structured',action='store_true');a=p.parse_args()
X=baseline(a.dimension);rng=np.random.default_rng(1192100+a.dimension);holes=[];start=time.monotonic()
starts=[]
if a.structured and a.dimension==21:
    z=np.load(DATA/'baseline21.npz');missing=np.flatnonzero(z['labels']<0)
    for _ in range(a.restarts):
        w=int(z['words'][rng.choice(missing)]);profile=z['profiles'][rng.integers(21)]
        starts.append(profile*np.array([1-2*((w>>i)&1) for i in range(21)]))
elif a.structured and a.dimension==19:
    z=np.load(DATA/'baseline19.npz')
    for _ in range(a.restarts):
        v=z['axis'][rng.integers(len(z['axis']))].astype(float);v[-3:]*=-1;starts.append(v)
for restart in range(a.restarts):
    v=starts[restart].astype(float) if starts else rng.normal(size=a.dimension);v/=np.linalg.norm(v)
    for iteration in range(18):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            r=linprog(-v,A_ub=X,b_ub=np.full(len(X),.5),bounds=[(-2,2)]*a.dimension,
                      method='highs',options={'threads':1})
        if not r.success:break
        d=r.x/np.linalg.norm(r.x)
        if np.linalg.norm(d-v)<1e-8:break
        v=d
    if r.success and float((X@r.x).max())<=.5+1e-9:
        holes.append((float((X@d).max()),d.copy()))
    if restart%10==9:
        print(json.dumps(dict(restarts=restart+1,best_max_dot=min(h[0] for h in holes),seconds=time.monotonic()-start)),flush=True)
    if holes and holes[-1][0]<=.5-1e-10:
        print('POTENTIAL FREE POINT; requires exact certification',flush=True);break
holes.sort(key=lambda h:h[0]);out=np.asarray([h[1] for h in holes]);suffix='_structured' if a.structured else ''
np.save(DATA/f'holes{a.dimension}{suffix}.npy',out)
report=dict(dimension=a.dimension,directions=len(out),best_max_dot=holes[0][0],seconds=time.monotonic()-start,
            scope='Numerical search for insertion directions, not a maximality proof.')
(DATA/f'holes{a.dimension}{suffix}_result.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report),flush=True)
