#!/usr/bin/env python3
"""Insertion directions for subsequent collective deformation, not a maximality proof."""
import argparse,json,warnings,time
import numpy as np
from scipy.optimize import linprog
from continuous import baseline
from orbits17 import DATA

p=argparse.ArgumentParser();p.add_argument('--dimension',type=int,required=True)
p.add_argument('--restarts',type=int,default=80);a=p.parse_args()
X=baseline(a.dimension);rng=np.random.default_rng(10000+a.dimension);holes=[];start=time.monotonic()
for restart in range(a.restarts):
    v=rng.normal(size=a.dimension);v/=np.linalg.norm(v)
    for iteration in range(16):
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
holes.sort(key=lambda h:h[0]);out=np.asarray([h[1] for h in holes])
np.save(DATA/f'holes{a.dimension}.npy',out)
print(json.dumps(dict(dimension=a.dimension,directions=len(out),best_max_dot=holes[0][0],seconds=time.monotonic()-start)),flush=True)
