"""Finite multi-start search for unrestricted-coordinate normal-sector holes.

For points of squared norm32 the inequalities are A x <=16. A point of this
polytope with norm>=sqrt(32) can be radially contracted to an insertion.
Maximizing a sequence of supporting linear forms is only a heuristic for the
nonconvex norm maximum. It neither proves maximality nor certifies a bound.
"""
import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import linprog


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input',type=Path,default=Path(__file__).with_name('normal26.npz'))
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--starts',type=int,default=512)
    ap.add_argument('--seed',type=int,default=20261113)
    ap.add_argument('--remove',type=int,nargs='*',default=[])
    a=ap.parse_args()
    z=np.load(a.input); y=z['Y'];b=z['B']; keep=np.ones(len(y),bool);keep[a.remove]=False
    A=np.vstack([b,y[keep]]);d=A.shape[1]
    rng=np.random.default_rng(a.seed);rows=[];best=0.;start=time.monotonic()
    for i in range(a.starts):
        c=rng.normal(size=d);last=0.;p=None
        for it in range(30):
            r=linprog(-c,A_ub=A,b_ub=np.full(len(A),16.),bounds=[(None,None)]*d,
                      method='highs',options={'threads':2})
            if not r.success:raise RuntimeError(r.message)
            p=r.x;norm2=float(p@p)
            if norm2-last<1e-10:break
            last=norm2;c=p
        row={'start':i,'norm2':norm2,'steps':it+1,'seconds':time.monotonic()-start}
        rows.append(row)
        if norm2>best+1e-10:
            best=norm2
            active=np.flatnonzero(np.abs(A@p-16)<1e-7)
            np.savez_compressed(a.out.with_suffix('.npz'),point=p,active=active,
                               constraints=A[active],indices=z['indices'])
            print(json.dumps(dict(row,best=True)),flush=True)
        a.out.write_text(json.dumps({'scope':'Finite numerical vertex search; no global upper bound',
                         'source':str(a.input),'seed':a.seed,'removed':a.remove,'best_norm2':best,
                         'required_norm2':32,'starts':rows},indent=2)+'\n')
        if norm2>32+1e-9:
            print('POTENTIAL INSERTION NEEDS EXACT REPLAY',flush=True);break


if __name__=='__main__':main()
