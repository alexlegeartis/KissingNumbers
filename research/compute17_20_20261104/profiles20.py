#!/usr/bin/env python3
"""Search nonuniform dense magnitudes over the complete fixed 17400-point base.

All decisions used to admit profiles and construct conflict graphs use integers.
A magnitude vector p gives norm-eight points sqrt(8/Q)*p*(-1)^c, Q=sum p_i^2.
"""
import argparse
import functools
import itertools
import json
import math
import multiprocessing as mp
from pathlib import Path
import sys
from fractions import Fraction
import warnings
import numpy as np
from scipy.optimize import linprog
from orbits17 import DATA,write_graph
from prepare20 import basis_of,reduction
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'common'))
from golay import golay24

_constraints = None


def constraints():
    G=golay24()
    blocks=G[(G[:,:4].sum(1)==0)&(G.sum(1)==8),4:]
    positions=[np.flatnonzero(b) for b in blocks]
    A=[];b=[]
    for i,j in itertools.combinations(range(20),2):
        row=np.zeros(20);row[i]=row[j]=1;A.append(row);b.append(2)
    for block in positions:
        for i in block:
            row=np.zeros(20);row[block]=1;row[i]=-1;A.append(row);b.append(4)
    return np.asarray(A),np.asarray(b),positions


def sample(task):
    global _constraints
    if _constraints is None:_constraints=constraints()
    A,b,positions=_constraints
    seed,sigma=task;rng=np.random.default_rng(seed)
    direction=np.ones(20)+sigma*rng.normal(size=20)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        result=linprog(-direction,A_ub=A,b_ub=b,bounds=(0,2),method='highs',options={'threads':1})
    if not result.success:return None
    fractions=[Fraction(float(x)).limit_denominator(10000) for x in result.x]
    # The uniform polytope vertex is deeper than the norm-eight sphere. Edges
    # towards shallower vertices still supply nonuniform feasible directions.
    v=np.asarray([float(x) for x in fractions]);u=np.full(20,2./3.)
    delta=v-u;aa=delta@delta;bb=2.*(u@delta);cc=u@u-8.
    if aa>1e-12:
        disc=bb*bb-4*aa*cc
        if disc>0:
            roots=[t for t in ((-bb-np.sqrt(disc))/(2*aa),(-bb+np.sqrt(disc))/(2*aa)) if 0<t<1]
            if roots:
                t=Fraction(max(0,int(min(roots)*4096)-1),4096)
                fractions=[(1-t)*Fraction(2,3)+t*x for x in fractions]
    denominator=math.lcm(*(x.denominator for x in fractions))
    if denominator>1000000:return None
    p=[x.numerator*(denominator//x.denominator) for x in fractions]
    gcd=math.gcd(*p);p=[x//gcd for x in p]
    Q=sum(x*x for x in p)
    if max(p)>1000000:return None
    if any(2*(p[i]+p[j])**2>Q for i,j in itertools.combinations(range(20),2)):return None
    if any((sum(p[i] for i in block)-2*min(p[i] for i in block))**2>2*Q for block in positions):return None
    return tuple(p)


def generate(trials,workers):
    tasks=[(100000+i, [.1,.3,.5,.8,1.2,2.,4.][i%7]) for i in range(trials)]
    profiles={tuple([1]*20)};found=0
    with mp.Pool(workers) as pool:
        for p in pool.imap_unordered(sample,tasks,chunksize=4):
            found+=1
            if p is not None:profiles.add(p)
            if found%200==0:print(json.dumps(dict(event='profile_sampling',trials=found,profiles=len(profiles))),flush=True)
    profiles=sorted(profiles,key=lambda p:(len(set(p)),sum(p),p))
    (DATA/'profiles20.json').write_text(json.dumps(dict(trials=trials,profiles=profiles),indent=2)+'\n')
    print(json.dumps(dict(event='profile_sampling_complete',trials=trials,profiles=len(profiles))),flush=True)


def make(rank,seed,limit):
    saved=json.loads((DATA/'profiles20.json').read_text())
    profiles=saved['profiles']
    uniform=[1]*20
    profiles=[uniform]+[p for p in profiles if p!=uniform]
    if limit and len(profiles)>limit:
        rng=np.random.default_rng(seed)
        G0=golay24()[:,4:].astype(np.int64)
        low=G0[G0.sum(1)==4]
        interesting=[i for i,p in enumerate(profiles) if i and
                     any(4*sum(int(p[j])**2 for j in np.flatnonzero(w))>=sum(int(x)**2 for x in p) for w in low)]
        chosen=interesting[:limit-1]
        remaining=[i for i in range(1,len(profiles)) if i not in chosen]
        chosen+=list(map(int,rng.choice(remaining,size=limit-1-len(chosen),replace=False)))
        profiles=[uniform]+[profiles[i] for i in sorted(chosen)]
    P=np.asarray(profiles,dtype=np.int64);Q=(P*P).sum(1)
    G=golay24();signs=1-2*G[:,4:].astype(np.int64)
    shortened_ids=np.flatnonzero(G[:,:4].sum(1)==0)
    rng=np.random.default_rng(seed)
    hbase=basis_of(rng.permutation(shortened_ids))[:rank]
    H=[0]
    for b in hbase:H += [x^b for x in H]
    _,coordinate,nonpivots=reduction(hbase,12)
    nc=2**(12-rank);ids=np.arange(nc,dtype=np.uint32)
    reps=np.zeros(nc,dtype=np.uint32)
    for i,p in enumerate(nonpivots):reps|=((ids>>i)&1)<<p
    rows=[set() for _ in range(len(P)*nc)]
    allowed=[]
    for a in range(len(P)):
        for b in range(a,len(P)):
            products=P[a]*P[b]
            dots=signs@products
            # Long double separates distant values; all possible ties are resolved
            # with arbitrary-precision Python integer squares.
            lhs=4*dots.astype(np.longdouble)**2
            rhs=np.longdouble(Q[a])*np.longdouble(Q[b])
            bad=(dots>0)&(lhs>rhs)
            near=np.flatnonzero(np.abs(lhs-rhs)<1e-15*max(1,rhs))
            for j in near:
                bad[j]=int(dots[j])>0 and 4*int(dots[j])**2>int(Q[a])*int(Q[b])
            if a==b:
                bad[0]=False
                if bad[H].any():raise RuntimeError('subgroup orbit has internal conflicts')
            connection=np.unique(coordinate(np.flatnonzero(bad)))
            for i in range(nc):
                u=a*nc+i
                for c in connection:
                    v=b*nc+(i^int(c))
                    if u!=v:rows[u].add(v);rows[v].add(u)
    rows=[np.asarray(sorted(row),dtype=np.uint32) for row in rows]
    name=f'profiles20_r{rank}_s{seed}_p{len(P)}'
    edges=write_graph(DATA/(name+'.graph'),rows)
    selected=np.unique(coordinate(np.flatnonzero(G[:,0]==0)))
    assert len(selected)*2**rank==2048
    chosen=np.zeros(len(rows),dtype=bool);chosen[selected]=True
    assert all(not chosen[rows[i]].any() for i in selected)
    np.savez_compressed(DATA/(name+'.npz'),profiles=P,norms=Q,weights=np.full(len(rows),2**rank),
                        seed=selected,subgroup=H,word_representatives=reps,
                        fixed=17400,dimension=20,profile_model=True)
    report=dict(name=name,vertices=len(rows),profiles=len(P),directed_edges=edges,
                fixed_base=17400,seed_total=19448,group_rank=rank)
    (DATA/(name+'.json')).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--trials',type=int,default=1400)
    p.add_argument('--workers',type=int,default=8);p.add_argument('--build',action='store_true')
    p.add_argument('--rank',type=int,default=4);p.add_argument('--seed',type=int,default=1)
    p.add_argument('--limit',type=int,default=32);a=p.parse_args()
    if a.build:make(a.rank,a.seed,a.limit)
    else:generate(a.trials,a.workers)
