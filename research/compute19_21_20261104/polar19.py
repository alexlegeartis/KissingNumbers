#!/usr/bin/env python3
"""Explore unrestricted rational directions in the polar of the fixed 19D base."""
import argparse,json,math,multiprocessing as mp,sys,time,warnings
from fractions import Fraction
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from baselines import DATA
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'compute17_20_20261104'))
from orbits17 import write_graph

_base=None


def sample(seed):
    global _base
    if _base is None:_base=np.load(DATA/'baseline19.npz')['base'].astype(float)
    rng=np.random.default_rng(seed);v=rng.normal(size=19);v/=np.linalg.norm(v)
    for _ in range(8):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            r=linprog(-v,A_ub=_base,b_ub=np.full(len(_base),4.),bounds=[(-5,5)]*19,
                      method='highs',options={'threads':1})
        if not r.success:return None
        norm=np.linalg.norm(r.x);direction=r.x/norm
        if np.linalg.norm(direction-v)<1e-8:break
        v=direction
    if norm**2<8-1e-8:return None
    fractions=[Fraction(float(x)).limit_denominator(1000) for x in r.x]
    den=math.lcm(*(x.denominator for x in fractions))
    if den>10000:return None
    p=[x.numerator*(den//x.denominator) for x in fractions]
    gcd=math.gcd(*p);p=[x//gcd for x in p]
    if max(map(abs,p))>1000:return None
    Q=sum(x*x for x in p)
    dots=np.rint(_base@np.asarray(p,dtype=float)).astype(np.int64)
    if max(0,int(dots.max()))**2>2*Q:return None
    return tuple(p)


def generate(trials,workers,limit):
    started=time.monotonic();unique=set();accepted=0
    with mp.Pool(workers) as pool:
        for i,p in enumerate(pool.imap_unordered(sample,range(190000,190000+trials),chunksize=4)):
            if p is not None:unique.add(p);accepted+=1
            if i%256==255:print(json.dumps(dict(event='polar_sampling',trials=i+1,unique=len(unique),seconds=time.monotonic()-started)),flush=True)
    z=np.load(DATA/'baseline19.npz');T=list(map(int,z['tetrad_positions']));H=list(map(int,z['head_positions']))
    tetrads=list(map(int,z['tetrads']));U=[0]
    for t in tetrads:U += [u^t for u in U]
    U0=[u for u in U if sum(bool(u&t) for t in tetrads)%2==0]
    transformations=[]
    for u in U0:
        flips=[-1 if u>>p&1 else 1 for p in H]+[-1 if u>>T[0]&1 else 1]*3
        transformations.append(tuple(flips))
    points=set()
    for p in sorted(unique):
        for s in transformations:points.add(tuple(a*b for a,b in zip(p,s)))
    rng=np.random.default_rng(191919);points=sorted(points)
    if len(points)>limit:points=[points[i] for i in sorted(rng.choice(len(points),limit,replace=False))]
    P=np.asarray(points,dtype=np.int64)
    np.save(DATA/'polar19_rays.npy',P)
    report=dict(trials=trials,accepted=accepted,unique_lp_rays=len(unique),group_order=len(U0),
                candidates=len(P),seconds=time.monotonic()-started,
                exact_base_admission=True)
    (DATA/'polar19_sampling.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)


def build(exact=False):
    z=np.load(DATA/'baseline19.npz');P=np.load(DATA/'polar19_rays.npy')
    A=z['axis'].astype(np.int64);D=z['diagR'].astype(np.int64);S=z['diagS'].astype(np.int64)
    old=np.concatenate([A/np.sqrt(4889129.),(D+np.sqrt(6.)*S)/np.sqrt(162.)])
    Q=(P*P).sum(1)
    X=np.concatenate([old,P/np.sqrt(Q)[:,None]])
    # A generous outward margin includes every true conflict. Old-old pairs are
    # restored from the already exhaustive exact verification of the baseline.
    rows=[];margin=1e-11
    if exact:
        from prepare19 import le0
        R=np.concatenate([A,D,P]);SS=np.concatenate([np.zeros_like(A),S,np.zeros_like(P)])
        norms=(R*R+6*SS*SS).sum(1)
        RF=R.astype(float);SF=SS.astype(float)
    for a in range(0,len(X),256):
        if exact:
            aa=RF[a:a+256]@RF.T+6*SF[a:a+256]@SF.T
            bb=RF[a:a+256]@SF.T+SF[a:a+256]@RF.T
            assert np.array_equal(aa,np.rint(aa)) and np.array_equal(bb,np.rint(bb))
            aa=aa.astype(np.int64);bb=bb.astype(np.int64)
            C=4*aa*aa+24*bb*bb-norms[a:a+256,None]*norms[None,:];Dfield=8*aa*bb
            invalid=~(le0(aa,bb)|le0(C,Dfield))
            for k,bad in enumerate(invalid):
                i=a+k
                if i<1792:assert not bad[:1792].any() or (bad[:1792].sum()==1 and bad[i])
                bad[i]=False;rows.append(np.flatnonzero(bad).astype(np.uint32))
        else:
            G=X[a:a+256]@X.T
            for k,row in enumerate(G):
                i=a+k;bad=row>.5-margin
                if i<1792:bad[:1792]=False
                bad[i]=False;rows.append(np.flatnonzero(bad).astype(np.uint32))
    name='polar19_exact' if exact else 'polar19'
    edges=write_graph(DATA/(name+'.graph'),rows)
    with (DATA/'polar19.search').open('w') as f:
        f.write(f'{len(X)} 10668 19\n'+' '.join(['1']*len(X))+'\n')
        f.write('1792\n'+' '.join(map(str,range(1792)))+'\n')
    np.savez_compressed(DATA/'polar19_pool.npz',rational_rays=P,units=X,seed=np.arange(1792))
    report=dict(vertices=len(X),new_rational_directions=len(P),directed_edges=edges,
                baseline=12460,scope='Fixed complete base, arbitrary sampled rational polar directions, and the published auxiliary layer.',
                conservative_pair_margin=None if exact else margin,exact_field_arithmetic=exact)
    (DATA/(name+'_graph.json')).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--trials',type=int,default=2048);p.add_argument('--workers',type=int,default=16)
    p.add_argument('--limit',type=int,default=12000);p.add_argument('--build',action='store_true');p.add_argument('--exact',action='store_true');a=p.parse_args()
    if a.build:build(a.exact)
    else:generate(a.trials,a.workers,a.limit)
