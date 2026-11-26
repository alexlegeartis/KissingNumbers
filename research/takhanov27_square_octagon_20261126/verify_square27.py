#!/usr/bin/env python3
"""Exact square-carrier gates for Takhanov's published 210108-point 27D code.
Usage: OPENBLAS_NUM_THREADS=1 python verify_square27.py /path/to/d27_certificate.npz
No optimization or tolerance; all Gram products of these short integer vectors are exact.
"""
import sys, hashlib, itertools
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np

def verify(path):
    assert hashlib.sha256(path.read_bytes()).hexdigest() == "8b704bb49503e3a4f6ad6b567690a7a2a9abfb30600f453313eb0510e6d29a94"
    with np.load(path, allow_pickle=False) as f:
        X = np.asarray(f["added_integer"], dtype=np.int32)
        mids = f["midpoint_rows"]
        E = f["midpoint_endpoints"]
    assert X.shape == (26796,27)
    # Pin either signed row to its actual source certificate row.
    ids = defaultdict(list)
    for row,(r,s) in zip(mids,E):
        ids[tuple(sorted((int(r),int(s))))].append(int(row))
    edges=list(ids)
    endpoints=np.array(edges)
    erows=np.array([ids[k] for k in edges])
    assert erows.shape==(13248,2)
    assert np.array_equal(X[erows[:,0],:24], X[erows[:,1],:24])
    assert np.array_equal(X[erows[:,0],24:], -X[erows[:,1],24:])
    inc=defaultdict(list)
    for i,(r,s) in enumerate(edges):
        inc[r].append(i); inc[s].append(i)
    assert len(inc)==13248 and {len(v) for v in inc.values()}=={2}
    visited=set(); squares=[]; lens=Counter()
    for r in sorted(inc):
        if r in visited: continue
        v=r; last=-1; cycle=[]
        while v not in visited:
            visited.add(v)
            e=next(x for x in inc[v] if x!=last)
            cycle.append(e)
            v=sum(edges[e])-v;last=e
        assert v==r
        lens[len(cycle)]+=1
        if len(cycle)==4:squares.append(cycle)
    assert lens=={4:3248,8:32}
    squares=np.array(squares)
    M=X[erows[squares,0],:24]
    assert np.all((M[:,0]+M[:,2])%2==0)
    O=(M[:,0]+M[:,2])//2;A=M[:,0]+M[:,3]-2*O;B=M[:,0]+M[:,1]-2*O
    assert np.array_equal(O*2,M[:,1]+M[:,3])
    for x,y in [(O,A),(O,B),(A,B)]:
        assert np.all(np.sum(x*y,axis=1)==0)
    for x in (O,A,B):assert np.all(np.sum(x*x,axis=1)==64)
    assert np.unique(O,axis=0).shape==O.shape
    tsq=np.sum(X[:,24:]**2,axis=1)
    poles=X[tsq==128]
    assert len(poles)==12 and np.all(poles[:,:24]==0)
    assert set(map(tuple,poles[:,24:]))=={tuple(0 if k==i else d[(k-i)%3-1] for k in range(3)) for i in range(3) for d in itertools.product([-8,8],repeat=2)}
    for i in range(3):
        for s in [-8,8]:
            h=X[(tsq==64)&(X[:,24+i]==s),:24]; hh=h@h.T
            assert hh.shape==(48,48) and np.all(np.diag(hh)==64)
            assert np.all(np.sum(hh==64,axis=1)==1) and np.all(np.sum(hh==-64,axis=1)==1) and np.all(np.isin(hh,[-64,0,64]))
    # Guarantee all the guards really occur AS PAIRED rows.
    hids=np.where(tsq==64)[0]; byh={tuple(X[i]):i for i in hids}
    hpair=[]
    for i in hids:
        t=X[i,24:]
        if t[np.flatnonzero(t)[0]]>0:
            hpair.append([i,byh[tuple(np.r_[X[i,:24],-t])]])
    pairs=np.concatenate((erows,np.asarray(hpair)),axis=0)
    assert len(pairs)==13248+144
    source=X[pairs[:,0]]
    t2=(source[:,24:]**2).sum(1)
    slopes=np.r_[np.ones(len(erows),int)*32, np.ones(len(hpair),int)*64]
    assert np.all(slopes==8*np.max(np.abs(source[:,24:]),axis=1))
    axis_miss=np.argmin(np.abs(X[erows[squares[:,0],0],24:]),axis=1)
    guard=[None]*len(squares); tilted=[None]*len(squares); special=[]; midsF=source[:,:24].T.astype(np.float32)
    for start in range(0,len(squares),32):
        sl=slice(start,start+32)
        proj=[(q[sl].astype(np.float32)@midsF).astype(np.int16) for q in [O,A,B]]
        for off,(o,pa,pb) in enumerate(zip(*proj)):
            i=start+off
            ix=np.where(o>=40)[0] # all claimed paired witnesses are in this cutoff
            h=64-o[ix].astype(np.int64);p=pa[ix].astype(np.int64);q=pb[ix].astype(np.int64)
            b2=p*p+q*q;d=64*t2[ix]-5*h*h-b2
            direct=ix[(b2==0)&(d>0)]; tilt=ix[(d>0)&(d*d>16*h*h*b2)]
            assert not (set(squares[i])&set(tilt))
            guard[i]=tuple(sorted(map(int,direct)))
            tilted[i]=tuple(sorted(map(int,tilt)))
            if not len(tilt):special.append(i)
    hist=Counter(map(len,guard))
    assert hist=={0:224,1:16,2:256,4:1408,6:512,8:384,10:256,12:128,16:64}
    assert len(special)==160
    quad=defaultdict(list)
    for i,gg in enumerate(guard):
        if len(gg)==2: quad[gg].append(i)
    assert len(quad)==64 and set(map(len,quad.values()))=={4}
    assert all(len({edge for row in S for edge in squares[row]})==16 for S in quad.values())
    assert len({g[0] for g in guard if len(g)==1})==16
    assert all(g[0]>=len(erows) for g in guard if len(g)==1)
    assert Counter(Counter(k for gg in guard for k in gg if k<len(erows)).values())=={34:512}
    quadruple={i:(g,S) for g,S in quad.items() for i in S}
    pX=X[:,:24].T.astype(np.float32)
    status=Counter();N=len(X)
    for i in special+list(quadruple):
        exposed=i in quadruple
        if exposed:
            g,S=quadruple[i]
            drop=set(erows[squares[S]].ravel())|set(pairs[list(g)].ravel())
        else: drop=set(erows[squares[i]].ravel())
        o,pa,pb=[(v[i].astype(np.float32)@pX).astype(int) for v in [O,A,B]]
        eligible=np.ones(N,dtype=bool);eligible[list(drop)]=False
        m=int(axis_miss[i]); other=[j for j in range(3) if j!=m]
        gates=Counter(); caps=Counter()
        for k in np.where(eligible & (o>=32) & ((tsq==32)|(tsq==64)))[0]:
            z=X[k,24:]
            if tsq[k]==64 and pa[k]==pb[k]==0:
                j=int(np.nonzero(z)[0][0]);caps[(j,int(np.sign(z[j])))]+=1
            if tsq[k]==32 and z[m]!=0:
                j=next(j for j in other if z[j]!=0)
                for y,(p,q) in enumerate(((pa[k],pb[k]),(pb[k],pa[k]))):
                    if abs(p)==32 and q==0:
                        gates[(y,int(np.sign(p)),int(np.sign(z[m])),j,int(np.sign(z[j])))]+=1
        for y,ss,mm,g0,g1 in itertools.product(range(2),[-1,1],[-1,1],[-1,1],[-1,1]):
            a0=gates[y,ss,mm,other[0],g0]; a1=gates[y,ss,mm,other[1],g1]
            if exposed:
                assert min(a0,a1,caps[other[0],g0],caps[other[1],g1])>=4
            else:
                assert (a0 and a1) or (a0 and caps[other[1],g1]) or (a1 and caps[other[0],g0])
        status['four-square carriers' if exposed else 'last fixed carriers']+=1
    print('source=210108; cycles=',dict(lens),'; pure guard pairs per square=',dict(sorted(hist.items())))
    print('tilted extra=',224-len(special),'; fixed final=',len(special),'; independent four-square certificates=',len(quad),'; slots checked=',dict(status))
if __name__=='__main__':
    verify(Path(sys.argv[1]))

