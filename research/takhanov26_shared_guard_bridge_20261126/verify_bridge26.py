#!/usr/bin/env python3
"""Takhanov v0.2, 26D: exact seven-coordinate shared-guard replay (no new bound)."""
from pathlib import Path
from itertools import combinations, product, permutations
import argparse
import hashlib
import json
import numpy as np

PUBLISHED=201212
EXPECTED={"d26_blocks.json":"5018425664146ee2ace418021537b8b13863fc57c7dc64376ab6f9e811aa5162","d26_certificate.npz":"35fd8143cf0c7e6b0146719485e00c86952cbd5e404c7881c5f6b92fa25b21ec"}
A=[0,1,2,3,4,8,9]
D=[4,5,6,7,0,8,9]
B=[0,3,5,7,8]
C=[1,2,4,6,8]
CUT=(0,4,8,9)
MAD=(0,1,2,3,6)
keep={4,10,14,16}
even=[s for s in product((1,-1),repeat=3) if s.count(-1)%2==0]

def dot(x,y):
    return sum(int(a)*int(b) for a,b in zip(x,y))

def embed(v,S,n=10):
    w=[0]*n
    for k,j in enumerate(S):
        w[j]=int(v[k])
    return tuple(w)

def leech_shell():
    g=[sum(1<<(j+i)%23 for j in (0,2,4,5,6,10,11)) | (1<<23) for i in range(12)]
    words=[0]
    for x in g: words += [w^x for w in words]
    assert len(set(words))==4096
    assert sorted(set(w.bit_count() for w in words))==[0,8,12,16,24]
    z=np.zeros((196560,24), dtype=np.int16)
    i=0
    for a,b in combinations(range(24),2):
        for h,k in product((-8,8),repeat=2):
            z[i,a]=h; z[i,b]=k; i+=1
    signs=np.array([s for s in product((-4,4),repeat=8) if s.count(-4)%2==0])
    for w in words:
        if w.bit_count()!=8: continue
        js=[i for i in range(24) if w>>i&1]
        z[i:i+128,js]=signs
        i+=128
    for w in words:
        pat=np.array([2 if w>>j&1 else -2 for j in range(24)])
        z[i:i+24,:]=pat
        z[i+np.arange(24),np.arange(24)]*=-3
        i+=24
    assert i==196560 and ((z.astype('i4')**2).sum(1)==128).all()
    return z[np.lexsort(z.T[::-1])]

def exterior(cert,d):
    x=leech_shell();K=d["kernel_coordinate_indices"];J=[i for i in range(24) if i not in K]
    with np.load(cert/"d26_certificate.npz",allow_pickle=False) as f:
        ix=np.ones(len(x),dtype=bool);ix[f["deleted_leech_indices"]]=False
        v=f["added_integer"].astype(np.int64)
    lx=np.pad(x[ix],((0,0),(0,2)))
    lx=lx[np.any(lx[:,J],axis=1)]
    v=v[np.any(v[:,J],axis=1)]
    assert (len(lx),len(v))==(191904,8960)
    f=np.array(d["frame_rows"],dtype=np.int64)
    assert np.array_equal(f@f.T,128*np.eye(10,dtype=np.int64))
    projected=np.concatenate([lx,v])[:,K].astype(np.int64)@f.T
    assert np.all(projected%32==0)
    P=np.unique(projected//32,axis=0)
    assert P.shape==(2601,10)
    return P

def all_dot_legal(x,y):
    w=dot(x,y)
    return w<=0 or 4*w*w<=dot(x,x)*dot(y,y)

def orient(r,order,mask):
    return tuple((-1 if mask>>i&1 else 1)*r[j] for i,j in enumerate(order))

def checking(cert, path=None):
    ddir=cert/"data"
    for nm,h in EXPECTED.items():
        assert hashlib.sha256((ddir/nm).read_bytes()).hexdigest()==h, nm
    d=json.loads((ddir/"d26_blocks.json").read_text());P=exterior(ddir,d)
    T=d["template_numerators_and_denominators"]
    for r in T:
        assert dot(r[:5],r[:5])==r[5]**2 and 2*max(map(abs,r[:5]))<=r[5]
    # Signed four-support rigidity with all four original chambers retained.
    old_boxes=d["five_coordinate_supports"]
    for k in (3,4):
        scores=[sum(sum(s*r[j] for j,s in zip(S,sg))>r[5] for r in T)
                for S in combinations(range(5),k) for sg in product((1,-1),repeat=k)]
        assert (min(scores),max(scores))==((2,4) if k==3 else (2,8))
    safe=[]
    for S in combinations(range(10),4):
        V=np.array(list(product((1,-1),repeat=4)),dtype=np.int64)
        count=int(np.sum(np.max(P[:,S]@V.T,axis=0)<=4))
        assert count in (0,16)
        if count==16: safe.append(S)
    same_skeleton=[S for S in safe if all(len(set(S)&set(H))<=2 for H in old_boxes)]
    assert len(safe)==70 and sorted(same_skeleton)==sorted(map(tuple,d["four_coordinate_supports"]))
    sk=[]
    for S in d["four_coordinate_supports"]:
        for sig in product((1,-1),repeat=4):
            if tuple(S)==CUT and sig==(1,1,1,1): continue
            sk.append(embed(sig,S))
    assert len(sk)==223
    aa=[embed((s,),[i]) for i,s in product(range(10),(1,-1))]
    old=( [embed(r[:5],B) for r in T] + [embed(r[:5],C) for r in T])
    if path is None:
        vv= [orient(r,[1,0,2,3,4],17) for i,r in enumerate(T) if i not in keep]
        xx=[q[:4]+(0,0,q[4]) for q in vv] + [(5,4*s[0],4*s[1],4*s[2],1,1,5) for s in even]
        sqs=[dot(r,r) for r in xx]
        for q in xx[-4:]:
            assert dot(q,q)==100
            assert sum(q[j] for j in (0,4,5,6))==12
        # Verify the reported four and only four removed rows for this template orientation.
        hits={i for i,r in enumerate(T) if any(not all_dot_legal(q,orient(r,[1,0,2,3,4],17)[:4]+(0,0,orient(r,[1,0,2,3,4],17)[4])) for q in xx[-4:])}
        assert hits==keep
        # Complete stock-orientation audit for this particular bridge quartet.
        toll=26
        for o in permutations(range(5)):
            for m in range(32):
                k=0
                for r in T:
                    z=orient(r,o,m); zz=z[:4]+(0,0,z[4])
                    if any(not all_dot_legal(q,zz) for q in xx[-4:]): k+=1
                    if k>=toll: break
                toll=min(toll,k)
        assert toll==4
    else:
        raw=np.loadtxt(path,skiprows=1)
        assert raw.shape==(27,7)
        # Candidate rays, with exact normalization; do not compare floats to the geometric limit.
        xx=[tuple(int(round(y*10**8)) for y in r) for r in raw]
        sqs=[dot(r,r) for r in xx]
    for q,n in zip(xx,sqs):
        assert n>0 and 4*max(v*v for v in q)<=n
        assert 25*max(q[4]**2,q[5]**2)<=n
    ray=old+aa+sk+[embed(q,S) for S in (A,D) for q in xx]
    for i,x in enumerate(ray):
        l=dot(x,x);assert l
        gp=P@np.array(x,dtype=np.int64)
        assert max(map(abs,gp),default=0)<2**31
        mx=int(gp.max())
        assert mx<=0 or mx*mx<=4*l, ("exterior",i)
        assert all(all_dot_legal(x,y) for y in ray[:i]), ("kernel pair",i)
    N=200864+len(ray)
    assert N==(201211 if path is None else 201213)
    print(json.dumps({"same_source_exterior":200864,"signed_axis":20,"signed_tetrad":223,"old_chambers":52,"mirror_size_each":len(xx),"total":N,"new_bound":N>PUBLISHED,"exact":True}))

def main():
    if not __debug__: raise RuntimeError('Run without python -O')
    ap=argparse.ArgumentParser();ap.add_argument("certificates",type=Path)
    ap.add_argument("--candidate",type=Path)
    a=ap.parse_args();checking(a.certificates,a.candidate)
if __name__=="__main__": main()
