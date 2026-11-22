#!/usr/bin/env python3
"""Exact alternate 63/64D clone-owner construction, retaining the cyclic equator."""
from fractions import Fraction as Q
from pathlib import Path
from hashlib import sha1,sha256
from itertools import combinations
import importlib.util,json
import numpy as np
H=Path(__file__).resolve().parent
V=H.parent
ONE=(1<<64)-1

def pin(path,digest):
    p=V/path
    a=p.read_bytes()
    if len(digest)==40:
        d=sha1(f"blob {len(a)}\0".encode()+a).hexdigest()
    else:
        d=sha256(a).hexdigest()
    assert d==digest,p
    return p

def load(path,digest):
    return json.loads(pin(path,digest).read_text())

def bits(m):
    while m:
        z=m&-m
        yield z.bit_length()-1
        m^=z

def stock():
    pars=pin("dim61-67-compatible-sign-families/verify_family.py","222c034a630bfab0db9df1249797d0e633436bfc")
    spec=importlib.util.spec_from_file_location("old",pars)
    old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    D,I=old.rref(old.generator(),64)
    assert len(D)==len(I)==28
    raw=pin("dim61-68-95-96-code-layers/source/dim61_96/code64_supports.npz","cdbd5e5a7bbda9ddb84ee86052e2d4d00ccc331f0188de958a69ed83b0c529f6")
    O=set(map(int,np.load(raw,allow_pickle=False)["words"]))
    assert len(O)==30828 and all(w.bit_count()==16 for w in O)
    # All unaltered W16 signs are even. They may be reordered in the 63D
    # section, but their total parity is preserved by the coordinate order.
    assert all((w&r).bit_count()%2==0 for w in O for r in D)
    g=load("dim61-67-compatible-sign-families/data/global64.json","d33d84363248951c8a0b348d8fb8a9c8e667d081")
    j=load("dim63-65-affine-reconstruction/certificate.json","f78b62b28cfb1f28495b6f50de52ba48d43568b8")
    assert g["odd_total_parity"] is False
    R=g["old"]+j["modified_parent_classes"]+g["new"]+j["new_classes"]
    assert len(R)==len({r["support"] for r in R})==3846
    assert all(r["support"] in O for r in g["old"]+j["modified_parent_classes"])
    assert all(r["support"] not in O for r in g["new"]+j["new_classes"])
    bymark=[[] for i in range(64)]
    for row in R:
        J=list(bits(row["support"]))
        assert len(J)==16
        cs=row["checks"];t=row["shift"]
        h=[sum(1<<J[k] for k in bits(c)) for c in cs]
        b=[(t&c).bit_count()%2 for c in cs]
        for a in J:bymark[a].append(tuple(zip(h,b)))
    return D,I,bymark

def verify():
    assert __debug__
    D,I,R=stock()
    section=load("dim63-65-affine-reconstruction/section63.json","cdc01424d1789116ac83c5bb6509d88269bfc16b")
    Z=section["octad"]
    assert Z.bit_count()==8 and section["terms"]["poles"]==112
    a63=load("dim63-balanced-nine-refill/claims.json","e85ed9d87813b767a188823fbd144a6a7b355d1f")
    a64=load("dim64-65-partial-nonary-signs/claims.json","96a743ef12f41a0c144c7ca9c86cf2c52749f573")
    oldbounds={63:a63["dimensions"][0]["lower_bound"],64:a64["dimensions"][0]["lower_bound"]}
    assert a63["dimensions"][0]["dimension"]==63 and a64["dimensions"][0]["dimension"]==64
    # Spike coordinates are ±7 on the mark and ±1 elsewhere, then all
    # multiplied by b=2/sqrt(7). Squared norm is 64.
    b2=Q(4,7)
    assert b2*(49+63)==64
    # Against unmodified W4 and W9 sparse parts: 4(7+3)b and (8/3)(7+8)b.
    assert 40**2*b2<32**2 and 8**2*b2<32**2
    # W16 guard: all signs agreeing is the only possibility >32 for a
    # support with the mark. Otherwise one mismatch decreases raw dot by 4b.
    assert 44**2*b2>32**2 and 40**2*b2<32**2 and 32**2*b2<32**2
    # With any other cyclic dense word, distinct owners are >=16 apart
    # (inherited verified cyclic-code minimum). On a matching heavy sign
    # one disagreement has been spent at the mark, and the unscaled dot
    # is at most (7+63-2*15)b=40b <32. The heavy mismatch case is smaller.
    assert (7+63-2*15)**2*b2<32**2 and -7+63-2*16<40
    out=[]
    for d in (63,64):
        W=set();axes=0;spikes=[]
        A=Z if d==63 else 0
        for ln in (H/f"win{d}.txt").read_text().splitlines():
            w,M=map(int,ln.split())
            assert 0<=w<=ONE and 0<M<=ONE and w not in W and M& A==0
            assert d==64 or (w&A).bit_count()==4
            v=w
            for r,j in zip(D,I):
                if v>>j&1:v^=r
            assert v==0
            for m in bits(M):
                s=w^(1<<m)
                # Ordinary unrepaired default even classes cannot block:
                # parent word is even on each original support, and flips
                # the heavy sign once. Test every nondefault phase exactly.
                assert all(any((x&s).bit_count()%2!=b for x,b in rows) for rows in R[m])
                assert d==64 or (s&A).bit_count()==4
                spikes.append((m,s))
            W.add(w)
            axes|=M
        assert len(spikes)==len(set(spikes))
        mx=-112
        for (a,x),(b,y) in combinations(spikes,2):
            z=x^y
            value=64-2*z.bit_count()
            if a==b:value+=48*(1-2*(z>>a&1))
            else:value+=6*(2-2*((z>>a&1)+(z>>b&1)))
            mx=max(mx,value)
        assert mx<=56
        removed=2*axes.bit_count()
        out.append(dict(dimension=d,previous=oldbounds[d],lower_bound=oldbounds[d]+len(spikes)-len(W)-removed,
                        gain=len(spikes)-len(W)-removed,inserted=len(spikes),removed_dense=len(W),
                        removed_poles=removed,max_scaled_dot=mx))
    later=load("dim63-67-gf8-norm-poles/claims.json","45a826ff2f34261540adf71e299c4c4c58a36c27")
    old={x['dimension']:x['lower_bound'] for x in later['dimensions']}
    for x in out:
        d=x['dimension']
        x.update(construction_base=x['previous'],construction_gain=x['gain'],
                 previous=old[d],gain=x['lower_bound']-old[d])
    assert out==json.loads((H/'claims.json').read_text())['dimensions']
    return dict(dimensions=out)

if __name__=="__main__":print(json.dumps(verify(),indent=2))
