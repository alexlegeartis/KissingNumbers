#!/usr/bin/env python3
"""Verify signed Netto quadruple trades against pinned 66/67D stocks."""
from pathlib import Path
from hashlib import sha1
import json, gzip, importlib.util
from itertools import combinations
from collections import defaultdict
from fractions import Fraction as Q
HERE=Path(__file__).resolve().parent
V=HERE.parent

def pin(rel, h):
    p=V / rel
    raw=p.read_bytes()
    assert sha1(f"blob {len(raw)}\0".encode()+raw).hexdigest()==h, rel
    return p

def load(rel,h):
    return json.loads(pin(rel,h).read_text())

def project_signs(mask, coords):
    return sum((mask>>j & 1)<<k for k,j in enumerate(coords))

def positions(m):
    return tuple(i for i in range(67) if m>>i&1)

def signed_triples(t, s):
    return [tuple((t[i],s>>i&1) for i in J) for J in combinations(range(4),3)]

def nine(d, T):
    p=f'dim{d}-{"tail-signed-nine" if d==66 else "signed-nine"}-refill/'
    hf={66:("a0ef843c1dc40f4c866d1a83c194dbeec47ec0f1","b1bb30b93a74b1a0d403782599140f244cd1e2d3"),
        67:("7b989108d5515e0c91bc17b7cecf76a49e6fc60d","1c85a14bd4c5d9ee1e227f23b7d22cb9d6a4252e")}
    S=load(p+"supports.json",hf[d][0])
    if d==66:
        n=load("dim61-66-nonary-affine-caps/certificate.json","366aab714f0ba186d57e19f784be52e27aee9f85")
        S=n["supports"]["66"]+S
    extra=json.loads(gzip.decompress(pin(p+"partial.json.gz",hf[d][1]).read_bytes()))
    phase=set()
    for seed in ('01112','03110','21310','21132'):
        for k in range(5):
            for a in (-1,1):
                p9=sum((0,1,3,2)[(a*int(seed[(j+k)%5]))%4]<<(2*j) for j in range(5))>>1
                phase.add(p9)
    assert len(phase)==40
    for m in S:
        J=positions(m)
        sponsor=[t for t in combinations(J,4) if t in T]
        assert len(sponsor)<=1
        if sponsor:
            t=sponsor[0]
            L=tuple(i for i in J if i not in t)
            J=(L[0],t[0],t[1],L[1],L[2],t[2],t[3],L[3],L[4])
        yield m, tuple(sum((a>>i&1)<<j for i,j in enumerate(J)) for a in phase)
    for m, P in extra:
        J=positions(m)
        yield m, tuple(sum((a>>i&1)<<j for i,j in enumerate(J)) for a in P)

def verify(d):
    modpath=pin("dim61-68-95-96-code-layers/source/continuation49_72/netto.py","604812c9d06ef60a851b587f4497c910c30200d9")
    spec=importlib.util.spec_from_file_location("old_netto",modpath)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    T=[t for t in mod.blocks68() if max(t)<d]
    assert len(T)==(11088 if d==66 else 11792)
    S=set(T)
    W=[]
    for line in (HERE/f"win{d}.txt").read_text().splitlines():
        a,b,c,e,s,*prices=map(int,line.split())
        t=(a,b,c,e)
        assert 0<=a<b<c<e<d and 0<=s<16
        assert t not in S
        W.append((t,s,tuple(prices)))
    assert len({(t,s) for t,s,b in W})==len(W)
    # Project all current nine-sign phases onto the old and changed quads only.
    lookup={sum(1<<i for i in t):t for t in S|{t for t,s,b in W}}
    forbidden=defaultdict(set)
    hits=set(); counts=0
    for m, signs in nine(d,S):
        assert m not in hits and all(a&m==a for a in signs)
        hits.add(m);counts+=len(signs)
        for t in combinations(positions(m),4):
            q=sum(1<<i for i in t)
            if q in lookup:forbidden[t].update(project_signs(a,t) for a in signs)
    assert counts==(86349 if d==66 else 82485)
    assert sum(len(forbidden[t]) for t in S)==(4823 if d==66 else 7110)
    triples=defaultdict(list)
    for i,t in enumerate(T):
        for s in range(16):
            if s in forbidden[t]:continue
            for j in signed_triples(t,s):triples[j].append(16*i+s)
    removed=set(); used={}
    for t,s,prices in W:
        assert s not in forbidden[t]
        U=signed_triples(t,s)
        for j in U:
            assert used.setdefault(j,t)==t
        actual={i for j in U for i in triples[j]}
        assert actual==set(prices)
        removed.update(actual)
    # Squares/norms in exact arithmetic; unmodified dense, 16-sites, poles,
    # and cap directions use the envelope checked by the pinned parent.
    assert 4*4*4==64 and 4*4*2==32 and 4*4<32 and 8*4==32
    assert 3*4*Q(8,3)==32 and 4*4*Q(8,3)>32
    assert 4*4*(3-0)==48 and 4*4*(3-2)==16
    # A cap can share at most one nonzero tail, replacing one ordinary entry.
    s,r=Q(61,80),Q(21,4)
    assert s*s>Q(37,64) and r*r>27
    assert max(16*s,12*s+4*r)<32
    base=load("dim65-67-bch-exterior-caps/claims.json","1dafcee005375381c2dad801cf50d942cfde1afc")
    prior=next(x["lower_bound"] for x in base["dimensions"] if x["dimension"]==d)
    gain=len(W)-len(removed)
    return dict(dimension=d,previous=prior,lower_bound=prior+gain,gain=gain,
                inserted=len(W),removed=len(removed))

if __name__=="__main__":
    assert __debug__
    print(json.dumps(dict(dimensions=[verify(d) for d in (66,67)]),indent=2))
