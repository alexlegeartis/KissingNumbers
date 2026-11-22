#!/usr/bin/env python3
"""Exact omitted-folded-fibre spikes in 61D and 62D."""
from fractions import Fraction as Q
from pathlib import Path
from hashlib import sha1,sha256
from itertools import combinations
import importlib.util,json
import numpy as np
HERE=Path(__file__).resolve().parent
V=HERE.parent
REM=(24,42,53,55)
mask=sum(1<<j for j in REM)
body=[i for i in range(64) if i not in REM]
BODY=((1<<64)-1)^mask

def pin(rel,h):
    p=V/rel;data=p.read_bytes()
    d=sha1(f"blob {len(data)}\0".encode()+data).hexdigest() if len(h)==40 else sha256(data).hexdigest()
    assert d==h,rel
    return p

def read(rel,h):
    return json.loads(pin(rel,h).read_text())

def span(rows):
    S=[0]
    for w in rows:S.extend([u^w for u in S])
    return S

def prefix(w):
    return sum((w>>i&1)<<j for j,i in enumerate(REM))

def split(rows):
    P={};K=[]
    for w in rows:
        while j:=w&mask:
            k=j.bit_length()-1
            if k not in P:
                P[k]=w
                break
            w^=P[k]
        else:
            assert w
            K.append(w)
    assert len(P)==4 and len(K)==24
    return P,K

def representative(P,j):
    z=sum((j>>k&1)<<i for k,i in enumerate(REM))
    y=0
    for k,w in sorted(P.items(),reverse=True):
        if (z^y)>>k&1:y^=w
    assert prefix(y)==j
    return y

def choose_rows(K,z):
    u=0
    for i,k in enumerate(K):
        if z>>i&1:u^=k
    return u

def verify():
    assert __debug__
    src=pin("dim61-67-compatible-sign-families/verify_family.py","222c034a630bfab0db9df1249797d0e633436bfc")
    s=importlib.util.spec_from_file_location("old",src)
    old=importlib.util.module_from_spec(s);s.loader.exec_module(old)
    G=old.generator();assert len(G)==28
    P,K=split(G)
    # Verify against every old unsigned W16 support; all associated phases
    # have odd parity. Body signs are from old D0 and are not flipped.
    p=pin("dim61-68-95-96-code-layers/source/dim61_96/code64_supports.npz","cdbd5e5a7bbda9ddb84ee86052e2d4d00ccc331f0188de958a69ed83b0c529f6")
    base={int(w) for w in np.load(p,allow_pickle=False)["words"] if int(w)&mask==0}
    assert len(base)==9673
    g=read("dim61-67-compatible-sign-families/data/short61_62.json","953a659d602217016254ec1e60e2ac46cf863eb1")
    h=read("dim61-62-66-affine-reconstruction/dim61-62/certificate.json","6fd52d4aae0b32d613369499c7c1e680f4b3c46e")
    assert g["removed_coordinates"]==list(REM) and g["odd_total_parity"] is True
    assert all(x["shift"].bit_count()%2==1 and 65535 in span(x["checks"])
               for x in g["new"]+g["old"]+h["new_classes"]+h["modified_parent_classes"])
    W=base|{x["support"] for x in g["new"]+h["new_classes"]}
    assert len(W)==9673+937+25
    assert all(x["support"] in base for x in g["old"]+h["modified_parent_classes"])
    assert all(w.bit_count()==16 and not (w&mask) and (w&a).bit_count()%2==0 for w in W for a in G)
    marks=read("dim61-62-marked-caps/certificate.json","ac51c49e4c555cc864be045a91ef34dad244cf7f")
    parent=read("dim61-62-marked-caps/claims.json","ed18047c48d93edd60bb58cee3777b221b0d2267")
    witnesses={
      61:dict(g=(5688674,10221373,4888685,14304727,14003591,9115396,
                 14968292,4482293,8771124,6591576),minimum=24, label=6,
              r=Q(2,5),a=Q(11,2), bmax=Q(847,1000)),
      62:dict(g=(4061247,13630160,6265978,6806694,6310578,16012736,12298798,5833322,11818628,9950434,11630021,10183805,15204052,2505,1258),minimum=20,label=2,
              r=Q(0),a=Q(9,2),bmax=Q(899,1000))}
    out=[]
    for d,dat in witnesses.items():
        # A sparse array of affine subcode generators, all in the actual
        # shortened kernel K. Verify all distances by enumerating one span.
        Gs=dat["g"]
        J=span(Gs)
        assert len(J)==len(set(J))==1<<len(Gs)
        Z=span([choose_rows(K,z) for z in Gs])
        assert len(Z)==len(set(Z))==len(J)
        dist=min(w.bit_count() for w in Z[1:])
        assert dist>=dat["minimum"] and all(prefix(w)==0 for w in Z)
        u=dat["label"]
        L=representative(P,u)
        occ=(0,1,2,3) if d==61 else (0,1,3,7,6,4)
        assert u not in occ
        banned={int(t,16).bit_length()-1 for t,_ in marks[str(d)]}
        assert len(marks[str(d)])==35 and len(banned)==27
        use=sorted(set(range(60))-banned)
        # Choose one distinct affine coset per retained mark in K/H.
        chosen=[]
        reserved=set()
        for x in range(1<<24):
            if x in reserved:continue
            chosen.append(x)
            reserved.update(x^y for y in J)
            if len(chosen)==len(use):break
        assert len(chosen)==len(use)==33 and len(reserved)==len(chosen)*len(J)
        assert all(0<=m<60 for m in use)
        # A point copies the signs of L + K[x] + Z on the physical
        # 60-coordinate body (no sign flips), with heavy factor a versus
        # regular 1. Its head is (r) in 61D or (0,0) in 62D.
        a,r,B=dat["a"],dat["r"],dat["bmax"]
        r2=r*r
        b2=(64-r2)/(59+a*a)
        assert b2<B*B and b2*a*a>16 and 59*b2+a*a*b2+r2==64
        # Old even-versus-odd weight sixteen ensures one flipped sign on
        # every overlap with the mark, and primary checked body has no tail.
        assert 2*(a+13)*B<32 and 32*B<32
        # Old retained signed W4, W9 and individual axes.
        assert 4*max((a+3)*B,(a+2)*B+r)<32
        assert Q(8,3)*max((a+8)*B,(a+7)*B+r)<32
        assert 8*max(B,r)<32
        # The marked capsules omit our mark as a heavy coordinate; even
        # if the mark is one of their twelve light coordinates, dot <= this.
        assert (26+2*a)*B<32
        # Distinct dense old words: original binary distance at least 16.
        # If all deleted differences number k, the body loses precisely k.
        # The two cases (agree, disagree on the heavy sign) are bounded by
        # the first, as a>1. A signed matching body cross has sum<=a+27+2k.
        raw=[a+27+2*(v^u).bit_count() for v in occ]
        if d==61:
            c=Q(33,32);assert c*c>Q(17,16)
            assert raw==[Q(73,2),Q(77,2),Q(69,2),Q(73,2)]
            assert max(raw[0]*B+2*r,raw[1]*B-2*r,
                       raw[2]*B*c+r/2,raw[3]*B*c-r/2)<32
        else:
            assert set(raw)=={Q(67,2),Q(71,2)}
            assert max(raw)*B<32
        # Exact new-pair inequalities. On the same mark the span's min
        # distance dist gives >=2 dist loss of scaled body product; if
        # heavy sign differs the loss is even greater. Different marks
        # have disjoint owner cosets, hence body distance >=16, since
        # every owner has the same deleted label. The heavy factors can
        # add at most 2(a-1) to the ordinary signed sum.
        assert (a*a+59-2*dist)*b2+r2<=32
        assert (60-2*16+2*(a-1))*b2+r2<=32
        n0=next(r0["lower_bound"] for r0 in parent if r0["dimension"]==d)
        n=len(Z)*len(use);gone=2*len(use)
        out.append(dict(dimension=d,previous=n0,lower_bound=n0+n-gone,
                        gain=n-gone,inserted=n,removed_poles=gone,
                        shortened_rank=len(Gs),minimum_hamming=dist,
                        omitted_label=u,weight_four_verification="amplitude envelope"))
    return dict(dimensions=out)

if __name__=="__main__":print(json.dumps(verify(),indent=2))
