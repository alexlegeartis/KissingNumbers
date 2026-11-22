#!/usr/bin/env python3
"""Exact omitted-folded-fibre spikes in 61D and 62D."""
from fractions import Fraction as Q
from pathlib import Path
from hashlib import sha1,sha256
from itertools import combinations
from functools import reduce
from operator import xor
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

def rem(w,p):
    for i in sorted(p,reverse=True):
        if w>>i&1:w^=p[i]
    return w

def ech(rows):
    p={}
    for w in rows:
        if v:=rem(w,p):p[v.bit_length()-1]=v
    return p

def dg61(G,K):
    p=pin("dim65-dg-dense-owners/verify.py","094053bd8c4ddf34546ae6c7deda4075d07fad90")
    a=importlib.util.spec_from_file_location("published_dg",p)
    dg=importlib.util.module_from_spec(a);a.loader.exec_module(dg)
    # A D0 coordinate labeling found by aligning the known dual-octad
    # geometry with affine F8^2. The inclusion is verified literally below.
    X=[23, 55, 61, 15, 31, 1, 37, 26, 4, 36, 9, 60, 62, 39, 17, 21, 29, 7, 53, 6, 2, 19, 32, 49, 20, 10, 24, 13, 5, 3, 18, 63, 42, 45, 35, 48, 47, 54, 16, 40, 25, 14, 38, 34, 41, 50, 43, 56, 27, 11, 59, 22, 30, 8, 33, 58, 51, 0, 52, 57, 44, 12, 28, 46]
    assert len(X)==64 and len(set(X))==64 and set(X)==set(range(64))
    A=[(1<<64)-1]+[sum(((x>>j)&1)<<i for i,x in enumerate(X)) for j in range(6)]
    B=[A[i+1]&A[j+1] for i,j in combinations(range(6),2)]
    R=ech(A+B)
    D=ech(G)
    assert len(R)==22 and all(not rem(x,D) for x in A+B)
    assert len(ech(prefix(x)for x in A))==3 and len(ech(prefix(x)for x in A+B))==4
    H=span(A)
    assert len(H)==128 and min(x.bit_count() for x in H[1:])==32
    forms={dg.form(a,b) for a in range(32) for b in range(32)}
    assert len(forms)==1024
    differences={a^b for a,b in combinations(forms,2)}
    assert len(differences)==15903
    def ev(f):return reduce(xor,(x for i,x in enumerate(B) if f>>i&1),0)
    assert min((ev(f)^t).bit_count() for f in differences for t in H)==24
    Z=[ev(f)^t for f in forms for t in H]
    assert len(Z)==len(set(Z))==131072
    labels=(4,5,6,7)
    Z=[x for x in Z if prefix(x) in labels]
    assert all(sum(prefix(x)==u for x in Z)==8192 for u in labels)
    # Prefix-zero representatives of K / (K intersect RM2). Lift each
    # residual using an actual K word so that translation preserves labels.
    extras=[]
    for k in K:
        if k2:=rem(k,R):
            R[k2.bit_length()-1]=k2
            extras.append(k)
    assert len(extras)==6
    shifts=span(extras)[:60]
    R=ech(A+B)
    assert len({rem(x,R) for x in shifts})==60 and all(prefix(x)==0 for x in shifts)
    return len(Z),24,shifts

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
    nonary=read("dim61-66-nonary-affine-caps/certificate.json","d2f7743669ca5ef39a6f7a58fdffc859bbc84550c815fcd89a6d958cb8b3f6de")
    parent=read("dim61-62-marked-caps/claims.json","ed18047c48d93edd60bb58cee3777b221b0d2267")
    witnesses={
      61:dict(g=(5688674,10221373,4888685,14304727,14003591,9115396,
                 14968292,4482293,8771124,6591576),minimum=24, label=6,
              r=Q(2,5),a=Q(11,2), bmax=Q(847,1000)),
      62:dict(g=(4061247,13630160,6265978,6806694,6310578,16012736,12298798,5833322,11818628,9950434,11630021,10183805,15204052,2505,1258),minimum=20,label=2,
              r=Q(0),a=Q(9,2),bmax=Q(899,1000))}
    out=[]
    for d,dat in witnesses.items():

        u=dat["label"]
        L=representative(P,u)
        occ=(0,1,2,3) if d==61 else (0,1,3,7,6,4)
        assert u not in occ
        # We use the nonary stock before the old marked capsules. All
        # 60 ordinary coordinates are consequently available.
        use=list(range(60))
        if d==61:
            per_axis,dist,chosen=dg61(G,K)
            assert not set((4,5,6,7)).intersection(occ)
            assert len(chosen)==60
        else:
            # An explicit binary rank-15 distance-20 subcode.
            Gs=dat["g"]
            J=span(Gs)
            assert len(J)==len(set(J))==1<<len(Gs)
            Z=span([choose_rows(K,z) for z in Gs])
            assert len(Z)==len(set(Z))==len(J)
            dist=min(w.bit_count() for w in Z[1:])
            assert dist>=dat["minimum"] and all(prefix(w)==0 for w in Z)
            per_axis=len(Z)
            chosen=[]
            reserved=set()
            for x in range(1<<24):
                if x in reserved:continue
                chosen.append(x)
                reserved.update(x^y for y in J)
                if len(chosen)==len(use):break
            assert len(chosen)==len(use)==60 and len(reserved)==len(chosen)*len(J)
        # Points copy the signs of one omitted parent codeword. Each
        # mark has its own quotient coset; each deleted label specifies
        # an exact latitude and its signed tail.
        if d==61:
            profiles=[(4,Q(23,4),Q(3,2),Q(819,1000)),
                      (5,Q(23,4),Q(-3,2),Q(819,1000)),
                      (6,Q(11,2),Q(2,5),Q(847,1000)),
                      (7,Q(11,2),Q(-2,5),Q(847,1000))]
            old=[(0,Q(1),Q(2)),(1,Q(1),Q(-2)),
                 (2,Q(17,16),Q(1,2)),(3,Q(17,16),Q(-1,2))]
        else:
            profiles=[(u,dat["a"],dat["r"],dat["bmax"])]
            old=[(v,Q(1),Q(0)) for v in occ]
        coords=[]
        def cap(scale2, coeff, tail):
            # If the body envelope is positive, remove its radical by
            # squaring against the nonnegative budget.
            if coeff>0:
                assert 32-tail>=0 and scale2*coeff*coeff<=(32-tail)**2
            else:
                assert tail<=32
        for u,a,r,B in profiles:
            assert u not in occ
            b2=(64-r*r)/(59+a*a)
            assert b2<B*B and b2*a*a>16 and b2*(59+a*a)+r*r==64
            # Every old W16 is supported only in ordinary coordinates,
            # and at least one sign disagrees by the dual-code argument.
            cap(b2,2*(a+13),0);cap(b2,32,0)
            # Old W4, W9 and retained axes. A support containing the
            # head and mark spends only the stated number of other slots.
            cap(b2,4*(a+3),0);cap(b2,4*(a+2),4*abs(r))
            cap(b2,Q(8,3)*(a+8),0);cap(b2,Q(8,3)*(a+7),Q(8,3)*abs(r))
            cap(b2,8,0);assert 8*abs(r)<=32
            # Dense owners of each old latitude have squared entry
            # magnitude c2; distance in the ordinary body is at least
            # 16 minus the difference on the four deleted positions.
            for v,c2,t in old:
                cap(b2*c2,a+27+2*(v^u).bit_count(),r*t)
            coords.append((u,a,r,b2))
        for i,(u,a,r,b2) in enumerate(coords):
            for v,c,t,f2 in coords[i:]:
                k=(u^v).bit_count()
                cap(b2*f2,a*c+59-2*(dist-k),r*t)
                cap(b2*f2,58+a+c-2*(16-k),r*t)
        # Reconstruct the certified stock immediately before the
        # distinguished-capsule exchanges, including its original
        # append-only nonary layer (every such point satisfies the above).
        prior=next(r0 for r0 in parent if r0["dimension"]==d)
        nb=len(nonary["supports"][str(d)])
        assert nb==(2331 if d==61 else 2891)
        aff=sum(x["count"] for x in h["new_classes"])
        assert aff==4416
        dense=1<<26 if d==61 else 6*(1<<24)
        four=8555 if d==61 else 9145
        stock_total=dense+(9673+937)*2048+aff+16*four+2*d+40*nb
        assert prior["previous"]==stock_total
        n=per_axis*len(use);gone=2*len(use)
        total=stock_total+n-gone
        out.append(dict(dimension=d,previous=prior["lower_bound"],construction_base=stock_total,
                        lower_bound=total,gain=total-prior["lower_bound"],inserted=n,
                        removed_poles=gone,points_per_axis=per_axis,minimum_hamming=dist,
                        omitted_labels=[4,5,6,7] if d==61 else [2]))
    return dict(dimensions=out)

if __name__=="__main__":print(json.dumps(verify(),indent=2))
