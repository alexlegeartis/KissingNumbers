from fractions import Fraction as F
from functools import reduce
from operator import xor
from itertools import combinations
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path

V = Path(__file__).resolve().parent.parent

def pin(p, digest):
    data = (V / p).read_bytes()
    assert hashlib.sha1(f"blob {len(data)}\0".encode()+data).hexdigest() == digest, p
    return V / p

def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

assert __debug__
pinned=load(pin("dim64-65-signed-parity/construction.py","5657a1072801fd2ded5e9fa02124976904e44e3e"),"pinned")
pinned_dg=load(pin("dim65-dg-dense-owners/verify.py","094053bd8c4ddf34546ae6c7deda4075d07fad90"),"pinned_dg")
base=json.loads(pin("dim64-65-partial-nonary-signs/claims.json","96a743ef12f41a0c144c7ca9c86cf2c52749f573").read_text())
prior65=json.loads(pin("dim65-dg-dense-owners/claims.json","2b42d9f2e2826f65becb9fefadc0d0cccee00f8b").read_text())
prior66_67=json.loads(pin("dim66-67-kerdock-axis-owners/claims.json","5ae3124788637df245a9ad497d1927ed4f59ea35").read_text())
assert base["dimensions"][0]["dimension"]==64 and prior65["dimension"]==65
assert prior65["sparse_owners"]==2688 and prior65["lower_bound"]==335929116

def elim(rows):
    d={}
    for r in rows:
        while r:
            j=r.bit_length()-1
            if j in d: r^=d[j]
            else: d[j]=r;break
    return d
def rem(v,d):
    for j in sorted(d,reverse=True):
        if v>>j&1:v^=d[j]
    return v
def span(rows):
    xs=[0]
    for r in rows: xs.extend([x^r for x in xs])
    return xs
a=[1]
for i in range(62):
    v=a[-1]<<1
    a.append(v^67 if v&64 else v)
assert len(set(a))==63 and set(a)==set(range(1,64)) and a[-1]==33
def ker(rows,n):
    d=elim(rows)
    out=[]
    for j in range(n):
        if j in d:continue
        v=1<<j
        for k,r in sorted(d.items()):
            if (v&r).bit_count()&1:v|=1<<k
        out.append(v)
    assert all((w&r).bit_count()%2==0 for r in rows for w in out)
    return out
def check_bch(n):
    return [(1<<64)-1]+[sum(((a[(i*j)%63]>>k)&1)<<j for j in range(63)) for i in range(1,n) for k in range(6)]
pinnedP,g=pinned.extended_bch(13)
P=ker(check_bch(13),64)
E=ker(check_bch(21),64)
Pb=elim(pinnedP)
assert len(P)==len(Pb)==30 and all(rem(w,Pb)==0 for w in P)
assert len(E)==18 and all(rem(w,Pb)==0 for w in E)
D=ker(check_bch(13)+P,64)
Db=elim(D)
assert len(D)==28 and all(w.bit_count()%4==0 and all((w&z).bit_count()%2==0 for z in D) for w in D)
shift=0x800000037cd0eb67
assert rem(shift,Pb)==0 and rem(shift,Db)!=0
Es=span(E)
h=Counter(w.bit_count() for w in Es)
Q=elim([rem(w,Db) for w in E])
overlap=1<<(len(E)-len(Q))
assert len(Es)==262144 and h[0]==1 and min(k for k in h if k)==22 and overlap==65536
# Independently rebuild signed-owner binary words from pinned quadratic mask definition.
coords=a+[0]
A=[(1<<64)-1]+[sum(((x>>i)&1)<<j for j,x in enumerate(coords)) for i in range(6)]
B=[A[i+1]&A[j+1] for i,j in combinations(range(6),2)]
assert all(rem(w,Db)==0 for w in A+B)
forms={pinned_dg.form(a,0) for a in range(32)}
lin=span(A)
sets=[]
for b in range(3):
    owners={shift^l^reduce(xor,(w for j,w in enumerate(B) if (t^b)>>j&1),0) for l in lin for t in forms}
    assert len(owners)==4096 and all(rem(w,Pb)==0 and rem(w,Db)!=0 for w in owners)
    sets.append(owners)
assert len(set.union(*sets))==3*4096
diff={a^b for a,b in combinations(forms,2)}
assert min((l^reduce(xor,(w for j,w in enumerate(B) if t>>j&1),0)).bit_count() for l in lin for t in diff)==28
# Height 4; min-distance-22 same-tail norm, and all cross constraints.
assert 36**2*F(3,4)<32**2
assert F(3,4)*(64-44)+16<32
assert F(3,4)*64-16==32 and 8*4==32
assert max(16*2,9*F(8,3),4*4,8)==32
assert F(3,4)*16+16<32
# Kerdock: norm-squared 64, same axis, diagonal or opposite, distinct axes.
s,r=F(61,80),F(21,4)
assert s*s>F(37,64) and r*r>27
assert F(37,64)*36**2<32**2
assert F(37,64)*36<32
assert 27+F(37,64)*8<32 and 37-27<32
assert max(32*s,28*s+2*r,16*s,12*s+4*r,F(8,3)*9*s,F(8,3)*(8*s+r),8*s)<32
res={65:base["dimensions"][0]["lower_bound"]+prior65["sparse_owners"]+2+len(Es)*2-overlap}
for row in prior66_67["dimensions"]:
    d=row["dimension"]
    restored=sum(len(s) for s in sets[:d-64])
    assert restored==row["owner_words"]
    res[d]=row["lower_bound"]+restored
assert res=={65:336256796,66:342389558,67:355195391}
claims=json.loads((Path(__file__).resolve().parent / "claims.json").read_text())["dimensions"]
previous={65:prior65["lower_bound"], **{row["dimension"]: row["lower_bound"] for row in prior66_67["dimensions"]}}
assert claims==[dict(dimension=d,previous=previous[d],lower_bound=v,gain=v-previous[d]) for d,v in sorted(res.items())]
print(res, sorted(h.items()),"overlap",overlap,"pinned generator",hex(g))
